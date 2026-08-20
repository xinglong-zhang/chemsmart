from __future__ import annotations

import socket
import ssl
import threading
import time
from contextlib import contextmanager
from http.client import BadStatusLine
from urllib.error import HTTPError
from urllib.request import Request

import pytest

from chemsmart.agent.runtime.alibaba import AlibabaTokenPlanHttpsTransport
from chemsmart.agent.runtime.deepseek import (
    DeepSeekHttpsTransport,
    DeepSeekTransportError,
)
from chemsmart.agent.runtime.transport import (
    ProviderDeadlineExceeded,
    ProviderTurnDeadline,
    ProviderTurnDeadlinesV1,
    _MonotonicHTTPSConnection,
    _perform_bounded_tls_handshake,
    open_bounded_https_response,
)


def test_response_header_drip_cannot_reset_monotonic_absolute_cap(monkeypatch):
    client, peer = socket.socketpair()
    peer_closed = threading.Event()

    def attach_socket(self):
        self.sock = client

    monkeypatch.setattr(_MonotonicHTTPSConnection, "connect", attach_socket)

    def drip_headers() -> None:
        try:
            request = b""
            while b"\r\n\r\n" not in request:
                request += peer.recv(4096)
            response = (
                b"HTTP/1.1 200 OK\r\n"
                + b"X-Drip: "
                + (b"a" * 100)
                + b"\r\nContent-Length: 0\r\nConnection: close\r\n\r\n"
            )
            for octet in response:
                peer.send(bytes((octet,)))
                time.sleep(0.025)
        except (BrokenPipeError, ConnectionResetError, OSError):
            pass
        finally:
            peer.close()
            peer_closed.set()

    worker = threading.Thread(
        target=drip_headers, name="synthetic-header-peer"
    )
    worker.start()
    deadline = ProviderTurnDeadline(
        ProviderTurnDeadlinesV1(
            connect_seconds=0.2,
            first_event_seconds=0.2,
            inter_event_seconds=0.2,
            absolute_seconds=0.2,
        ),
        turn_limit_seconds=0.1,
    )
    started = time.monotonic()

    with pytest.raises(ProviderDeadlineExceeded) as observed:
        with open_bounded_https_response(
            Request("https://fixed.example/chat/completions"),
            deadline=deadline,
        ):
            pytest.fail("dripped headers must not be accepted")

    elapsed = time.monotonic() - started
    worker.join(timeout=1.0)
    assert observed.value.phase == "absolute"
    assert observed.value.timeout_seconds == 0.1
    assert elapsed < 0.35
    assert peer_closed.is_set()
    assert not worker.is_alive()
    assert client.fileno() == -1


def test_resolved_tcp_addresses_share_one_monotonic_connect_cap(monkeypatch):
    import chemsmart.agent.runtime.transport as transport

    addresses = (
        (
            socket.AF_INET,
            socket.SOCK_STREAM,
            socket.IPPROTO_TCP,
            "",
            ("192.0.2.1", 443),
        ),
        (
            socket.AF_INET,
            socket.SOCK_STREAM,
            socket.IPPROTO_TCP,
            "",
            ("198.51.100.1", 443),
        ),
    )
    created = []

    class BlackHoledSocket:
        def __init__(self, *_args):
            self.timeout = 0.0
            self.closed = False
            self.close_calls = 0
            self.bound_address = None
            created.append(self)

        def settimeout(self, value):
            self.timeout = float(value)

        def connect(self, _address):
            # Each address remains responsive enough to evade a per-address
            # idle timeout. Only a recomputed shared allowance stops address 2.
            time.sleep(min(0.07, self.timeout))
            raise socket.timeout

        def bind(self, address):
            self.bound_address = address

        def close(self):
            self.closed = True
            self.close_calls += 1

    monkeypatch.setattr(
        transport.socket,
        "getaddrinfo",
        lambda *_args, **_kwargs: addresses,
    )
    monkeypatch.setattr(transport.socket, "socket", BlackHoledSocket)
    deadline = ProviderTurnDeadline(
        ProviderTurnDeadlinesV1(
            connect_seconds=0.1,
            first_event_seconds=0.1,
            inter_event_seconds=0.1,
            absolute_seconds=0.2,
        )
    )
    connection = _MonotonicHTTPSConnection(
        "fixed.example",
        port=443,
        deadline=deadline,
        source_address=("127.0.0.1", 0),
    )
    started = time.monotonic()

    with pytest.raises(ProviderDeadlineExceeded) as observed:
        connection.connect()

    elapsed = time.monotonic() - started
    assert observed.value.phase == "connect"
    assert elapsed < 0.13
    assert len(created) == 2
    assert created[1].timeout < created[0].timeout
    assert all(item.bound_address == ("127.0.0.1", 0) for item in created)
    assert all(item.closed and item.close_calls == 1 for item in created)
    assert connection.sock is None


def test_effective_turn_limit_is_reported_as_absolute_limiting_phase():
    class Clock:
        value = 0.0

        def __call__(self):
            return self.value

    clock = Clock()
    deadline = ProviderTurnDeadline(
        ProviderTurnDeadlinesV1(
            connect_seconds=5.0,
            first_event_seconds=5.0,
            inter_event_seconds=5.0,
            absolute_seconds=10.0,
        ),
        turn_limit_seconds=2.0,
        clock=clock,
    )
    clock.value = 0.1
    deadline.response_acquired()
    deadline.before_read(object())

    error = deadline.timeout_error()
    assert error.phase == "absolute"
    assert error.timeout_seconds == 2.0


def test_tls_handshake_wait_is_bounded_without_a_worker_thread():
    client, peer = socket.socketpair()

    class NeverCompletesHandshake:
        def setblocking(self, value):
            client.setblocking(value)

        def settimeout(self, value):
            client.settimeout(value)

        def fileno(self):
            return client.fileno()

        def do_handshake(self):
            raise ssl.SSLWantReadError()

    deadline = ProviderTurnDeadline(
        ProviderTurnDeadlinesV1(
            connect_seconds=0.05,
            first_event_seconds=0.05,
            inter_event_seconds=0.05,
            absolute_seconds=0.1,
        )
    )
    started = time.monotonic()
    try:
        with pytest.raises(ProviderDeadlineExceeded) as observed:
            _perform_bounded_tls_handshake(NeverCompletesHandshake(), deadline)
    finally:
        client.close()
        peer.close()

    assert observed.value.phase == "connect"
    assert time.monotonic() - started < 0.25
    assert not any(
        thread.name.startswith("provider") for thread in threading.enumerate()
    )


@pytest.mark.parametrize(
    ("transport_factory", "payload"),
    (
        (
            lambda: AlibabaTokenPlanHttpsTransport(api_key="sk-sp-test"),
            {"model": "deepseek-v4-flash-0731", "messages": []},
        ),
        (
            lambda: DeepSeekHttpsTransport(api_key="synthetic-lease"),
            {"model": "deepseek-v4-flash", "messages": []},
        ),
    ),
)
@pytest.mark.parametrize(
    "raw_error",
    (
        ConnectionResetError("raw peer authorization secret"),
        BadStatusLine("raw status authorization secret"),
        ssl.SSLError("raw TLS authorization secret"),
    ),
)
def test_raw_connection_failures_are_sanitized(
    monkeypatch, transport_factory, payload, raw_error
):
    import chemsmart.agent.runtime.alibaba as alibaba
    import chemsmart.agent.runtime.deepseek as deepseek

    def fail_open(*_args, **_kwargs):
        raise raw_error

    monkeypatch.setattr(alibaba, "open_bounded_https_response", fail_open)
    monkeypatch.setattr(deepseek, "open_bounded_https_response", fail_open)

    with pytest.raises(DeepSeekTransportError) as observed:
        transport_factory()(payload)

    assert observed.value.error_class == "transport"
    assert str(observed.value) == "transport"
    assert "authorization" not in str(observed.value)


def test_alibaba_429_is_classified_without_reading_dripped_error_body(
    monkeypatch,
):
    import chemsmart.agent.runtime.alibaba as alibaba

    class UnreadHTTPError(HTTPError):
        def __init__(self):
            super().__init__(
                "https://fixed.example", 429, "synthetic", {}, None
            )
            self.read_calls = 0

        def read(self, *_args, **_kwargs):
            self.read_calls += 1
            raise AssertionError("429 response body must remain unread")

    error = UnreadHTTPError()

    @contextmanager
    def fail_open(*_args, **_kwargs):
        raise error
        yield  # pragma: no cover

    monkeypatch.setattr(alibaba, "open_bounded_https_response", fail_open)
    transport = AlibabaTokenPlanHttpsTransport(api_key="sk-sp-test")

    with pytest.raises(DeepSeekTransportError) as observed:
        transport({"model": "deepseek-v4-flash-0731", "messages": []})

    assert observed.value.error_class == "rate_limited"
    assert observed.value.http_status == 429
    assert error.read_calls == 0
