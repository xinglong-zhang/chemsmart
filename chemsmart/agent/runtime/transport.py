"""Bounded synchronous reads for provider HTTP responses.

Socket timeouts are idle-operation bounds: a server that sends periodic
heartbeats can keep resetting them forever.  Runtime V2 additionally needs a
monotonic turn deadline that no byte, SSE comment, or reasoning token can
extend.  This module keeps that policy in the provider transport layer and
uses the calling thread only, so closing the response also closes the only
owned socket and no blocked worker can outlive the attempt.
"""

from __future__ import annotations

import errno
import http.client
import io
import math
import select
import socket
import ssl
import time
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Callable
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request

from chemsmart.agent._contracts import ContractError


@dataclass(frozen=True)
class ProviderTurnDeadlinesV1:
    """Independent connect, first-data, idle-data, and absolute turn bounds."""

    connect_seconds: float = 15.0
    first_event_seconds: float = 90.0
    inter_event_seconds: float = 90.0
    absolute_seconds: float = 300.0

    def __post_init__(self) -> None:
        values = (
            self.connect_seconds,
            self.first_event_seconds,
            self.inter_event_seconds,
            self.absolute_seconds,
        )
        if any(not math.isfinite(value) or value <= 0 for value in values):
            raise ContractError(
                "provider turn deadlines must be finite and positive"
            )
        if any(value > self.absolute_seconds for value in values[:-1]):
            raise ContractError(
                "provider connect and event deadlines cannot exceed the absolute turn deadline"
            )

    def public_record(
        self, *, effective_absolute_turn_seconds: float | None = None
    ) -> dict[str, float]:
        """Return non-secret transport policy for Runtime evidence."""

        effective = (
            self.absolute_seconds
            if effective_absolute_turn_seconds is None
            else min(
                self.absolute_seconds,
                max(0.0, float(effective_absolute_turn_seconds)),
            )
        )
        return {
            "connect_seconds": float(self.connect_seconds),
            "first_event_seconds": float(self.first_event_seconds),
            "inter_event_seconds": float(self.inter_event_seconds),
            "absolute_turn_seconds": float(self.absolute_seconds),
            "effective_absolute_turn_seconds": float(effective),
        }

    def configuration_record(self) -> dict[str, float]:
        """Return the four immutable values encoded by ``agent.yaml``."""

        return {
            "connect_seconds": float(self.connect_seconds),
            "first_event_seconds": float(self.first_event_seconds),
            "inter_event_seconds": float(self.inter_event_seconds),
            "absolute_turn_seconds": float(self.absolute_seconds),
        }


class ProviderDeadlineExceeded(TimeoutError):
    """Private controller signal converted to a sanitized transport error."""

    def __init__(self, phase: str, *, timeout_seconds: float) -> None:
        super().__init__(phase)
        self.phase = str(phase)
        self.timeout_seconds = max(0.0, float(timeout_seconds))


class ProviderTurnDeadline:
    """One monotonic deadline controller owned by one synchronous HTTP call."""

    def __init__(
        self,
        policy: ProviderTurnDeadlinesV1,
        *,
        turn_limit_seconds: float | None = None,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        limit = (
            policy.absolute_seconds
            if turn_limit_seconds is None
            else min(policy.absolute_seconds, float(turn_limit_seconds))
        )
        if not math.isfinite(limit) or limit <= 0:
            raise ContractError(
                "provider turn allowance must be finite and positive"
            )
        self.policy = policy
        self.clock = clock
        self.started_at = clock()
        self.absolute_seconds = limit
        self.absolute_deadline = self.started_at + limit
        self.connected_at: float | None = None
        self.last_event_at: float | None = None
        self.current_phase = "connect"
        self.current_timeout_seconds = min(
            self.policy.connect_seconds, self.absolute_seconds
        )

    def public_record(self) -> dict[str, float]:
        return self.policy.public_record(
            effective_absolute_turn_seconds=self.absolute_seconds
        )

    @property
    def connect_timeout_seconds(self) -> float:
        return self._remaining_for(
            "connect",
            self.started_at
            + min(self.policy.connect_seconds, self.absolute_seconds),
        )

    def before_response_io(self) -> float:
        """Return the remaining connect/TLS/status/header allowance.

        This is intentionally called before *each* underlying send/receive,
        rather than once before ``HTTPResponse.begin``.  A peer therefore
        cannot reset the hard cap by dripping header bytes.
        """

        if self.connected_at is not None:
            raise ContractError("provider response is already acquired")
        return self._remaining_for(
            "connect",
            self.started_at + self.policy.connect_seconds,
        )

    def connected(self) -> None:
        self._check_absolute()
        self.connected_at = self.clock()
        self.current_phase = "first_event"

    def response_acquired(self) -> None:
        """Mark connect/TLS/status/headers complete under the same hard cap."""

        self.connected()

    def before_read(self, response: Any) -> float:
        now = self.clock()
        self._check_absolute(now)
        if self.connected_at is None:
            raise ContractError("provider response was not marked connected")
        if self.last_event_at is None:
            phase = "first_event"
            phase_deadline = (
                self.connected_at + self.policy.first_event_seconds
            )
        else:
            phase = "inter_event"
            phase_deadline = (
                self.last_event_at + self.policy.inter_event_seconds
            )
        timeout = self._remaining_for(phase, phase_deadline, now=now)
        _set_response_socket_timeout(response, timeout)
        return timeout

    def event_observed(self) -> None:
        self._check_absolute()
        self.last_event_at = self.clock()
        self.current_phase = "inter_event"

    def check(self) -> None:
        self._check_absolute()

    def timeout_error(self) -> ProviderDeadlineExceeded:
        phase = self.current_phase
        return ProviderDeadlineExceeded(
            phase, timeout_seconds=self.current_timeout_seconds
        )

    def _remaining_for(
        self, phase: str, phase_deadline: float, *, now: float | None = None
    ) -> float:
        observed = self.clock() if now is None else now
        absolute_remaining = self.absolute_deadline - observed
        phase_remaining = phase_deadline - observed
        remaining = min(absolute_remaining, phase_remaining)
        selected_phase = (
            "absolute" if absolute_remaining <= phase_remaining else phase
        )
        selected_timeout = (
            self.absolute_seconds
            if selected_phase == "absolute"
            else {
                "connect": self.policy.connect_seconds,
                "first_event": self.policy.first_event_seconds,
                "inter_event": self.policy.inter_event_seconds,
            }[phase]
        )
        self.current_phase = selected_phase
        self.current_timeout_seconds = selected_timeout
        if remaining <= 0:
            raise ProviderDeadlineExceeded(
                selected_phase, timeout_seconds=selected_timeout
            )
        return remaining

    def _check_absolute(self, now: float | None = None) -> None:
        observed = self.clock() if now is None else now
        if observed >= self.absolute_deadline:
            self.current_phase = "absolute"
            self.current_timeout_seconds = self.absolute_seconds
            raise ProviderDeadlineExceeded(
                "absolute", timeout_seconds=self.absolute_seconds
            )


@contextmanager
def open_bounded_https_response(
    request: Request,
    *,
    deadline: ProviderTurnDeadline,
):
    """Open one direct HTTPS request with a bounded status/header parser.

    ``urllib`` passes one socket timeout into ``HTTPResponse.begin``. A peer
    can drip one header byte before each socket timeout and extend that method
    indefinitely. This owner-layer opener recomputes the monotonic allowance
    before every TLS operation and raw header receive. No worker thread or
    second request control plane is created; closing the response owns and
    closes the sole connection.

    The supported provider profiles bind fixed official HTTPS endpoints. This
    direct connection intentionally does not interpret process proxy settings;
    proxy traversal would require a separately bounded CONNECT/TLS authority.
    Platform DNS resolution itself is synchronous and cannot be cancelled by
    CPython without a worker; all socket work after resolution is hard-bounded.
    """

    parsed = urlsplit(request.full_url)
    if parsed.scheme.lower() != "https" or not parsed.hostname:
        raise ContractError("bounded provider transport requires HTTPS")
    port = parsed.port or 443
    path = parsed.path or "/"
    if parsed.query:
        path += "?" + parsed.query
    headers = {key: value for key, value in request.header_items()}
    headers.setdefault("Connection", "close")
    connection = _MonotonicHTTPSConnection(
        parsed.hostname,
        port=port,
        deadline=deadline,
    )
    response = None
    try:
        connection.request(
            request.get_method(),
            path,
            body=request.data,
            headers=headers,
        )
        response = connection.getresponse()
        deadline.response_acquired()
        if response.status >= 400:
            raise HTTPError(
                request.full_url,
                response.status,
                response.reason,
                response.headers,
                response,
            )
        yield response
    except ProviderDeadlineExceeded:
        raise
    except (socket.timeout, TimeoutError) as exc:
        raise deadline.timeout_error() from exc
    except HTTPError:
        raise
    except http.client.IncompleteRead:
        # DeepSeek gives this stable no-retry classification; Alibaba maps it
        # through its generic HTTPException boundary.
        raise
    except (
        http.client.HTTPException,
        ssl.SSLError,
        ConnectionError,
        OSError,
    ) as exc:
        raise URLError(OSError("provider transport failed")) from exc
    finally:
        if response is not None:
            response.close()
        connection.close()


class _MonotonicHTTPSConnection(http.client.HTTPSConnection):
    """HTTPSConnection whose TLS and response-header reads share one cap."""

    def __init__(
        self,
        host: str,
        *,
        port: int,
        deadline: ProviderTurnDeadline,
        source_address: tuple[str, int] | None = None,
    ) -> None:
        self._deadline = deadline
        super().__init__(
            host,
            port=port,
            timeout=deadline.connect_timeout_seconds,
            source_address=source_address,
        )
        self.response_class = (
            lambda sock, *args, **kwargs: _DeadlineHTTPResponse(
                sock,
                *args,
                deadline=self._deadline,
                **kwargs,
            )
        )

    def connect(self) -> None:
        if self._tunnel_host is not None:
            # The fixed provider opener never configures a proxy tunnel. A
            # future CONNECT path needs its own bounded header/TLS authority.
            raise OSError("bounded provider proxy tunnel is unsupported")
        self.sock = _open_bounded_tcp_socket(
            host=self.host,
            port=self.port,
            deadline=self._deadline,
            source_address=self.source_address,
        )
        server_hostname = self.host
        self.sock = self._context.wrap_socket(
            self.sock,
            server_hostname=server_hostname,
            do_handshake_on_connect=False,
        )
        _perform_bounded_tls_handshake(self.sock, self._deadline)

    def send(self, data: Any) -> None:
        if self.sock is None:
            self.connect()
        if self.sock is None:  # pragma: no cover - connect invariant
            raise ConnectionError("provider socket is unavailable")
        if not isinstance(data, (bytes, bytearray, memoryview)):
            raise TypeError("bounded provider request body must be bytes")
        _send_bounded(self.sock, data, self._deadline)

    def getresponse(self) -> http.client.HTTPResponse:
        if self.sock is None:
            raise http.client.ResponseNotReady()
        return super().getresponse()


class _DeadlineHTTPResponse(http.client.HTTPResponse):
    """HTTP response whose status/header file checks the cap per raw recv."""

    def __init__(
        self,
        sock: socket.socket,
        *args: Any,
        deadline: ProviderTurnDeadline,
        **kwargs: Any,
    ) -> None:
        super().__init__(
            _DeadlineSocketProxy(sock, deadline),
            *args,
            **kwargs,
        )


def _open_bounded_tcp_socket(
    *,
    host: str,
    port: int,
    deadline: ProviderTurnDeadline,
    source_address: tuple[str, int] | None = None,
) -> socket.socket:
    """Connect one resolved address at a time under one monotonic cap.

    ``socket.create_connection`` gives every address returned by
    ``getaddrinfo`` the same timeout. Multiple black-holed IPv4/IPv6 results
    can therefore multiply a nominal connect bound. Resolution remains the
    documented synchronous platform limitation; after it returns, each socket
    receives only the remaining shared allowance and every failure is closed.
    """

    deadline.before_response_io()
    addresses = socket.getaddrinfo(host, port, 0, socket.SOCK_STREAM)
    if not addresses:
        raise OSError("provider DNS returned no TCP addresses")
    for family, socktype, protocol, _canonical_name, address in addresses:
        candidate = None
        try:
            timeout = deadline.before_response_io()
            candidate = socket.socket(family, socktype, protocol)
            candidate.settimeout(max(timeout, 1e-6))
            if source_address is not None:
                candidate.bind(source_address)
            candidate.connect(address)
            try:
                candidate.setsockopt(
                    socket.IPPROTO_TCP,
                    socket.TCP_NODELAY,
                    1,
                )
            except OSError as exc:
                if exc.errno != errno.ENOPROTOOPT:
                    raise
            return candidate
        except ProviderDeadlineExceeded:
            if candidate is not None:
                candidate.close()
            raise
        except OSError:
            if candidate is not None:
                candidate.close()
            # An immediate refusal may leave enough time for another resolved
            # address. A timeout at the shared cap raises here before another
            # socket can be created.
            deadline.before_response_io()
        except BaseException:
            if candidate is not None:
                candidate.close()
            raise
    raise OSError("provider TCP connection failed")


class _DeadlineSocketProxy:
    """Delegate socket operations but own a deadline-aware ``makefile``."""

    def __init__(
        self, sock: socket.socket, deadline: ProviderTurnDeadline
    ) -> None:
        self._sock = sock
        self._deadline = deadline

    def makefile(self, mode: str, buffering: int = -1, **kwargs: Any):
        if mode != "rb":
            return self._sock.makefile(mode, buffering, **kwargs)
        raw = self._sock.makefile("rb", buffering=0)
        bounded = _DeadlineRawIO(raw, self._deadline)
        buffer_size = (
            io.DEFAULT_BUFFER_SIZE if buffering in {-1, 0} else buffering
        )
        return io.BufferedReader(bounded, buffer_size=buffer_size)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._sock, name)


class _DeadlineRawIO(io.RawIOBase):
    """Recompute the hard response-acquisition cap before every receive."""

    def __init__(self, raw: Any, deadline: ProviderTurnDeadline) -> None:
        super().__init__()
        self._raw = raw
        self._sock = getattr(raw, "_sock", None)
        self._deadline = deadline

    def readable(self) -> bool:
        return True

    def fileno(self) -> int:
        return self._raw.fileno()

    def readinto(self, target: Any) -> int | None:
        if self._deadline.connected_at is None:
            timeout = self._deadline.before_response_io()
            if self._sock is not None:
                self._sock.settimeout(max(timeout, 1e-6))
        return self._raw.readinto(target)

    def settimeout(self, value: float) -> None:
        if self._sock is not None:
            self._sock.settimeout(value)

    def close(self) -> None:
        if not self.closed:
            self._raw.close()
        super().close()


def _perform_bounded_tls_handshake(
    sock: ssl.SSLSocket,
    deadline: ProviderTurnDeadline,
) -> None:
    sock.setblocking(False)
    while True:
        deadline.before_response_io()
        try:
            sock.do_handshake()
            break
        except ssl.SSLWantReadError:
            _wait_bounded(sock, deadline, readable=True)
        except ssl.SSLWantWriteError:
            _wait_bounded(sock, deadline, readable=False)
    sock.settimeout(max(deadline.before_response_io(), 1e-6))


def _send_bounded(
    sock: socket.socket,
    data: bytes | bytearray | memoryview,
    deadline: ProviderTurnDeadline,
) -> None:
    view = memoryview(data)
    sock.setblocking(False)
    while view:
        deadline.before_response_io()
        try:
            sent = sock.send(view)
            if sent <= 0:
                raise ConnectionError("provider socket send failed")
            view = view[sent:]
        except (ssl.SSLWantWriteError, BlockingIOError):
            _wait_bounded(sock, deadline, readable=False)
        except ssl.SSLWantReadError:
            _wait_bounded(sock, deadline, readable=True)
    sock.settimeout(max(deadline.before_response_io(), 1e-6))


def _wait_bounded(
    sock: socket.socket,
    deadline: ProviderTurnDeadline,
    *,
    readable: bool,
) -> None:
    timeout = deadline.before_response_io()
    readers = [sock] if readable else []
    writers = [] if readable else [sock]
    ready_read, ready_write, _ = select.select(readers, writers, [], timeout)
    if not ready_read and not ready_write:
        # Recompute so the raised phase reflects the actual limiting clock.
        deadline.before_response_io()


def read_bounded_response_body(
    response: Any,
    *,
    deadline: ProviderTurnDeadline,
    chunk_size: int = 64 * 1024,
) -> bytes:
    """Read one non-streaming body without allowing progress to extend time."""

    chunks: list[bytes] = []
    reader = getattr(response, "read1", None)
    if reader is None:
        reader = response.read
    while True:
        deadline.before_read(response)
        try:
            chunk = reader(chunk_size)
        except TypeError:
            # Small in-memory test doubles commonly expose ``read()`` only.
            # Production urllib responses use ``read1(size)`` or ``read(size)``.
            chunk = reader()
        deadline.check()
        if not chunk:
            break
        if not isinstance(chunk, bytes):
            raise TypeError("provider response body must be bytes")
        chunks.append(chunk)
        deadline.event_observed()
    return b"".join(chunks)


def iter_bounded_response_lines(
    response: Any,
    *,
    deadline: ProviderTurnDeadline,
    chunk_size: int = 64 * 1024,
):
    """Yield complete byte lines while checking the hard cap per socket read.

    Using ``HTTPResponse.readline`` is insufficient for a hard deadline: a
    peer can drip bytes without a newline and keep each underlying receive
    below the socket idle timeout.  ``read1`` returns after one buffered/raw
    read, allowing the monotonic cap to be checked even for a partial line.
    Semantic event deadlines are advanced by the SSE parser, not by this byte
    iterator, so comments and incomplete data cannot count as model progress.
    """

    pending = bytearray()
    reader = getattr(response, "read1", None)
    if reader is None:
        reader = response.read
    while True:
        deadline.before_read(response)
        try:
            chunk = reader(chunk_size)
        except TypeError:
            chunk = reader()
        deadline.check()
        if not chunk:
            if pending:
                yield bytes(pending)
            return
        if not isinstance(chunk, bytes):
            raise TypeError("provider response body must be bytes")
        pending.extend(chunk)
        while True:
            newline = pending.find(b"\n")
            if newline < 0:
                break
            line = bytes(pending[: newline + 1])
            del pending[: newline + 1]
            yield line


def _set_response_socket_timeout(response: Any, timeout: float) -> None:
    """Set the current TLS socket read timeout on a urllib HTTP response.

    CPython's ``HTTPResponse`` exposes the owned socket through
    ``fp.raw._sock``.  The additional candidates keep this usable with focused
    transport doubles and compatible response wrappers.  An in-memory object
    may have no socket; its read is immediate and the monotonic post-read check
    still enforces synthetic tests.
    """

    candidates = [response]
    fp = getattr(response, "fp", None)
    raw = getattr(fp, "raw", None)
    candidates.extend((fp, raw, getattr(response, "raw", None)))
    candidates.extend(
        getattr(item, "_sock", None) for item in tuple(candidates) if item
    )
    for candidate in candidates:
        setter = getattr(candidate, "settimeout", None)
        if setter is not None:
            setter(max(float(timeout), 1e-6))
            return


def is_socket_timeout(error: BaseException) -> bool:
    """Return whether an exception denotes one bounded socket wait."""

    return isinstance(error, (TimeoutError, socket.timeout))


__all__ = [
    "ProviderDeadlineExceeded",
    "ProviderTurnDeadline",
    "ProviderTurnDeadlinesV1",
    "is_socket_timeout",
    "iter_bounded_response_lines",
    "open_bounded_https_response",
    "read_bounded_response_body",
]
