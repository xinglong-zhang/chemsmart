from __future__ import annotations

from chemsmart.agent import providers


def _write_api_env(path, api_key: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"ai_api_key={api_key}\n", encoding="utf-8")
    return str(path)


def _configure_openai(monkeypatch):
    class DummyOpenAIProvider:
        def __init__(self, api_key: str) -> None:
            self.api_key = api_key

    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.delenv("ai_api_key", raising=False)
    monkeypatch.delenv("CHEMSMART_API_ENV", raising=False)
    monkeypatch.setattr(providers, "OpenAIProvider", DummyOpenAIProvider)
    return DummyOpenAIProvider


def test_get_provider_prefers_explicit_api_env_path(monkeypatch, tmp_path):
    dummy_provider = _configure_openai(monkeypatch)
    home_dir = tmp_path / "home"
    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.chdir(cwd_dir)

    explicit_env = _write_api_env(tmp_path / "explicit.env", "explicit-key")
    env_var_env = _write_api_env(tmp_path / "env-var.env", "env-var-key")
    _write_api_env(
        home_dir / ".chemsmart" / "api.env",
        "home-key",
    )
    _write_api_env(cwd_dir / "api.env", "cwd-key")
    monkeypatch.setenv("CHEMSMART_API_ENV", env_var_env)

    provider = providers.get_provider(env_path=explicit_env)

    assert isinstance(provider, dummy_provider)
    assert provider.api_key == "explicit-key"


def test_get_provider_prefers_chemsmart_api_env(monkeypatch, tmp_path):
    dummy_provider = _configure_openai(monkeypatch)
    home_dir = tmp_path / "home"
    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.chdir(cwd_dir)

    env_var_env = _write_api_env(tmp_path / "env-var.env", "env-var-key")
    _write_api_env(home_dir / ".chemsmart" / "api.env", "home-key")
    _write_api_env(cwd_dir / "api.env", "cwd-key")
    monkeypatch.setenv("CHEMSMART_API_ENV", env_var_env)

    provider = providers.get_provider()

    assert isinstance(provider, dummy_provider)
    assert provider.api_key == "env-var-key"


def test_get_provider_prefers_home_api_env_over_cwd(monkeypatch, tmp_path):
    dummy_provider = _configure_openai(monkeypatch)
    home_dir = tmp_path / "home"
    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.chdir(cwd_dir)

    _write_api_env(home_dir / ".chemsmart" / "api.env", "home-key")
    _write_api_env(cwd_dir / "api.env", "cwd-key")

    provider = providers.get_provider()

    assert isinstance(provider, dummy_provider)
    assert provider.api_key == "home-key"


def test_get_provider_uses_cwd_api_env_when_home_missing(
    monkeypatch, tmp_path
):
    dummy_provider = _configure_openai(monkeypatch)
    home_dir = tmp_path / "home"
    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.chdir(cwd_dir)

    _write_api_env(cwd_dir / "api.env", "cwd-key")

    provider = providers.get_provider()

    assert isinstance(provider, dummy_provider)
    assert provider.api_key == "cwd-key"


def test_get_provider_skips_dotenv_when_no_api_env_file_found(
    monkeypatch, tmp_path
):
    dummy_provider = _configure_openai(monkeypatch)
    home_dir = tmp_path / "home"
    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.chdir(cwd_dir)
    monkeypatch.setenv("ai_api_key", "shell-key")

    def fail_load_dotenv(*args, **kwargs):
        raise AssertionError(
            "load_dotenv should be skipped when no file exists"
        )

    monkeypatch.setattr(providers, "load_dotenv", fail_load_dotenv)

    provider = providers.get_provider()

    assert isinstance(provider, dummy_provider)
    assert provider.api_key == "shell-key"


def test_resolve_api_env_path_returns_none_when_no_candidate_exists(
    monkeypatch, tmp_path
):
    home_dir = tmp_path / "home"
    cwd_dir = tmp_path / "cwd"
    cwd_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.chdir(cwd_dir)
    monkeypatch.delenv("CHEMSMART_API_ENV", raising=False)

    assert providers._resolve_api_env_path(None) is None
