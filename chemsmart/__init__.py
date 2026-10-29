from pathlib import Path

# Define __version__ from the VERSION file shipped with the package
try:
    _version_file = Path(__file__).with_name("VERSION")
    __version__ = _version_file.read_text(encoding="utf-8").strip()
except (
    FileNotFoundError,
    PermissionError,
    UnicodeDecodeError,
):  # pragma: no cover
    __version__ = "0.0.0"
