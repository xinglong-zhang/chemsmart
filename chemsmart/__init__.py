from pathlib import Path

# Read version from VERSION file
with open(Path(__file__).parent / "VERSION", "r") as f:
    __version__ = f.read().strip()
