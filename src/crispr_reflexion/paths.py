"""Where crispr-reflexion keeps downloaded data and cached API responses.

Two situations, two answers:

* **Installed with pip.** Data (a 396 MB atlas) and the API-response cache must not land in the
  Python install directory or in whatever folder the user happens to be in. They go to a per-user
  cache directory, overridable with `CRISPR_REFLEXION_HOME`.
* **Source checkout** (a `pyproject.toml` sits above the package). Keep using the project folder,
  so a developer's existing `data/` and `.ag2_cache/` are found and nothing is re-downloaded.

Reports and figures are different: they are the user's output, so they go to `./outputs` by default.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().parent
_CHECKOUT_ROOT = _PACKAGE_DIR.parents[1]  # <repo>/ when running from <repo>/src/crispr_reflexion


def is_source_checkout() -> bool:
    return (_CHECKOUT_ROOT / "pyproject.toml").exists() and (_CHECKOUT_ROOT / "src" / "crispr_reflexion").is_dir()


def _user_cache_root() -> Path:
    if sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = str(Path.home() / "Library" / "Caches")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "crispr-reflexion"


def home() -> Path:
    """Root directory for downloaded data and cached API responses."""
    override = os.environ.get("CRISPR_REFLEXION_HOME", "").strip()
    if override:
        return Path(override).expanduser()
    return _CHECKOUT_ROOT if is_source_checkout() else _user_cache_root()


def data_dir() -> Path:
    return home() / "data"


def api_cache_dir() -> Path:
    return home() / (".ag2_cache" if home() == _CHECKOUT_ROOT else "ag2_cache")
