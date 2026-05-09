"""Local path helpers for project-side Python modules."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"


def ensure_external_paths() -> None:
    """Make the editable external package importable from project code."""

    if str(ADAPT_ENV_ROOT) not in sys.path:
        sys.path.insert(0, str(ADAPT_ENV_ROOT))
