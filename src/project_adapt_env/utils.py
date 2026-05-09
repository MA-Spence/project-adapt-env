"""Small project-local utilities."""

from __future__ import annotations

import hashlib
import importlib.metadata
from pathlib import Path
from typing import Any

import numpy as np


def hash_mod(text: str, modulo: int) -> int:
    digest = hashlib.sha1(text.encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % modulo


def flatten_mapping(prefix: str, payload: dict[str, Any]) -> dict[str, Any]:
    return {f"{prefix}__{key}": value for key, value in payload.items()}


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def to_builtin(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, dict):
        return {str(key): to_builtin(item) for key, item in value.items()}
    if isinstance(value, list):
        return [to_builtin(item) for item in value]
    return value


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantile: float,
) -> float:
    x = np.asarray(values, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if x.size == 0:
        raise ValueError("weighted_quantile requires at least one value")
    if x.shape != w.shape:
        raise ValueError("values and weights must have the same shape")
    if not 0.0 <= quantile <= 1.0:
        raise ValueError("quantile must lie in [0, 1]")
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    total = float(np.sum(w))
    if total <= 0:
        raise ValueError("weights must sum to a positive value")
    cumulative = np.cumsum(w) / total
    return float(np.interp(quantile, cumulative, x))
