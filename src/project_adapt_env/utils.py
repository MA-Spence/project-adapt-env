"""Small project-local utilities."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


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


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def atomic_write_text(path: Path, text: str, *, encoding: str = "utf-8") -> Path:
    ensure_parent(path)
    tmp_path = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    tmp_path.write_text(text, encoding=encoding)
    tmp_path.replace(path)
    return path


def atomic_write_json(path: Path, payload: Any, *, sort_keys: bool = True) -> Path:
    return atomic_write_text(
        path,
        json.dumps(to_builtin(payload), indent=2, sort_keys=sort_keys) + "\n",
        encoding="utf-8",
    )


def atomic_write_dataframe_csv(path: Path, frame: pd.DataFrame, *, index: bool = False) -> Path:
    ensure_parent(path)
    tmp_path = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    frame.to_csv(tmp_path, index=index)
    tmp_path.replace(path)
    return path


def format_duration(seconds: float | None) -> str | None:
    if seconds is None or not np.isfinite(seconds):
        return None
    total_seconds = max(int(round(float(seconds))), 0)
    hours, remainder = divmod(total_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def detect_worker_count(default: int = 1) -> int:
    for key in ("LABPROJ_CPUS", "SLURM_CPUS_PER_TASK"):
        value = os.environ.get(key)
        if value:
            try:
                parsed = int(value)
            except ValueError:
                continue
            if parsed > 0:
                return parsed
    cpu_count = os.cpu_count()
    if cpu_count is not None and cpu_count > 0:
        return int(cpu_count)
    return int(default)


def build_logger(name: str, log_path: Path) -> logging.Logger:
    ensure_parent(log_path)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter(
        fmt="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    return logger


@dataclass
class ProgressTracker:
    """Persist lightweight progress state with ETA estimates."""

    progress_path: Path
    logger: logging.Logger
    run_label: str = "run"
    run_started_at: float = field(default_factory=time.time)

    def write(
        self,
        *,
        stage: str,
        completed: int,
        total: int,
        message: str,
        details: dict[str, Any] | None = None,
        stage_started_at: float | None = None,
    ) -> None:
        details = details or {}
        elapsed_seconds = max(time.time() - float(stage_started_at or self.run_started_at), 0.0)
        rate = float(completed / elapsed_seconds) if elapsed_seconds > 0 and completed > 0 else None
        remaining = max(int(total) - int(completed), 0)
        eta_seconds = float(remaining / rate) if rate and rate > 0 else None
        payload = {
            "run_label": self.run_label,
            "stage": stage,
            "completed": int(completed),
            "total": int(total),
            "fraction_complete": float(completed / total) if total > 0 else 1.0,
            "elapsed_seconds": elapsed_seconds,
            "elapsed_hms": format_duration(elapsed_seconds),
            "eta_seconds": eta_seconds,
            "eta_hms": format_duration(eta_seconds),
            "updated_at_epoch": time.time(),
            "message": message,
            "details": to_builtin(details),
        }
        atomic_write_json(self.progress_path, payload)
        percent = 100.0 * payload["fraction_complete"]
        eta_text = payload["eta_hms"] or "unknown"
        self.logger.info(
            "[%s] %s: %d/%d (%.1f%%) ETA %s | %s",
            self.run_label,
            stage,
            int(completed),
            int(total),
            percent,
            eta_text,
            message,
        )
