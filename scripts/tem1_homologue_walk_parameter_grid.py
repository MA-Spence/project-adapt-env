#!/usr/bin/env python3
"""Run TEM-1 homologue adaptive/random walk parameter-grid diagnostics."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
SCRIPT_DIR = Path(__file__).resolve().parent
for path in (ADAPT_ENV_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig, MMseqsServerClient  # noqa: E402
from adaptenv.amino_acids import AA_ORDER  # noqa: E402

from adaptive_walk_peak_epistasis_baselines import (  # noqa: E402
    evaluate_batch,
    finite_mean,
    json_default,
    merge_dict,
    single_mutant_scan,
    spearman,
    threshold_for,
    write_csv,
    write_json,
)


TEM1_SEQUENCE = (
    "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRP"
    "EERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVREL"
    "CSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTM"
    "PAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGS"
    "RGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
)


LITERATURE_BASELINES = [
    {
        "id": "tem1_stability_robustness_axis",
        "expectation": (
            "TEM-1 beta-lactamase and related marginally stable proteins should "
            "show a robustness/stability axis: homologues closer to the family "
            "consensus should be more buffered, while more consensus-distant "
            "backgrounds should expose more deleterious effects and stronger "
            "constraints."
        ),
        "sources": [
            "Stiffler et al. 2015, Cell, doi:10.1016/j.cell.2015.01.035",
            "Bershtein et al. 2006, Nature, doi:10.1038/nature04806",
            "Bloom et al. 2006, PNAS, doi:10.1073/pnas.0510098103",
        ],
        "metrics": [
            "alignment-consensus distance versus stability margin",
            "alignment-consensus distance versus neutral/deleterious fractions",
            "alignment-consensus distance versus negative and sign epistasis",
        ],
    },
    {
        "id": "adaptive_walk_diminishing_returns",
        "expectation": (
            "Adaptive walks should show early available gains and then decelerate "
            "or stall as local beneficial moves are exhausted."
        ),
        "sources": [
            "Wiser et al. 2013, Science, doi:10.1126/science.1243357",
            "Kryazhimskiy et al. 2014, Science, doi:10.1126/science.1250939",
        ],
        "metrics": [
            "adaptive-walk total gain and length by consensus distance",
            "late/early adaptive gain ratio",
            "stall probability and available beneficial fraction along adaptive walks",
        ],
    },
    {
        "id": "random_walk_ruggedness_autocorrelation",
        "expectation": (
            "Random walks on smoother landscapes should have higher fitness "
            "autocorrelation and smaller step variance; stronger epistasis should "
            "shorten correlation length and increase ruggedness proxies."
        ),
        "sources": [
            "Kauffman and Levin 1987, Journal of Theoretical Biology, doi:10.1016/S0022-5193(87)80029-2",
            "Stadler 1996, Journal of Mathematical Chemistry, doi:10.1007/BF01165157",
        ],
        "metrics": [
            "random-walk lag-1 and lag-5 fitness autocorrelation",
            "random-walk step-effect variance",
            "beneficial, neutral, and deleterious random-step fractions",
        ],
    },
    {
        "id": "directional_parameter_knobs",
        "expectation": (
            "Stability margin and epistasis strength should move their target "
            "diagnostics monotonically over realistic ranges."
        ),
        "sources": [
            "Martin and Lenormand 2008, Genetics, doi:10.1534/genetics.108.087122",
            "de Visser and Krug 2014, Nature Reviews Genetics, doi:10.1038/nrg3744",
        ],
        "metrics": [
            "stability_margin versus neutral/deleterious fractions",
            "epistasis_strength versus absolute, negative, and sign epistasis",
            "epistasis_strength versus random-walk autocorrelation and step variance",
        ],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def load_config(path: Path, *, quick: bool) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if quick:
        merged = merge_dict(raw, raw.get("quick", {}))
        merged.pop("quick", None)
        return merged
    raw.pop("quick", None)
    return raw


def seq_to_string(seq: np.ndarray) -> str:
    return "".join(AA_ORDER[int(aa)] for aa in np.asarray(seq, dtype=np.int64))


def write_alignment_fasta(path: Path, sequences: list[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for idx, sequence in enumerate(sequences):
            handle.write(f">tem1_mmseqs_{idx}\n")
            handle.write(f"{sequence}\n")


def synthetic_smoke_alignment(reference: str, *, n_sequences: int, seed: int) -> list[str]:
    """Small deterministic alignment for local smoke tests only."""
    rng = np.random.RandomState(seed)
    aa = np.asarray(list("ACDEFGHIKLMNPQRSTVWY"))
    rows = [reference]
    ref = np.asarray(list(reference))
    for idx in range(max(0, n_sequences - 1)):
        seq = ref.copy()
        rate = 0.01 + 0.18 * (idx / max(n_sequences - 1, 1))
        mask = rng.random_sample(ref.size) < rate
        for pos in np.flatnonzero(mask):
            choices = aa[aa != seq[pos]]
            seq[pos] = choices[int(rng.randint(choices.size))]
        rows.append("".join(seq.tolist()))
    return rows


def resolve_path(project_root: Path, config_path: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    config_relative = config_path.parent / path
    if config_relative.exists():
        return config_relative
    return project_root / path


def build_alignment_profile(
    cfg: dict[str, Any],
    *,
    project_root: Path,
    config_path: Path,
    output_dir: Path,
) -> AlignmentProfile:
    alignment_cfg = cfg["alignment"]
    reference = str(cfg["reference"]["sequence"])
    fasta_path = resolve_path(project_root, config_path, alignment_cfg.get("fasta_path"))
    max_sequences = int(alignment_cfg.get("max_sequences", 1000))
    if bool(alignment_cfg.get("synthetic_smoke_sequences", False)):
        sequences = synthetic_smoke_alignment(
            reference,
            n_sequences=max_sequences,
            seed=int(alignment_cfg.get("synthetic_smoke_seed", 12345)),
        )
        profile = AlignmentProfile.from_sequences(sequences, max_sequences=max_sequences)
    elif fasta_path is not None and fasta_path.exists():
        profile = AlignmentProfile.from_alignment(fasta_path, max_sequences=max_sequences)
    else:
        cache_dir = resolve_path(project_root, config_path, alignment_cfg.get("mmseqs_cache_dir"))
        client = MMseqsServerClient(
            cache_dir=cache_dir,
            poll_interval=float(alignment_cfg.get("poll_interval_seconds", 5.0)),
            timeout=float(alignment_cfg.get("timeout_seconds", 10.0)),
            max_retries=int(alignment_cfg.get("max_retries", 5)),
        )
        profile = AlignmentProfile.from_query_sequence(
            reference,
            mmseqs_client=client,
            max_sequences=max_sequences,
        )
    write_alignment_fasta(output_dir / "tem1_mmseqs_alignment.fasta", list(profile._raw_seqs))
    write_json(
        output_dir / "alignment_profile_summary.json",
        profile.to_dict(include_sequences=False, include_statistics=True),
    )
    return profile


def profile_log_odds_bits(profile: AlignmentProfile) -> np.ndarray:
    freq = np.maximum(np.asarray(profile.frequencies, dtype=np.float64), 1e-12)
    background = 1.0 / 20.0
    return np.log2(freq / background)


def consensus_metrics(
    seq: np.ndarray,
    *,
    consensus: np.ndarray,
    log_odds_bits: np.ndarray,
    consensus_score: float,
) -> dict[str, float | int]:
    seq = np.asarray(seq, dtype=np.int64)
    positions = np.arange(seq.size, dtype=np.int64)
    score = float(np.sum(log_odds_bits[positions, seq]))
    hamming = int(np.sum(seq != consensus))
    return {
        "consensus_bitscore": score,
        "consensus_distance_bits": float(consensus_score - score),
        "consensus_hamming_distance": hamming,
        "consensus_identity": float(1.0 - hamming / max(seq.size, 1)),
    }


def filled_unique_homologues(profile: AlignmentProfile) -> tuple[np.ndarray, list[int]]:
    matrix = np.asarray(profile.seq_matrix, dtype=np.int64)
    consensus = np.asarray(profile.consensus_idx, dtype=np.int64)
    filled = matrix.copy()
    gap_rows, gap_cols = np.where(filled < 0)
    if gap_rows.size:
        filled[gap_rows, gap_cols] = consensus[gap_cols]

    unique: list[np.ndarray] = []
    source_indices: list[int] = []
    seen: set[tuple[int, ...]] = set()
    for idx, row in enumerate(filled):
        key = tuple(int(x) for x in row)
        if key in seen:
            continue
        seen.add(key)
        unique.append(row.astype(np.int64, copy=True))
        source_indices.append(int(idx))
    if not unique:
        raise RuntimeError("No usable TEM-1 homologues were available after alignment filtering")
    return np.asarray(unique, dtype=np.int64), source_indices


def choose_evenly(indices: np.ndarray, count: int) -> np.ndarray:
    if indices.size <= count:
        return indices
    positions = np.linspace(0, indices.size - 1, count)
    chosen = np.unique(np.rint(positions).astype(np.int64))
    if chosen.size < count:
        missing = [i for i in range(indices.size) if i not in set(chosen)]
        chosen = np.asarray(list(chosen) + missing[: count - chosen.size], dtype=np.int64)
    return indices[np.sort(chosen[:count])]


def select_homologue_starts(
    profile: AlignmentProfile,
    cfg: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    seqs, source_indices = filled_unique_homologues(profile)
    consensus = np.asarray(profile.consensus_idx, dtype=np.int64)
    log_odds = profile_log_odds_bits(profile)
    positions = np.arange(profile.L, dtype=np.int64)
    consensus_score = float(np.sum(log_odds[positions, consensus]))

    rows: list[dict[str, Any]] = []
    for unique_idx, seq in enumerate(seqs):
        metrics = consensus_metrics(
            seq,
            consensus=consensus,
            log_odds_bits=log_odds,
            consensus_score=consensus_score,
        )
        rows.append({
            "unique_homologue_index": int(unique_idx),
            "source_alignment_index": int(source_indices[unique_idx]),
            "sequence": seq_to_string(seq),
            **metrics,
        })

    distances = np.asarray([row["consensus_distance_bits"] for row in rows], dtype=np.float64)
    order = np.argsort(distances, kind="mergesort")
    selected_rows: list[dict[str, Any]] = []
    starts_per_stratum = int(cfg["sampling"]["starts_per_stratum"])
    for stratum in cfg["sampling"]["start_strata"]:
        q_min = float(stratum["quantile_min"])
        q_max = float(stratum["quantile_max"])
        low = float(np.quantile(distances, q_min))
        high = float(np.quantile(distances, q_max))
        candidates = np.asarray(
            [idx for idx in order if low <= distances[idx] <= high],
            dtype=np.int64,
        )
        if candidates.size == 0:
            midpoint = 0.5 * (low + high)
            candidates = np.asarray(
                [int(np.argmin(np.abs(distances - midpoint)))],
                dtype=np.int64,
            )
        chosen = choose_evenly(candidates, starts_per_stratum)
        for replicate, idx in enumerate(chosen):
            row = dict(rows[int(idx)])
            row["start_id"] = f"{stratum['name']}__rep-{replicate}"
            row["start_stratum"] = str(stratum["name"])
            row["stratum_quantile_min"] = q_min
            row["stratum_quantile_max"] = q_max
            row["start_replicate"] = int(replicate)
            selected_rows.append(row)

    metadata = {
        "n_alignment_sequences": int(profile.N),
        "n_unique_homologues": int(len(rows)),
        "n_selected_starts": int(len(selected_rows)),
        "consensus_sequence": profile.consensus_seq,
        "consensus_bitscore": consensus_score,
        "distance_bits_min": float(np.min(distances)),
        "distance_bits_median": float(np.median(distances)),
        "distance_bits_max": float(np.max(distances)),
    }
    return selected_rows, metadata


def expand_parameter_grid(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    if "parameter_sets" in cfg:
        defaults = dict(cfg.get("parameter_defaults", {}))
        rows = []
        for idx, raw_row in enumerate(cfg["parameter_sets"]):
            row = dict(defaults)
            row.update(dict(raw_row))
            row.setdefault("parameter_set_id", f"pset-{idx:03d}")
            rows.append(add_ruggedness_score(row, cfg))
        return rows

    grid = cfg["parameter_grid"]
    keys = sorted(grid)
    rows = []
    for idx, values in enumerate(itertools.product(*(grid[key] for key in keys))):
        row = {key: values[pos] for pos, key in enumerate(keys)}
        row["parameter_set_id"] = f"pset-{idx:03d}"
        rows.append(add_ruggedness_score(row, cfg))
    return rows


def parameter_fields(parameter_grid: list[dict[str, Any]]) -> list[str]:
    preferred = [
        "stability_margin",
        "epistasis_strength",
        "n_epistatic_terms",
        "epistasis_covariation_bias",
        "pairwise_epistasis_rank",
        "n_higher_order_epistatic_terms",
        "epistasis_order_max",
        "epistasis_order_decay",
        "empirical_pairwise_strength",
        "generic_epistasis_target",
        "ruggedness_parameter_score",
    ]
    keys: set[str] = set()
    for row in parameter_grid:
        keys.update(str(key) for key in row if key != "parameter_set_id")
    ordered = [key for key in preferred if key in keys]
    ordered.extend(sorted(key for key in keys if key not in set(ordered)))
    return ordered


def parameter_value(row: dict[str, Any], cfg: dict[str, Any], key: str, default: Any = 0.0) -> Any:
    if key in row:
        return row[key]
    if key in cfg.get("adapt_env", {}).get("base_config", {}):
        return cfg["adapt_env"]["base_config"][key]
    return default


def add_ruggedness_score(row: dict[str, Any], cfg: dict[str, Any]) -> dict[str, Any]:
    """Expected perturbation scale used only to order ruggedness diagnostics."""
    out = dict(row)
    strength = float(parameter_value(out, cfg, "epistasis_strength", 0.0) or 0.0)
    n_pairs = max(float(parameter_value(out, cfg, "n_epistatic_terms", 0.0) or 0.0), 0.0)
    rank = max(float(parameter_value(out, cfg, "pairwise_epistasis_rank", 1.0) or 1.0), 1.0)
    covar = max(float(parameter_value(out, cfg, "epistasis_covariation_bias", 0.0) or 0.0), 0.0)
    higher = max(float(parameter_value(out, cfg, "n_higher_order_epistatic_terms", 0.0) or 0.0), 0.0)
    pairwise_scale = strength * math.sqrt(n_pairs / 320.0) * math.sqrt(rank / 2.0) * (1.0 + 0.5 * covar)
    higher_order_scale = 0.25 * strength * math.sqrt(higher / 1000.0) if higher > 0 else 0.0
    out["ruggedness_parameter_score"] = float(pairwise_scale + higher_order_scale)
    return out


def build_landscape(
    cfg: dict[str, Any],
    *,
    seed: int,
    profile: AlignmentProfile,
    parameter_set: dict[str, Any],
) -> FitnessLandscapeEnv:
    base = dict(cfg["adapt_env"]["base_config"])
    base.update({
        key: value
        for key, value in parameter_set.items()
        if key not in {"parameter_set_id", "ruggedness_parameter_score"}
    })
    base["L"] = len(str(cfg["reference"]["sequence"]))
    base["seed"] = int(seed)
    landscape_config = LandscapeConfig(**base)
    return FitnessLandscapeEnv(
        landscape_config,
        alignment_profile=profile,
        reference_sequence=str(cfg["reference"]["sequence"]),
    )


def components_for_seq(landscape: FitnessLandscapeEnv, seq: np.ndarray) -> dict[str, float]:
    components = landscape.evaluate(seq, return_components=True)
    return {
        "stability": float(components["stability"]),
        "effective_stability": float(components["effective_stability"]),
        "stability_margin_used": float(components["stability_margin_used"]),
        "effective_stability_margin_used": float(components["effective_stability_margin_used"]),
        "stability_gate": float(components["stability_gate"]),
        "epistatic_score": float(components["epistatic_score"]),
        "generic_pairwise_score": float(components["generic_pairwise_score"]),
        "empirical_pairwise_score": float(components["empirical_pairwise_score"]),
        "latent_readout": float(components["latent_readout"]),
    }


def classify_delta(delta: float, parent_fitness: float, sampling_cfg: dict[str, Any]) -> str:
    beneficial_threshold = threshold_for(parent_fitness, sampling_cfg, "beneficial_threshold")
    neutral_threshold = threshold_for(parent_fitness, sampling_cfg, "neutral_threshold")
    threshold = max(beneficial_threshold, neutral_threshold)
    if delta > threshold:
        return "beneficial"
    if delta < -threshold:
        return "deleterious"
    return "neutral"


def run_adaptive_walk(
    landscape: FitnessLandscapeEnv,
    start_seq: np.ndarray,
    *,
    mode: str,
    rng: np.random.RandomState,
    sampling_cfg: dict[str, Any],
    consensus: np.ndarray,
    log_odds_bits: np.ndarray,
    consensus_score: float,
) -> dict[str, Any]:
    seq = np.asarray(start_seq, dtype=np.int64).copy()
    trajectory: list[dict[str, Any]] = []
    previous_fitness: float | None = None
    stop_reason = "max_steps"
    max_steps = int(sampling_cfg["max_adaptive_walk_steps"])

    for step in range(max_steps + 1):
        scan = single_mutant_scan(landscape, seq, sampling_cfg=sampling_cfg)
        current_fitness = float(scan["parent_fitness"])
        beneficial_idx = np.flatnonzero(scan["beneficial"])
        components = components_for_seq(landscape, seq)
        row = {
            "step": int(step),
            "fitness": current_fitness,
            "step_gain": 0.0 if previous_fitness is None else current_fitness - previous_fitness,
            "available_beneficial_fraction": float(scan["fraction_beneficial"]),
            "available_neutral_fraction": float(scan["fraction_neutral"]),
            "available_deleterious_fraction": float(scan["fraction_deleterious"]),
            "available_mean_beneficial_effect": float(scan["mean_beneficial_effect"]),
            "best_available_effect": float(scan["best_delta"]),
            "chosen_effect": math.nan,
            **components,
            **consensus_metrics(
                seq,
                consensus=consensus,
                log_odds_bits=log_odds_bits,
                consensus_score=consensus_score,
            ),
        }
        trajectory.append(row)
        if step == max_steps:
            break
        if beneficial_idx.size == 0:
            stop_reason = "no_beneficial_single_mutant"
            break
        if mode == "greedy_adaptive":
            chosen = int(beneficial_idx[np.argmax(scan["deltas"][beneficial_idx])])
        elif mode == "random_beneficial_adaptive":
            chosen = int(beneficial_idx[int(rng.randint(beneficial_idx.size))])
        else:
            raise ValueError(f"Unsupported adaptive mode: {mode}")
        trajectory[-1]["chosen_effect"] = float(scan["deltas"][chosen])
        seq = scan["mutants"][chosen].copy()
        previous_fitness = current_fitness

    fitnesses = np.asarray([row["fitness"] for row in trajectory], dtype=np.float64)
    if fitnesses.size >= 4:
        midpoint = max(1, int((fitnesses.size - 1) // 2))
        early_gain = float(fitnesses[midpoint] - fitnesses[0])
        late_gain = float(fitnesses[-1] - fitnesses[midpoint])
        deceleration_ratio = late_gain / (early_gain + 1e-12) if early_gain > 0 else math.nan
    else:
        early_gain = math.nan
        late_gain = math.nan
        deceleration_ratio = math.nan
    return {
        "trajectory": trajectory,
        "final_sequence": seq,
        "summary": {
            "mode": mode,
            "n_steps": int(max(0, len(trajectory) - 1)),
            "start_fitness": float(fitnesses[0]),
            "final_fitness": float(fitnesses[-1]),
            "total_gain": float(fitnesses[-1] - fitnesses[0]),
            "early_gain": early_gain,
            "late_gain": late_gain,
            "deceleration_ratio": deceleration_ratio,
            "stop_reason": stop_reason,
        },
    }


def autocorrelation(values: np.ndarray, lag: int) -> float | None:
    if lag <= 0 or values.size <= lag + 2:
        return None
    x = values[:-lag]
    y = values[lag:]
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def run_random_walk(
    landscape: FitnessLandscapeEnv,
    start_seq: np.ndarray,
    *,
    rng: np.random.RandomState,
    sampling_cfg: dict[str, Any],
    consensus: np.ndarray,
    log_odds_bits: np.ndarray,
    consensus_score: float,
) -> dict[str, Any]:
    seq = np.asarray(start_seq, dtype=np.int64).copy()
    rows: list[dict[str, Any]] = []
    previous_fitness: float | None = None
    previous_parent_fitness: float | None = None
    step_classes: list[str] = []

    for step in range(int(sampling_cfg["random_walk_steps"]) + 1):
        fitness = float(landscape.evaluate(seq))
        step_gain = 0.0 if previous_fitness is None else fitness - previous_fitness
        step_class = (
            "start"
            if previous_parent_fitness is None
            else classify_delta(step_gain, previous_parent_fitness, sampling_cfg)
        )
        if step_class != "start":
            step_classes.append(step_class)
        components = components_for_seq(landscape, seq)
        rows.append({
            "step": int(step),
            "fitness": fitness,
            "step_gain": step_gain,
            "step_class": step_class,
            "available_beneficial_fraction": math.nan,
            "available_neutral_fraction": math.nan,
            "available_deleterious_fraction": math.nan,
            "available_mean_beneficial_effect": math.nan,
            "best_available_effect": math.nan,
            "chosen_effect": step_gain if step > 0 else math.nan,
            **components,
            **consensus_metrics(
                seq,
                consensus=consensus,
                log_odds_bits=log_odds_bits,
                consensus_score=consensus_score,
            ),
        })
        if step == int(sampling_cfg["random_walk_steps"]):
            break
        previous_fitness = fitness
        previous_parent_fitness = fitness
        pos = int(rng.randint(landscape.L))
        aa = int(rng.randint(landscape.A - 1))
        if aa >= int(seq[pos]):
            aa += 1
        seq[pos] = aa

    fitnesses = np.asarray([row["fitness"] for row in rows], dtype=np.float64)
    deltas = np.diff(fitnesses)
    lag1 = autocorrelation(fitnesses, 1)
    lag5 = autocorrelation(fitnesses, 5)
    correlation_length = (
        float(-1.0 / math.log(lag1))
        if lag1 is not None and 0.0 < lag1 < 1.0
        else math.nan
    )
    counts = {name: step_classes.count(name) for name in ("beneficial", "neutral", "deleterious")}
    denom = max(len(step_classes), 1)
    return {
        "trajectory": rows,
        "final_sequence": seq,
        "summary": {
            "mode": "unconditional_random",
            "n_steps": int(max(0, len(rows) - 1)),
            "start_fitness": float(fitnesses[0]),
            "final_fitness": float(fitnesses[-1]),
            "total_gain": float(fitnesses[-1] - fitnesses[0]),
            "early_gain": math.nan,
            "late_gain": math.nan,
            "deceleration_ratio": math.nan,
            "stop_reason": "completed_random_walk",
            "random_step_mean_delta": finite_mean(deltas),
            "random_step_sd_delta": float(np.std(deltas)) if deltas.size else math.nan,
            "random_lag1_autocorrelation": lag1,
            "random_lag5_autocorrelation": lag5,
            "random_correlation_length": correlation_length,
            "random_beneficial_step_fraction": float(counts["beneficial"] / denom),
            "random_neutral_step_fraction": float(counts["neutral"] / denom),
            "random_deleterious_step_fraction": float(counts["deleterious"] / denom),
        },
    }


def sample_epistasis_pairs(
    landscape: FitnessLandscapeEnv,
    background_seq: np.ndarray,
    *,
    background_id: str,
    row_context: dict[str, Any],
    rng: np.random.RandomState,
    sampling_cfg: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n_pairs = int(sampling_cfg["epistasis_pairs_per_background"])
    background_seq = np.asarray(background_seq, dtype=np.int64)
    background_fitness = float(landscape.evaluate(background_seq))
    threshold = max(
        threshold_for(background_fitness, sampling_cfg, "beneficial_threshold"),
        threshold_for(background_fitness, sampling_cfg, "neutral_threshold"),
        float(sampling_cfg.get("min_epistasis_effect_threshold", 0.0)),
    )

    pos_i = rng.randint(0, landscape.L, size=n_pairs)
    pos_j = rng.randint(0, landscape.L - 1, size=n_pairs)
    pos_j = pos_j + (pos_j >= pos_i)
    aa_i = rng.randint(0, landscape.A - 1, size=n_pairs)
    aa_j = rng.randint(0, landscape.A - 1, size=n_pairs)
    aa_i = aa_i + (aa_i >= background_seq[pos_i])
    aa_j = aa_j + (aa_j >= background_seq[pos_j])

    single_i = np.repeat(background_seq[None, :], n_pairs, axis=0)
    single_j = np.repeat(background_seq[None, :], n_pairs, axis=0)
    double = np.repeat(background_seq[None, :], n_pairs, axis=0)
    idx = np.arange(n_pairs, dtype=np.int64)
    single_i[idx, pos_i] = aa_i
    single_j[idx, pos_j] = aa_j
    double[idx, pos_i] = aa_i
    double[idx, pos_j] = aa_j

    fitness = evaluate_batch(
        landscape,
        np.vstack([single_i, single_j, double]),
        batch_size=int(sampling_cfg["batch_size"]),
    )
    fit_i = fitness[:n_pairs]
    fit_j = fitness[n_pairs:2 * n_pairs]
    fit_double = fitness[2 * n_pairs:]
    effect_i = fit_i - background_fitness
    effect_j = fit_j - background_fitness
    effect_i_on_j = fit_double - fit_j
    effect_j_on_i = fit_double - fit_i
    epistasis = fit_double - fit_i - fit_j + background_fitness

    rows: list[dict[str, Any]] = []
    flags: dict[str, list[bool]] = defaultdict(list)
    nonzero_eps = float(sampling_cfg.get("nonzero_epistasis_epsilon", 1e-9))
    for k in range(n_pairs):
        sign_i = 1 if effect_i[k] > threshold else (-1 if effect_i[k] < -threshold else 0)
        sign_j = 1 if effect_j[k] > threshold else (-1 if effect_j[k] < -threshold else 0)
        sign_i_on_j = 1 if effect_i_on_j[k] > threshold else (-1 if effect_i_on_j[k] < -threshold else 0)
        sign_j_on_i = 1 if effect_j_on_i[k] > threshold else (-1 if effect_j_on_i[k] < -threshold else 0)
        flip_i = sign_i != 0 and sign_i_on_j != 0 and sign_i != sign_i_on_j
        flip_j = sign_j != 0 and sign_j_on_i != 0 and sign_j != sign_j_on_i
        sign_epistasis = bool(flip_i or flip_j)
        reciprocal_sign = bool(flip_i and flip_j)
        negative = bool(epistasis[k] < -threshold)
        positive = bool(epistasis[k] > threshold)
        diminishing = bool((effect_i[k] > threshold or effect_j[k] > threshold) and negative)
        nonzero = bool(abs(epistasis[k]) > nonzero_eps)
        strong = bool(abs(epistasis[k]) > threshold)
        flags["sign"].append(sign_epistasis)
        flags["reciprocal"].append(reciprocal_sign)
        flags["negative"].append(negative)
        flags["positive"].append(positive)
        flags["diminishing"].append(diminishing)
        flags["nonzero"].append(nonzero)
        flags["strong"].append(strong)
        rows.append({
            **row_context,
            "background_id": background_id,
            "pair_index": int(k),
            "pos_i": int(pos_i[k]),
            "aa_i": int(aa_i[k]),
            "pos_j": int(pos_j[k]),
            "aa_j": int(aa_j[k]),
            "background_fitness": background_fitness,
            "fitness_i": float(fit_i[k]),
            "fitness_j": float(fit_j[k]),
            "fitness_double": float(fit_double[k]),
            "effect_i": float(effect_i[k]),
            "effect_j": float(effect_j[k]),
            "effect_i_on_j": float(effect_i_on_j[k]),
            "effect_j_on_i": float(effect_j_on_i[k]),
            "epistasis": float(epistasis[k]),
            "abs_epistasis": float(abs(epistasis[k])),
            "negative_epistasis": int(negative),
            "positive_epistasis": int(positive),
            "nonzero_epistasis": int(nonzero),
            "strong_epistasis": int(strong),
            "sign_epistasis": int(sign_epistasis),
            "reciprocal_sign_epistasis": int(reciprocal_sign),
            "diminishing_returns_epistasis": int(diminishing),
        })
    summary = {
        "epistasis_pair_count": int(n_pairs),
        "epistasis_mean": finite_mean(epistasis),
        "epistasis_abs_mean": finite_mean(np.abs(epistasis)),
        "epistasis_variance": float(np.var(epistasis)) if epistasis.size else math.nan,
        "negative_epistasis_fraction": float(np.mean(flags["negative"])) if flags["negative"] else math.nan,
        "positive_epistasis_fraction": float(np.mean(flags["positive"])) if flags["positive"] else math.nan,
        "nonzero_epistasis_fraction": float(np.mean(flags["nonzero"])) if flags["nonzero"] else math.nan,
        "strong_epistasis_fraction": float(np.mean(flags["strong"])) if flags["strong"] else math.nan,
        "sign_epistasis_fraction": float(np.mean(flags["sign"])) if flags["sign"] else math.nan,
        "reciprocal_sign_epistasis_fraction": float(np.mean(flags["reciprocal"])) if flags["reciprocal"] else math.nan,
        "diminishing_returns_epistasis_fraction": float(np.mean(flags["diminishing"])) if flags["diminishing"] else math.nan,
    }
    return rows, summary


def relation_test(
    rows: list[dict[str, Any]],
    *,
    x_key: str,
    y_key: str,
    direction: str,
    min_abs_spearman: float,
) -> dict[str, Any]:
    def as_float(value: Any) -> float:
        try:
            if value is None:
                return math.nan
            return float(value)
        except (TypeError, ValueError):
            return math.nan

    x = np.asarray([as_float(row.get(x_key)) for row in rows], dtype=np.float64)
    y = np.asarray([as_float(row.get(y_key)) for row in rows], dtype=np.float64)
    rho = spearman(x, y)
    mask = np.isfinite(x) & np.isfinite(y)
    low = high = None
    if np.count_nonzero(mask) >= 4:
        low_cut = float(np.quantile(x[mask], 1 / 3))
        high_cut = float(np.quantile(x[mask], 2 / 3))
        low = finite_mean(y[mask & (x <= low_cut)])
        high = finite_mean(y[mask & (x >= high_cut)])
    if direction == "increase":
        passed = rho is not None and rho >= min_abs_spearman and low is not None and high is not None and high > low
    elif direction == "decrease":
        passed = rho is not None and rho <= -min_abs_spearman and low is not None and high is not None and high < low
    else:
        raise ValueError(direction)
    return {
        "x": x_key,
        "metric": y_key,
        "direction": direction,
        "spearman": rho,
        "low_x_mean": low,
        "high_x_mean": high,
        "min_abs_spearman": float(min_abs_spearman),
        "pass": bool(passed),
    }


def group_means(rows: list[dict[str, Any]], group_keys: list[str], metric_keys: list[str]) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    group_meta: dict[tuple[Any, ...], dict[str, Any]] = {}
    for row in rows:
        key = tuple(row[k] for k in group_keys)
        group_meta.setdefault(key, {k: row[k] for k in group_keys})
        for metric in metric_keys:
            value = row.get(metric)
            try:
                value_float = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(value_float):
                grouped[key][metric].append(value_float)
    out: list[dict[str, Any]] = []
    for key, values_by_metric in grouped.items():
        row = dict(group_meta[key])
        row["n"] = int(max((len(v) for v in values_by_metric.values()), default=0))
        for metric in metric_keys:
            row[metric] = finite_mean(values_by_metric.get(metric, []))
        out.append(row)
    return out


def summarize_validation(
    background_rows: list[dict[str, Any]],
    walk_summaries: list[dict[str, Any]],
    parameter_summary_rows: list[dict[str, Any]],
    cfg: dict[str, Any],
) -> dict[str, Any]:
    min_abs = float(cfg["validation_thresholds"]["min_abs_spearman"])
    adaptive = [row for row in walk_summaries if row["mode"] == "greedy_adaptive"]
    random_rows = [row for row in walk_summaries if row["mode"] == "unconditional_random"]
    decel = np.asarray(
        [row["deceleration_ratio"] for row in adaptive if np.isfinite(row["deceleration_ratio"])],
        dtype=np.float64,
    )
    median_decel = float(np.median(decel)) if decel.size else None
    tests = {
        "fitness_declines_with_consensus_distance": relation_test(
            background_rows,
            x_key="consensus_distance_bits",
            y_key="background_fitness",
            direction="decrease",
            min_abs_spearman=min_abs,
        ),
        "beneficial_fraction_increases_with_consensus_distance": relation_test(
            background_rows,
            x_key="consensus_distance_bits",
            y_key="fraction_beneficial",
            direction="increase",
            min_abs_spearman=min_abs,
        ),
        "adaptive_total_gain_increases_with_consensus_distance": relation_test(
            adaptive,
            x_key="start_consensus_distance_bits",
            y_key="total_gain",
            direction="increase",
            min_abs_spearman=min_abs,
        ),
        "adaptive_walks_decelerate": {
            "metric": "deceleration_ratio",
            "median": median_decel,
            "expected_max": float(cfg["validation_thresholds"]["max_median_deceleration_ratio"]),
            "pass": bool(
                median_decel is not None
                and median_decel < float(cfg["validation_thresholds"]["max_median_deceleration_ratio"])
            ),
        },
        "random_walk_autocorrelation_declines_with_consensus_distance": relation_test(
            random_rows,
            x_key="start_consensus_distance_bits",
            y_key="random_lag1_autocorrelation",
            direction="decrease",
            min_abs_spearman=float(cfg["validation_thresholds"]["min_random_walk_autocorrelation_spearman"]),
        ),
    }

    def varied(key: str) -> bool:
        values = []
        for row in parameter_summary_rows:
            try:
                value = float(row.get(key))
            except (TypeError, ValueError):
                continue
            if np.isfinite(value):
                values.append(value)
        return len(set(values)) > 1

    parameter_tests: dict[str, dict[str, Any]] = {}
    if varied("stability_margin"):
        parameter_tests.update({
            "stability_margin_increases_neutral_fraction": relation_test(
                parameter_summary_rows,
                x_key="stability_margin",
                y_key="fraction_neutral",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
            "stability_margin_decreases_deleterious_fraction": relation_test(
                parameter_summary_rows,
                x_key="stability_margin",
                y_key="fraction_deleterious",
                direction="decrease",
                min_abs_spearman=min_abs,
            ),
        })
    if varied("epistasis_strength"):
        parameter_tests.update({
            "epistasis_strength_increases_abs_epistasis": relation_test(
                parameter_summary_rows,
                x_key="epistasis_strength",
                y_key="epistasis_abs_mean",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
            "epistasis_strength_increases_random_step_variance": relation_test(
                parameter_summary_rows,
                x_key="epistasis_strength",
                y_key="random_step_sd_delta",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
        })
    if varied("ruggedness_parameter_score"):
        parameter_tests.update({
            "ruggedness_score_increases_abs_epistasis": relation_test(
                parameter_summary_rows,
                x_key="ruggedness_parameter_score",
                y_key="epistasis_abs_mean",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
            "ruggedness_score_increases_nonzero_epistasis": relation_test(
                parameter_summary_rows,
                x_key="ruggedness_parameter_score",
                y_key="nonzero_epistasis_fraction",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
            "ruggedness_score_increases_epistasis_variance": relation_test(
                parameter_summary_rows,
                x_key="ruggedness_parameter_score",
                y_key="epistasis_variance",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
            "ruggedness_score_increases_random_step_variance": relation_test(
                parameter_summary_rows,
                x_key="ruggedness_parameter_score",
                y_key="random_step_sd_delta",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
            "ruggedness_score_decreases_random_lag1_autocorrelation": relation_test(
                parameter_summary_rows,
                x_key="ruggedness_parameter_score",
                y_key="random_lag1_autocorrelation",
                direction="decrease",
                min_abs_spearman=float(cfg["validation_thresholds"]["min_random_walk_autocorrelation_spearman"]),
            ),
            "ruggedness_score_increases_sign_epistasis": relation_test(
                parameter_summary_rows,
                x_key="ruggedness_parameter_score",
                y_key="sign_epistasis_fraction",
                direction="increase",
                min_abs_spearman=min_abs,
            ),
        })
    all_tests = {**tests, **parameter_tests}
    pass_count = sum(1 for test in all_tests.values() if bool(test.get("pass")))
    return {
        "tests": all_tests,
        "pass_count": int(pass_count),
        "test_count": int(len(all_tests)),
        "all_pass": bool(pass_count == len(all_tests)),
    }


def run_experiment(cfg: dict[str, Any], *, config_path: Path, project_root: Path, output_dir: Path) -> dict[str, Any]:
    profile = build_alignment_profile(
        cfg,
        project_root=project_root,
        config_path=config_path,
        output_dir=output_dir,
    )
    starts, start_metadata = select_homologue_starts(profile, cfg)
    consensus = np.asarray(profile.consensus_idx, dtype=np.int64)
    log_odds = profile_log_odds_bits(profile)
    consensus_score = float(np.sum(log_odds[np.arange(profile.L), consensus]))
    parameter_grid = expand_parameter_grid(cfg)
    param_fields = parameter_fields(parameter_grid)

    background_rows: list[dict[str, Any]] = []
    epistasis_rows: list[dict[str, Any]] = []
    walk_rows: list[dict[str, Any]] = []
    walk_summaries: list[dict[str, Any]] = []

    sampling_cfg = cfg["sampling"]
    for parameter_idx, parameter_set in enumerate(parameter_grid):
        for seed in sampling_cfg["seeds"]:
            landscape_seed = int(seed) + parameter_idx * 100_000
            landscape = build_landscape(
                cfg,
                seed=landscape_seed,
                profile=profile,
                parameter_set=parameter_set,
            )
            for start in starts:
                seq = np.asarray([AA_ORDER.index(aa) for aa in start["sequence"]], dtype=np.int64)
                background_id = f"{parameter_set['parameter_set_id']}__seed-{seed}__{start['start_id']}"
                row_context = {
                    "experiment_id": cfg["experiment_id"],
                    "parameter_set_id": parameter_set["parameter_set_id"],
                    **{field: parameter_set.get(field) for field in param_fields},
                    "seed": int(seed),
                    "landscape_seed": int(landscape_seed),
                    "start_id": start["start_id"],
                    "start_stratum": start["start_stratum"],
                    "start_replicate": int(start["start_replicate"]),
                }
                scan = single_mutant_scan(landscape, seq, sampling_cfg=sampling_cfg)
                components = components_for_seq(landscape, seq)
                pair_rows, pair_summary = sample_epistasis_pairs(
                    landscape,
                    seq,
                    background_id=background_id,
                    row_context=row_context,
                    rng=np.random.RandomState(landscape_seed + int(start["source_alignment_index"])),
                    sampling_cfg=sampling_cfg,
                )
                epistasis_rows.extend(pair_rows)
                background_rows.append({
                    **row_context,
                    "background_id": background_id,
                    "source_alignment_index": int(start["source_alignment_index"]),
                    "unique_homologue_index": int(start["unique_homologue_index"]),
                    "sequence": start["sequence"],
                    "background_fitness": float(scan["parent_fitness"]),
                    "fraction_beneficial": float(scan["fraction_beneficial"]),
                    "fraction_neutral": float(scan["fraction_neutral"]),
                    "fraction_deleterious": float(scan["fraction_deleterious"]),
                    "mean_delta": float(scan["mean_delta"]),
                    "best_delta": float(scan["best_delta"]),
                    "mean_beneficial_effect": float(scan["mean_beneficial_effect"]),
                    **{key: start[key] for key in (
                        "consensus_bitscore",
                        "consensus_distance_bits",
                        "consensus_hamming_distance",
                        "consensus_identity",
                    )},
                    **components,
                    **pair_summary,
                })

                modes = list(sampling_cfg["adaptive_walk_modes"]) + ["unconditional_random"]
                for mode_idx, mode in enumerate(modes):
                    walk_rng = np.random.RandomState(
                        landscape_seed + int(start["source_alignment_index"]) * 1009 + mode_idx * 7919
                    )
                    if mode == "unconditional_random":
                        walk = run_random_walk(
                            landscape,
                            seq,
                            rng=walk_rng,
                            sampling_cfg=sampling_cfg,
                            consensus=consensus,
                            log_odds_bits=log_odds,
                            consensus_score=consensus_score,
                        )
                    else:
                        walk = run_adaptive_walk(
                            landscape,
                            seq,
                            mode=str(mode),
                            rng=walk_rng,
                            sampling_cfg=sampling_cfg,
                            consensus=consensus,
                            log_odds_bits=log_odds,
                            consensus_score=consensus_score,
                        )
                    walk_summary = {
                        **row_context,
                        "background_id": background_id,
                        "mode": str(mode),
                        "start_consensus_bitscore": float(start["consensus_bitscore"]),
                        "start_consensus_distance_bits": float(start["consensus_distance_bits"]),
                        "start_consensus_hamming_distance": int(start["consensus_hamming_distance"]),
                        "start_consensus_identity": float(start["consensus_identity"]),
                        **walk["summary"],
                    }
                    walk_summaries.append(walk_summary)
                    start_fitness = float(walk["trajectory"][0]["fitness"])
                    for row in walk["trajectory"]:
                        walk_rows.append({
                            **row_context,
                            "background_id": background_id,
                            "mode": str(mode),
                            "start_consensus_distance_bits": float(start["consensus_distance_bits"]),
                            "start_consensus_hamming_distance": int(start["consensus_hamming_distance"]),
                            "step": int(row["step"]),
                            "fitness": float(row["fitness"]),
                            "cumulative_gain": float(row["fitness"] - start_fitness),
                            **row,
                        })

    parameter_summary_rows = group_means(
        background_rows + [
            row for row in walk_summaries if row.get("mode") == "unconditional_random"
        ],
        group_keys=["parameter_set_id", *param_fields],
        metric_keys=[
            "fraction_beneficial",
            "fraction_neutral",
            "fraction_deleterious",
            "mean_beneficial_effect",
            "epistasis_mean",
            "epistasis_abs_mean",
            "epistasis_variance",
            "negative_epistasis_fraction",
            "positive_epistasis_fraction",
            "nonzero_epistasis_fraction",
            "strong_epistasis_fraction",
            "sign_epistasis_fraction",
            "reciprocal_sign_epistasis_fraction",
            "diminishing_returns_epistasis_fraction",
            "random_step_mean_delta",
            "random_step_sd_delta",
            "random_lag1_autocorrelation",
            "random_lag5_autocorrelation",
            "random_correlation_length",
            "random_beneficial_step_fraction",
            "random_neutral_step_fraction",
            "random_deleterious_step_fraction",
        ],
    )
    validation = summarize_validation(background_rows, walk_summaries, parameter_summary_rows, cfg)
    return {
        "start_rows": starts,
        "start_metadata": start_metadata,
        "parameter_grid": parameter_grid,
        "background_rows": background_rows,
        "epistasis_rows": epistasis_rows,
        "walk_rows": walk_rows,
        "walk_summaries": walk_summaries,
        "parameter_summary_rows": parameter_summary_rows,
        "parameter_fields": param_fields,
        "summary": {
            "experiment_id": cfg["experiment_id"],
            "hypothesis_ids": ["HYP-004", "HYP-005"],
            "title": cfg["title"],
            "reference": cfg["reference"],
            "alignment": start_metadata,
            "parameter_fields": param_fields,
            "parameter_grid": parameter_grid,
            "n_parameter_sets": len(parameter_grid),
            "n_backgrounds": len(background_rows),
            "n_walk_summaries": len(walk_summaries),
            "n_walk_trajectory_rows": len(walk_rows),
            "n_epistasis_pairs": len(epistasis_rows),
            "literature_baselines": LITERATURE_BASELINES,
            "validation": validation,
        },
    }


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config_path, quick=bool(args.quick))
    payload = run_experiment(
        cfg,
        config_path=config_path,
        project_root=project_root,
        output_dir=output_dir,
    )

    write_json(output_dir / "literature_baselines.json", LITERATURE_BASELINES)
    write_json(output_dir / "summary.json", payload["summary"])
    write_json(output_dir / "walk_summaries.json", payload["walk_summaries"])
    write_csv(
        output_dir / "homologue_starts.csv",
        payload["start_rows"],
        [
            "start_id",
            "start_stratum",
            "start_replicate",
            "source_alignment_index",
            "unique_homologue_index",
            "consensus_bitscore",
            "consensus_distance_bits",
            "consensus_hamming_distance",
            "consensus_identity",
            "sequence",
        ],
    )
    write_csv(
        output_dir / "parameter_grid_summary.csv",
        payload["parameter_summary_rows"],
        [
            "parameter_set_id",
            *payload["parameter_fields"],
            "n",
            "fraction_beneficial",
            "fraction_neutral",
            "fraction_deleterious",
            "mean_beneficial_effect",
            "epistasis_mean",
            "epistasis_abs_mean",
            "epistasis_variance",
            "negative_epistasis_fraction",
            "positive_epistasis_fraction",
            "nonzero_epistasis_fraction",
            "strong_epistasis_fraction",
            "sign_epistasis_fraction",
            "reciprocal_sign_epistasis_fraction",
            "diminishing_returns_epistasis_fraction",
            "random_step_mean_delta",
            "random_step_sd_delta",
            "random_lag1_autocorrelation",
            "random_lag5_autocorrelation",
            "random_correlation_length",
            "random_beneficial_step_fraction",
            "random_neutral_step_fraction",
            "random_deleterious_step_fraction",
        ],
    )
    context_fields = [
        "experiment_id",
        "parameter_set_id",
        *payload["parameter_fields"],
        "seed",
        "landscape_seed",
        "start_id",
        "start_stratum",
        "start_replicate",
    ]
    background_fields = [
        *context_fields,
        "background_id",
        "source_alignment_index",
        "unique_homologue_index",
        "background_fitness",
        "consensus_bitscore",
        "consensus_distance_bits",
        "consensus_hamming_distance",
        "consensus_identity",
        "stability_margin_used",
        "effective_stability_margin_used",
        "stability_gate",
        "epistatic_score",
        "fraction_beneficial",
        "fraction_neutral",
        "fraction_deleterious",
        "mean_delta",
        "best_delta",
        "mean_beneficial_effect",
        "epistasis_pair_count",
        "epistasis_mean",
        "epistasis_abs_mean",
        "epistasis_variance",
        "negative_epistasis_fraction",
        "positive_epistasis_fraction",
        "nonzero_epistasis_fraction",
        "strong_epistasis_fraction",
        "sign_epistasis_fraction",
        "reciprocal_sign_epistasis_fraction",
        "diminishing_returns_epistasis_fraction",
    ]
    write_csv(output_dir / "background_neighborhoods.csv", payload["background_rows"], background_fields)
    walk_fields = [
        *context_fields,
        "background_id",
        "mode",
        "start_consensus_distance_bits",
        "start_consensus_hamming_distance",
        "step",
        "fitness",
        "cumulative_gain",
        "step_gain",
        "step_class",
        "available_beneficial_fraction",
        "available_neutral_fraction",
        "available_deleterious_fraction",
        "available_mean_beneficial_effect",
        "best_available_effect",
        "chosen_effect",
        "consensus_bitscore",
        "consensus_distance_bits",
        "consensus_hamming_distance",
        "consensus_identity",
        "stability_margin_used",
        "effective_stability_margin_used",
        "stability_gate",
        "epistatic_score",
    ]
    write_csv(output_dir / "walk_trajectories.csv", payload["walk_rows"], walk_fields)
    summary_fields = [
        *context_fields,
        "background_id",
        "mode",
        "start_consensus_bitscore",
        "start_consensus_distance_bits",
        "start_consensus_hamming_distance",
        "start_consensus_identity",
        "n_steps",
        "start_fitness",
        "final_fitness",
        "total_gain",
        "early_gain",
        "late_gain",
        "deceleration_ratio",
        "stop_reason",
        "random_step_mean_delta",
        "random_step_sd_delta",
        "random_lag1_autocorrelation",
        "random_lag5_autocorrelation",
        "random_correlation_length",
        "random_beneficial_step_fraction",
        "random_neutral_step_fraction",
        "random_deleterious_step_fraction",
    ]
    write_csv(output_dir / "walk_summaries.csv", payload["walk_summaries"], summary_fields)
    write_csv(
        output_dir / "epistasis_pairs.csv",
        payload["epistasis_rows"],
        [
            *context_fields,
            "background_id",
            "pair_index",
            "pos_i",
            "aa_i",
            "pos_j",
            "aa_j",
            "background_fitness",
            "fitness_i",
            "fitness_j",
            "fitness_double",
            "effect_i",
            "effect_j",
            "effect_i_on_j",
            "effect_j_on_i",
            "epistasis",
            "abs_epistasis",
            "negative_epistasis",
            "positive_epistasis",
            "nonzero_epistasis",
            "strong_epistasis",
            "sign_epistasis",
            "reciprocal_sign_epistasis",
            "diminishing_returns_epistasis",
        ],
    )
    print(json.dumps(payload["summary"]["validation"], indent=2, sort_keys=True, default=json_default))


if __name__ == "__main__":
    main()
