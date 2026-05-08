#!/usr/bin/env python3
"""Run EXP-001: ProteinGym DMS distributional realism panel."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib import request

import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
if str(ADAPT_ENV_ROOT) not in sys.path:
    sys.path.insert(0, str(ADAPT_ENV_ROOT))

from adaptenv import AlignmentProfile, FitnessLandscapeEnv, LandscapeConfig  # noqa: E402
from adaptenv.calibration import (  # noqa: E402
    CalibrationOptions,
    calibrate_synthetic_landscape,
    summarize_empirical_landscape,
)
from adaptenv.mmseqs import MMseqsServerClient  # noqa: E402


AA_SET = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class EmpiricalSequence:
    """Minimal sequence object compatible with adaptenv calibration helpers."""

    residues: tuple[str, ...]
    alphabet: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")
    name: str | None = None
    id: str | None = None

    @classmethod
    def from_string(
        cls,
        value: str,
        *,
        alphabet: tuple[str, ...] | list[str] | None = None,
        name: str | None = None,
        id: str | None = None,
    ) -> "EmpiricalSequence":
        return cls(
            residues=tuple(str(value)),
            alphabet=tuple(alphabet) if alphabet is not None else tuple("ACDEFGHIKLMNPQRSTVWY"),
            name=name,
            id=id,
        )

    @classmethod
    def from_iterable(
        cls,
        value: Any,
        *,
        alphabet: tuple[str, ...] | list[str] | None = None,
        name: str | None = None,
        id: str | None = None,
    ) -> "EmpiricalSequence":
        return cls(
            residues=tuple(str(symbol) for symbol in value),
            alphabet=tuple(alphabet) if alphabet is not None else tuple("ACDEFGHIKLMNPQRSTVWY"),
            name=name,
            id=id,
        )

    def to_array(self) -> np.ndarray:
        return np.asarray(self.residues, dtype=object)

    def to_str(self) -> str:
        return "".join(self.residues)

    def mutate(self, positions: Any, values: Any) -> "EmpiricalSequence":
        residues = list(self.residues)
        if np.isscalar(positions):
            pos_list = [int(positions)]
            value_list = [str(values)]
        else:
            pos_list = [int(pos) for pos in positions]
            value_list = [str(value) for value in values]
        for pos, value in zip(pos_list, value_list):
            residues[pos] = value
        return EmpiricalSequence(
            residues=tuple(residues),
            alphabet=self.alphabet,
        )


class ProteinGymEmpiricalLandscape:
    """Adapter exposing a ProteinGym assay to adaptenv calibration code."""

    def __init__(
        self,
        *,
        dms_id: str,
        wildtype_sequence: str,
        assay_frame: pd.DataFrame,
        score_column: str,
    ) -> None:
        self.name = dms_id
        self.wildtype = EmpiricalSequence.from_string(
            wildtype_sequence,
            name="wt",
            id="wt",
        )
        self.wildtype_sequence = self.wildtype
        self.sequences: list[EmpiricalSequence] = [self.wildtype]
        self._fitness_by_sequence: dict[str, float] = {wildtype_sequence: 0.0}

        for row in assay_frame.itertuples(index=False):
            sequence = str(row.mutated_sequence)
            score = float(getattr(row, score_column))
            seq_obj = EmpiricalSequence.from_string(sequence)
            self.sequences.append(seq_obj)
            self._fitness_by_sequence[sequence] = score

    def get_fitness(self, sequence: Any) -> float:
        if hasattr(sequence, "to_str"):
            key = str(sequence.to_str())
        else:
            key = "".join(str(symbol) for symbol in np.asarray(sequence, dtype=object))
        if key not in self._fitness_by_sequence:
            raise KeyError(f"Sequence not present in empirical assay: {key}")
        return float(self._fitness_by_sequence[key])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    return parser.parse_args()


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def download_file(url: str, destination: Path) -> Path:
    ensure_parent(destination)
    if destination.is_file():
        return destination
    tmp_path = destination.with_suffix(destination.suffix + ".tmp")
    with request.urlopen(url, timeout=120) as response:
        tmp_path.write_bytes(response.read())
    tmp_path.replace(destination)
    return destination


def canonical_sequence(sequence: str) -> str:
    cleaned = "".join(str(sequence).split()).upper()
    if not cleaned:
        raise ValueError("Encountered empty sequence")
    bad = sorted(set(cleaned) - AA_SET)
    if bad:
        raise ValueError(f"Sequence contains unsupported residues: {''.join(bad)}")
    return cleaned


def normalize_includes_multiple(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"true", "1", "yes"}


def select_assay_panel(reference_df: pd.DataFrame, panel_cfg: dict[str, Any]) -> pd.DataFrame:
    override_ids = [str(item) for item in panel_cfg.get("assay_ids", [])]
    if override_ids:
        selected = reference_df[reference_df["DMS_id"].isin(override_ids)].copy()
        if selected.empty:
            raise ValueError("No assay_ids from config matched the ProteinGym reference table")
        return selected.reset_index(drop=True)

    allowed_taxa = [str(item) for item in panel_cfg["allowed_taxa"]]
    assays_per_taxon = int(panel_cfg["assays_per_taxon"])
    panel_size = int(panel_cfg["panel_size"])
    min_single_mutants = int(panel_cfg["min_single_mutants"])
    max_sequence_length = int(panel_cfg["max_sequence_length"])

    filtered = reference_df.copy()
    filtered = filtered[filtered["taxon"].isin(allowed_taxa)]
    filtered = filtered[~filtered["includes_multiple_mutants"].map(normalize_includes_multiple)]
    filtered = filtered[
        filtered["DMS_total_number_mutants"].astype(int)
        == filtered["DMS_number_single_mutants"].astype(int)
    ]
    filtered = filtered[filtered["DMS_number_single_mutants"].astype(int) >= min_single_mutants]
    filtered = filtered[filtered["seq_len"].astype(int) <= max_sequence_length]
    filtered = filtered[filtered["target_seq"].notna()]
    filtered = filtered[filtered["target_seq"].map(lambda seq: set(str(seq).upper()) <= AA_SET)]
    filtered = filtered.sort_values(
        by=["taxon", "MSA_N_eff", "DMS_number_single_mutants", "seq_len"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )

    chosen_frames: list[pd.DataFrame] = []
    chosen_ids: set[str] = set()
    for taxon in allowed_taxa:
        subset = filtered[filtered["taxon"] == taxon]
        if subset.empty:
            continue
        take = subset.head(assays_per_taxon).copy()
        chosen_frames.append(take)
        chosen_ids.update(str(item) for item in take["DMS_id"])

    if chosen_frames:
        panel = pd.concat(chosen_frames, ignore_index=True)
    else:
        panel = filtered.iloc[0:0].copy()

    if len(panel) < panel_size:
        extras = filtered[~filtered["DMS_id"].isin(chosen_ids)].head(panel_size - len(panel))
        if not extras.empty:
            panel = pd.concat([panel, extras], ignore_index=True)

    panel = panel.head(panel_size).reset_index(drop=True)
    if panel.empty:
        raise ValueError("Panel selection produced no assays")
    return panel


def zscore(values: pd.Series) -> pd.Series:
    mean = float(values.mean())
    std = float(values.std(ddof=0))
    if not np.isfinite(std) or std <= 1e-12:
        return pd.Series(np.zeros(len(values), dtype=np.float64), index=values.index)
    return (values - mean) / std


def mutation_positions(wildtype: str, mutated_sequences: pd.Series) -> list[int]:
    positions: list[int] = []
    wt_chars = np.asarray(list(wildtype), dtype=object)
    for sequence in mutated_sequences.astype(str):
        seq_chars = np.asarray(list(sequence), dtype=object)
        diff = np.flatnonzero(seq_chars != wt_chars)
        if diff.size != 1:
            raise ValueError(
                f"Expected a single-substitution assay, but observed {diff.size} differences"
            )
        positions.append(int(diff[0]))
    return positions


def write_fasta(path: Path, sequences: list[str], *, prefix: str) -> None:
    ensure_parent(path)
    with path.open("w", encoding="utf-8") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f">{prefix}_{index:04d}\n")
            handle.write(f"{sequence}\n")


def build_alignment_fasta(
    *,
    dms_id: str,
    wildtype_sequence: str,
    fasta_path: Path,
    mmseqs_cfg: dict[str, Any],
) -> Path:
    if fasta_path.is_file():
        return fasta_path
    client = MMseqsServerClient(
        cache_dir=Path(mmseqs_cfg["cache_dir"]),
        timeout=float(mmseqs_cfg["timeout_seconds"]),
        poll_interval=float(mmseqs_cfg["poll_interval_seconds"]),
        max_retries=int(mmseqs_cfg["max_retries"]),
    )
    sequences = client.fetch_alignment(wildtype_sequence)[: int(mmseqs_cfg["max_alignment_sequences"])]
    write_fasta(fasta_path, sequences, prefix=dms_id)
    return fasta_path


def assay_score_summary(scores: np.ndarray) -> dict[str, float]:
    return {
        "count": float(scores.size),
        "mean": float(np.mean(scores)),
        "std": float(np.std(scores)),
        "skewness": float(pd.Series(scores).skew()),
        "q05": float(np.quantile(scores, 0.05)),
        "q25": float(np.quantile(scores, 0.25)),
        "q50": float(np.quantile(scores, 0.50)),
        "q75": float(np.quantile(scores, 0.75)),
        "q95": float(np.quantile(scores, 0.95)),
        "fraction_gt_1sd": float(np.mean(scores > 1.0)),
        "fraction_lt_minus_1sd": float(np.mean(scores < -1.0)),
    }


def empirical_position_sensitivity(
    assay_frame: pd.DataFrame,
    *,
    sequence_length: int,
    score_column: str,
) -> np.ndarray:
    by_position = np.full(sequence_length, np.nan, dtype=np.float64)
    for position, subset in assay_frame.groupby("mutation_position"):
        by_position[int(position)] = float(np.mean(np.abs(subset[score_column].to_numpy(dtype=np.float64))))
    return by_position


def correlation_or_none(x: np.ndarray, y: np.ndarray) -> float | None:
    mask = np.isfinite(x) & np.isfinite(y)
    if int(np.sum(mask)) < 2:
        return None
    x_masked = x[mask]
    y_masked = y[mask]
    if np.ptp(x_masked) <= 1e-12 or np.ptp(y_masked) <= 1e-12:
        return None
    return float(np.corrcoef(x_masked, y_masked)[0, 1])


def compute_empirical_metrics(
    assay_frame: pd.DataFrame,
    *,
    wildtype_sequence: str,
    alignment_profile: AlignmentProfile,
    score_column: str,
) -> dict[str, Any]:
    scores = assay_frame[score_column].to_numpy(dtype=np.float64)
    position_sensitivity = empirical_position_sensitivity(
        assay_frame,
        sequence_length=len(wildtype_sequence),
        score_column=score_column,
    )
    conservation_corr = correlation_or_none(alignment_profile.conservation, position_sensitivity)
    return {
        "score_summary": assay_score_summary(scores),
        "fraction_positive_bin": float(np.mean(assay_frame["DMS_score_bin"].to_numpy(dtype=np.float64) > 0.5))
            if "DMS_score_bin" in assay_frame
            else None,
        "conservation_sensitivity_correlation": conservation_corr,
        "mutation_coverage_fraction": float(len(assay_frame) / (19 * len(wildtype_sequence))),
    }


def sequence_array_from_strings(sequences: pd.Series) -> np.ndarray:
    return np.asarray([[symbol for symbol in sequence] for sequence in sequences.astype(str)], dtype=object)


def sequences_to_indices(sequences: pd.Series) -> np.ndarray:
    arrays = []
    for sequence in sequences.astype(str):
        arrays.append([AA_TO_INT[symbol] for symbol in sequence])
    return np.asarray(arrays, dtype=np.int64)


AA_TO_INT = {aa: idx for idx, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")}


def compute_synthetic_metrics(
    assay_frame: pd.DataFrame,
    *,
    wildtype_sequence: str,
    alignment_profile: AlignmentProfile,
    base_config: dict[str, Any],
    fitted_config: dict[str, Any] | None,
) -> dict[str, Any]:
    config_kwargs = dict(base_config)
    if fitted_config is not None:
        config_kwargs.update(fitted_config)
    config_kwargs["L"] = len(wildtype_sequence)
    landscape = FitnessLandscapeEnv(
        LandscapeConfig(**config_kwargs),
        alignment_profile=alignment_profile,
        reference_sequence=wildtype_sequence,
    )
    seqs = sequences_to_indices(assay_frame["mutated_sequence"])
    fitness = landscape.evaluate_batch(seqs)
    scores = zscore(pd.Series(fitness, dtype=np.float64)).to_numpy(dtype=np.float64)
    assay_with_scores = assay_frame.copy()
    assay_with_scores["synthetic_score"] = scores
    position_sensitivity = empirical_position_sensitivity(
        assay_with_scores,
        sequence_length=len(wildtype_sequence),
        score_column="synthetic_score",
    )
    conservation_corr = correlation_or_none(alignment_profile.conservation, position_sensitivity)
    return {
        "score_summary": assay_score_summary(scores),
        "conservation_sensitivity_correlation": conservation_corr,
        "reference_distance_to_peak": int(
            landscape.hamming_distance(landscape.reference, landscape.peak_sequence)
        ),
        "reference_fraction_of_peak": float(landscape.evaluate(landscape.reference) / landscape.peak_fitness)
            if landscape.peak_fitness != 0
            else None,
    }


def compare_metrics(empirical: dict[str, Any], synthetic: dict[str, Any]) -> dict[str, Any]:
    emp = empirical["score_summary"]
    syn = synthetic["score_summary"]
    return {
        "score_skewness_abs_diff": abs(float(emp["skewness"]) - float(syn["skewness"])),
        "q05_abs_diff": abs(float(emp["q05"]) - float(syn["q05"])),
        "q95_abs_diff": abs(float(emp["q95"]) - float(syn["q95"])),
        "fraction_gt_1sd_abs_diff": abs(float(emp["fraction_gt_1sd"]) - float(syn["fraction_gt_1sd"])),
        "fraction_lt_minus_1sd_abs_diff": abs(
            float(emp["fraction_lt_minus_1sd"]) - float(syn["fraction_lt_minus_1sd"])
        ),
        "conservation_sensitivity_abs_diff": None
        if empirical["conservation_sensitivity_correlation"] is None
        or synthetic["conservation_sensitivity_correlation"] is None
        else abs(
            float(empirical["conservation_sensitivity_correlation"])
            - float(synthetic["conservation_sensitivity_correlation"])
        ),
    }


def build_empirical_landscape(
    *,
    dms_id: str,
    assay_frame: pd.DataFrame,
    wildtype_sequence: str,
) -> ProteinGymEmpiricalLandscape:
    return ProteinGymEmpiricalLandscape(
        dms_id=dms_id,
        wildtype_sequence=wildtype_sequence,
        assay_frame=assay_frame,
        score_column="DMS_score_z",
    )


def calibration_options_from_config(config: dict[str, Any]) -> CalibrationOptions:
    return CalibrationOptions(**config)


def to_builtin(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, dict):
        return {str(key): to_builtin(item) for key, item in value.items()}
    if isinstance(value, list):
        return [to_builtin(item) for item in value]
    return value


def installed_proteingym_version() -> str | None:
    try:
        return importlib.metadata.version("proteingym")
    except importlib.metadata.PackageNotFoundError:
        return None


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config))

    proteingym_cfg = config["proteingym"]
    panel_cfg = config["panel"]
    mmseqs_cfg = dict(config["mmseqs"])
    base_config = dict(config["adapt_env"]["base_config"])

    reference_csv_path = project_root / proteingym_cfg["reference_csv_path"]
    substitutions_parquet_path = project_root / proteingym_cfg["substitutions_parquet_path"]
    alignment_dir = project_root / mmseqs_cfg["alignment_dir"]
    mmseqs_cfg["cache_dir"] = str(project_root / mmseqs_cfg["cache_dir"])

    download_file(proteingym_cfg["reference_csv_url"], reference_csv_path)
    download_file(proteingym_cfg["substitutions_parquet_url"], substitutions_parquet_path)

    reference_df = pd.read_csv(reference_csv_path)
    panel_df = select_assay_panel(reference_df, panel_cfg)
    panel_ids = panel_df["DMS_id"].tolist()
    assay_df = pd.read_parquet(
        substitutions_parquet_path,
        columns=["DMS_id", "mutant", "mutated_sequence", "DMS_score", "DMS_score_bin"],
        filters=[("DMS_id", "in", panel_ids)],
    )
    assay_df["mutated_sequence"] = assay_df["mutated_sequence"].map(canonical_sequence)

    raw_assay_dir = project_root / proteingym_cfg["assay_output_dir"]
    raw_assay_dir.mkdir(parents=True, exist_ok=True)

    empirical_landscapes = []
    alignment_paths = []
    wildtypes = []
    per_assay_payload = []

    for row in panel_df.itertuples(index=False):
        dms_id = str(row.DMS_id)
        wildtype_sequence = canonical_sequence(str(row.target_seq))
        subset = assay_df[assay_df["DMS_id"] == dms_id].copy()
        if subset.empty:
            raise ValueError(f"No ProteinGym rows found for selected assay {dms_id}")
        subset["mutation_position"] = mutation_positions(wildtype_sequence, subset["mutated_sequence"])
        subset["DMS_score_z"] = zscore(subset["DMS_score"].astype(np.float64))

        assay_output_path = raw_assay_dir / f"{dms_id}.csv"
        subset.to_csv(assay_output_path, index=False)

        alignment_path = alignment_dir / f"{dms_id}.fasta"
        build_alignment_fasta(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            fasta_path=alignment_path,
            mmseqs_cfg=mmseqs_cfg,
        )
        alignment_profile = AlignmentProfile.from_alignment(
            alignment_path,
            max_sequences=int(mmseqs_cfg["max_alignment_sequences"]),
        )

        empirical_metrics = compute_empirical_metrics(
            subset,
            wildtype_sequence=wildtype_sequence,
            alignment_profile=alignment_profile,
            score_column="DMS_score_z",
        )
        empirical_landscape = build_empirical_landscape(
            dms_id=dms_id,
            assay_frame=subset,
            wildtype_sequence=wildtype_sequence,
        )
        empirical_summary = summarize_empirical_landscape(
            empirical_landscape,
            kind="functional",
            alignment_profile=alignment_path,
            wildtype=wildtype_sequence,
            options=calibration_options_from_config(config["calibration"]["options"]),
            index=len(empirical_landscapes),
        )

        empirical_landscapes.append(empirical_landscape)
        alignment_paths.append(alignment_path)
        wildtypes.append(wildtype_sequence)
        per_assay_payload.append(
            {
                "dms_id": dms_id,
                "uniprot_id": str(row.UniProt_ID),
                "taxon": str(row.taxon),
                "selection_type": str(row.selection_type),
                "sequence_length": int(row.seq_len),
                "n_single_mutants": int(row.DMS_number_single_mutants),
                "alignment_path": str(alignment_path),
                "assay_csv_path": str(assay_output_path),
                "empirical_metrics": to_builtin(empirical_metrics),
                "empirical_effect_summary": to_builtin(empirical_summary.effect_summary),
                "empirical_conservation_correlation": empirical_summary.conservation_sensitivity_correlation,
            }
        )

    calibration_result = calibrate_synthetic_landscape(
        functional_landscapes=empirical_landscapes,
        functional_alignment_profiles=alignment_paths,
        functional_wildtypes=wildtypes,
        base_config=LandscapeConfig(**base_config),
        options=calibration_options_from_config(config["calibration"]["options"]),
    )

    fitted_config = to_builtin(calibration_result.fitted_parameters)
    for assay_payload, row in zip(per_assay_payload, panel_df.itertuples(index=False)):
        dms_id = str(row.DMS_id)
        subset = pd.read_csv(project_root / proteingym_cfg["assay_output_dir"] / f"{dms_id}.csv")
        alignment_profile = AlignmentProfile.from_alignment(assay_payload["alignment_path"])
        synthetic_metrics = compute_synthetic_metrics(
            subset,
            wildtype_sequence=canonical_sequence(str(row.target_seq)),
            alignment_profile=alignment_profile,
            base_config=base_config,
            fitted_config=fitted_config,
        )
        assay_payload["synthetic_metrics"] = to_builtin(synthetic_metrics)
        assay_payload["comparison"] = to_builtin(
            compare_metrics(assay_payload["empirical_metrics"], synthetic_metrics)
        )

    panel_csv_path = output_dir / "selected_panel.csv"
    panel_df.to_csv(panel_csv_path, index=False)

    summary_payload = {
        "experiment_id": config["experiment"]["id"],
        "proteingym_version": proteingym_cfg["version"],
        "proteingym_package_version": installed_proteingym_version(),
        "reference_csv_path": str(reference_csv_path),
        "substitutions_parquet_path": str(substitutions_parquet_path),
        "alignment_dir": str(alignment_dir),
        "panel_csv_path": str(panel_csv_path),
        "fitted_parameters": fitted_config,
        "calibration_validation": to_builtin(calibration_result.validation),
        "calibration_objective_terms": to_builtin(calibration_result.objective_terms),
        "per_assay": per_assay_payload,
    }

    (output_dir / "summary.json").write_text(
        json.dumps(to_builtin(summary_payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output_dir / "summary.json")


if __name__ == "__main__":
    main()
