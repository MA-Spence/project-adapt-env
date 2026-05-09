"""Reusable ProteinGym panel preparation helpers."""

from __future__ import annotations

import importlib.metadata
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib import request

import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr

from ._compat import ensure_external_paths
from .utils import hash_mod, to_builtin

ensure_external_paths()

from adaptenv.calibration import CalibrationOptions  # noqa: E402
from adaptenv.mmseqs import MMseqsServerClient  # noqa: E402


AA_SET = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class EmpiricalSequence:
    """Minimal sequence object compatible with Adapt-Env calibration helpers."""

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
    """Expose a ProteinGym assay on an arbitrary coordinate system to Adapt-Env."""

    def __init__(
        self,
        *,
        dms_id: str,
        wildtype_sequence: str,
        wildtype_score: float,
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
        self._fitness_by_sequence: dict[str, float] = {wildtype_sequence: float(wildtype_score)}

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


@dataclass
class PreparedAssay:
    """One prepared ProteinGym assay plus derived calibration objects."""

    dms_id: str
    uniprot_id: str
    taxon: str
    selection_type: str
    sequence_length: int
    n_single_mutants_reference: int
    n_multiple_mutants_reference: int
    n_variants_loaded: int
    n_variants_for_calibration: int
    wildtype_sequence: str
    alignment_path: Path
    assay_csv_path: Path
    mavenn_metrics: dict[str, Any]
    raw_score_summary: dict[str, float]
    latent_phi_summary: dict[str, float]
    raw_landscape: ProteinGymEmpiricalLandscape
    latent_landscape: ProteinGymEmpiricalLandscape

    def summary_payload(self) -> dict[str, Any]:
        return {
            "dms_id": self.dms_id,
            "uniprot_id": self.uniprot_id,
            "taxon": self.taxon,
            "selection_type": self.selection_type,
            "sequence_length": self.sequence_length,
            "n_single_mutants_reference": self.n_single_mutants_reference,
            "n_multiple_mutants_reference": self.n_multiple_mutants_reference,
            "n_variants_loaded": self.n_variants_loaded,
            "n_variants_for_calibration": self.n_variants_for_calibration,
            "alignment_path": str(self.alignment_path),
            "assay_csv_path": str(self.assay_csv_path),
            "mavenn": to_builtin(self.mavenn_metrics),
            "raw_score_summary": self.raw_score_summary,
            "latent_phi_summary": self.latent_phi_summary,
        }


@dataclass
class PreparedPanel:
    """Prepared ProteinGym panel and derived empirical landscapes."""

    panel_df: pd.DataFrame
    assays: list[PreparedAssay]

    @property
    def raw_landscapes(self) -> list[ProteinGymEmpiricalLandscape]:
        return [assay.raw_landscape for assay in self.assays]

    @property
    def latent_landscapes(self) -> list[ProteinGymEmpiricalLandscape]:
        return [assay.latent_landscape for assay in self.assays]

    @property
    def alignment_paths(self) -> list[Path]:
        return [assay.alignment_path for assay in self.assays]

    @property
    def wildtypes(self) -> list[str]:
        return [assay.wildtype_sequence for assay in self.assays]

    def summary_payload(self) -> list[dict[str, Any]]:
        return [assay.summary_payload() for assay in self.assays]

    def mavenn_metrics_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "dms_id": assay.dms_id,
                    "taxon": assay.taxon,
                    "sequence_length": assay.sequence_length,
                    "fit_variant_count": assay.mavenn_metrics["fit_variant_count"],
                    "train_variant_count": assay.mavenn_metrics["train_variant_count"],
                    "test_variant_count": assay.mavenn_metrics["test_variant_count"],
                    "wildtype_phi": assay.mavenn_metrics["wildtype_phi"],
                    "wildtype_yhat": assay.mavenn_metrics["wildtype_yhat"],
                    "test_nrmse": assay.mavenn_metrics["test_metrics"]["nrmse"],
                    "test_spearman": assay.mavenn_metrics["test_metrics"]["spearman"],
                    "I_var_test_bits": assay.mavenn_metrics["I_var_test_bits"],
                    "I_pred_test_bits": assay.mavenn_metrics["I_pred_test_bits"],
                }
                for assay in self.assays
            ]
        )


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_dict(dict(merged[key]), value)
        else:
            merged[key] = value
    return merged


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

    allowed_taxa = [str(item) for item in panel_cfg.get("allowed_taxa", [])]
    assays_per_taxon = int(panel_cfg.get("assays_per_taxon", 0))
    panel_size = int(panel_cfg["panel_size"])
    min_single_mutants = int(panel_cfg["min_single_mutants"])
    min_multiple_mutants = int(panel_cfg.get("min_multiple_mutants", 0))
    max_sequence_length = int(panel_cfg["max_sequence_length"])
    coarse_selection_type = panel_cfg.get("coarse_selection_type")
    selection_type = panel_cfg.get("selection_type")
    require_multiple_mutants = bool(panel_cfg.get("require_multiple_mutants", False))

    filtered = reference_df.copy()
    filtered = filtered[filtered["target_seq"].notna()]
    filtered = filtered[filtered["target_seq"].map(lambda seq: set(str(seq).upper()) <= AA_SET)]
    if allowed_taxa:
        filtered = filtered[filtered["taxon"].isin(allowed_taxa)]
    if coarse_selection_type:
        filtered = filtered[filtered["coarse_selection_type"] == coarse_selection_type]
    if selection_type:
        filtered = filtered[filtered["selection_type"] == selection_type]
    filtered = filtered[filtered["DMS_number_single_mutants"].astype(int) >= min_single_mutants]
    filtered = filtered[
        filtered["DMS_number_multiple_mutants"].fillna(0).astype(int) >= min_multiple_mutants
    ]
    filtered = filtered[filtered["seq_len"].astype(int) <= max_sequence_length]
    if require_multiple_mutants:
        filtered = filtered[filtered["includes_multiple_mutants"].map(normalize_includes_multiple)]
    filtered = filtered.sort_values(
        by=["taxon", "DMS_number_multiple_mutants", "DMS_number_single_mutants", "seq_len"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )

    chosen_frames: list[pd.DataFrame] = []
    chosen_ids: set[str] = set()
    if allowed_taxa and assays_per_taxon > 0:
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


def mutation_counts(wildtype: str, mutated_sequences: pd.Series) -> np.ndarray:
    wt_chars = np.asarray(list(wildtype), dtype=object)
    counts: list[int] = []
    for sequence in mutated_sequences.astype(str):
        seq_chars = np.asarray(list(sequence), dtype=object)
        if seq_chars.shape[0] != wt_chars.shape[0]:
            counts.append(-1)
            continue
        counts.append(int(np.count_nonzero(seq_chars != wt_chars)))
    return np.asarray(counts, dtype=np.int64)


def mutation_positions(wildtype: str, mutated_sequences: pd.Series) -> list[int]:
    positions: list[int] = []
    wt_chars = np.asarray(list(wildtype), dtype=object)
    for sequence in mutated_sequences.astype(str):
        seq_chars = np.asarray(list(sequence), dtype=object)
        diff = np.flatnonzero(seq_chars != wt_chars)
        if diff.size != 1:
            raise ValueError(
                f"Expected a single-substitution sequence, but observed {diff.size} differences"
            )
        positions.append(int(diff[0]))
    return positions


def prediction_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float | None]:
    obs = np.asarray(observed, dtype=np.float64)
    pred = np.asarray(predicted, dtype=np.float64)
    mask = np.isfinite(obs) & np.isfinite(pred)
    obs = obs[mask]
    pred = pred[mask]
    if obs.size == 0:
        return {
            "count": 0.0,
            "rmse": None,
            "nrmse": None,
            "pearson": None,
            "spearman": None,
        }
    rmse = float(np.sqrt(np.mean((pred - obs) ** 2)))
    scale = max(float(np.std(obs)), 1e-8)
    pearson = float(np.corrcoef(obs, pred)[0, 1]) if obs.size >= 2 and np.std(pred) > 0 else None
    spearman = float(spearmanr(obs, pred).statistic) if obs.size >= 2 else None
    return {
        "count": float(obs.size),
        "rmse": rmse,
        "nrmse": float(rmse / scale),
        "pearson": pearson,
        "spearman": spearman,
    }


def _coerce_info_metric(value: Any) -> float | None:
    if isinstance(value, tuple):
        return _coerce_info_metric(value[0])
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def assay_score_summary(scores: np.ndarray) -> dict[str, float]:
    series = pd.Series(scores, dtype=np.float64)
    return {
        "count": float(series.size),
        "mean": float(series.mean()),
        "std": float(series.std(ddof=0)),
        "skewness": float(series.skew()),
        "q05": float(np.quantile(series, 0.05)),
        "q25": float(np.quantile(series, 0.25)),
        "q50": float(np.quantile(series, 0.50)),
        "q75": float(np.quantile(series, 0.75)),
        "q95": float(np.quantile(series, 0.95)),
    }


def fit_mavenn_model(
    *,
    assay_frame: pd.DataFrame,
    wildtype_sequence: str,
    cfg: dict[str, Any],
) -> tuple[Any, dict[str, Any]]:
    try:
        import mavenn
    except ModuleNotFoundError as exc:  # pragma: no cover - runtime guard
        raise RuntimeError(
            "mavenn is required for ProteinGym calibration experiments. "
            "Ensure the experiment overlay is applied."
        ) from exc

    filtered = assay_frame[assay_frame["mutation_count"] <= int(cfg["max_mutation_count"])].copy()
    filtered = filtered[filtered["mutation_count"] >= 1].copy()
    max_variants = int(cfg.get("max_variants_per_assay", 0))
    if max_variants > 0 and len(filtered) > max_variants:
        filtered = (
            filtered.sample(
                n=max_variants,
                random_state=int(cfg.get("sample_seed", 0)),
                replace=False,
            )
            .sort_index()
            .reset_index(drop=True)
        )

    if filtered.empty:
        raise ValueError("MAVE-NN fit subset is empty after mutation-count filtering")

    x_all = filtered["mutated_sequence"].to_numpy(dtype=object)
    y_all = filtered["DMS_score"].to_numpy(dtype=np.float64)
    holdout_modulo = max(int(cfg["test_holdout_modulo"]), 2)
    holdout_remainder = int(cfg["test_holdout_remainder"]) % holdout_modulo
    test_mask = np.asarray(
        [hash_mod(str(sequence), holdout_modulo) == holdout_remainder for sequence in x_all],
        dtype=bool,
    )
    if int(np.sum(~test_mask)) < 50 or int(np.sum(test_mask)) < 25:
        test_mask = np.zeros_like(test_mask, dtype=bool)
        step = max(len(test_mask) // 10, 1)
        test_mask[::step] = True
        if int(np.sum(~test_mask)) == 0:
            test_mask[0] = False

    train_x = x_all[~test_mask]
    train_y = y_all[~test_mask]
    test_x = x_all[test_mask]
    test_y = y_all[test_mask]

    model = mavenn.Model(
        L=len(wildtype_sequence),
        alphabet=str(cfg.get("alphabet", "protein")),
        regression_type=str(cfg.get("regression_type", "GE")),
        gpmap_type=str(cfg.get("gpmap_type", "additive")),
        ge_nonlinearity_type=str(cfg.get("ge_nonlinearity_type", "nonlinear")),
        ge_nonlinearity_monotonic=bool(cfg.get("ge_nonlinearity_monotonic", True)),
        ge_noise_model_type=str(cfg.get("ge_noise_model_type", "SkewedT")),
        ge_heteroskedasticity_order=int(cfg.get("ge_heteroskedasticity_order", 0)),
        theta_regularization=float(cfg.get("theta_regularization", 0.01)),
        eta_regularization=float(cfg.get("eta_regularization", 0.001)),
    )
    model.set_data(
        x=train_x,
        y=train_y,
        validation_frac=float(cfg.get("validation_fraction", 0.1)),
        shuffle=True,
        verbose=False,
    )
    model.fit(
        epochs=int(cfg.get("epochs", 50)),
        learning_rate=float(cfg.get("learning_rate", 0.005)),
        early_stopping=True,
        early_stopping_patience=int(cfg.get("early_stopping_patience", 10)),
        batch_size=int(cfg.get("batch_size", 128)),
        verbose=False,
        try_tqdm=False,
    )

    phi_train = np.asarray(model.x_to_phi(train_x), dtype=np.float64).reshape(-1)
    yhat_train = np.asarray(model.phi_to_yhat(phi_train), dtype=np.float64).reshape(-1)
    train_metrics = prediction_metrics(train_y, yhat_train)

    report_information_metrics = bool(cfg.get("report_information_metrics", False))
    if test_x.size:
        phi_test = np.asarray(model.x_to_phi(test_x), dtype=np.float64).reshape(-1)
        yhat_test = np.asarray(model.phi_to_yhat(phi_test), dtype=np.float64).reshape(-1)
        test_metrics = prediction_metrics(test_y, yhat_test)
        if report_information_metrics:
            i_var = _coerce_info_metric(model.I_variational(test_x, test_y, uncertainty=False))
            i_pred = _coerce_info_metric(
                model.I_predictive(test_x, test_y, uncertainty=False, verbose=False)
            )
        else:
            i_var = None
            i_pred = None
    else:
        test_metrics = prediction_metrics(np.empty(0), np.empty(0))
        i_var = None
        i_pred = None

    wt_phi = float(np.asarray(model.x_to_phi(np.asarray([wildtype_sequence], dtype=object))).reshape(-1)[0])
    wt_yhat = float(np.asarray(model.phi_to_yhat(np.asarray([wt_phi], dtype=np.float64))).reshape(-1)[0])

    return model, {
        "fit_variant_count": int(len(filtered)),
        "train_variant_count": int(train_x.size),
        "test_variant_count": int(test_x.size),
        "wildtype_phi": wt_phi,
        "wildtype_yhat": wt_yhat,
        "train_metrics": train_metrics,
        "test_metrics": test_metrics,
        "I_var_test_bits": i_var,
        "I_pred_test_bits": i_pred,
    }


def build_empirical_landscape(
    *,
    dms_id: str,
    wildtype_sequence: str,
    wildtype_score: float,
    assay_frame: pd.DataFrame,
    score_column: str,
) -> ProteinGymEmpiricalLandscape:
    return ProteinGymEmpiricalLandscape(
        dms_id=dms_id,
        wildtype_sequence=wildtype_sequence,
        wildtype_score=wildtype_score,
        assay_frame=assay_frame,
        score_column=score_column,
    )


def calibration_options_from_config(config: dict[str, Any]) -> CalibrationOptions:
    return CalibrationOptions(**config)


def installed_proteingym_version() -> str | None:
    try:
        return importlib.metadata.version("proteingym")
    except importlib.metadata.PackageNotFoundError:
        return None


def installed_mavenn_version() -> str | None:
    try:
        return importlib.metadata.version("mavenn")
    except importlib.metadata.PackageNotFoundError:
        return None


def prepare_proteingym_panel(
    *,
    project_root: Path,
    proteingym_cfg: dict[str, Any],
    panel_cfg: dict[str, Any],
    mmseqs_cfg: dict[str, Any],
    mavenn_cfg: dict[str, Any],
    calibration_max_mutation_count: int,
) -> PreparedPanel:
    reference_csv_path = project_root / proteingym_cfg["reference_csv_path"]
    substitutions_parquet_path = project_root / proteingym_cfg["substitutions_parquet_path"]
    alignment_dir = project_root / mmseqs_cfg["alignment_dir"]

    download_file(proteingym_cfg["reference_csv_url"], reference_csv_path)
    download_file(proteingym_cfg["substitutions_parquet_url"], substitutions_parquet_path)

    reference_df = pd.read_csv(reference_csv_path)
    panel_df = select_assay_panel(reference_df, panel_cfg)
    panel_ids = panel_df["DMS_id"].tolist()
    assay_df = pd.read_parquet(
        substitutions_parquet_path,
        columns=["DMS_id", "mutant", "mutated_sequence", "DMS_score"],
        filters=[("DMS_id", "in", panel_ids)],
    )
    assay_df["mutated_sequence"] = assay_df["mutated_sequence"].map(canonical_sequence)

    raw_assay_dir = project_root / proteingym_cfg["assay_output_dir"]
    raw_assay_dir.mkdir(parents=True, exist_ok=True)

    assays: list[PreparedAssay] = []
    for row in panel_df.itertuples(index=False):
        dms_id = str(row.DMS_id)
        wildtype_sequence = canonical_sequence(str(row.target_seq))
        subset = assay_df[assay_df["DMS_id"] == dms_id].copy()
        if subset.empty:
            raise ValueError(f"No ProteinGym rows found for selected assay {dms_id}")
        subset = subset[subset["mutated_sequence"].map(len) == len(wildtype_sequence)].copy()
        subset["mutation_count"] = mutation_counts(wildtype_sequence, subset["mutated_sequence"])
        subset = subset[subset["mutation_count"] >= 1].copy()
        if subset.empty:
            raise ValueError(f"No valid ProteinGym rows remained for selected assay {dms_id}")
        singles = subset["mutation_count"] == 1
        subset.loc[singles, "mutation_position"] = mutation_positions(
            wildtype_sequence,
            subset.loc[singles, "mutated_sequence"],
        )

        assay_output_path = raw_assay_dir / f"{dms_id}.csv"
        subset.to_csv(assay_output_path, index=False)

        alignment_path = alignment_dir / f"{dms_id}.fasta"
        build_alignment_fasta(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            fasta_path=alignment_path,
            mmseqs_cfg=mmseqs_cfg,
        )

        mavenn_model, mavenn_metrics = fit_mavenn_model(
            assay_frame=subset,
            wildtype_sequence=wildtype_sequence,
            cfg=mavenn_cfg,
        )
        subset["latent_phi"] = np.asarray(
            mavenn_model.x_to_phi(subset["mutated_sequence"].to_numpy(dtype=object)),
            dtype=np.float64,
        ).reshape(-1)
        subset["mavenn_yhat"] = np.asarray(
            mavenn_model.phi_to_yhat(subset["latent_phi"].to_numpy(dtype=np.float64)),
            dtype=np.float64,
        ).reshape(-1)

        calibration_subset = subset[subset["mutation_count"] <= int(calibration_max_mutation_count)].copy()
        if calibration_subset.empty:
            raise ValueError(f"No calibration subset remained for assay {dms_id}")

        raw_landscape = build_empirical_landscape(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            wildtype_score=float(mavenn_metrics["wildtype_yhat"]),
            assay_frame=calibration_subset,
            score_column="DMS_score",
        )
        latent_landscape = build_empirical_landscape(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            wildtype_score=float(mavenn_metrics["wildtype_phi"]),
            assay_frame=calibration_subset,
            score_column="latent_phi",
        )
        assays.append(
            PreparedAssay(
                dms_id=dms_id,
                uniprot_id=str(row.UniProt_ID),
                taxon=str(row.taxon),
                selection_type=str(row.selection_type),
                sequence_length=int(row.seq_len),
                n_single_mutants_reference=int(row.DMS_number_single_mutants),
                n_multiple_mutants_reference=int(row.DMS_number_multiple_mutants),
                n_variants_loaded=int(len(subset)),
                n_variants_for_calibration=int(len(calibration_subset)),
                wildtype_sequence=wildtype_sequence,
                alignment_path=alignment_path,
                assay_csv_path=assay_output_path,
                mavenn_metrics=mavenn_metrics,
                raw_score_summary=assay_score_summary(subset["DMS_score"].to_numpy(dtype=np.float64)),
                latent_phi_summary=assay_score_summary(subset["latent_phi"].to_numpy(dtype=np.float64)),
                raw_landscape=raw_landscape,
                latent_landscape=latent_landscape,
            )
        )
    return PreparedPanel(panel_df=panel_df, assays=assays)
