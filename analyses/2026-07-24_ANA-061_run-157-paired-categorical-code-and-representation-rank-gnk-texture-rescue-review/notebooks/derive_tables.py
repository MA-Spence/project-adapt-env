from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


NOTEBOOK_DIR = Path(__file__).resolve().parent
ANALYSIS_DIR = NOTEBOOK_DIR.parent
PROJECT_ROOT = ANALYSIS_DIR.parents[1]
TABLE_DIR = ANALYSIS_DIR / "tables"
SOURCE_DIR = ANALYSIS_DIR / "source_tables"
RUN_DIR = (
    PROJECT_ROOT
    / "data/processed/paired-categorical-code-and-representation-rank-gnk-texture-rescue/RUN-157"
)
DVC_POINTER = (
    PROJECT_ROOT
    / "data/processed/paired-categorical-code-and-representation-rank-gnk-texture-rescue/RUN-157.dvc"
)
CONFIG_PATH = (
    PROJECT_ROOT
    / "experiments/2026-07-23_EXP-063_paired-categorical-code-and-representation-rank-gnk-texture-rescue/config.yaml"
)
RUN_RECORD = (
    PROJECT_ROOT
    / "experiments/2026-07-23_EXP-063_paired-categorical-code-and-representation-rank-gnk-texture-rescue/runs/RUN-157.yaml"
)
REMOTE_AUDIT = SOURCE_DIR / "run-157_remote_audit.json"

ARMS = ("binary_r1", "dense_q20_r1", "binary_r8", "dense_q20_r8")
RESCUE_ARMS = ARMS[1:]
CONTEXT_KEYS = ["reference_id", "length_bin", "length", "core_seed", "gnk_seed"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if pd.isna(value):
        return None
    return value


def bool_value(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def rank_correlation(values: list[float] | np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    if array.size < 2 or not np.isfinite(array).all():
        return float("nan")
    ranks = pd.Series(array).rank(method="average").to_numpy(dtype=float)
    if float(np.std(ranks)) == 0.0:
        return 0.0
    return float(np.corrcoef(np.arange(array.size, dtype=float), ranks)[0, 1])


def simple_bootstrap_interval(
    values: np.ndarray, *, replicates: int, seed: int
) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    rng = np.random.default_rng(seed)
    boot = np.asarray(
        [
            np.median(rng.choice(array, size=array.size, replace=True))
            for _ in range(replicates)
        ],
        dtype=float,
    )
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def crossed_bootstrap_lower(
    frame: pd.DataFrame, *, metric: str, replicates: int, seed: int
) -> float:
    selected = frame[np.isfinite(frame[metric])].copy().sort_values(
        ["reference_id", "core_seed", "gnk_seed"], kind="mergesort"
    )
    references = sorted(str(value) for value in selected["reference_id"].unique())
    matrices: dict[str, np.ndarray] = {}
    for reference_id in references:
        block = selected[selected["reference_id"] == reference_id]
        cores = sorted(int(value) for value in block["core_seed"].unique())
        gnks = sorted(int(value) for value in block["gnk_seed"].unique())
        matrix = (
            block.pivot(index="core_seed", columns="gnk_seed", values=metric)
            .reindex(index=cores, columns=gnks)
            .to_numpy(dtype=float)
        )
        if matrix.shape != (len(cores), len(gnks)) or not np.isfinite(matrix).all():
            raise ValueError(f"incomplete crossed design for {reference_id}")
        matrices[reference_id] = matrix
    rng = np.random.default_rng(seed)
    boot = np.empty(replicates, dtype=float)
    for replicate in range(replicates):
        sampled_references = rng.choice(
            references, size=len(references), replace=True
        )
        reference_values: list[float] = []
        for reference_id in sampled_references:
            matrix = matrices[str(reference_id)]
            sampled_cores = rng.integers(0, matrix.shape[0], size=matrix.shape[0])
            sampled_gnks = rng.integers(0, matrix.shape[1], size=matrix.shape[1])
            reference_values.append(
                float(np.median(matrix[np.ix_(sampled_cores, sampled_gnks)]))
            )
        boot[replicate] = float(np.median(reference_values))
    return float(np.quantile(boot, 0.025))


def minimum_bin_fraction(
    frame: pd.DataFrame, pass_column: str
) -> tuple[float, dict[str, float]]:
    fractions = frame.groupby("length_bin")[pass_column].mean().sort_index()
    if fractions.empty:
        return 0.0, {}
    return float(fractions.min()), {
        str(key): float(value) for key, value in fractions.items()
    }


def recompute_pairing(
    representation: pd.DataFrame, config: dict[str, Any]
) -> tuple[bool, int]:
    technical = representation[
        representation["row_type"] == "technical_context_audit"
    ].copy()
    tolerance = config["acceptance"]["technical"]
    manipulation = config["acceptance"]["manipulation"]
    passes: list[bool] = []
    for _, group in technical.groupby(CONTEXT_KEYS, sort=True):
        arms = sorted(str(value) for value in group["representation_arm"])
        dense = group[group["representation_arm"].str.startswith("dense_q20")]
        binary = group[group["representation_arm"].str.startswith("binary")]
        rank_pass = all(
            int(row["minimum_sampled_operator_gram_rank"])
            >= (
                1
                if str(row["representation_arm"]).endswith("_r1")
                else int(manipulation["rank8_effective_rank_min"])
            )
            for _, row in group.iterrows()
        )
        passes.append(
            bool(
                arms == sorted(ARMS)
                and group["support_fingerprint"].nunique() == 1
                and group["pair_pool_hash"].nunique() == 1
                and group["component_zero_base_coefficient_hash"].nunique() == 1
                and group["pair_design_hash"].nunique() == 1
                and group["square_design_hash"].nunique() == 1
                and group["all_outputs_finite"].map(bool_value).all()
                and group["maximum_operator_energy_relative_error"].max()
                <= float(tolerance["operator_energy_relative_error_max"])
                and group["maximum_formula_parity_abs_difference"].max()
                <= float(tolerance["formula_parity_atol"])
                and group["candidate_reference_center_abs"].max()
                <= float(tolerance["reference_center_atol"])
                and group["maximum_sampled_contrast_mean_abs"].max()
                <= float(tolerance["contrast_mean_atol"])
                and group["maximum_sampled_contrast_second_moment_error"].max()
                <= float(tolerance["contrast_second_moment_atol"])
                and group["surrogate_latent_epsilon_max_abs"].max()
                <= float(tolerance["additive_surrogate_epistasis_atol"])
                and group["maximum_zero_amplitude_residual_rms"].max()
                <= float(tolerance["zero_amplitude_atol"])
                and group.loc[
                    group["representation_arm"] == "binary_r1",
                    "binary_rank1_legacy_max_abs_difference",
                ].max()
                <= float(tolerance["binary_rank1_score_parity_atol"])
                and dense["sampled_nonself_tie_fraction"].max()
                <= float(manipulation["dense_q20_nonself_tie_fraction_max"])
                and rank_pass
                and binary["sampled_binary_10_10_balance_pass"]
                .map(bool_value)
                .all()
            )
        )
    return bool(len(passes) == 360 and all(passes)), len(passes)


def recompute_directionality(
    texture: pd.DataFrame, arm: str, config: dict[str, Any]
) -> tuple[bool, dict[str, Any], pd.DataFrame]:
    gates = config["acceptance"]["dose_response"]
    doses = [float(value) for value in config["intervention"]["directionality_amplitudes"]]
    selected = texture[
        (texture["representation_arm"] == arm)
        & (texture["additive_core_mode"] == "actual_additive_core")
        & texture["pooled_supported_reference"].map(bool_value)
        & texture["amplitude"].isin(doses)
    ].copy()

    def curve(group: pd.DataFrame) -> dict[str, Any]:
        medians = (
            group.groupby("amplitude")[
                ["direct_tau_y_rms", "direct_tau_practical_resolution_fraction"]
            ]
            .median()
            .reindex(doses)
        )
        rms = medians["direct_tau_y_rms"].to_numpy(dtype=float)
        practical = medians[
            "direct_tau_practical_resolution_fraction"
        ].to_numpy(dtype=float)
        rms_rho = rank_correlation(rms)
        practical_rho = rank_correlation(practical)
        rms_delta = float(rms[-1] - rms[0])
        practical_delta = float(practical[-1] - practical[0])
        passed = bool(
            rms_rho + float(gates["rho_comparison_atol"])
            >= float(gates["spearman_min"])
            and practical_rho + float(gates["rho_comparison_atol"])
            >= float(gates["spearman_min"])
            and rms_delta >= float(gates["rms_growth_floor"])
            and practical_delta >= float(gates["practical_fraction_growth_floor"])
        )
        return {
            "tau_rms_medians": json.dumps([float(value) for value in rms]),
            "practical_fraction_medians": json.dumps(
                [float(value) for value in practical]
            ),
            "tau_rms_spearman": rms_rho,
            "practical_fraction_spearman": practical_rho,
            "tau_rms_a2_minus_a0.25": rms_delta,
            "practical_fraction_a2_minus_a0.25": practical_delta,
            "directionality_pass": passed,
        }

    rows: list[dict[str, Any]] = []
    for keys, group in selected.groupby(CONTEXT_KEYS, sort=True):
        rows.append(
            {
                "aggregation_unit": "protein_core_gnk_cell",
                **dict(zip(CONTEXT_KEYS, keys)),
                **curve(group),
            }
        )
    for keys, group in selected.groupby(
        ["reference_id", "length_bin", "length"], sort=True
    ):
        rows.append(
            {
                "aggregation_unit": "protein_median_over_cells",
                **dict(zip(["reference_id", "length_bin", "length"], keys)),
                "core_seed": None,
                "gnk_seed": None,
                **curve(group),
            }
        )
    frame = pd.DataFrame(rows)
    cells = frame[frame["aggregation_unit"] == "protein_core_gnk_cell"].copy()
    proteins = frame[frame["aggregation_unit"] == "protein_median_over_cells"].copy()
    cell_min, cell_bins = minimum_bin_fraction(cells, "directionality_pass")
    protein_min, protein_bins = minimum_bin_fraction(
        proteins, "directionality_pass"
    )
    rms_lower, _ = simple_bootstrap_interval(
        proteins["tau_rms_a2_minus_a0.25"].to_numpy(dtype=float),
        replicates=int(config["acceptance"]["bootstrap_replicates"]),
        seed=int(config["acceptance"]["bootstrap_seed"]) + 200,
    )
    practical_lower, _ = simple_bootstrap_interval(
        proteins["practical_fraction_a2_minus_a0.25"].to_numpy(dtype=float),
        replicates=int(config["acceptance"]["bootstrap_replicates"]),
        seed=int(config["acceptance"]["bootstrap_seed"]) + 201,
    )
    protein_fraction = float(proteins["directionality_pass"].mean())
    passed = bool(
        protein_fraction >= float(gates["minimum_protein_pass_fraction"])
        and protein_min
        >= float(gates["minimum_protein_pass_fraction_per_length_bin"])
        and cell_min >= float(gates["minimum_cell_pass_fraction_per_length_bin"])
        and rms_lower > 0.0
        and practical_lower > 0.0
    )
    return (
        passed,
        {
            "protein_count": int(len(proteins)),
            "cell_count": int(len(cells)),
            "protein_pass_fraction": protein_fraction,
            "minimum_protein_pass_fraction_per_length_bin": protein_min,
            "protein_pass_fraction_per_length_bin": protein_bins,
            "minimum_cell_pass_fraction_per_length_bin": cell_min,
            "cell_pass_fraction_per_length_bin": cell_bins,
            "tau_rms_growth_bootstrap_lower": rms_lower,
            "practical_fraction_growth_bootstrap_lower": practical_lower,
        },
        frame,
    )


def recompute_length(
    texture: pd.DataFrame, arm: str, config: dict[str, Any]
) -> tuple[bool, dict[str, Any]]:
    selected = texture[
        (texture["representation_arm"] == arm)
        & (texture["additive_core_mode"] == "actual_additive_core")
        & texture["pooled_supported_reference"].map(bool_value)
        & np.isclose(texture["amplitude"], 1.0)
    ]
    protein = (
        selected.groupby(["reference_id", "length_bin", "length"], as_index=False)[
            "direct_tau_y_rms"
        ]
        .median()
        .query("direct_tau_y_rms > 0")
    )
    x = np.log(protein["length"].to_numpy(dtype=float))
    y = np.log(protein["direct_tau_y_rms"].to_numpy(dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    rng = np.random.default_rng(int(config["acceptance"]["bootstrap_seed"]) + 400)
    boot: list[float] = []
    bins = sorted(str(value) for value in protein["length_bin"].unique())
    for _ in range(int(config["acceptance"]["bootstrap_replicates"])):
        pieces = []
        for length_bin in bins:
            block = protein[protein["length_bin"] == length_bin]
            indices = rng.integers(0, len(block), size=len(block))
            pieces.append(block.iloc[indices])
        replicate = pd.concat(pieces, ignore_index=True)
        bx = np.log(replicate["length"].to_numpy(dtype=float))
        by = np.log(replicate["direct_tau_y_rms"].to_numpy(dtype=float))
        if float(np.std(bx)) > 0:
            boot.append(float(np.polyfit(bx, by, 1)[0]))
    low, high = np.quantile(np.asarray(boot), [0.025, 0.975])
    gates = config["acceptance"]["length_stability"]
    passed = bool(
        len(bins) == int(config["acceptance"]["expected_length_bins"])
        and abs(float(slope)) <= float(gates["slope_abs_max"])
        and float(low) >= -float(gates["slope_ci_abs_max"])
        and float(high) <= float(gates["slope_ci_abs_max"])
    )
    return passed, {
        "representation_arm": arm,
        "n_proteins": int(len(protein)),
        "n_length_bins": len(bins),
        "slope": float(slope),
        "intercept": float(intercept),
        "slope_ci_low": float(low),
        "slope_ci_high": float(high),
        "bootstrap_replicates_realized": len(boot),
    }


def classify_decision(
    *,
    pairing_pass: bool,
    core_pass: bool,
    arm_flags: dict[str, dict[str, bool]],
    superiority: dict[str, bool],
) -> str:
    if not pairing_pass:
        return "EXECUTION_PAIRING_OR_MANIPULATION_INVALID"
    if not core_pass:
        return "FROZEN_EXP062_CORE_OPERATING_POINT_NOT_REPRODUCED"
    if arm_flags["binary_r1"]["package_pass_before_paired_superiority"]:
        return "INCUMBENT_PREVALENCE_FAILURE_NOT_REPLICATED"
    eligible = {
        arm
        for arm in RESCUE_ARMS
        if arm_flags[arm]["package_pass_before_paired_superiority"]
        and superiority[arm]
    }
    if "dense_q20_r1" in eligible and "binary_r8" in eligible:
        return "ALTERNATIVE_SINGLE_FACTOR_REPAIRS_BOTH_SUFFICIENT"
    if "dense_q20_r1" in eligible:
        return "DENSE_Q20_RANK1_SUFFICIENT"
    if "binary_r8" in eligible:
        return "BINARY_RANK8_SUFFICIENT"
    if "dense_q20_r8" in eligible:
        return "COMBINED_DENSE_Q20_RANK8_REQUIRED"
    if any(
        arm_flags[arm]["package_pass_before_paired_superiority"]
        for arm in RESCUE_ARMS
    ):
        return "APPARENT_RESCUE_NOT_PAIRED_SUPERIOR_TO_BINARY_R1"
    if not any(arm_flags[arm]["midpoint_control_pass"] for arm in RESCUE_ARMS):
        return "NO_INTERVENTION_PASSES_MIDPOINT_CONTROL"
    if not any(arm_flags[arm]["edge_pass"] for arm in RESCUE_ARMS):
        return "NO_INTERVENTION_PASSES_EDGE_TRANSFER"
    if not any(arm_flags[arm]["actual_texture_pass"] for arm in RESCUE_ARMS):
        return "NO_INTERVENTION_RESCUES_TEXTURE_PREVALENCE"
    if not any(arm_flags[arm]["directionality_pass"] for arm in RESCUE_ARMS):
        return "REPRESENTATION_RESCUE_NOT_DIRECTIONAL"
    if not any(arm_flags[arm]["length_stability_pass"] for arm in RESCUE_ARMS):
        return "REPRESENTATION_RESCUE_NOT_LENGTH_STABLE"
    return "REPRESENTATION_REPAIR_INSUFFICIENT_AT_FIXED_SIGMOID"


TABLE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
summary = json.loads((RUN_DIR / "summary.json").read_text(encoding="utf-8"))
preflight = json.loads((RUN_DIR / "preflight.json").read_text(encoding="utf-8"))
run_record = yaml.safe_load(RUN_RECORD.read_text(encoding="utf-8"))
core = pd.read_csv(RUN_DIR / "core_guardrail_metrics.csv")
representation = pd.read_csv(RUN_DIR / "representation_audit.csv")
latent = pd.read_csv(RUN_DIR / "latent_geometry_metrics.csv")
edge = pd.read_csv(RUN_DIR / "validation_edge_metrics.csv")
texture = pd.read_csv(RUN_DIR / "validation_texture_metrics.csv")
recorded_checks = pd.read_csv(RUN_DIR / "acceptance_checks.csv")
recorded_decisions = pd.read_csv(RUN_DIR / "decision_table.csv")
recorded_paired = pd.read_csv(RUN_DIR / "paired_arm_effects.csv")
recorded_factorial = pd.read_csv(RUN_DIR / "factorial_effects.csv")
recorded_fits = pd.read_csv(RUN_DIR / "scale_law_fits.csv")

replicates = int(config["acceptance"]["bootstrap_replicates"])
base_seed = int(config["acceptance"]["bootstrap_seed"])
transmission = config["acceptance"]["transmission"]

# Core and pairing gates.
core_min, core_bins = minimum_bin_fraction(core, "core_pass")
core_pass = bool(
    core_min
    >= float(config["acceptance"]["core"]["minimum_cell_pass_fraction_per_length_bin"])
)
pairing_pass, pairing_contexts = recompute_pairing(representation, config)

# Arm gates, summaries, length strata, and dose curves.
arm_flags: dict[str, dict[str, bool]] = {}
arm_values: dict[str, dict[str, Any]] = {}
arm_summary_rows: list[dict[str, Any]] = []
arm_length_rows: list[dict[str, Any]] = []
dose_curve_rows: list[dict[str, Any]] = []
direction_rows: list[pd.DataFrame] = []
length_rows: list[dict[str, Any]] = []

for arm_index, arm in enumerate(ARMS):
    edge_selected = edge[
        (edge["representation_arm"] == arm)
        & np.isclose(
            edge["amplitude"], float(transmission["primary_amplitude"])
        )
        & (edge["sample_id"] == "reference_shell_h1_complete")
    ].copy()
    edge_selected["cell_pass"] = (
        edge_selected["direct_practical_resolution_fraction"]
        >= float(transmission["edge_practical_fraction_min"])
    ) & (
        edge_selected["direct_correlation"]
        >= float(transmission["correlation_min"])
    )
    edge_min, edge_bins = minimum_bin_fraction(edge_selected, "cell_pass")
    edge_pass = bool(
        edge_min
        >= float(transmission["minimum_cell_pass_fraction_per_length_bin"])
    )

    texture_values: dict[str, dict[str, Any]] = {}
    texture_passes: dict[str, bool] = {}
    selected_by_mode: dict[str, pd.DataFrame] = {}
    for core_mode, label, seed_offset in (
        ("actual_additive_core", "actual", 10),
        ("frozen_midpoint_occupancy_0.5", "midpoint", 20),
    ):
        selected = texture[
            (texture["representation_arm"] == arm)
            & (texture["additive_core_mode"] == core_mode)
            & texture["pooled_supported_reference"].map(bool_value)
            & np.isclose(
                texture["amplitude"], float(transmission["primary_amplitude"])
            )
        ].copy()
        selected["cell_pass"] = (
            selected["direct_tau_practical_resolution_fraction"]
            >= float(transmission["texture_practical_fraction_min"])
        ) & (
            selected["direct_tau_correlation"]
            >= float(transmission["correlation_min"])
        ) & (
            selected["direct_tau_projection"] > 0.0
        ) & (
            selected["one_sided_pvalue"]
            <= float(config["acceptance"]["texture_permutation_alpha"])
        )
        minimum, bins = minimum_bin_fraction(selected, "cell_pass")
        lower = crossed_bootstrap_lower(
            selected,
            metric="observed_projection",
            replicates=replicates,
            seed=base_seed + seed_offset + 100 * arm_index,
        )
        passed = bool(
            minimum
            >= float(transmission["minimum_cell_pass_fraction_per_length_bin"])
            and lower
            > float(transmission["projection_bootstrap_lower_min"])
        )
        selected_by_mode[label] = selected
        texture_passes[label] = passed
        texture_values[label] = {
            "cell_count": int(len(selected)),
            "cell_pass_count": int(selected["cell_pass"].sum()),
            "minimum_length_bin_pass_fraction": minimum,
            "per_length_bin_pass_fraction": bins,
            "projection_bootstrap_lower": lower,
        }

    direction_pass, direction_values, arm_direction = recompute_directionality(
        texture, arm, config
    )
    arm_direction["representation_arm"] = arm
    direction_rows.append(arm_direction)
    length_pass, length_values = recompute_length(texture, arm, config)
    length_rows.append(length_values)

    flags = {
        "technical_pass": pairing_pass,
        "edge_pass": edge_pass,
        "actual_texture_pass": texture_passes["actual"],
        "midpoint_control_pass": texture_passes["midpoint"],
        "directionality_pass": direction_pass,
        "length_stability_pass": length_pass,
    }
    flags["package_pass_before_paired_superiority"] = all(flags.values())
    arm_flags[arm] = flags
    arm_values[arm] = {
        "edge": {"minimum": edge_min, "bins": edge_bins},
        "actual": texture_values["actual"],
        "midpoint": texture_values["midpoint"],
        "directionality": direction_values,
        "length": length_values,
    }

    actual = selected_by_mode["actual"]
    midpoint = selected_by_mode["midpoint"]
    latent_selected = latent[
        (latent["representation_arm"] == arm)
        & (latent["additive_core_mode"] == "actual_additive_core")
        & np.isclose(latent["amplitude"], 1.0)
    ]
    arm_summary_rows.append(
        {
            "representation_arm": arm,
            "actual_cells": len(actual),
            "actual_cell_pass_count": int(actual["cell_pass"].sum()),
            "actual_cell_pass_fraction": float(actual["cell_pass"].mean()),
            "actual_practical_fraction_median": float(
                actual["direct_tau_practical_resolution_fraction"].median()
            ),
            "actual_tau_rms_median": float(actual["direct_tau_y_rms"].median()),
            "actual_correlation_median": float(
                actual["direct_tau_correlation"].median()
            ),
            "actual_projection_median": float(
                actual["direct_tau_projection"].median()
            ),
            "actual_projection_bootstrap_lower": texture_values["actual"][
                "projection_bootstrap_lower"
            ],
            "midpoint_cell_pass_count": int(midpoint["cell_pass"].sum()),
            "midpoint_cell_pass_fraction": float(midpoint["cell_pass"].mean()),
            "midpoint_practical_fraction_median": float(
                midpoint["direct_tau_practical_resolution_fraction"].median()
            ),
            "edge_cell_pass_count": int(edge_selected["cell_pass"].sum()),
            "edge_cell_pass_fraction": float(edge_selected["cell_pass"].mean()),
            "edge_practical_fraction_median": float(
                edge_selected["direct_practical_resolution_fraction"].median()
            ),
            "edge_tau_rms_median": float(edge_selected["direct_y_rms"].median()),
            "latent_rms_median": float(latent_selected["latent_epsilon_rms"].median()),
            "latent_abs_median_median": float(
                latent_selected["latent_epsilon_abs_median"].median()
            ),
            "latent_near_zero_fraction_median": float(
                latent_selected["latent_epsilon_near_zero_fraction"].median()
            ),
            "latent_exact_zero_fraction_median": float(
                latent_selected["latent_epsilon_exact_zero_fraction"].median()
            ),
            "latent_excess_kurtosis_median": float(
                latent_selected["latent_epsilon_excess_kurtosis"].median()
            ),
            **flags,
        }
    )

    for length_bin, group in actual.groupby("length_bin", sort=True):
        midpoint_bin = midpoint[midpoint["length_bin"] == length_bin]
        edge_bin = edge_selected[edge_selected["length_bin"] == length_bin]
        arm_length_rows.append(
            {
                "representation_arm": arm,
                "length_bin": length_bin,
                "cell_count": len(group),
                "actual_cell_pass_fraction": float(group["cell_pass"].mean()),
                "actual_practical_fraction_median": float(
                    group["direct_tau_practical_resolution_fraction"].median()
                ),
                "actual_tau_rms_median": float(group["direct_tau_y_rms"].median()),
                "actual_correlation_median": float(
                    group["direct_tau_correlation"].median()
                ),
                "actual_projection_median": float(
                    group["direct_tau_projection"].median()
                ),
                "midpoint_cell_pass_fraction": float(
                    midpoint_bin["cell_pass"].mean()
                ),
                "edge_cell_pass_fraction": float(edge_bin["cell_pass"].mean()),
            }
        )

    curve = texture[
        (texture["representation_arm"] == arm)
        & (texture["additive_core_mode"] == "actual_additive_core")
        & texture["pooled_supported_reference"].map(bool_value)
    ]
    for amplitude, group in curve.groupby("amplitude", sort=True):
        finite_correlations = group["direct_tau_correlation"].dropna()
        dose_curve_rows.append(
            {
                "representation_arm": arm,
                "amplitude": float(amplitude),
                "cell_count": len(group),
                "tau_rms_median": float(group["direct_tau_y_rms"].median()),
                "practical_fraction_median": float(
                    group["direct_tau_practical_resolution_fraction"].median()
                ),
                "correlation_median": (
                    float(finite_correlations.median())
                    if not finite_correlations.empty
                    else None
                ),
            }
        )

# Paired rescue estimates and factorial effects.
selected = texture[
    texture["representation_arm"].isin(ARMS)
    & (texture["additive_core_mode"] == "actual_additive_core")
    & texture["pooled_supported_reference"].map(bool_value)
    & np.isclose(texture["amplitude"], 1.0)
][
    CONTEXT_KEYS
    + ["representation_arm", "direct_tau_practical_resolution_fraction"]
].copy()
pivot = selected.pivot(
    index=CONTEXT_KEYS,
    columns="representation_arm",
    values="direct_tau_practical_resolution_fraction",
).reset_index()

paired_rows: list[dict[str, Any]] = []
superiority: dict[str, bool] = {}
for arm in RESCUE_ARMS:
    frame = pivot[CONTEXT_KEYS].copy()
    frame["difference"] = pivot[arm] - pivot["binary_r1"]
    protein = frame.groupby(
        ["reference_id", "length_bin", "length"], as_index=False
    )["difference"].median()
    lower = crossed_bootstrap_lower(
        frame,
        metric="difference",
        replicates=replicates,
        seed=base_seed + 600 + ARMS.index(arm),
    )
    median_gain = float(protein["difference"].median())
    passed = bool(
        median_gain
        >= float(
            config["acceptance"]["manipulation"][
                "median_practical_fraction_gain_min"
            ]
        )
        and lower
        > float(
            config["acceptance"]["manipulation"][
                "candidate_minus_incumbent_practical_fraction_bootstrap_lower_min"
            ]
        )
    )
    superiority[arm] = passed
    recorded = recorded_paired[
        (recorded_paired["row_type"] == "comparison_summary")
        & (recorded_paired["comparison_arm"] == arm)
    ].iloc[0]
    paired_rows.append(
        {
            "comparison_arm": arm,
            "protein_count": len(protein),
            "median_gain": median_gain,
            "bootstrap_lower": lower,
            "paired_superiority_pass": passed,
            "recorded_median_gain": float(recorded["median_gain"]),
            "recorded_bootstrap_lower": float(recorded["bootstrap_lower"]),
            "recorded_pass": bool_value(recorded["paired_superiority_pass"]),
        }
    )

factorial = pd.DataFrame(
    {
        **{key: pivot[key] for key in CONTEXT_KEYS},
        "dense_main_effect": 0.5
        * (
            (pivot["dense_q20_r1"] - pivot["binary_r1"])
            + (pivot["dense_q20_r8"] - pivot["binary_r8"])
        ),
        "rank_main_effect": 0.5
        * (
            (pivot["binary_r8"] - pivot["binary_r1"])
            + (pivot["dense_q20_r8"] - pivot["dense_q20_r1"])
        ),
        "dense_by_rank_interaction": (
            (pivot["dense_q20_r8"] - pivot["binary_r8"])
            - (pivot["dense_q20_r1"] - pivot["binary_r1"])
        ),
    }
)
factorial_rows: list[dict[str, Any]] = []
for effect_index, effect in enumerate(
    ("dense_main_effect", "rank_main_effect", "dense_by_rank_interaction")
):
    protein = factorial.groupby(
        ["reference_id", "length_bin", "length"], as_index=False
    )[effect].median()
    low, high = simple_bootstrap_interval(
        protein[effect].to_numpy(dtype=float),
        replicates=replicates,
        seed=base_seed + 700 + effect_index,
    )
    recorded = recorded_factorial[
        (recorded_factorial["row_type"] == "factorial_effect_summary")
        & (recorded_factorial["effect"] == effect)
    ].iloc[0]
    factorial_rows.append(
        {
            "effect": effect,
            "protein_count": len(protein),
            "median": float(protein[effect].median()),
            "bootstrap_low": low,
            "bootstrap_high": high,
            "recorded_median": float(recorded["median"]),
            "recorded_bootstrap_low": float(recorded["bootstrap_lower"]),
            "recorded_bootstrap_high": float(recorded["bootstrap_upper"]),
        }
    )

recomputed_decision = classify_decision(
    pairing_pass=pairing_pass,
    core_pass=core_pass,
    arm_flags=arm_flags,
    superiority=superiority,
)

# Compare all recorded acceptance rows to independent reaggregation.
recomputed_checks: dict[str, bool] = {
    "shared_frozen_core_guardrail": core_pass,
    "cross_arm_pairing_and_manipulation": pairing_pass,
    "registered_terminal_decision_resolved": recomputed_decision
    in set(config["decision"]["precedence"]),
}
for arm in ARMS:
    recomputed_checks.update(
        {
            f"{arm}_technical_and_pairing": arm_flags[arm]["technical_pass"],
            f"{arm}_complete_h1_edge_transfer": arm_flags[arm]["edge_pass"],
            f"{arm}_actual_core_texture": arm_flags[arm][
                "actual_texture_pass"
            ],
            f"{arm}_midpoint_frozen_positive_control": arm_flags[arm][
                "midpoint_control_pass"
            ],
            f"{arm}_dose_directionality": arm_flags[arm]["directionality_pass"],
            f"{arm}_length_stability": arm_flags[arm]["length_stability_pass"],
        }
    )
for arm in RESCUE_ARMS:
    recomputed_checks[f"{arm}_paired_superiority"] = superiority[arm]

acceptance_rows: list[dict[str, Any]] = []
for _, row in recorded_checks.iterrows():
    check_id = str(row["check_id"])
    recomputed = recomputed_checks[check_id]
    recorded = bool_value(row["passed"])
    acceptance_rows.append(
        {
            "check_id": check_id,
            "representation_arm": row["representation_arm"],
            "family": row["family"],
            "required": bool_value(row["required"]),
            "recorded_pass": recorded,
            "recomputed_pass": recomputed,
            "matches": recorded == recomputed,
        }
    )

# Technical summaries.
technical = representation[
    representation["row_type"] == "technical_context_audit"
].copy()
arm_audit = representation[
    representation["row_type"] == "arm_context_audit"
].copy()
technical_rows: list[dict[str, Any]] = []
for arm in ARMS:
    tech = technical[technical["representation_arm"] == arm]
    audit = arm_audit[arm_audit["representation_arm"] == arm]
    technical_rows.append(
        {
            "representation_arm": arm,
            "context_count": len(tech),
            "representation_rank": int(audit["representation_rank"].median()),
            "component_count_median": float(audit["component_count"].median()),
            "sampled_nonself_tie_fraction_median": float(
                audit["sampled_nonself_tie_fraction"].median()
            ),
            "minimum_operator_gram_rank": int(
                audit["minimum_sampled_operator_gram_rank"].min()
            ),
            "operator_gram_participation_ratio_median": float(
                audit["median_sampled_operator_gram_participation_ratio"].median()
            ),
            "operator_energy_relative_error_max": float(
                audit["maximum_operator_energy_relative_error"].max()
            ),
            "complete_h1_rms_ratio_min": float(
                tech["complete_h1_gnk_mutation_rms_ratio"].min()
            ),
            "complete_h1_rms_ratio_median": float(
                tech["complete_h1_gnk_mutation_rms_ratio"].median()
            ),
            "complete_h1_rms_ratio_max": float(
                tech["complete_h1_gnk_mutation_rms_ratio"].max()
            ),
            "formula_parity_max_abs": float(
                tech["maximum_formula_parity_abs_difference"].max()
            ),
            "zero_amplitude_residual_rms_max": float(
                tech["maximum_zero_amplitude_residual_rms"].max()
            ),
            "surrogate_latent_epsilon_max_abs": float(
                tech["surrogate_latent_epsilon_max_abs"].max()
            ),
        }
    )

# Source, protocol, strict JSON, and content-addressed collection audits.
source_rows: list[dict[str, Any]] = []
source_rows.append(
    {
        "path": str(CONFIG_PATH.relative_to(PROJECT_ROOT)),
        "recorded_sha256": preflight["config_file_sha256"],
        "current_sha256": sha256_file(CONFIG_PATH),
        "matches": preflight["config_file_sha256"] == sha256_file(CONFIG_PATH),
    }
)
reference_catalog_path = PROJECT_ROOT / str(config["references"]["catalog"])
reference_catalog_hash = sha256_file(reference_catalog_path)
source_rows.append(
    {
        "path": str(reference_catalog_path.relative_to(PROJECT_ROOT)),
        "recorded_sha256": preflight["reference_catalog_sha256"],
        "current_sha256": reference_catalog_hash,
        "matches": reference_catalog_hash
        == preflight["reference_catalog_sha256"],
    }
)
for relative, recorded_hash in sorted(
    preflight["implementation_hash_components"].items()
):
    path = PROJECT_ROOT / relative
    current_hash = sha256_file(path)
    source_rows.append(
        {
            "path": relative,
            "recorded_sha256": recorded_hash,
            "current_sha256": current_hash,
            "matches": current_hash == recorded_hash,
        }
    )

expected_outputs = [str(value) for value in config["outputs"]["expected"]]
protocol_rows: list[dict[str, Any]] = []
for name in expected_outputs:
    path = RUN_DIR / name
    row_count = len(pd.read_csv(path)) if path.suffix == ".csv" else None
    protocol_rows.append(
        {
            "artifact": name,
            "expected": True,
            "observed": path.is_file(),
            "row_count": row_count,
        }
    )
checkpoints = sorted((RUN_DIR / "checkpoints").glob("*.json"))
logs = sorted((RUN_DIR / "logs").glob("*"))
protocol_rows.extend(
    [
        {
            "artifact": "checkpoints/*.json",
            "expected": 360,
            "observed": len(checkpoints),
            "row_count": None,
        },
        {
            "artifact": "logs/*",
            "expected": 3,
            "observed": len(logs),
            "row_count": None,
        },
        {
            "artifact": "logs/labproj-exit-code.txt",
            "expected": "0",
            "observed": (RUN_DIR / "logs/labproj-exit-code.txt")
            .read_text(encoding="utf-8")
            .strip(),
            "row_count": None,
        },
        {
            "artifact": "logs/stderr.txt bytes",
            "expected": 0,
            "observed": (RUN_DIR / "logs/stderr.txt").stat().st_size,
            "row_count": None,
        },
    ]
)

strict_json_failures: list[str] = []
nonstandard_tokens = 0
for path in sorted(RUN_DIR.rglob("*.json")):
    text = path.read_text(encoding="utf-8")
    nonstandard_tokens += sum(
        text.count(token) for token in ("NaN", "Infinity", "-Infinity")
    )
    try:
        json.loads(
            text,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"non-standard JSON constant: {value}")
            ),
        )
    except (json.JSONDecodeError, ValueError):
        strict_json_failures.append(str(path.relative_to(RUN_DIR)))

pointer = yaml.safe_load(DVC_POINTER.read_text(encoding="utf-8"))["outs"][0]
dvc_hash = str(pointer["md5"])
dvc_manifest_path = (
    PROJECT_ROOT / ".dvc/cache/files/md5" / dvc_hash[:2] / dvc_hash[2:]
)
dvc_entries = json.loads(dvc_manifest_path.read_text(encoding="utf-8"))
local_files = sorted(path for path in RUN_DIR.rglob("*") if path.is_file())
local_bytes = sum(path.stat().st_size for path in local_files)
local_md5_match = True
cache_md5_match = True
for entry in dvc_entries:
    relative = Path(entry["relpath"])
    expected_md5 = str(entry["md5"])
    local_path = RUN_DIR / relative
    cache_path = (
        PROJECT_ROOT
        / ".dvc/cache/files/md5"
        / expected_md5[:2]
        / expected_md5[2:]
    )
    local_md5_match &= (
        local_path.is_file() and md5_file(local_path) == expected_md5
    )
    cache_md5_match &= (
        cache_path.is_file() and md5_file(cache_path) == expected_md5
    )

manifest_lines: list[str] = []
for path in local_files:
    relative = path.relative_to(RUN_DIR).as_posix()
    manifest_lines.append(f"{sha256_file(path)}  {path.stat().st_size}  {relative}")
local_manifest_sha256 = hashlib.sha256(
    ("\n".join(manifest_lines) + "\n").encode("utf-8")
).hexdigest()

remote = (
    json.loads(REMOTE_AUDIT.read_text(encoding="utf-8"))
    if REMOTE_AUDIT.exists()
    else {}
)
collection_rows = [
    {
        "check": "dvc_pointer_nfiles",
        "expected": int(pointer["nfiles"]),
        "observed": len(local_files),
        "passed": len(local_files) == int(pointer["nfiles"]),
        "note": "",
    },
    {
        "check": "dvc_pointer_size_bytes",
        "expected": int(pointer["size"]),
        "observed": local_bytes,
        "passed": local_bytes == int(pointer["size"]),
        "note": "",
    },
    {
        "check": "dvc_local_md5",
        "expected": True,
        "observed": local_md5_match,
        "passed": local_md5_match,
        "note": "",
    },
    {
        "check": "dvc_cache_md5",
        "expected": True,
        "observed": cache_md5_match,
        "passed": cache_md5_match,
        "note": "",
    },
    {
        "check": "strict_json",
        "expected": 0,
        "observed": len(strict_json_failures),
        "passed": not strict_json_failures,
        "note": f"non-standard token count={nonstandard_tokens}",
    },
    {
        "check": "remote_runtime_sha256",
        "expected": True,
        "observed": remote.get("all_sha256_and_size_match"),
        "passed": bool(remote.get("all_sha256_and_size_match", False)),
        "note": (
            f"remote files={remote.get('file_count')}; "
            f"canonical manifest={remote.get('manifest_sha256')}"
        ),
    },
    {
        "check": "configured_dvc_remote_content_objects",
        "expected": len(dvc_entries),
        "observed": remote.get("dvc_remote_present_content_objects"),
        "passed": remote.get("dvc_remote_present_content_objects")
        == len(dvc_entries),
        "note": (
            f"missing={remote.get('dvc_remote_missing_content_objects')}; "
            f"directory object present={remote.get('dvc_remote_dir_object_present')}"
        ),
    },
    {
        "check": "configured_dvc_backup",
        "expected": True,
        "observed": False,
        "passed": False,
        "note": "labproj collection and local DVC tracking succeeded; redundant push to /Volumes/DutchZFS failed with permission denied",
    },
]

# Persist compact source copies needed to rerun the interpretation.
(SOURCE_DIR / "run-157_summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
(SOURCE_DIR / "run-157_preflight.json").write_text(
    json.dumps(preflight, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)

pd.DataFrame(acceptance_rows).to_csv(
    TABLE_DIR / "run-157_acceptance_audit.csv", index=False
)
pd.DataFrame(arm_summary_rows).to_csv(
    TABLE_DIR / "run-157_arm_summary.csv", index=False
)
pd.DataFrame(arm_length_rows).to_csv(
    TABLE_DIR / "run-157_arm_by_length.csv", index=False
)
pd.DataFrame(paired_rows).to_csv(
    TABLE_DIR / "run-157_paired_effects.csv", index=False
)
pd.DataFrame(factorial_rows).to_csv(
    TABLE_DIR / "run-157_factorial_effects.csv", index=False
)
pd.DataFrame(dose_curve_rows).to_csv(
    TABLE_DIR / "run-157_dose_curve.csv", index=False
)
pd.DataFrame(length_rows).to_csv(
    TABLE_DIR / "run-157_scale_law_summary.csv", index=False
)
pd.DataFrame(technical_rows).to_csv(
    TABLE_DIR / "run-157_representation_audit.csv", index=False
)
pd.DataFrame(protocol_rows).to_csv(
    TABLE_DIR / "run-157_protocol_output_audit.csv", index=False
)
pd.DataFrame(source_rows).to_csv(
    TABLE_DIR / "run-157_source_hash_audit.csv", index=False
)
pd.DataFrame(collection_rows).to_csv(
    TABLE_DIR / "run-157_collection_audit.csv", index=False
)

arm_summary = {
    row["representation_arm"]: clean(row) for row in arm_summary_rows
}
metrics = {
    "analysis_id": "ANA-061",
    "experiment_id": "EXP-063",
    "run_id": "RUN-157",
    "scheduler_job_id": run_record["scheduler_job_id"],
    "scheduler_exit_code": run_record["scheduler_exit_code"],
    "registered_decision": summary["decision"],
    "recomputed_decision": recomputed_decision,
    "decision_matches": summary["decision"] == recomputed_decision,
    "scientific_pass": summary["scientific_pass"],
    "required_checks": {
        "passed": summary["required_checks_passed"],
        "total": summary["required_check_count"],
        "all_recorded_states_reproduced": all(
            row["matches"] for row in acceptance_rows
        ),
    },
    "shared_core": {
        "cells": len(core),
        "passes": int(core["core_pass"].sum()),
        "minimum_length_bin_pass_fraction": core_min,
        "per_length_bin_pass_fraction": core_bins,
    },
    "pairing": {
        "context_count": pairing_contexts,
        "pass": pairing_pass,
    },
    "arms": arm_summary,
    "paired_effects": {
        row["comparison_arm"]: clean(row) for row in paired_rows
    },
    "factorial_effects": {
        row["effect"]: clean(row) for row in factorial_rows
    },
    "collection": {
        "dvc_dir_hash": dvc_hash,
        "file_count": len(local_files),
        "size_bytes": local_bytes,
        "local_and_cache_md5_match": local_md5_match and cache_md5_match,
        "strict_json_failures": strict_json_failures,
        "nonstandard_json_token_count": nonstandard_tokens,
        "local_manifest_sha256": local_manifest_sha256,
        "remote_manifest_sha256": remote.get("manifest_sha256"),
        "remote_sha256_match": bool(
            remote.get("all_sha256_and_size_match", False)
        ),
        "dvc_remote_present_content_objects": remote.get(
            "dvc_remote_present_content_objects"
        ),
        "dvc_remote_missing_content_objects": remote.get(
            "dvc_remote_missing_content_objects"
        ),
        "dvc_remote_dir_object_present": remote.get(
            "dvc_remote_dir_object_present"
        ),
        "configured_backup_succeeded": False,
    },
    "hypothesis_effect": {
        "supports": ["HYP-005"],
        "weakens": [],
        "refutes": [],
        "scope": "narrow experiment-local observable texture-rescue evidence on a reused 12-protein challenge panel",
    },
    "candidate_is_experiment_local": True,
    "public_model_changed": False,
    "independent_protein_generalization_tested": False,
    "empirical_realism_tested": False,
    "optimizer_or_accessibility_tested": False,
    "stage_4c_4d_4e_hard_stop_retained": True,
}
(TABLE_DIR / "metrics.json").write_text(
    json.dumps(clean(metrics), indent=2, sort_keys=True) + "\n",
    encoding="utf-8",
)

print(json.dumps(clean(metrics), indent=2, sort_keys=True))
