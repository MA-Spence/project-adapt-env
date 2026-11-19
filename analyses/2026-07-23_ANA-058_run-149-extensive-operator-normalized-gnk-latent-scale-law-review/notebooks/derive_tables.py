from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
RUN_DIR = (
    PROJECT_ROOT
    / "data/processed/extensive-operator-normalized-gnk-latent-scale-law-validation/RUN-149"
)
TABLE_DIR = ANALYSIS_DIR / "tables"
SOURCE_DIR = ANALYSIS_DIR / "source_tables"

OUTPUTS = [
    "preflight.json",
    "reference_catalog.csv",
    "pairing_audit.csv",
    "support_audit.csv",
    "contrast_audit.csv",
    "latent_scale_metrics.csv",
    "order_component_metrics.csv",
    "scale_law_fits.csv",
    "acceptance_checks.csv",
    "summary.json",
]
LOGS = ["logs/labproj-exit-code.txt", "logs/stderr.txt", "logs/stdout.txt"]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def log_slope(frame: pd.DataFrame, metric: str) -> float:
    aggregate = frame.groupby("length", sort=True)[metric].median()
    return float(
        np.polyfit(
            np.log(aggregate.index.to_numpy(dtype=float)),
            np.log(aggregate.to_numpy(dtype=float)),
            1,
        )[0]
    )


def select(
    frame: pd.DataFrame,
    *,
    arm: str,
    regime: str,
    distance: str,
    weight: float = 1.0,
) -> pd.DataFrame:
    return frame[
        (frame["arm"] == arm)
        & (frame["regime"] == regime)
        & (frame["distance_label"] == distance)
        & np.isclose(frame["weight"], weight)
    ].copy()


TABLE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

summary = json.loads((RUN_DIR / "summary.json").read_text())
preflight = json.loads((RUN_DIR / "preflight.json").read_text())
checks = pd.read_csv(RUN_DIR / "acceptance_checks.csv")
fits = pd.read_csv(RUN_DIR / "scale_law_fits.csv")
metrics = pd.read_csv(RUN_DIR / "latent_scale_metrics.csv")
orders = pd.read_csv(RUN_DIR / "order_component_metrics.csv")
supports = pd.read_csv(RUN_DIR / "support_audit.csv")
contrasts = pd.read_csv(RUN_DIR / "contrast_audit.csv")
pairing = pd.read_csv(RUN_DIR / "pairing_audit.csv")
references = pd.read_csv(RUN_DIR / "reference_catalog.csv")

# Collection and execution audit.
collection = pd.DataFrame(
    [
        {
            "relative_path": relative,
            "exists": (RUN_DIR / relative).is_file(),
            "size_bytes": (RUN_DIR / relative).stat().st_size,
            "sha256": file_sha256(RUN_DIR / relative),
        }
        for relative in OUTPUTS + LOGS
    ]
)
collection.to_csv(TABLE_DIR / "run-149_collection_audit.csv", index=False)

checkpoint_count = len(list((RUN_DIR / "checkpoints").glob("*.json")))
assert collection["exists"].all()
assert (RUN_DIR / "logs/labproj-exit-code.txt").read_text().strip() == "0"
assert (RUN_DIR / "logs/stderr.txt").stat().st_size == 0
assert checkpoint_count == 180
assert summary["quick"] is False and summary["scientific_decision_enabled"] is True
assert summary["task_count"] == 180
assert len(metrics) == summary["metric_row_count"] == 32_400
assert len(pairing) == 2_160
assert len(supports) == 1_620
assert len(contrasts) == 1_620
assert len(orders) == 540
assert references["length"].tolist() == [48, 93, 118, 220, 238, 448]
assert not summary["observable_fitness_metrics"]
assert not summary["gnk_transmitted_through_sigmoid"]
assert not summary["accessibility_or_optimizer_metrics"]
assert not summary["public_model_changed"]

for name in ("summary.json", "preflight.json"):
    shutil.copy2(RUN_DIR / name, SOURCE_DIR / f"run-149_{name}")

# Registered gate audit.
checks["passed"] = checks["passed"].astype(bool)
checks["required"] = checks["required"].astype(bool)
required = checks[checks["required"]]
advisory = checks[~checks["required"]]
assert len(required) == 26 and required["passed"].all()
assert len(advisory) == 2 and advisory["passed"].all()
assert summary["scientific_pass"] is True
assert summary["decision"] == "LATENT_GNK_EXTENSIVE_SCALE_LAW_VALIDATED"
checks.to_csv(TABLE_DIR / "run-149_acceptance_audit.csv", index=False)

# Local mutation scale.
local_all = metrics[
    (metrics["regime"] == "uniform_background")
    & (metrics["distance_label"] == "fixed_1")
    & np.isclose(metrics["weight"], 1.0)
].copy()
local_ratio_rows: list[dict[str, object]] = []
for (arm, length), group in local_all.groupby(["arm", "length"], sort=True):
    ratio = group["gnk_to_additive_rms_ratio"].astype(float)
    local_ratio_rows.append(
        {
            "arm": arm,
            "length": int(length),
            "n_seeds": int(group["seed"].nunique()),
            "median_ratio": float(ratio.median()),
            "minimum_ratio": float(ratio.min()),
            "maximum_ratio": float(ratio.max()),
            "fraction_in_0.70_1.30": float(((ratio >= 0.70) & (ratio <= 1.30)).mean()),
        }
    )
local_ratio_table = pd.DataFrame(local_ratio_rows)
local_ratio_table.to_csv(TABLE_DIR / "run-149_local_ratio_by_length.csv", index=False)

fit_index = fits.set_index("fit_id")
candidate_local = select(
    metrics,
    arm="extensive_balanced_bounded",
    regime="uniform_background",
    distance="fixed_1",
)
normalized_local = select(
    metrics,
    arm="legacy_support_local_normalized",
    regime="uniform_background",
    distance="fixed_1",
)
legacy_local = select(
    metrics,
    arm="legacy_public_global_fixed",
    regime="uniform_background",
    distance="fixed_1",
)
candidate_local_slope = log_slope(candidate_local, "gnk_delta_rms")
normalized_local_slope = log_slope(normalized_local, "gnk_delta_rms")
legacy_local_slope = log_slope(legacy_local, "gnk_delta_rms")
assert np.isclose(candidate_local_slope, fit_index.loc["extensive_balanced_bounded_local", "slope"])
assert np.isclose(normalized_local_slope, fit_index.loc["legacy_support_local_normalized_local", "slope"])
assert np.isclose(legacy_local_slope, fit_index.loc["legacy_public_global_fixed_local", "slope"])

candidate_ratio = candidate_local["gnk_to_additive_rms_ratio"].astype(float)
candidate_ratio_medians = candidate_local.groupby("length")["gnk_to_additive_rms_ratio"].median()
assert ((candidate_ratio_medians >= 0.85) & (candidate_ratio_medians <= 1.15)).all()
assert ((candidate_ratio >= 0.70) & (candidate_ratio <= 1.30)).all()

# Scalar weight linearity, reconstructed relative to weight one.
linearity_errors: list[float] = []
for _, group in metrics.groupby(
    ["reference_id", "seed", "arm", "regime", "distance_label"], sort=False
):
    unit = float(group.loc[np.isclose(group["weight"], 1.0), "gnk_delta_rms"].iloc[0])
    for row in group.itertuples(index=False):
        expected = unit * float(row.weight)
        error = abs(float(row.gnk_delta_rms) - expected)
        linearity_errors.append(error if expected == 0.0 else error / abs(expected))
max_linearity_error = float(max(linearity_errors))
assert max_linearity_error < 1.0e-12

# Candidate and normalized-control scale laws.
field_rows = candidate_local.assign(
    variance_per_length=candidate_local["expected_field_variance"] / candidate_local["length"]
)
variance_per_length = field_rows.groupby("length")["variance_per_length"].median()
candidate_field_slope = float(fit_index.loc["candidate_field_variance", "slope"])
variance_per_length_cv = float(np.std(variance_per_length) / np.mean(variance_per_length))

candidate_fractional_rms_slopes = {
    label: float(fit_index.loc[f"candidate_reference_shell_delta_{label}", "slope"])
    for label in ("fraction_0.1", "fraction_0.5")
}
candidate_fractional_seed_slope_summary: dict[str, dict[str, float | int]] = {}
for label in ("fraction_0.1", "fraction_0.5"):
    fractional_rows = select(
        metrics,
        arm="extensive_balanced_bounded",
        regime="reference_shell",
        distance=label,
    )
    seed_slopes = []
    for _, seed_rows in fractional_rows.groupby("seed", sort=True):
        seed_slopes.append(log_slope(seed_rows, "gnk_delta_rms"))
    seed_slopes_array = np.asarray(seed_slopes, dtype=float)
    in_band = (seed_slopes_array >= 0.35) & (seed_slopes_array <= 0.65)
    candidate_fractional_seed_slope_summary[label] = {
        "in_registered_band": int(in_band.sum()),
        "total": int(seed_slopes_array.size),
        "fraction_in_registered_band": float(in_band.mean()),
        "minimum": float(seed_slopes_array.min()),
        "median": float(np.median(seed_slopes_array)),
        "maximum": float(seed_slopes_array.max()),
    }
candidate_fractional_sd_slopes = {
    label: log_slope(
        select(
            metrics,
            arm="extensive_balanced_bounded",
            regime="reference_shell",
            distance=label,
        ),
        "gnk_delta_sd",
    )
    for label in ("fraction_0.1", "fraction_0.5")
}
background_ratio_slopes = {
    label: float(fit_index.loc[f"candidate_ratio_{label}", "slope"])
    for label in ("fraction_0.1", "fraction_0.5")
}
fixed_slopes = {
    "fixed_1": log_slope(
        select(
            metrics,
            arm="extensive_balanced_bounded",
            regime="reference_shell",
            distance="fixed_1",
        ),
        "gnk_delta_rms",
    ),
    **{
        label: float(fit_index.loc[f"candidate_reference_shell_delta_{label}", "slope"])
        for label in ("fixed_2", "fixed_4", "fixed_8")
    },
}
normalized_control_scaling = {
    "empirical_field_variance_slope": log_slope(normalized_local, "source_field_variance"),
    **{
        f"uniform_background_delta_{label}_slope": log_slope(
            select(
                metrics,
                arm="legacy_support_local_normalized",
                regime="uniform_background",
                distance=label,
            ),
            "gnk_delta_rms",
        )
        for label in ("fraction_0.1", "fraction_0.5")
    },
    **{
        f"reference_shell_delta_{label}_slope": log_slope(
            select(
                metrics,
                arm="legacy_support_local_normalized",
                regime="reference_shell",
                distance=label,
            ),
            "gnk_delta_rms",
        )
        for label in ("fraction_0.1", "fraction_0.5")
    },
}

assert abs(candidate_local_slope) <= 0.10
assert -0.70 <= legacy_local_slope <= -0.30
assert 0.80 <= candidate_field_slope <= 1.20 and variance_per_length_cv <= 0.20
assert all(0.35 <= value <= 0.65 for value in candidate_fractional_rms_slopes.values())
assert candidate_fractional_seed_slope_summary["fraction_0.1"]["in_registered_band"] == 30
assert candidate_fractional_seed_slope_summary["fraction_0.5"]["in_registered_band"] == 24
assert all(abs(value) <= 0.15 for value in background_ratio_slopes.values())
assert all(abs(value) <= 0.15 for value in fixed_slopes.values())

selected_fit_ids = [
    "legacy_public_global_fixed_local",
    "legacy_support_local_normalized_local",
    "extensive_balanced_bounded_local",
    "candidate_field_variance",
    "candidate_empirical_field_variance",
    "candidate_reference_shell_delta_fixed_2",
    "candidate_reference_shell_delta_fixed_4",
    "candidate_reference_shell_delta_fixed_8",
    "candidate_ratio_fraction_0.1",
    "candidate_ratio_fraction_0.5",
    "candidate_reference_shell_delta_fraction_0.1",
    "candidate_reference_shell_delta_fraction_0.5",
]
fits[fits["fit_id"].isin(selected_fit_ids)].to_csv(
    TABLE_DIR / "run-149_scale_law_summary.csv", index=False
)

# Candidate support and contrast audit.
candidate_support = supports[supports["arm"] == "extensive_balanced_bounded"].copy()
candidate_support["degree_range"] = candidate_support["degree_max"] - candidate_support["degree_min"]
support_table = (
    candidate_support.groupby(["length", "order"], as_index=False)
    .agg(
        term_count=("term_count", "first"),
        maximum_degree_range=("degree_range", "max"),
        maximum_degree_cv=("degree_cv", "max"),
        maximum_duplicate_fraction=("duplicate_support_fraction", "max"),
        uncovered_sites=("uncovered_sites", "sum"),
        n_seeds=("seed", "nunique"),
    )
)
support_table.to_csv(TABLE_DIR / "run-149_candidate_support_by_length_order.csv", index=False)

max_degree_range = int(candidate_support["degree_range"].max())
max_degree_cv = float(candidate_support["degree_cv"].max())
max_duplicate = float(candidate_support["duplicate_support_fraction"].max())
uncovered = int(candidate_support["uncovered_sites"].sum())
assert max_degree_range == 1 and max_degree_cv <= 0.10
assert max_duplicate == 0.0 and uncovered == 0

candidate_contrast = contrasts[contrasts["arm"] == "extensive_balanced_bounded"]
contrast_mean_error = float(candidate_contrast["max_abs_mean"].max())
contrast_second_error = float(np.max(np.abs(candidate_contrast["mean_second_moment"] - 1.0)))
contrast_max_abs = float(candidate_contrast["max_abs_value"].max())
contrast_substitution_error = float(
    np.max(
        np.abs(
            candidate_contrast["mean_substitution_energy"]
            - candidate_contrast["expected_q20_substitution_energy"]
        )
    )
)
assert contrast_mean_error <= 1.0e-12
assert contrast_second_error <= 1.0e-12
assert contrast_max_abs <= 1.0
assert contrast_substitution_error <= 1.0e-12

# Despite the legacy output column names, these are mutation-energy shares.
order_table = (
    orders.groupby(["length", "order"], as_index=False)
    .agg(
        n_seeds=("seed", "nunique"),
        target_mutation_energy_fraction=("target_variance_fraction", "first"),
        median_sample_mutation_energy_fraction=("sample_variance_fraction", "median"),
        median_curvature_rms=("curvature_rms", "median"),
        minimum_curvature_rms=("curvature_rms", "min"),
    )
)
order_table["absolute_mutation_energy_share_error"] = np.abs(
    order_table["median_sample_mutation_energy_fraction"]
    - order_table["target_mutation_energy_fraction"]
)
curvature_slopes: dict[str, float] = {}
for order, group in orders.groupby("order", sort=True):
    slope = log_slope(group, "curvature_rms")
    curvature_slopes[str(int(order))] = slope
    order_table.loc[order_table["order"] == order, "curvature_log_length_slope"] = slope
order_table.to_csv(TABLE_DIR / "run-149_order_allocation_and_curvature.csv", index=False)

max_order_energy_share_error = float(order_table["absolute_mutation_energy_share_error"].max())
minimum_curvature_rms = float(orders["curvature_rms"].min())
assert max_order_energy_share_error <= 0.10
assert minimum_curvature_rms > 0.0
assert all(abs(value) <= 0.15 for value in curvature_slopes.values())

target_energy = {2: 0.35, 3: 0.45, 4: 0.20}
field_weights = {order: share / order for order, share in target_energy.items()}
field_weight_total = sum(field_weights.values())
actual_field_variance_shares = {
    str(order): value / field_weight_total for order, value in field_weights.items()
}

candidate_tail = candidate_local["gnk_delta_pearson_kurtosis"].astype(float)
tail_by_arm = {
    str(arm): {
        "median": float(group["gnk_delta_pearson_kurtosis"].median()),
        "p95": float(group["gnk_delta_pearson_kurtosis"].quantile(0.95)),
        "maximum": float(group["gnk_delta_pearson_kurtosis"].max()),
    }
    for arm, group in local_all.groupby("arm", sort=True)
}
legacy_l448_uncovered = {
    str(int(order)): float(value)
    for order, value in (
        supports[
            (supports["arm"] == "legacy_public_global_fixed")
            & (supports["length"] == 448)
        ]
        .groupby("order")["uncovered_sites"]
        .median()
        .items()
    )
}

assert not pairing["arm_or_weight_in_sampling_key"].astype(bool).any()
assert pairing["n_arms"].eq(3).all() and pairing["n_pairs"].eq(512).all()
assert pairing["pair_hash"].nunique() == len(pairing)

structural = pd.DataFrame(
    [
        ("candidate_support_max_degree_range", max_degree_range),
        ("candidate_support_max_degree_cv", max_degree_cv),
        ("candidate_support_max_duplicate_fraction", max_duplicate),
        ("candidate_support_uncovered_sites", uncovered),
        ("candidate_contrast_max_abs_mean", contrast_mean_error),
        ("candidate_contrast_second_moment_error", contrast_second_error),
        ("candidate_contrast_max_abs_value", contrast_max_abs),
        ("candidate_contrast_substitution_energy_error", contrast_substitution_error),
        ("maximum_order_mutation_energy_share_error", max_order_energy_share_error),
        ("minimum_supported_pair_curvature_rms", minimum_curvature_rms),
        ("candidate_local_effect_pearson_kurtosis_median", float(candidate_tail.median())),
        ("maximum_weight_linearity_relative_error", max_linearity_error),
    ],
    columns=["metric", "value"],
)
structural.to_csv(TABLE_DIR / "run-149_structural_audit.csv", index=False)

audited = {
    "analysis_id": "ANA-058",
    "experiment_id": "EXP-060",
    "run_id": "RUN-149",
    "decision": summary["decision"],
    "scientific_pass": True,
    "checks": {
        "required": {"passed": 26, "total": 26},
        "advisory": {"passed": 2, "total": 2},
    },
    "execution": {
        "references": 6,
        "seeds": 30,
        "tasks": 180,
        "metric_rows": 32_400,
        "pair_designs": 2_160,
        "checkpoints": checkpoint_count,
        "stderr_bytes": 0,
        "exit_code": 0,
    },
    "local_scale": {
        "candidate_slope": candidate_local_slope,
        "candidate_bootstrap_ci": [
            float(fit_index.loc["extensive_balanced_bounded_local", "bootstrap_ci_low"]),
            float(fit_index.loc["extensive_balanced_bounded_local", "bootstrap_ci_high"]),
        ],
        "candidate_length_median_ratios": {
            str(int(length)): float(value) for length, value in candidate_ratio_medians.items()
        },
        "candidate_seed_fraction_in_0.70_1.30": 1.0,
        "legacy_public_slope": legacy_local_slope,
        "legacy_support_local_normalized_slope": normalized_local_slope,
    },
    "candidate_field_and_distance_scaling": {
        "analytic_field_variance_slope": candidate_field_slope,
        "variance_per_length_cv": variance_per_length_cv,
        "reference_shell_fractional_rms_slopes": candidate_fractional_rms_slopes,
        "reference_shell_fractional_seed_slope_summary": (
            candidate_fractional_seed_slope_summary
        ),
        "reference_shell_fractional_sd_slopes": candidate_fractional_sd_slopes,
        "uniform_background_fractional_ratio_slopes": background_ratio_slopes,
        "reference_shell_fixed_distance_slopes": fixed_slopes,
    },
    "diagnostic_comparator": {
        "legacy_support_local_normalized": normalized_control_scaling,
        "tail_by_arm": tail_by_arm,
        "legacy_L448_median_uncovered_sites_by_order": legacy_l448_uncovered,
    },
    "structure": {
        "support_max_degree_range": max_degree_range,
        "support_max_degree_cv": max_degree_cv,
        "support_duplicate_fraction": max_duplicate,
        "support_uncovered_sites": uncovered,
        "contrast_max_abs": contrast_max_abs,
        "contrast_substitution_energy_error": contrast_substitution_error,
        "maximum_order_mutation_energy_share_error": max_order_energy_share_error,
        "actual_uniform_background_field_variance_shares_by_order": actual_field_variance_shares,
        "curvature_slopes_by_order": curvature_slopes,
        "minimum_supported_pair_curvature_rms": minimum_curvature_rms,
    },
    "scope_and_deviations": {
        "latent_only": True,
        "sigmoid_transmission_tested": False,
        "observable_fitness_tested": False,
        "public_model_changed": False,
        "protein_identity_confounded_with_length": True,
        "component_necessity_ablated": False,
        "fixed_distance_h1_omitted_from_runner_gate": True,
        "fixed_distance_h1_independent_slope": fixed_slopes["fixed_1"],
        "deterministic_replay_is_same_platform": True,
        "implementation_hash_excludes_adaptenv_dependency": True,
        "order_share_columns_are_mutation_energy_not_field_variance": True,
    },
}
(TABLE_DIR / "metrics.json").write_text(json.dumps(audited, indent=2, sort_keys=True) + "\n")
print(json.dumps(audited, indent=2, sort_keys=True))
