from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
RUN_DIR = (
    PROJECT_ROOT
    / "data/processed/operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/RUN-151"
)
OUTPUT_DIR = RUN_DIR / "outputs"
TABLE_DIR = ANALYSIS_DIR / "tables"
SOURCE_DIR = ANALYSIS_DIR / "source_tables"
CONFIG_PATH = (
    PROJECT_ROOT
    / "experiments/2026-07-23_EXP-061_operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/config.yaml"
)
REFERENCE_CONFIG_PATH = (
    PROJECT_ROOT
    / "experiments/2026-07-22_EXP-060_extensive-operator-normalized-gnk-latent-scale-law-validation/config.yaml"
)
DVC_POINTER = (
    PROJECT_ROOT
    / "data/processed/operator-normalized-gnk-sigmoid-transmission-paired-residual-validation/RUN-151.dvc"
)
DVC_DIR_HASH = "5aace1d2b1f33df736417bf314cfe5bb.dir"
DVC_MANIFEST = PROJECT_ROOT / ".dvc/cache/files/md5/5a/ace1d2b1f33df736417bf314cfe5bb.dir"
REMOTE_RUN_DIR = "/tank/data/labproj/runtime/project-adapt-env/runs/RUN-151"
REMOTE_HOST = "matthew-spence@raspi-1"

EXPECTED_OUTPUTS = [
    "preflight.json",
    "design_catalog.csv",
    "occupancy_calibration.csv",
    "paired_residual_metrics.csv",
    "mutation_transmission_metrics.csv",
    "epistasis_transmission_metrics.csv",
    "topology_audit.csv",
    "scale_law_fits.csv",
    "acceptance_checks.csv",
    "summary.json",
]


def file_hashes(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def remote_manifest() -> dict[str, Any]:
    manifest_path = SOURCE_DIR / "run-151_remote_sha256.json"
    if manifest_path.is_file():
        return json.loads(manifest_path.read_text(encoding="utf-8"))

    command = (
        f"cd {REMOTE_RUN_DIR} && "
        "find logs outputs -type f -print0 | sort -z | xargs -0 sha256sum"
    )
    result = subprocess.run(
        ["ssh", "-o", "BatchMode=yes", REMOTE_HOST, command],
        check=True,
        capture_output=True,
        text=True,
    )
    files: dict[str, str] = {}
    for line in result.stdout.splitlines():
        digest, relative_path = line.split(maxsplit=1)
        files[relative_path.lstrip(" *")] = digest
    payload = {
        "run_id": "RUN-151",
        "remote_host": REMOTE_HOST,
        "remote_run_dir": REMOTE_RUN_DIR,
        "algorithm": "sha256",
        "file_count": len(files),
        "files": dict(sorted(files.items())),
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def bootstrap_median_lower(values: np.ndarray, *, replicates: int, seed: int) -> float:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    boot = np.empty(replicates, dtype=np.float64)
    for index in range(replicates):
        boot[index] = float(np.median(rng.choice(values, size=values.size, replace=True)))
    return float(np.quantile(boot, 0.025))


def log_length_slope(frame: pd.DataFrame, metric: str) -> float:
    medians = frame.groupby("length", sort=True)[metric].median()
    medians = medians[np.isfinite(medians) & (medians > 0.0)]
    return float(
        np.polyfit(
            np.log(medians.index.to_numpy(dtype=float)),
            np.log(medians.to_numpy(dtype=float)),
            1,
        )[0]
    )


def json_value(checks: pd.DataFrame, check_id: str) -> Any:
    value = checks.loc[checks["check_id"] == check_id, "value"].iloc[0]
    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return value


TABLE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

summary = json.loads((OUTPUT_DIR / "summary.json").read_text(encoding="utf-8"))
preflight = json.loads((OUTPUT_DIR / "preflight.json").read_text(encoding="utf-8"))
config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
checks = pd.read_csv(OUTPUT_DIR / "acceptance_checks.csv")
fits = pd.read_csv(OUTPUT_DIR / "scale_law_fits.csv")
design = pd.read_csv(OUTPUT_DIR / "design_catalog.csv")
occupancy = pd.read_csv(OUTPUT_DIR / "occupancy_calibration.csv")
paired = pd.read_csv(OUTPUT_DIR / "paired_residual_metrics.csv")
mutation = pd.read_csv(OUTPUT_DIR / "mutation_transmission_metrics.csv")
epistasis = pd.read_csv(OUTPUT_DIR / "epistasis_transmission_metrics.csv")
topology = pd.read_csv(OUTPUT_DIR / "topology_audit.csv")

expected_lengths = [int(value) for value in config["references"].get("lengths", summary["lengths"])]
if expected_lengths != [48, 93, 118, 220, 238, 448]:
    expected_lengths = [int(value) for value in summary["lengths"]]
expected_seeds = [int(value) for value in config["replication"]["seeds"]]
expected_index = pd.MultiIndex.from_product(
    [expected_lengths, expected_seeds], names=["length", "generator_seed"]
)
gates = config["acceptance"]

# Collection, remote parity, and DVC content-address audit.
remote = remote_manifest()
dvc_entries = json.loads(DVC_MANIFEST.read_text(encoding="utf-8"))
dvc_by_path = {str(row["relpath"]): str(row["md5"]) for row in dvc_entries}
local_files = sorted(path for path in RUN_DIR.rglob("*") if path.is_file())
collection_rows: list[dict[str, Any]] = []
for path in local_files:
    relative = path.relative_to(RUN_DIR).as_posix()
    local_md5, local_sha256 = file_hashes(path)
    collection_rows.append(
        {
            "relative_path": relative,
            "size_bytes": path.stat().st_size,
            "local_md5": local_md5,
            "dvc_md5": dvc_by_path.get(relative),
            "dvc_md5_matches": local_md5 == dvc_by_path.get(relative),
            "local_sha256": local_sha256,
            "remote_sha256": remote["files"].get(relative),
            "remote_sha256_matches": local_sha256 == remote["files"].get(relative),
        }
    )
collection = pd.DataFrame(collection_rows)
collection.to_csv(TABLE_DIR / "run-151_collection_audit.csv", index=False)

assert len(collection) == len(dvc_by_path) == int(remote["file_count"]) == 193
assert set(collection["relative_path"]) == set(dvc_by_path) == set(remote["files"])
assert collection["dvc_md5_matches"].all()
assert collection["remote_sha256_matches"].all()
assert int(collection["size_bytes"].sum()) == 1_513_486_733
assert DVC_POINTER.is_file() and DVC_MANIFEST.is_file()
assert f"md5: {DVC_DIR_HASH}" in DVC_POINTER.read_text(encoding="utf-8")
assert (RUN_DIR / "logs/labproj-exit-code.txt").read_text(encoding="utf-8").strip() == "0"
assert (RUN_DIR / "logs/stderr.txt").stat().st_size == 0
assert len(list((OUTPUT_DIR / "checkpoints").glob("*.json"))) == 180
assert all((OUTPUT_DIR / name).is_file() for name in EXPECTED_OUTPUTS)

for name in ("summary.json", "preflight.json"):
    shutil.copy2(OUTPUT_DIR / name, SOURCE_DIR / f"run-151_{name}")

# Runtime source/config hashes are checked against the present registered inputs.
source_hash_rows = [
    {
        "relative_path": str(CONFIG_PATH.relative_to(PROJECT_ROOT)),
        "recorded_sha256": preflight["config_file_sha256"],
        "current_sha256": file_sha256(CONFIG_PATH),
        "hash_role": "experiment_config",
    },
    {
        "relative_path": str(REFERENCE_CONFIG_PATH.relative_to(PROJECT_ROOT)),
        "recorded_sha256": preflight["reference_source_sha256"],
        "current_sha256": file_sha256(REFERENCE_CONFIG_PATH),
        "hash_role": "reference_source_config",
    },
]
for relative, recorded in preflight["implementation_hash_components"].items():
    source_hash_rows.append(
        {
            "relative_path": relative,
            "recorded_sha256": recorded,
            "current_sha256": file_sha256(PROJECT_ROOT / relative),
            "hash_role": "implementation_component",
        }
    )
source_hash_audit = pd.DataFrame(source_hash_rows)
source_hash_audit["matches"] = (
    source_hash_audit["recorded_sha256"] == source_hash_audit["current_sha256"]
)
source_hash_audit.to_csv(TABLE_DIR / "run-151_source_hash_audit.csv", index=False)
assert source_hash_audit["matches"].all()

# Basic execution and row-count audit.
assert summary["completed_task_count"] == summary["task_count"] == 180
assert summary["reference_count"] == 6 and summary["generator_seed_count"] == 30
assert summary["design_row_count"] == len(design) == 3_240
assert summary["paired_residual_row_count"] == len(paired) == 46_080
assert summary["mutation_transmission_row_count"] == len(mutation) == 171_520
assert summary["epistasis_transmission_row_count"] == len(epistasis) == 146_880
assert summary["topology_audit_row_count"] == len(topology) == 112_320
assert design["task_signature"].nunique() == 180
assert design["sampling_key_excludes_treatment"].dropna().astype(bool).all()
assert design["design_replay_match"].dropna().astype(bool).all()
assert design["candidate_replay_match"].dropna().astype(bool).all()

# Pointwise calibration is a technical transfer-function assay, not one landscape.
point_rows: list[dict[str, Any]] = []
for operating_point, group in paired[
    (paired["arm"] == "extensive_balanced_bounded")
    & (paired["readout_scope"] == "pointwise_occupancy")
    & (paired["sample_id"] == "reference_shell_h1_complete")
    & (paired["weight_label"] == "a_0.5")
].groupby("operating_point", sort=True):
    for length, length_rows in group.groupby("length", sort=True):
        point_rows.append(
            {
                "operating_point": operating_point,
                "occupancy": float(length_rows["occupancy"].iloc[0]),
                "length": int(length),
                "n_seeds": int(length_rows["generator_seed"].nunique()),
                "median_gamma": float(length_rows["gamma_median"].median()),
                "median_practical_resolution_fraction": float(
                    length_rows["stable_practical_resolution_fraction"].median()
                ),
                "median_stable_y_rms": float(length_rows["stable_y_rms"].median()),
                "median_linear_correlation": float(length_rows["linear_correlation"].median()),
                "median_linear_nrmse": float(length_rows["linear_nrmse"].median()),
            }
        )
pointwise = pd.DataFrame(point_rows)
pointwise.to_csv(TABLE_DIR / "run-151_pointwise_transfer_by_length.csv", index=False)

point_primary = paired[
    (paired["arm"] == "extensive_balanced_bounded")
    & (paired["readout_scope"] == "pointwise_occupancy")
    & (paired["operating_point"] == "p_point_0.5")
    & (paired["sample_id"] == "reference_shell_h1_complete")
    & (paired["weight_label"] == "a_0.5")
]
point_slope = log_length_slope(point_primary, "stable_y_rms")
point_fit = fits.set_index("fit_id").loc["pointwise_p0.5_a0.5_h1"]
assert np.isclose(point_slope, float(point_fit["slope"]))
point_gamma_by_length = point_primary.groupby("length")["gamma_median"].median()
point_resolution_by_length = point_primary.groupby("length")[
    "stable_practical_resolution_fraction"
].median()
point_linear_correlation = float(point_primary["linear_correlation"].median())
point_linear_nrmse = float(point_primary["linear_nrmse"].median())
shoulder_primary = paired[
    (paired["arm"] == "extensive_balanced_bounded")
    & (paired["readout_scope"] == "pointwise_occupancy")
    & (paired["sample_id"] == "reference_shell_h1_complete")
    & (paired["operating_point"].isin(["p_point_0.1", "p_point_0.9"]))
    & (paired["weight_label"] == "a_0.5")
]
shoulder_gamma = float(
    shoulder_primary.groupby("operating_point")["gamma_median"].median().min()
)
pointwise_pass = bool(
    (point_gamma_by_length >= gates["pointwise_gamma_median_min"]).all()
    and (point_resolution_by_length >= gates["pointwise_detectable_fraction_min"]).all()
    and shoulder_gamma >= gates["shoulder_gamma_median_min"]
    and point_linear_correlation >= gates["linearization_correlation_min"]
    and point_linear_nrmse <= gates["linearization_nrmse_max"]
    and abs(float(point_fit["slope"])) <= gates["local_slope_abs_max"]
    and float(point_fit["slope_ci_low"]) >= -gates["local_slope_ci_abs_max"]
    and float(point_fit["slope_ci_high"]) <= gates["local_slope_ci_abs_max"]
)

# Common-threshold reference-anchored edges.
controlled_all = mutation[
    (mutation["arm"] == "extensive_balanced_bounded")
    & (mutation["readout_scope"] == "reference_anchored")
    & (mutation["operating_point"] == "p_ref_0.5")
    & (mutation["sample_id"] == "reference_shell_h1_complete")
    & (mutation["weight_label"] == "a_0.5")
    & (mutation["occupancy_stratum"] == "both_dynamic")
].copy()
controlled_eligible = controlled_all[
    controlled_all["n"] >= int(gates["minimum_stratum_observations"])
].copy()
all_edges = mutation[
    (mutation["arm"] == "extensive_balanced_bounded")
    & (mutation["readout_scope"] == "reference_anchored")
    & (mutation["operating_point"] == "p_ref_0.5")
    & (mutation["sample_id"] == "reference_shell_h1_complete")
    & (mutation["weight_label"] == "a_0.5")
    & (mutation["occupancy_stratum"] == "all")
].copy()
occupancy_half = occupancy[
    (occupancy["readout_scope"] == "reference_anchored")
    & (occupancy["operating_point"] == "p_ref_0.5")
    & (occupancy["sample_id"] == "reference_shell_h1_complete")
].copy()

edge_rows: list[dict[str, Any]] = []
for length in expected_lengths:
    raw = controlled_all[controlled_all["length"] == length]
    eligible = controlled_eligible[controlled_eligible["length"] == length]
    all_length = all_edges[all_edges["length"] == length]
    occ = occupancy_half[occupancy_half["length"] == length]
    edge_rows.append(
        {
            "length": length,
            "complete_panel_n": 19 * length,
            "n_seeds": int(raw["generator_seed"].nunique()),
            "median_both_dynamic_n": float(raw["n"].median()),
            "median_both_dynamic_fraction": float(occ["both_dynamic_fraction"].median()),
            "median_one_crossing_fraction": float(occ["one_crossing_fraction"].median()),
            "eligible_seed_count_n_ge_32": int(eligible["generator_seed"].nunique()),
            "eligible_seed_fraction_n_ge_32": float(
                eligible["generator_seed"].nunique() / len(expected_seeds)
            ),
            "eligible_median_stable_correlation": (
                float(eligible["stable_correlation"].median()) if len(eligible) else np.nan
            ),
            "eligible_median_practical_resolution_fraction": (
                float(eligible["stable_practical_resolution_fraction"].median())
                if len(eligible)
                else np.nan
            ),
            "eligible_median_stable_y_rms": (
                float(eligible["stable_y_rms"].median()) if len(eligible) else np.nan
            ),
            "all_edge_median_stable_correlation": float(all_length["stable_correlation"].median()),
            "all_edge_median_practical_resolution_fraction": float(
                all_length["stable_practical_resolution_fraction"].median()
            ),
            "all_edge_median_exact_zero_fraction": float(
                all_length["stable_exact_zero_fraction"].median()
            ),
            "median_across_seed_source_core_occupancy_mean": float(
                all_length["core_occupancy_source_mean"].median()
            ),
            "median_across_seed_target_core_occupancy_mean": float(
                all_length["core_occupancy_target_mean"].median()
            ),
            "median_panel_source_core_occupancy": float(
                occ["source_occupancy_median"].median()
            ),
            "median_panel_target_core_occupancy": float(
                occ["target_occupancy_median"].median()
            ),
        }
    )
edge_table = pd.DataFrame(edge_rows)
edge_table.to_csv(TABLE_DIR / "run-151_reference_anchored_edge_by_length.csv", index=False)

edge_fit = fits.set_index("fit_id").loc["reference_anchored_p0.5_a0.5_h1_both_dynamic"]
assert np.isclose(log_length_slope(controlled_eligible, "stable_y_rms"), float(edge_fit["slope"]))
controlled_edge_pass = bool(
    (edge_table["eligible_seed_fraction_n_ge_32"] >= gates["minimum_stratum_seed_coverage_fraction_per_length"]).all()
    and edge_table["eligible_median_stable_correlation"].dropna().min()
    >= gates["transmission_correlation_min"]
    and edge_table["eligible_median_practical_resolution_fraction"].dropna().min()
    >= gates["texture_detectable_fraction_min"]
    and abs(float(edge_fit["slope"])) <= gates["local_slope_abs_max"]
    and float(edge_fit["slope_ci_low"]) >= -gates["local_slope_ci_abs_max"]
    and float(edge_fit["slope_ci_high"]) <= gates["local_slope_ci_abs_max"]
)


def texture_gate_summary(
    core_mode: str,
    *,
    deployment: bool,
) -> tuple[bool, list[dict[str, Any]]]:
    frame = epistasis[
        (epistasis["arm"] == "extensive_balanced_bounded")
        & (epistasis["support_order"] > 0)
        & (epistasis["square_background"] == "reference")
        & (epistasis["additive_core_mode"] == core_mode)
    ].copy()
    if deployment:
        projection = frame[
            (frame["operating_point"] == "public_base_margin8")
            & (frame["weight_label"] == "fixed_0.6")
        ].copy()
        practical = projection
        prefix = "direct_tau_"
        bootstrap_offset = 11
        fit = None
    else:
        projection = frame[
            (frame["operating_point"] == "p_ref_0.5")
            & (frame["weight_label"] == "a_0.5")
        ].copy()
        practical = frame[
            (frame["operating_point"] == "p_ref_0.5")
            & (frame["weight_label"] == "a_1")
        ].copy()
        prefix = "stable_tau_"
        bootstrap_offset = 8 if core_mode == "frozen_additive_core" else 7
        fit_id = (
            "texture_frozen_p0.5_a0.5"
            if core_mode == "frozen_additive_core"
            else "texture_actual_p0.5_a0.5"
        )
        fit = fits.set_index("fit_id").loc[fit_id]

    seed_projection = (
        projection.groupby(["length", "generator_seed"])[f"{prefix}projection"]
        .median()
        .reindex(expected_index)
    )
    rows: list[dict[str, Any]] = []
    for length_index, length in enumerate(expected_lengths):
        projected_length = projection[projection["length"] == length]
        practical_length = practical[practical["length"] == length]
        seed_values = seed_projection.xs(length, level="length")
        pooled_by_seed = projected_length.groupby("generator_seed").apply(
            lambda group: float(group[f"{prefix}sum_xy"].sum() / group[f"{prefix}sum_x2"].sum()),
            include_groups=False,
        )
        predictor_energy = (
            projected_length[f"{prefix}sum_x2"] / projected_length[f"{prefix}n"].clip(lower=1)
        )
        permutation = projected_length[
            projected_length["permutation_required_primary"].eq(True)
        ]
        rows.append(
            {
                "core_mode": core_mode,
                "scope": "margin8_deployment" if deployment else "reference_anchored_p0.5",
                "length": length,
                "projection_weight_label": "fixed_0.6" if deployment else "a_0.5",
                "practical_weight_label": "fixed_0.6" if deployment else "a_1",
                "n_projection_rows": int(len(projected_length)),
                "seed_coverage_fraction": float(seed_values.notna().mean()),
                "positive_seed_projection_fraction_order_median": float((seed_values > 0.0).mean()),
                "positive_seed_projection_fraction_pooled": float((pooled_by_seed > 0.0).mean()),
                "median_projection_order_median": float(seed_values.median()),
                "median_projection_pooled": float(pooled_by_seed.median()),
                "projection_bootstrap_lower": bootstrap_median_lower(
                    seed_values.dropna().to_numpy(dtype=float),
                    replicates=int(gates["bootstrap_replicates"]),
                    seed=int(gates["bootstrap_seed"]) + bootstrap_offset + 101 * length_index,
                ),
                "median_practical_resolution_fraction": float(
                    practical_length[f"{prefix}practical_resolution_fraction"].median()
                ),
                "pooled_practical_resolution_fraction": float(
                    np.average(
                        practical_length[f"{prefix}practical_resolution_fraction"],
                        weights=practical_length[f"{prefix}n"],
                    )
                ),
                "permutation_pass_fraction": float(permutation["permutation_pass"].eq(True).mean()),
                "median_linear_correlation": (
                    np.nan if deployment else float(projected_length["linear_correlation"].median())
                ),
                "minimum_predictor_energy": float(predictor_energy.min()),
                "minimum_n": int(projected_length["n"].min()),
                "median_tau_y_rms": float(projected_length[f"{prefix}y_rms"].median()),
                "median_tau_exact_zero_fraction": float(
                    projected_length[f"{prefix}exact_zero_fraction"].median()
                ),
                "median_epsilon_R_rms": float(
                    projected_length[
                        "epsilon_R_direct_rms" if deployment else "epsilon_R_stable_rms"
                    ].median()
                ),
                "length_slope": np.nan if fit is None else float(fit["slope"]),
                "length_slope_ci_low": np.nan if fit is None else float(fit["slope_ci_low"]),
                "length_slope_ci_high": np.nan if fit is None else float(fit["slope_ci_high"]),
            }
        )

    table = pd.DataFrame(rows)
    projection_threshold = gates[
        "deployment_projection_positive_seed_fraction_min"
        if deployment
        else "transmission_projection_positive_seed_fraction_min"
    ]
    scale_pass = True if fit is None else bool(
        abs(float(fit["slope"])) <= gates["local_slope_abs_max"]
        and float(fit["slope_ci_low"]) >= -gates["local_slope_ci_abs_max"]
        and float(fit["slope_ci_high"]) <= gates["local_slope_ci_abs_max"]
    )
    passed = bool(
        (table["seed_coverage_fraction"] >= gates["minimum_stratum_seed_coverage_fraction_per_length"]).all()
        and (table["positive_seed_projection_fraction_order_median"] >= projection_threshold).all()
        and (table["projection_bootstrap_lower"] > gates["texture_projection_bootstrap_lower_min"]).all()
        and (table["median_practical_resolution_fraction"] >= gates["texture_detectable_fraction_min"]).all()
        and (table["permutation_pass_fraction"] >= gates["texture_permutation_pass_fraction_min"]).all()
        and (table["minimum_predictor_energy"] >= gates["minimum_texture_predictor_variance"]).all()
        and (table["minimum_n"] >= gates["minimum_nonzero_texture_count"]).all()
        and (
            deployment
            or (table["median_linear_correlation"] >= gates["transmission_correlation_min"]).all()
        )
        and scale_pass
    )
    return passed, rows


frozen_texture_pass, frozen_rows = texture_gate_summary(
    "frozen_additive_core", deployment=False
)
actual_texture_pass, actual_rows = texture_gate_summary(
    "actual_additive_core", deployment=False
)
margin8_texture_pass, margin8_texture_rows = texture_gate_summary(
    "actual_additive_core", deployment=True
)
texture_table = pd.DataFrame(frozen_rows + actual_rows + margin8_texture_rows)
texture_table.to_csv(TABLE_DIR / "run-151_texture_transmission_by_length.csv", index=False)

for core_mode, fit_id in (
    ("frozen_additive_core", "texture_frozen_p0.5_a0.5"),
    ("actual_additive_core", "texture_actual_p0.5_a0.5"),
):
    selected = epistasis[
        (epistasis["arm"] == "extensive_balanced_bounded")
        & (epistasis["support_order"] > 0)
        & (epistasis["square_background"] == "reference")
        & (epistasis["additive_core_mode"] == core_mode)
        & (epistasis["operating_point"] == "p_ref_0.5")
        & (epistasis["weight_label"] == "a_0.5")
    ]
    assert np.isclose(
        log_length_slope(selected, "stable_tau_y_rms"),
        float(fits.set_index("fit_id").loc[fit_id, "slope"]),
    )

# Margin 8 primary deployment and mandatory margin 10 comparator.
deployment_rows: list[dict[str, Any]] = []
for operating_point, margin in (
    ("public_base_margin8", 8.0),
    ("constructor_default_margin10", 10.0),
):
    edge = mutation[
        (mutation["arm"] == "extensive_balanced_bounded")
        & (mutation["operating_point"] == operating_point)
        & (mutation["weight_label"] == "fixed_0.6")
        & (mutation["sample_id"] == "reference_shell_h1_complete")
        & (mutation["occupancy_stratum"] == "all")
    ]
    square = epistasis[
        (epistasis["arm"] == "extensive_balanced_bounded")
        & (epistasis["support_order"] > 0)
        & (epistasis["square_background"] == "reference")
        & (epistasis["additive_core_mode"] == "actual_additive_core")
        & (epistasis["operating_point"] == operating_point)
        & (epistasis["weight_label"] == "fixed_0.6")
    ]
    for length in expected_lengths:
        edge_length = edge[edge["length"] == length]
        square_length = square[square["length"] == length]
        permutation = square_length[square_length["permutation_required_primary"].eq(True)]
        seed_projection = square_length.groupby("generator_seed")["direct_tau_projection"].median()
        deployment_rows.append(
            {
                "operating_point": operating_point,
                "deployment_margin": margin,
                "length": length,
                "n_seeds": int(edge_length["generator_seed"].nunique()),
                "edge_direct_practical_resolution_fraction": float(
                    edge_length["direct_practical_resolution_fraction"].median()
                ),
                "edge_stable_practical_resolution_fraction": float(
                    edge_length["stable_practical_resolution_fraction"].median()
                ),
                "edge_direct_exact_zero_fraction": float(
                    edge_length["direct_exact_zero_fraction"].median()
                ),
                "edge_stable_exact_zero_fraction": float(
                    edge_length["stable_exact_zero_fraction"].median()
                ),
                "edge_direct_y_rms": float(edge_length["direct_y_rms"].median()),
                "edge_stable_y_rms": float(edge_length["stable_y_rms"].median()),
                "texture_positive_seed_projection_fraction": float((seed_projection > 0.0).mean()),
                "texture_direct_practical_resolution_fraction": float(
                    square_length["direct_tau_practical_resolution_fraction"].median()
                ),
                "texture_direct_exact_zero_fraction": float(
                    square_length["direct_tau_exact_zero_fraction"].median()
                ),
                "texture_direct_tau_y_rms": float(square_length["direct_tau_y_rms"].median()),
                "texture_permutation_pass_fraction": float(
                    permutation["permutation_pass"].eq(True).mean()
                ),
            }
        )
deployment = pd.DataFrame(deployment_rows)
deployment.to_csv(
    TABLE_DIR / "run-151_deployment_margin_comparison_by_length.csv", index=False
)

margin8 = deployment[deployment["deployment_margin"] == 8.0]
margin8_edge_pass = bool(
    (margin8["edge_direct_practical_resolution_fraction"] >= gates["deployment_detectable_fraction_min"]).all()
)
margin8_public_pass = bool(margin8_edge_pass and margin8_texture_pass)

# Reproduce the registered decision and compare every pass state.
recorded_pass = checks.set_index("check_id")["passed"].astype(bool).to_dict()
parity_value = json_value(checks, "public_float64_formula_parity")
identity_value = float(json_value(checks, "identity_readout_algebra"))
finite_table_outputs = all(
    not np.isinf(frame.select_dtypes(include=[np.number]).to_numpy()).any()
    for frame in (design, occupancy, paired, mutation, epistasis, topology, fits)
)
reference_occupancy_error = float(
    np.max(
        np.abs(
            occupancy.loc[
                (occupancy["readout_scope"] == "reference_anchored")
                & (occupancy["sample_id"] == "reference_shell_h1_complete"),
                "source_occupancy_median",
            ]
            - occupancy.loc[
                (occupancy["readout_scope"] == "reference_anchored")
                & (occupancy["sample_id"] == "reference_shell_h1_complete"),
                "occupancy",
            ]
        )
    )
)
design_rows = design[design["row_type"] != "field"]
zero_max = float(
    max(
        paired.loc[np.isclose(paired["weight"], 0.0), ["direct_y_rms", "stable_y_rms"]]
        .abs()
        .to_numpy()
        .max(),
        mutation.loc[np.isclose(mutation["weight"], 0.0), ["direct_y_rms", "stable_y_rms"]]
        .abs()
        .to_numpy()
        .max(),
        epistasis.loc[
            np.isclose(epistasis["weight"], 0.0),
            [
                "direct_tau_y_rms",
                "stable_tau_y_rms",
                "epsilon_R_direct_rms",
                "epsilon_R_stable_rms",
            ],
        ]
        .abs()
        .to_numpy()
        .max(),
    )
)
surrogate_all = epistasis[epistasis["arm"] == "additive_surrogate"]
recomputed_pass = {
    "deterministic_field_replay": bool(
        design.loc[design["arm"] == "extensive_balanced_bounded", "candidate_replay_match"]
        .dropna()
        .astype(bool)
        .all()
    ),
    "reference_centering": bool(
        design.loc[
            (design["row_type"] == "field")
            & design["arm"].isin(
                [
                    "extensive_balanced_bounded",
                    "legacy_support_local_normalized",
                    "additive_surrogate",
                ]
            ),
            "reference_score",
        ].abs().max()
        <= gates["reference_center_atol"]
    ),
    "public_float64_formula_parity": bool(
        float(parity_value["core"]) <= gates["public_formula_parity_atol"]
        and float(parity_value["end_to_end_256_term"]) <= gates["public_formula_parity_atol"]
        and design["public_formula_parity_max_abs"].dropna().abs().max()
        <= gates["public_formula_parity_atol"]
    ),
    "finite_outputs": bool(finite_table_outputs and summary["completed_task_count"] == 180),
    "identity_readout_algebra": bool(identity_value <= gates["identity_atol"]),
    "reference_occupancy_calibration": bool(
        reference_occupancy_error <= gates["occupancy_atol"]
    ),
    "constructor_margin10_comparator_present": bool(
        occupancy.loc[
            occupancy["operating_point"] == "constructor_default_margin10",
            ["reference_id", "generator_seed"],
        ]
        .drop_duplicates()
        .shape[0]
        == 180
    ),
    "treatment_independent_sampling": bool(
        design_rows["sampling_key_excludes_treatment"].astype(bool).all()
        and design_rows["design_replay_match"].astype(bool).all()
        and len(design_rows) == 2_520
    ),
    "weight_zero_exact": bool(zero_max <= gates["zero_weight_atol"]),
    "additive_surrogate_latent_epistasis_zero": bool(
        surrogate_all["latent_epsilon_max_abs"].max()
        <= gates["additive_surrogate_latent_epistasis_atol"]
        and int((surrogate_all["epsilon_R_stable_rms"] > 0.0).sum()) == 12_960
    ),
    "strict_monotone_non_tie_topology": bool(
        topology["non_tie_rank_agreement"].min() >= gates["topology_non_tie_agreement_min"]
        and int(topology["non_tie_mismatch_count"].sum()) == 0
    ),
}
recomputed_pass["pointwise_transfer_calibration"] = pointwise_pass
recomputed_pass["reference_anchored_edge_transmission"] = controlled_edge_pass
recomputed_pass["frozen_core_texture_transmission"] = frozen_texture_pass
recomputed_pass["actual_core_texture_transmission"] = actual_texture_pass
recomputed_pass["margin8_direct_float64_deployment"] = margin8_public_pass

checks["recomputed_pass"] = checks["check_id"].map(recomputed_pass)
checks["pass_state_matches"] = checks["passed"].astype(bool) == checks["recomputed_pass"].astype(bool)
checks["audit_method"] = np.where(
    checks["check_id"].isin(
        [
            "pointwise_transfer_calibration",
            "reference_anchored_edge_transmission",
            "frozen_core_texture_transmission",
            "actual_core_texture_transmission",
            "margin8_direct_float64_deployment",
        ]
    ),
    "independent_reaggregation_from_sufficient_statistics",
    "recorded_task_or_aggregate_diagnostic_rechecked_for_consistency",
)
checks.to_csv(TABLE_DIR / "run-151_acceptance_audit.csv", index=False)
assert checks["pass_state_matches"].all()
assert int(checks["passed"].astype(bool).sum()) == 13
assert len(checks) == 16 and checks["required"].astype(bool).all()

technical_pass = bool(checks.loc[checks["family"] == "technical", "passed"].astype(bool).all())
decision_flags = {
    "technical": technical_pass,
    "controlled_edges": controlled_edge_pass,
    "frozen_texture": frozen_texture_pass,
    "actual_texture": actual_texture_pass,
    "public": margin8_public_pass,
}
if not technical_pass:
    decision = config["decision_rule"]["technical_fail"]
elif not controlled_edge_pass:
    decision = config["decision_rule"]["transfer_fail"]
elif frozen_texture_pass and not actual_texture_pass:
    decision = config["decision_rule"]["frozen_core_only"]
elif not actual_texture_pass:
    decision = config["decision_rule"]["texture_fail"]
elif margin8_public_pass:
    decision = config["decision_rule"]["public_pass"]
else:
    decision = config["decision_rule"]["transition_only"]
assert decision_flags == summary["decision_flags"]
assert decision == summary["decision"] == "SIGMOID_TRANSFER_FAILED"

# Controls and documented protocol/output caveats.
surrogate = epistasis[epistasis["arm"] == "additive_surrogate"]
surrogate_focus = surrogate[
    (surrogate["additive_core_mode"] == "actual_additive_core")
    & (surrogate["square_background"] == "reference")
    & (surrogate["operating_point"] == "p_ref_0.5")
    & (surrogate["weight_label"] == "a_0.5")
]
shoulder = pointwise[pointwise["occupancy"].isin([0.1, 0.9])]
protocol_audit = pd.DataFrame(
    [
        {
            "item": "edge_minimum_n_length_interaction",
            "classification": "registered_gate_caveat",
            "evidence": "Both-dynamic fractions remain about 1.4-1.5%, but complete panels scale as 19L; eligible seeds are 0,4,17,30,30,30 by length.",
            "decision_impact": "none; registered controlled-edge coverage gate still fails",
        },
        {
            "item": "square_endpoint_core_occupancies",
            "classification": "output_omission",
            "evidence": "epistasis_transmission_metrics.csv does not retain all square endpoint occupancies or p_ab_null promised by the prose.",
            "decision_impact": "none; actual-versus-frozen texture and edge occupancy still establish the failure, but double-endpoint saturation is inferred rather than directly tabulated",
        },
        {
            "item": "texture_projection_aggregation",
            "classification": "implementation_detail",
            "evidence": "The gate uses each seed's median of three support-order projections; pooled sum_xy/sum_x2 projections have the same positive-fraction pass states.",
            "decision_impact": "none",
        },
        {
            "item": "shoulder_linearization_scope",
            "classification": "prose_ambiguity",
            "evidence": f"The NRMSE <=0.25 gate is applied only at p=0.5; median shoulder NRMSE spans {shoulder['median_linear_nrmse'].min():.6f}-{shoulder['median_linear_nrmse'].max():.6f}.",
            "decision_impact": "none; shoulder acceptance is based on the separately registered gamma >=0.20 bound",
        },
        {
            "item": "finite_outputs_label_scope",
            "classification": "implementation_label_caveat",
            "evidence": "The task check covers raw arrays and selected RMS columns, not literally every CSV cell; structural NaNs remain in nonapplicable fields.",
            "decision_impact": "none; no consequential nonfinite scientific value was found",
        },
        {
            "item": "pointwise_landscape_scope",
            "classification": "scope_control",
            "evidence": "theta is recalibrated per genotype in pointwise rows, so the passing assay is scalar transfer calibration rather than one fitness landscape.",
            "decision_impact": "prevents a topology or public-method claim from the pointwise pass",
        },
        {
            "item": "additive_surrogate_global_epistasis",
            "classification": "required_control",
            "evidence": f"Maximum latent epsilon is {surrogate['latent_epsilon_max_abs'].max():.6g}, yet observable epsilon is nonzero in {(surrogate['epsilon_R_stable_rms'] > 0.0).sum()} rows; focused median epsilon_R RMS is {surrogate_focus['epsilon_R_stable_rms'].median():.6g}.",
            "decision_impact": "supports tau_texture rather than raw epsilon_R as the GNK-texture estimand",
        },
    ]
)
protocol_audit.to_csv(TABLE_DIR / "run-151_protocol_output_audit.csv", index=False)

fits.to_csv(TABLE_DIR / "run-151_scale_law_summary.csv", index=False)

controlled_check = json_value(checks, "reference_anchored_edge_transmission")
frozen_check = json_value(checks, "frozen_core_texture_transmission")
actual_check = json_value(checks, "actual_core_texture_transmission")
deployment_check = json_value(checks, "margin8_direct_float64_deployment")

metrics = {
    "analysis_id": "ANA-059",
    "experiment_id": "EXP-061",
    "run_id": "RUN-151",
    "decision": decision,
    "scientific_pass": False,
    "decision_flags": decision_flags,
    "required_checks": {"passed": 13, "total": 16},
    "technical_checks": {"passed": 12, "total": 12},
    "collection": {
        "file_count": 193,
        "data_bytes": 1_513_486_733,
        "remote_sha256_matches": 193,
        "dvc_md5_matches": 193,
        "dvc_dir_hash": DVC_DIR_HASH,
        "dvc_remote_backup_method": "manual content-addressed object transfer because the configured client mount was unavailable",
        "dvc_remote_object_bytes_including_dir_manifest": 1_513_511_029,
        "dvc_remote_checksum_dry_run_remaining_files": 0,
    },
    "execution": {
        "scheduler_job_id": "2041",
        "scheduler_exit_code": 0,
        "task_count": 180,
        "reference_count": 6,
        "generator_seed_count": 30,
        "lengths": expected_lengths,
        "stderr_empty": True,
        "source_hashes_match_current_files": bool(source_hash_audit["matches"].all()),
    },
    "pointwise": {
        "minimum_central_gamma": float(
            pointwise.loc[pointwise["occupancy"] == 0.5, "median_gamma"].min()
        ),
        "minimum_central_practical_fraction": float(
            pointwise.loc[
                pointwise["occupancy"] == 0.5,
                "median_practical_resolution_fraction",
            ].min()
        ),
        "minimum_shoulder_gamma": shoulder_gamma,
        "minimum_per_length_shoulder_gamma": float(
            pointwise.loc[pointwise["occupancy"].isin([0.1, 0.9]), "median_gamma"].min()
        ),
        "central_linear_correlation": point_linear_correlation,
        "central_linear_nrmse": point_linear_nrmse,
        "length_slope": float(point_fit["slope"]),
        "length_slope_ci": [float(point_fit["slope_ci_low"]), float(point_fit["slope_ci_high"])],
        "pass": pointwise_pass,
        "scientific_landscape_scope": False,
    },
    "reference_anchored_edges": {
        "eligible_seed_count_by_length": {
            str(row.length): int(row.eligible_seed_count_n_ge_32)
            for row in edge_table.itertuples(index=False)
        },
        "minimum_eligible_seed_fraction": float(
            edge_table["eligible_seed_fraction_n_ge_32"].min()
        ),
        "minimum_eligible_correlation": float(controlled_check["minimum_per_length_median_correlation"]),
        "minimum_eligible_practical_fraction": float(
            controlled_check["minimum_per_length_median_practical_fraction"]
        ),
        "all_edge_practical_fraction_range": [
            float(edge_table["all_edge_median_practical_resolution_fraction"].min()),
            float(edge_table["all_edge_median_practical_resolution_fraction"].max()),
        ],
        "length_slope": float(edge_fit["slope"]),
        "length_slope_ci": [float(edge_fit["slope_ci_low"]), float(edge_fit["slope_ci_high"])],
        "pass": controlled_edge_pass,
    },
    "frozen_core_texture": {**frozen_check, "pass": frozen_texture_pass},
    "actual_core_texture": {**actual_check, "pass": actual_texture_pass},
    "margin8_deployment": {**deployment_check, "pass": margin8_public_pass},
    "margin10_comparator": {
        "edge_practical_fraction_range": [
            float(
                deployment.loc[
                    deployment["deployment_margin"] == 10.0,
                    "edge_direct_practical_resolution_fraction",
                ].min()
            ),
            float(
                deployment.loc[
                    deployment["deployment_margin"] == 10.0,
                    "edge_direct_practical_resolution_fraction",
                ].max()
            ),
        ],
        "texture_practical_fraction_range": [
            float(
                deployment.loc[
                    deployment["deployment_margin"] == 10.0,
                    "texture_direct_practical_resolution_fraction",
                ].min()
            ),
            float(
                deployment.loc[
                    deployment["deployment_margin"] == 10.0,
                    "texture_direct_practical_resolution_fraction",
                ].max()
            ),
        ],
    },
    "scope": {
        "candidate_is_experiment_local": True,
        "candidate_replaces_public_generator": True,
        "observation_noise_enabled": False,
        "public_model_changed": False,
        "empirical_realism_tested": False,
        "optimizer_or_accessibility_tested": False,
        "stage_4c_4d_4e_hard_stop_retained": True,
    },
}
(TABLE_DIR / "metrics.json").write_text(
    json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)

print(json.dumps({
    "decision": decision,
    "decision_flags": decision_flags,
    "required_checks_passed": 13,
    "collection_files_verified": int(collection["remote_sha256_matches"].sum()),
}, indent=2, sort_keys=True))
