from __future__ import annotations

import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
RUN_DIR = (
    PROJECT_ROOT
    / "data/processed/dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-validation/RUN-153"
)
EXPERIMENT_DIR = (
    PROJECT_ROOT
    / "experiments/2026-07-23_EXP-062_dimensionless-additive-core-sigmoid-calibration-and-held-out-gnk-dose-response-validation"
)
CONFIG_PATH = EXPERIMENT_DIR / "config.yaml"
CATALOG_PATH = EXPERIMENT_DIR / "reference_catalog.yaml"
DVC_POINTER = RUN_DIR.with_suffix(".dvc")
REMOTE_RUN_DIR = Path(
    "/Volumes/DutchZFS/labproj/runtime/project-adapt-env/runs/RUN-153"
)
TABLE_DIR = ANALYSIS_DIR / "tables"
SOURCE_DIR = ANALYSIS_DIR / "source_tables"

GLOBAL = "combined_global"
EXACT = "exact_per_landscape_beta"
CANDIDATE = "extensive_balanced_bounded"
ACTUAL = "actual_additive_core"
FROZEN = "frozen_additive_core"
EXPECTED_OUTPUTS = [
    "preflight.json",
    "reference_catalog.csv",
    "core_calibration_metrics.csv",
    "calibration_selection.csv",
    "validation_edge_metrics.csv",
    "validation_texture_metrics.csv",
    "dose_response_metrics.csv",
    "scale_law_fits.csv",
    "acceptance_checks.csv",
    "summary.json",
]


def hashes(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def sha256(path: Path) -> str:
    return hashes(path)[1]


def remote_manifest() -> dict[str, Any]:
    persisted = SOURCE_DIR / "run-153_remote_sha256.json"
    if persisted.is_file():
        return json.loads(persisted.read_text(encoding="utf-8"))
    roots = [(REMOTE_RUN_DIR / "outputs", ""), (REMOTE_RUN_DIR / "logs", "logs/")]
    if not all(root.is_dir() for root, _ in roots):
        raise FileNotFoundError(
            "Shared RUN-153 source unavailable and no persisted manifest exists"
        )
    files: dict[str, dict[str, Any]] = {}
    for root, prefix in roots:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                relative = prefix + path.relative_to(root).as_posix()
                files[relative] = {
                    "sha256": sha256(path),
                    "size_bytes": int(path.stat().st_size),
                }
    payload = {
        "run_id": "RUN-153",
        "source": str(REMOTE_RUN_DIR),
        "algorithm": "sha256",
        "file_count": len(files),
        "files": dict(sorted(files.items())),
    }
    persisted.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return payload


def bin_fractions(frame: pd.DataFrame, column: str) -> tuple[float, dict[str, float]]:
    fractions = frame.groupby("length_bin", sort=True)[column].mean()
    return float(fractions.min()), {
        str(key): float(value) for key, value in fractions.items()
    }


def crossed_bootstrap_lower(
    frame: pd.DataFrame, metric: str, replicates: int, seed: int
) -> float:
    matrices: dict[str, np.ndarray] = {}
    for reference_id, group in frame.sort_values(
        ["reference_id", "core_seed", "gnk_seed"], kind="mergesort"
    ).groupby("reference_id", sort=True):
        matrix = group.pivot(
            index="core_seed", columns="gnk_seed", values=metric
        ).to_numpy(dtype=float)
        assert matrix.shape == (5, 6) and np.isfinite(matrix).all()
        matrices[str(reference_id)] = matrix
    references = sorted(matrices)
    rng = np.random.default_rng(seed)
    values = np.empty(replicates, dtype=float)
    for index in range(replicates):
        sampled_refs = rng.choice(references, len(references), replace=True)
        ref_values = []
        for reference_id in sampled_refs:
            matrix = matrices[str(reference_id)]
            core_index = rng.integers(0, 5, 5)
            gnk_index = rng.integers(0, 6, 6)
            ref_values.append(float(np.median(matrix[np.ix_(core_index, gnk_index)])))
        values[index] = float(np.median(ref_values))
    return float(np.quantile(values, 0.025))


def length_fit(
    texture: pd.DataFrame, mode: str, config: dict[str, Any]
) -> dict[str, Any]:
    selected = texture[
        (texture["scale_mode"] == mode)
        & (texture["additive_core_mode"] == ACTUAL)
        & (texture["field_arm"] == CANDIDATE)
        & texture["pooled_supported_reference"].eq(True)
        & np.isclose(texture["amplitude"], 1.0)
    ]
    protein = selected.groupby(
        ["reference_id", "length_bin", "length"], as_index=False
    )["direct_tau_y_rms"].median()
    x = np.log(protein["length"].to_numpy(dtype=float))
    y = np.log(protein["direct_tau_y_rms"].to_numpy(dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    rng = np.random.default_rng(
        int(config["acceptance"]["bootstrap_seed"]) + (400 if mode == GLOBAL else 500)
    )
    bins = sorted(protein["length_bin"].unique())
    boot: list[float] = []
    for _ in range(int(config["acceptance"]["bootstrap_replicates"])):
        pieces = []
        for length_bin in bins:
            block = protein[protein["length_bin"] == length_bin]
            pieces.append(block.iloc[rng.integers(0, len(block), len(block))])
        replicate = pd.concat(pieces)
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
    return {
        "fit_id": f"{mode}_actual_tau_a1",
        "metric": "protein_median_direct_tau_y_rms",
        "n_proteins": int(len(protein)),
        "n_length_bins": len(bins),
        "slope": float(slope),
        "intercept": float(intercept),
        "slope_ci_low": float(low),
        "slope_ci_high": float(high),
        "bootstrap_replicates_realized": len(boot),
        "recomputed_pass": passed,
    }


def clean(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    return value


TABLE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

summary = json.loads((RUN_DIR / "summary.json").read_text(encoding="utf-8"))
preflight = json.loads((RUN_DIR / "preflight.json").read_text(encoding="utf-8"))
config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
catalog = pd.read_csv(RUN_DIR / "reference_catalog.csv")
core = pd.read_csv(RUN_DIR / "core_calibration_metrics.csv")
selection = pd.read_csv(RUN_DIR / "calibration_selection.csv")
edge = pd.read_csv(RUN_DIR / "validation_edge_metrics.csv")
texture = pd.read_csv(RUN_DIR / "validation_texture_metrics.csv")
dose = pd.read_csv(RUN_DIR / "dose_response_metrics.csv")
fits = pd.read_csv(RUN_DIR / "scale_law_fits.csv")
checks = pd.read_csv(RUN_DIR / "acceptance_checks.csv")

# Collection parity and content-addressed DVC audit.
pointer = yaml.safe_load(DVC_POINTER.read_text(encoding="utf-8"))["outs"][0]
dvc_hash = str(pointer["md5"])
dvc_manifest_path = (
    PROJECT_ROOT / ".dvc/cache/files/md5" / dvc_hash[:2] / dvc_hash[2:]
)
dvc_entries = json.loads(dvc_manifest_path.read_text(encoding="utf-8"))
dvc_by_path = {str(row["relpath"]): str(row["md5"]) for row in dvc_entries}
remote = remote_manifest()
collection_rows = []
for path in sorted(item for item in RUN_DIR.rglob("*") if item.is_file()):
    relative = path.relative_to(RUN_DIR).as_posix()
    local_md5, local_sha256 = hashes(path)
    remote_entry = remote["files"].get(relative, {})
    collection_rows.append(
        {
            "relative_path": relative,
            "size_bytes": int(path.stat().st_size),
            "local_md5": local_md5,
            "dvc_md5": dvc_by_path.get(relative),
            "dvc_md5_matches": local_md5 == dvc_by_path.get(relative),
            "local_sha256": local_sha256,
            "remote_sha256": remote_entry.get("sha256"),
            "remote_sha256_matches": local_sha256 == remote_entry.get("sha256"),
            "remote_size_matches": int(path.stat().st_size)
            == remote_entry.get("size_bytes"),
        }
    )
collection = pd.DataFrame(collection_rows)
collection.to_csv(TABLE_DIR / "run-153_collection_audit.csv", index=False)
assert len(collection) == len(dvc_by_path) == remote["file_count"] == 373
assert set(collection["relative_path"]) == set(dvc_by_path) == set(remote["files"])
assert collection["dvc_md5_matches"].all()
assert collection["remote_sha256_matches"].all()
assert collection["remote_size_matches"].all()
assert int(collection["size_bytes"].sum()) == pointer["size"] == 332_322_278
assert pointer["nfiles"] == 373
assert (RUN_DIR / "logs/labproj-exit-code.txt").read_text().strip() == "0"
assert (RUN_DIR / "logs/stderr.txt").stat().st_size == 0
assert all((RUN_DIR / name).is_file() for name in EXPECTED_OUTPUTS)

shutil.copy2(RUN_DIR / "summary.json", SOURCE_DIR / "run-153_summary.json")
shutil.copy2(RUN_DIR / "preflight.json", SOURCE_DIR / "run-153_preflight.json")

# Registered source hashes.
source_rows = [
    {
        "relative_path": str(CONFIG_PATH.relative_to(PROJECT_ROOT)),
        "hash_role": "experiment_config",
        "recorded_sha256": preflight["config_file_sha256"],
        "current_sha256": sha256(CONFIG_PATH),
    },
    {
        "relative_path": str(CATALOG_PATH.relative_to(PROJECT_ROOT)),
        "hash_role": "reference_catalog",
        "recorded_sha256": preflight["reference_catalog_sha256"],
        "current_sha256": sha256(CATALOG_PATH),
    },
]
for relative, recorded in preflight["implementation_hash_components"].items():
    source_rows.append(
        {
            "relative_path": relative,
            "hash_role": "implementation_component",
            "recorded_sha256": recorded,
            "current_sha256": sha256(PROJECT_ROOT / relative),
        }
    )
source_audit = pd.DataFrame(source_rows)
source_audit["matches"] = (
    source_audit["recorded_sha256"] == source_audit["current_sha256"]
)
source_audit.to_csv(TABLE_DIR / "run-153_source_hash_audit.csv", index=False)
assert source_audit["matches"].all()

# Checkpoints: complete crossed design, aggregate counts, and technical controls.
technical_rows = []
checkpoint_edge_count = 0
checkpoint_texture_count = 0
nan_token_count = 0
for path in sorted((RUN_DIR / "checkpoints").glob("*.json")):
    text = path.read_text(encoding="utf-8")
    nan_token_count += text.count("NaN")
    payload = json.loads(text)
    technical_rows.append(payload["technical"])
    checkpoint_edge_count += len(payload["edge_rows"])
    checkpoint_texture_count += len(payload["texture_rows"])
technical = pd.DataFrame(technical_rows)
assert len(technical) == technical["task_signature"].nunique() == 360
assert technical.groupby("reference_id").size().eq(30).all()
assert technical.groupby("length_bin").size().eq(60).all()
assert checkpoint_edge_count == len(edge) == 17_280
assert checkpoint_texture_count == len(texture) == 78_120
assert nan_token_count > 0

# Independent core-only selection and held-out calibration.
scale = core[core["row_type"] == "scale_estimate"].copy()
grid = core[core["row_type"] == "core_grid_cell"].copy()
heldout = core[core["row_type"] == "heldout_core_confirmation"].copy()
assert (len(scale), len(grid), len(heldout)) == (60, 1_080, 120)
global_sigma_a = float(scale["sigma_a_complete_h1"].median())
assert np.isclose(global_sigma_a, summary["global_sigma_a_complete_h1_median"])
minimum_core = float(
    config["acceptance"]["core"]["minimum_cell_pass_fraction_per_length_bin"]
)
grid_rows = []
for keys, group in grid.groupby(
    ["beta", "mu", "reference_occupancy"], sort=True
):
    beta, mu, occupancy = keys
    fractions = group.groupby("length_bin")["core_pass"].mean()
    grid_rows.append(
        {
            "beta": float(beta),
            "mu": float(mu),
            "reference_occupancy": float(occupancy),
            "minimum_length_bin_pass_fraction": float(fractions.min()),
            "eligible": bool(float(fractions.min()) >= minimum_core),
        }
    )
grid_audit = pd.DataFrame(grid_rows)
selected = grid_audit[grid_audit["eligible"]].sort_values(
    ["beta", "mu"], ascending=False
).iloc[0]
assert np.isclose(selected["beta"], 4.0)
assert np.isclose(selected["mu"], math.log(99.0))
assert int(selection["selected"].sum()) == 1
assert not selection["gnk_outcomes_used"].astype(bool).any()

selected_grid = grid[
    np.isclose(grid["beta"], selected["beta"])
    & np.isclose(grid["mu"], selected["mu"])
]
core_calibration_rows = []
for length_bin, group in scale.groupby("length_bin", sort=True):
    chosen = selected_grid[selected_grid["length_bin"] == length_bin]
    core_calibration_rows.append(
        {
            "length_bin": length_bin,
            "scale_cell_count": int(len(group)),
            "sigma_a_complete_h1_mean": float(group["sigma_a_complete_h1"].mean()),
            "sigma_a_complete_h1_median": float(
                group["sigma_a_complete_h1"].median()
            ),
            "sigma_a_complete_h1_min": float(group["sigma_a_complete_h1"].min()),
            "sigma_a_complete_h1_max": float(group["sigma_a_complete_h1"].max()),
            "selected_grid_cell_count": int(len(chosen)),
            "selected_grid_pass_fraction": float(chosen["core_pass"].mean()),
            "reference_fitness_median": float(chosen["reference_fitness"].median()),
            "h1_dynamic_fraction_median": float(
                chosen["h1_target_dynamic_fraction"].median()
            ),
            "h1_core_effect_rms_median": float(
                chosen["h1_core_effect_rms"].median()
            ),
        }
    )
pd.DataFrame(core_calibration_rows).to_csv(
    TABLE_DIR / "run-153_core_calibration_by_length.csv", index=False
)

heldout_rows = []
for keys, group in heldout.groupby(["scale_mode", "length_bin"], sort=True):
    mode, length_bin = keys
    heldout_rows.append(
        {
            "scale_mode": mode,
            "length_bin": length_bin,
            "cell_count": int(len(group)),
            "core_pass_fraction": float(group["core_pass"].mean()),
            "realized_beta_median": float(group["realized_beta"].median()),
            "realized_beta_min": float(group["realized_beta"].min()),
            "realized_beta_max": float(group["realized_beta"].max()),
            "reference_fitness_median": float(group["reference_fitness"].median()),
            "h1_dynamic_fraction_median": float(
                group["h1_target_dynamic_fraction"].median()
            ),
            "h1_core_effect_rms_median": float(group["h1_core_effect_rms"].median()),
            "fraction_0.1_low_plateau_median": float(
                group["fraction_0.1_low_plateau_fraction"].median()
            ),
            "fraction_0.5_low_plateau_median": float(
                group["fraction_0.5_low_plateau_fraction"].median()
            ),
        }
    )
heldout_summary = pd.DataFrame(heldout_rows)
heldout_summary.to_csv(
    TABLE_DIR / "run-153_heldout_core_by_length.csv", index=False
)
global_core_pass = bool(
    heldout_summary[heldout_summary["scale_mode"] == GLOBAL][
        "core_pass_fraction"
    ].min()
    >= minimum_core
)
exact_core_pass = bool(
    heldout_summary[heldout_summary["scale_mode"] == EXACT][
        "core_pass_fraction"
    ].min()
    >= minimum_core
)

# Edge and texture gates reconstructed directly from the output tables.
transmission_gate = config["acceptance"]["transmission"]
transmission_frames: dict[tuple[str, str], tuple[pd.DataFrame, pd.DataFrame]] = {}
transmission_passes: dict[tuple[str, str], bool] = {}
transmission_values: dict[tuple[str, str], dict[str, Any]] = {}
for mode, core_mode in [(GLOBAL, ACTUAL), (GLOBAL, FROZEN), (EXACT, ACTUAL)]:
    edges = edge[
        (edge["scale_mode"] == mode)
        & (edge["field_arm"] == CANDIDATE)
        & (edge["additive_core_mode"] == ACTUAL)
        & (edge["sample_id"] == "reference_shell_h1_complete")
        & np.isclose(edge["amplitude"], 1.0)
    ].copy()
    edges["edge_pass"] = (
        edges["direct_practical_resolution_fraction"]
        >= float(transmission_gate["edge_practical_fraction_min"])
    ) & (
        edges["direct_correlation"] >= float(transmission_gate["correlation_min"])
    )
    cells = texture[
        (texture["scale_mode"] == mode)
        & (texture["additive_core_mode"] == core_mode)
        & (texture["field_arm"] == CANDIDATE)
        & texture["pooled_supported_reference"].eq(True)
        & np.isclose(texture["amplitude"], 1.0)
    ].copy()
    assert len(edges) == len(cells) == 360
    if core_mode == ACTUAL:
        frozen = texture[
            (texture["scale_mode"] == mode)
            & (texture["additive_core_mode"] == FROZEN)
            & (texture["field_arm"] == CANDIDATE)
            & texture["pooled_supported_reference"].eq(True)
            & np.isclose(texture["amplitude"], 1.0)
        ][["reference_id", "core_seed", "gnk_seed", "direct_tau_y_rms"]].rename(
            columns={"direct_tau_y_rms": "frozen_y_rms"}
        )
        cells = cells.merge(
            frozen,
            on=["reference_id", "core_seed", "gnk_seed"],
            validate="one_to_one",
        )
        cells["actual_to_frozen_rms_ratio"] = (
            cells["direct_tau_y_rms"] / cells["frozen_y_rms"]
        )
        ratio_pass = (
            cells["actual_to_frozen_rms_ratio"]
            >= float(transmission_gate["actual_to_frozen_rms_ratio_min"])
        )
    else:
        cells["actual_to_frozen_rms_ratio"] = np.nan
        ratio_pass = pd.Series(True, index=cells.index)
    cells["texture_pass"] = (
        (
            cells["direct_tau_practical_resolution_fraction"]
            >= float(transmission_gate["texture_practical_fraction_min"])
        )
        & (cells["direct_tau_projection"] > 0)
        & (
            cells["direct_tau_correlation"]
            >= float(transmission_gate["correlation_min"])
        )
        & cells["permutation_pass"].eq(True)
        & ratio_pass
    )
    if core_mode == ACTUAL:
        joint = cells.merge(
            edges[
                [
                    "reference_id",
                    "length_bin",
                    "core_seed",
                    "gnk_seed",
                    "edge_pass",
                ]
            ],
            on=["reference_id", "length_bin", "core_seed", "gnk_seed"],
            validate="one_to_one",
        )
        joint["joint_pass"] = joint["texture_pass"] & joint["edge_pass"]
        edge_minimum, edge_fractions = bin_fractions(edges, "edge_pass")
    else:
        joint = cells.copy()
        joint["joint_pass"] = joint["texture_pass"]
        edge_minimum, edge_fractions = 1.0, {"not_required": 1.0}
    texture_minimum, texture_fractions = bin_fractions(cells, "texture_pass")
    joint_minimum, joint_fractions = bin_fractions(joint, "joint_pass")
    bootstrap_lower = crossed_bootstrap_lower(
        cells,
        "direct_tau_projection",
        int(config["acceptance"]["bootstrap_replicates"]),
        int(config["acceptance"]["bootstrap_seed"])
        + (0 if mode == GLOBAL else 100)
        + (0 if core_mode == ACTUAL else 10),
    )
    cell_minimum = float(transmission_gate["minimum_cell_pass_fraction_per_length_bin"])
    passed = bool(
        edge_minimum >= cell_minimum
        and texture_minimum >= cell_minimum
        and joint_minimum >= cell_minimum
        and bootstrap_lower > float(transmission_gate["projection_bootstrap_lower_min"])
    )
    transmission_frames[(mode, core_mode)] = (cells, edges)
    transmission_passes[(mode, core_mode)] = passed
    transmission_values[(mode, core_mode)] = {
        "edge_minimum": edge_minimum,
        "edge_fractions": edge_fractions,
        "texture_minimum": texture_minimum,
        "texture_fractions": texture_fractions,
        "joint_minimum": joint_minimum,
        "joint_fractions": joint_fractions,
        "projection_bootstrap_lower": bootstrap_lower,
    }

transmission_rows = []
for (mode, core_mode), (cells, edges) in transmission_frames.items():
    for length_bin, group in cells.groupby("length_bin", sort=True):
        edge_group = edges[edges["length_bin"] == length_bin]
        transmission_rows.append(
            {
                "scale_mode": mode,
                "additive_core_mode": core_mode,
                "length_bin": length_bin,
                "texture_cell_count": int(len(group)),
                "texture_gate_pass_count": int(group["texture_pass"].sum()),
                "texture_gate_pass_fraction": float(group["texture_pass"].mean()),
                "texture_practical_pass_count": int(
                    (
                        group["direct_tau_practical_resolution_fraction"]
                        >= float(transmission_gate["texture_practical_fraction_min"])
                    ).sum()
                ),
                "texture_practical_fraction_median": float(
                    group["direct_tau_practical_resolution_fraction"].median()
                ),
                "texture_correlation_median": float(
                    group["direct_tau_correlation"].median()
                ),
                "texture_projection_median": float(
                    group["direct_tau_projection"].median()
                ),
                "texture_rms_median": float(group["direct_tau_y_rms"].median()),
                "texture_permutation_pass_fraction": float(
                    group["permutation_pass"].mean()
                ),
                "edge_gate_pass_fraction": (
                    float(edge_group["edge_pass"].mean())
                    if core_mode == ACTUAL
                    else np.nan
                ),
                "edge_practical_fraction_median": (
                    float(edge_group["direct_practical_resolution_fraction"].median())
                    if core_mode == ACTUAL
                    else np.nan
                ),
                "edge_correlation_median": (
                    float(edge_group["direct_correlation"].median())
                    if core_mode == ACTUAL
                    else np.nan
                ),
                "edge_rms_median": (
                    float(edge_group["direct_y_rms"].median())
                    if core_mode == ACTUAL
                    else np.nan
                ),
                "actual_to_frozen_rms_ratio_min": (
                    float(group["actual_to_frozen_rms_ratio"].min())
                    if core_mode == ACTUAL
                    else np.nan
                ),
                "actual_to_frozen_rms_ratio_median": (
                    float(group["actual_to_frozen_rms_ratio"].median())
                    if core_mode == ACTUAL
                    else np.nan
                ),
            }
        )
pd.DataFrame(transmission_rows).to_csv(
    TABLE_DIR / "run-153_transmission_by_length.csv", index=False
)

# Registered directionality plus a tolerance-aware boundary sensitivity.
dose_cells = dose[dose["aggregation_unit"] == "protein_core_gnk_cell"].copy()
dose_proteins = dose[dose["aggregation_unit"] == "protein_median_over_cells"].copy()
dose_gate = config["acceptance"]["dose_response"]
def tolerant_direction(frame: pd.DataFrame, tolerance: float) -> pd.Series:
    return (
        (frame["tau_rms_spearman"] >= float(dose_gate["spearman_min"]) - tolerance)
        & (
            frame["practical_fraction_spearman"]
            >= float(dose_gate["spearman_min"]) - tolerance
        )
        & (
            frame["tau_rms_a2_minus_a0.25"]
            >= float(dose_gate["rms_growth_floor"]) - tolerance
        )
        & (
            frame["practical_fraction_a2_minus_a0.25"]
            >= float(dose_gate["practical_fraction_growth_floor"]) - tolerance
        )
    )

dose_cells["strict_pass"] = tolerant_direction(dose_cells, 0.0)
dose_cells["tolerance_1e-12_pass"] = tolerant_direction(dose_cells, 1e-12)
dose_proteins["strict_pass"] = tolerant_direction(dose_proteins, 0.0)
dose_proteins["tolerance_1e-12_pass"] = tolerant_direction(
    dose_proteins, 1e-12
)
assert (
    dose["directionality_pass"].astype(bool)
    == pd.concat([dose_cells, dose_proteins])
    .sort_index()["strict_pass"]
    .sort_index()
).all()

direction_rows = []
direction_passes: dict[tuple[str, str], bool] = {}
for mode in [GLOBAL, EXACT]:
    cells = dose_cells[dose_cells["scale_mode"] == mode]
    proteins = dose_proteins[dose_proteins["scale_mode"] == mode]
    for adjudication, column in [
        ("strict", "strict_pass"),
        ("tolerance_1e-12", "tolerance_1e-12_pass"),
    ]:
        cell_minimum, cell_fractions = bin_fractions(cells, column)
        protein_minimum, protein_fractions = bin_fractions(proteins, column)
        protein_fraction = float(proteins[column].mean())
        passed = bool(
            cell_minimum
            >= float(dose_gate["minimum_cell_pass_fraction_per_length_bin"])
            and protein_minimum
            >= float(dose_gate["minimum_protein_pass_fraction_per_length_bin"])
            and protein_fraction >= float(dose_gate["minimum_protein_pass_fraction"])
        )
        direction_passes[(mode, adjudication)] = passed
        for length_bin, group in cells.groupby("length_bin", sort=True):
            protein_group = proteins[proteins["length_bin"] == length_bin]
            direction_rows.append(
                {
                    "scale_mode": mode,
                    "adjudication": adjudication,
                    "length_bin": length_bin,
                    "cell_count": int(len(group)),
                    "cell_pass_fraction": float(group[column].mean()),
                    "protein_count": int(len(protein_group)),
                    "protein_pass_fraction": float(protein_group[column].mean()),
                    "cell_tau_rms_spearman_median": float(
                        group["tau_rms_spearman"].median()
                    ),
                    "cell_practical_spearman_median": float(
                        group["practical_fraction_spearman"].median()
                    ),
                    "cell_tau_rms_growth_median": float(
                        group["tau_rms_a2_minus_a0.25"].median()
                    ),
                    "cell_practical_growth_median": float(
                        group["practical_fraction_a2_minus_a0.25"].median()
                    ),
                    "minimum_cell_fraction_all_bins": cell_minimum,
                    "minimum_protein_fraction_all_bins": protein_minimum,
                    "protein_fraction_overall": protein_fraction,
                    "cell_fractions_json": json.dumps(cell_fractions, sort_keys=True),
                    "protein_fractions_json": json.dumps(
                        protein_fractions, sort_keys=True
                    ),
                }
            )
pd.DataFrame(direction_rows).to_csv(
    TABLE_DIR / "run-153_directionality_by_length.csv", index=False
)

# Direct dose curve and saturation-stress summary.
dose_curve_rows = []
for keys, group in texture[
    (texture["additive_core_mode"] == ACTUAL)
    & (texture["field_arm"] == CANDIDATE)
    & texture["pooled_supported_reference"].eq(True)
].groupby(["scale_mode", "amplitude"], sort=True):
    mode, amplitude = keys
    dose_curve_rows.append(
        {
            "scale_mode": mode,
            "amplitude": float(amplitude),
            "cell_count": int(len(group)),
            "texture_rms_median": float(group["direct_tau_y_rms"].median()),
            "texture_practical_fraction_median": float(
                group["direct_tau_practical_resolution_fraction"].median()
            ),
            "texture_correlation_median": (
                float(group["direct_tau_correlation"].median())
                if group["direct_tau_correlation"].notna().any()
                else np.nan
            ),
        }
    )
pd.DataFrame(dose_curve_rows).to_csv(
    TABLE_DIR / "run-153_dose_curve.csv", index=False
)

global_fit = length_fit(texture, GLOBAL, config)
exact_fit = length_fit(texture, EXACT, config)
scale_audit = pd.DataFrame([global_fit, exact_fit])
recorded_fit = fits.set_index("fit_id")
scale_audit["recorded_slope"] = [
    float(recorded_fit.loc[row["fit_id"], "slope"])
    for row in [global_fit, exact_fit]
]
scale_audit["slope_matches"] = np.isclose(
    scale_audit["slope"], scale_audit["recorded_slope"]
)
scale_audit.to_csv(TABLE_DIR / "run-153_scale_law_summary.csv", index=False)
assert scale_audit["slope_matches"].all()

# Edge-only ablation reporting.
ablation_rows = []
ablation_modes = [
    "steepness_only_margin8",
    "threshold_only_k20",
    "candidate_raw_margin8_w0.6",
    "candidate_raw_margin10_w0.6",
]
for mode in ablation_modes:
    selected_edges = edge[
        (edge["scale_mode"] == mode)
        & (edge["sample_id"] == "reference_shell_h1_complete")
    ]
    for length_bin, group in selected_edges.groupby("length_bin", sort=True):
        ablation_rows.append(
            {
                "scale_mode": mode,
                "length_bin": length_bin,
                "cell_count": int(len(group)),
                "practical_fraction_median": float(
                    group["direct_practical_resolution_fraction"].median()
                ),
                "correlation_median": float(group["direct_correlation"].median()),
                "output_rms_median": float(group["direct_y_rms"].median()),
            }
        )
    ablation_rows.append(
        {
            "scale_mode": mode,
            "length_bin": "ALL",
            "cell_count": int(len(selected_edges)),
            "practical_fraction_median": float(
                selected_edges["direct_practical_resolution_fraction"].median()
            ),
            "correlation_median": float(
                selected_edges["direct_correlation"].median()
            ),
            "output_rms_median": float(selected_edges["direct_y_rms"].median()),
        }
    )
ablation = pd.DataFrame(ablation_rows)
ablation.to_csv(TABLE_DIR / "run-153_ablation_by_length.csv", index=False)

# Independent strict acceptance reconstruction.
calibration_hashes = set(
    catalog.loc[catalog["split"] == "calibration", "sequence_sha256"]
)
confirmation_hashes = set(
    catalog.loc[catalog["split"] == "confirmation", "sequence_sha256"]
)
split_disjoint = calibration_hashes.isdisjoint(confirmation_hashes)
scale_finite = bool(
    scale["finite"].all()
    and np.isfinite(scale["sigma_a_complete_h1"]).all()
    and (scale["sigma_a_complete_h1"] > 0).all()
)
technical_limits = config["acceptance"]["technical"]
zero_edge = edge[
    edge["scale_mode"].isin([GLOBAL, EXACT])
    & np.isclose(edge["amplitude"], 0.0)
]
zero_texture = texture[
    (texture["field_arm"] == CANDIDATE)
    & np.isclose(texture["amplitude"], 0.0)
]
technical_values = {
    "candidate_replay_all": bool(technical["candidate_replay_match"].all()),
    "pair_replay_all": bool(technical["pair_design_replay_match"].all()),
    "square_replay_all": bool(technical["square_design_replay_match"].all()),
    "finite_all": bool(technical["all_outputs_finite"].all()),
    "max_reference_center_abs": float(
        technical["candidate_reference_center_abs"].max()
    ),
    "max_formula_parity_abs": float(technical["formula_parity_max_abs"].max()),
    "max_surrogate_latent_epsilon_abs": float(
        technical["surrogate_latent_epsilon_max_abs"].max()
    ),
    "max_zero_amplitude_residual_rms": max(
        float(zero_edge["direct_y_rms"].max()),
        float(zero_texture["direct_tau_y_rms"].max()),
    ),
}
technical_pass = bool(
    technical_values["candidate_replay_all"]
    and technical_values["pair_replay_all"]
    and technical_values["square_replay_all"]
    and technical_values["finite_all"]
    and technical_values["max_reference_center_abs"]
    <= float(technical_limits["reference_center_atol"])
    and technical_values["max_formula_parity_abs"]
    <= float(technical_limits["formula_parity_atol"])
    and technical_values["max_surrogate_latent_epsilon_abs"]
    <= float(technical_limits["additive_surrogate_epistasis_atol"])
    and technical_values["max_zero_amplitude_residual_rms"]
    <= float(technical_limits["zero_amplitude_atol"])
)
strict_exact_full = bool(
    exact_core_pass
    and transmission_passes[(EXACT, ACTUAL)]
    and direction_passes[(EXACT, "strict")]
    and exact_fit["recomputed_pass"]
    and transmission_passes[(GLOBAL, FROZEN)]
)
tolerant_exact_full = bool(
    exact_core_pass
    and transmission_passes[(EXACT, ACTUAL)]
    and direction_passes[(EXACT, "tolerance_1e-12")]
    and exact_fit["recomputed_pass"]
    and transmission_passes[(GLOBAL, FROZEN)]
)
strict_map = {
    "reference_split_disjoint": split_disjoint,
    "calibration_scale_finite": scale_finite,
    "validation_execution_and_algebra": technical_pass,
    "core_operating_point_selected": True,
    "global_core_calibration_heldout": global_core_pass,
    "global_actual_core_gnk_transmission": transmission_passes[(GLOBAL, ACTUAL)],
    "global_actual_core_gnk_directionality": direction_passes[(GLOBAL, "strict")],
    "global_actual_core_gnk_length_stability": global_fit["recomputed_pass"],
    "frozen_core_texture_positive_control": transmission_passes[(GLOBAL, FROZEN)],
    "exact_per_landscape_beta_diagnostic": strict_exact_full,
    "steepness_threshold_margin_ablation_reported": all(
        len(
            edge[
                (edge["scale_mode"] == mode)
                & (edge["sample_id"] == "reference_shell_h1_complete")
            ]
        )
        == 360
        for mode in ablation_modes
    ),
}
tolerant_map = dict(strict_map)
tolerant_map["global_actual_core_gnk_directionality"] = direction_passes[
    (GLOBAL, "tolerance_1e-12")
]
tolerant_map["exact_per_landscape_beta_diagnostic"] = tolerant_exact_full
acceptance_rows = []
for row in checks.itertuples(index=False):
    acceptance_rows.append(
        {
            "check_id": row.check_id,
            "family": row.family,
            "required": bool(row.required),
            "recorded_pass": bool(row.passed),
            "recomputed_strict_pass": bool(strict_map[row.check_id]),
            "recomputed_tolerance_1e-12_pass": bool(tolerant_map[row.check_id]),
            "strict_matches_recorded": bool(
                strict_map[row.check_id] == bool(row.passed)
            ),
            "sensitivity_note": (
                "rho mathematically equal to 0.80 is stored as "
                "0.7999999999999999 and fails a literal comparison"
                if row.check_id == "global_actual_core_gnk_directionality"
                else ""
            ),
        }
    )
acceptance = pd.DataFrame(acceptance_rows)
acceptance.to_csv(TABLE_DIR / "run-153_acceptance_audit.csv", index=False)
assert acceptance["strict_matches_recorded"].all()
required_ids = list(acceptance.loc[acceptance["required"], "check_id"])
strict_required = sum(strict_map[check_id] for check_id in required_ids)
tolerant_required = sum(tolerant_map[check_id] for check_id in required_ids)
assert (strict_required, len(required_ids), tolerant_required) == (6, 9, 7)
assert summary["required_checks_passed"] == 6
assert summary["decision"] == "CALIBRATED_SIGMOID_TRANSMISSION_FAILED"

protocol_rows = [
    {
        "artifact": name,
        "expected": True,
        "observed": (RUN_DIR / name).is_file(),
        "row_count": (
            len(pd.read_csv(RUN_DIR / name)) if name.endswith(".csv") else np.nan
        ),
    }
    for name in EXPECTED_OUTPUTS
]
protocol_rows += [
    {
        "artifact": "checkpoints/*.json",
        "expected": 360,
        "observed": len(list((RUN_DIR / "checkpoints").glob("*.json"))),
        "row_count": np.nan,
    },
    {
        "artifact": "checkpoint non-standard NaN tokens",
        "expected": "strict JSON would contain 0",
        "observed": nan_token_count,
        "row_count": np.nan,
    },
    {
        "artifact": "logs/labproj-exit-code.txt",
        "expected": "0",
        "observed": (RUN_DIR / "logs/labproj-exit-code.txt").read_text().strip(),
        "row_count": np.nan,
    },
    {
        "artifact": "logs/stderr.txt bytes",
        "expected": 0,
        "observed": (RUN_DIR / "logs/stderr.txt").stat().st_size,
        "row_count": np.nan,
    },
]
pd.DataFrame(protocol_rows).to_csv(
    TABLE_DIR / "run-153_protocol_output_audit.csv", index=False
)

global_actual, global_edges = transmission_frames[(GLOBAL, ACTUAL)]
global_frozen, _ = transmission_frames[(GLOBAL, FROZEN)]
exact_actual, _ = transmission_frames[(EXACT, ACTUAL)]
practical_threshold = float(transmission_gate["texture_practical_fraction_min"])
metrics = {
    "analysis_id": "ANA-060",
    "experiment_id": "EXP-062",
    "run_id": "RUN-153",
    "registered_decision": summary["decision"],
    "registered_required_checks": {"passed": strict_required, "total": 9},
    "independent_tolerance_1e-12_required_checks": {
        "passed": tolerant_required,
        "total": 9,
        "overall_decision_unchanged": True,
    },
    "collection": {
        "file_count": len(collection),
        "size_bytes": int(collection["size_bytes"].sum()),
        "dvc_dir_hash": dvc_hash,
        "remote_sha256_all_match": bool(collection["remote_sha256_matches"].all()),
        "dvc_md5_all_match": bool(collection["dvc_md5_matches"].all()),
        "checkpoint_nonstandard_nan_token_count": nan_token_count,
    },
    "calibration": {
        "global_sigma_a": global_sigma_a,
        "global_sigma_a_cv": float(
            scale["sigma_a_complete_h1"].std(ddof=1)
            / scale["sigma_a_complete_h1"].mean()
        ),
        "selected_beta": float(selected["beta"]),
        "selected_mu": float(selected["mu"]),
        "selected_reference_occupancy": float(selected["reference_occupancy"]),
        "heldout_global_core_pass_cells": int(
            heldout[heldout["scale_mode"] == GLOBAL]["core_pass"].sum()
        ),
        "heldout_global_core_cells": 60,
    },
    "global_edges_a1": {
        "pass_cells": int(global_edges["edge_pass"].sum()),
        "cells": len(global_edges),
        "median_practical_fraction": float(
            global_edges["direct_practical_resolution_fraction"].median()
        ),
        "median_correlation": float(global_edges["direct_correlation"].median()),
        "median_rms": float(global_edges["direct_y_rms"].median()),
    },
    "global_actual_texture_a1": {
        "texture_gate_pass_cells": int(global_actual["texture_pass"].sum()),
        "practical_pass_cells": int(
            (
                global_actual["direct_tau_practical_resolution_fraction"]
                >= practical_threshold
            ).sum()
        ),
        "cells": len(global_actual),
        "median_practical_fraction": float(
            global_actual["direct_tau_practical_resolution_fraction"].median()
        ),
        "median_correlation": float(
            global_actual["direct_tau_correlation"].median()
        ),
        "median_projection": float(
            global_actual["direct_tau_projection"].median()
        ),
        "median_rms": float(global_actual["direct_tau_y_rms"].median()),
        "projection_bootstrap_lower": transmission_values[(GLOBAL, ACTUAL)][
            "projection_bootstrap_lower"
        ],
        "minimum_actual_to_frozen_rms_ratio": float(
            global_actual["actual_to_frozen_rms_ratio"].min()
        ),
        "median_actual_to_frozen_rms_ratio": float(
            global_actual["actual_to_frozen_rms_ratio"].median()
        ),
    },
    "global_frozen_texture_a1": {
        "texture_gate_pass_cells": int(global_frozen["texture_pass"].sum()),
        "practical_pass_cells": int(
            (
                global_frozen["direct_tau_practical_resolution_fraction"]
                >= practical_threshold
            ).sum()
        ),
        "cells": len(global_frozen),
        "median_practical_fraction": float(
            global_frozen["direct_tau_practical_resolution_fraction"].median()
        ),
        "median_correlation": float(
            global_frozen["direct_tau_correlation"].median()
        ),
        "median_projection": float(
            global_frozen["direct_tau_projection"].median()
        ),
        "median_rms": float(global_frozen["direct_tau_y_rms"].median()),
        "projection_bootstrap_lower": transmission_values[(GLOBAL, FROZEN)][
            "projection_bootstrap_lower"
        ],
    },
    "exact_actual_texture_a1": {
        "practical_pass_cells": int(
            (
                exact_actual["direct_tau_practical_resolution_fraction"]
                >= practical_threshold
            ).sum()
        ),
        "cells": len(exact_actual),
        "median_practical_fraction": float(
            exact_actual["direct_tau_practical_resolution_fraction"].median()
        ),
        "median_correlation": float(
            exact_actual["direct_tau_correlation"].median()
        ),
        "median_projection": float(
            exact_actual["direct_tau_projection"].median()
        ),
        "median_rms": float(exact_actual["direct_tau_y_rms"].median()),
    },
    "directionality": {
        "registered_global_pass": direction_passes[(GLOBAL, "strict")],
        "tolerance_1e-12_global_pass": direction_passes[
            (GLOBAL, "tolerance_1e-12")
        ],
        "registered_global_cell_passes": int(
            dose_cells[dose_cells["scale_mode"] == GLOBAL]["strict_pass"].sum()
        ),
        "tolerance_1e-12_global_cell_passes": int(
            dose_cells[dose_cells["scale_mode"] == GLOBAL][
                "tolerance_1e-12_pass"
            ].sum()
        ),
        "global_cell_count": 360,
    },
    "length_stability": {
        "global_slope": global_fit["slope"],
        "global_slope_ci": [
            global_fit["slope_ci_low"],
            global_fit["slope_ci_high"],
        ],
        "exact_slope": exact_fit["slope"],
        "exact_slope_ci": [exact_fit["slope_ci_low"], exact_fit["slope_ci_high"]],
    },
    "hypothesis_effect": {
        "supports": [],
        "weakens": ["HYP-005"],
        "refutes": [],
        "narrow_positive_evidence": [
            "global additive-core calibration",
            "complete-h1 GNK edge transmission",
            "tolerance-aware dose directionality",
            "registered length stability",
        ],
    },
    "public_model_changed": False,
    "stage_4c_4d_4e_hard_stop_retained": True,
}
(TABLE_DIR / "metrics.json").write_text(
    json.dumps(clean(metrics), indent=2, sort_keys=True, allow_nan=False) + "\n",
    encoding="utf-8",
)

print(
    json.dumps(
        {
            "registered_decision": summary["decision"],
            "registered_required_checks": "6/9",
            "tolerance_1e-12_required_checks": "7/9",
            "global_actual_texture_practical_passes": metrics[
                "global_actual_texture_a1"
            ]["practical_pass_cells"],
            "frozen_texture_practical_passes": metrics[
                "global_frozen_texture_a1"
            ]["practical_pass_cells"],
            "dvc_dir_hash": dvc_hash,
        },
        indent=2,
        sort_keys=True,
    )
)
