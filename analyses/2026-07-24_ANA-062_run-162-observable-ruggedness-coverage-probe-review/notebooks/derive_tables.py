from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ANALYSIS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = ANALYSIS_DIR.parents[1]
TABLE_DIR = ANALYSIS_DIR / "tables"
RUN_ROOT = (
    REPO_ROOT
    / "data/processed/observable-ruggedness-coverage-probe/RUN-162"
)
EXP_DIR = (
    REPO_ROOT
    / "experiments/2026-07-24_EXP-064_observable-ruggedness-coverage-probe"
)
RUN_RECORD = EXP_DIR / "runs/RUN-162.yaml"
POINTER = RUN_ROOT.parent / "RUN-162.dvc"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def strict_json(path: Path) -> object:
    def reject_constant(token: str) -> None:
        raise ValueError(f"Non-finite JSON token {token} in {path}")

    return json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=reject_constant,
    )


def finite(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def maximum_numeric_difference(
    left: pd.DataFrame,
    right: pd.DataFrame,
    keys: list[str],
    columns: list[str],
) -> float:
    merged = left.merge(right, on=keys, suffixes=("_independent", "_recorded"))
    differences: list[float] = []
    for column in columns:
        a = pd.to_numeric(merged[f"{column}_independent"], errors="coerce")
        b = pd.to_numeric(merged[f"{column}_recorded"], errors="coerce")
        valid = a.notna() & b.notna()
        differences.extend(np.abs(a[valid] - b[valid]).tolist())
        if bool((a.isna() != b.isna()).any()):
            return math.inf
    return float(max(differences, default=0.0))


def write_csv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(TABLE_DIR / name, index=False)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load((EXP_DIR / "config.yaml").read_text(encoding="utf-8"))
    run_record = yaml.safe_load(RUN_RECORD.read_text(encoding="utf-8"))
    pointer = yaml.safe_load(POINTER.read_text(encoding="utf-8"))["outs"][0]
    summary = strict_json(RUN_ROOT / "summary.json")
    preflight = strict_json(RUN_ROOT / "preflight.json")

    generated = pd.read_csv(RUN_ROOT / "generated_descriptor_metrics.csv")
    bands = pd.read_csv(RUN_ROOT / "empirical_bands_resolved.csv")
    recorded_coverage = pd.read_csv(RUN_ROOT / "coverage_matrix.csv")
    recorded_reach = pd.read_csv(RUN_ROOT / "rugged_end_reachability.csv")
    recorded_specificity = pd.read_csv(RUN_ROOT / "specificity_table.csv")
    recorded_operating = pd.read_csv(RUN_ROOT / "selected_operating_region.csv")
    technical = pd.read_csv(RUN_ROOT / "technical_checks.csv")

    descriptors = list(config["descriptors"]["rugged_direction"])
    rugged_templates = sorted(
        bands.loc[bands["rugged_end_template"].astype(bool), "template"].unique()
    )
    primary = generated[
        generated["representation"].eq("binary_r8")
        & generated["support_mode"].eq("structured")
        & generated["amplitude"].gt(0)
    ].copy()

    coverage_rows: list[dict[str, object]] = []
    for band in bands.to_dict(orient="records"):
        values = pd.to_numeric(
            primary.loc[primary["template"].eq(band["template"]), band["descriptor"]],
            errors="coerce",
        ).dropna()
        reachable_min = float(values.min()) if len(values) else math.nan
        reachable_max = float(values.max()) if len(values) else math.nan
        median = float(band["empirical_median"])
        rugged_target = float(band["empirical_rugged_q90"])
        contains_median = bool(
            finite(reachable_min)
            and finite(median)
            and reachable_min <= median + 1e-12
            and reachable_max >= median - 1e-12
        )
        if band["rugged_direction"] == "increase":
            reaches_rugged = bool(
                finite(reachable_max)
                and finite(rugged_target)
                and reachable_max >= rugged_target - 1e-12
            )
            shortfall = (
                max(0.0, rugged_target - reachable_max)
                if finite(rugged_target) and finite(reachable_max)
                else math.nan
            )
        else:
            reaches_rugged = bool(
                finite(reachable_min)
                and finite(rugged_target)
                and reachable_min <= rugged_target + 1e-12
            )
            shortfall = (
                max(0.0, reachable_min - rugged_target)
                if finite(rugged_target) and finite(reachable_min)
                else math.nan
            )
        coverage_rows.append(
            {
                **band,
                "reachable_min": reachable_min,
                "reachable_max": reachable_max,
                "contains_empirical_median": contains_median,
                "reaches_empirical_rugged_q90": reaches_rugged,
                "cover": bool(contains_median and reaches_rugged),
                "rugged_shortfall": shortfall,
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    coverage_numeric_difference = maximum_numeric_difference(
        coverage,
        recorded_coverage,
        ["template", "descriptor"],
        ["reachable_min", "reachable_max", "rugged_shortfall"],
    )
    coverage_boolean_match = bool(
        coverage.merge(
            recorded_coverage,
            on=["template", "descriptor"],
            suffixes=("_independent", "_recorded"),
        )
        .apply(
            lambda row: all(
                bool(row[f"{column}_independent"])
                == bool(row[f"{column}_recorded"])
                for column in (
                    "contains_empirical_median",
                    "reaches_empirical_rugged_q90",
                    "cover",
                )
            ),
            axis=1,
        )
        .all()
    )

    coverage_by_template = (
        coverage.groupby(["rugged_end_template", "template"], as_index=False)
        .agg(
            n_descriptors=("descriptor", "size"),
            n_available=("descriptor_available", "sum"),
            n_cover=("cover", "sum"),
            n_contains_median=("contains_empirical_median", "sum"),
            n_reaches_rugged_q90=("reaches_empirical_rugged_q90", "sum"),
        )
    )
    coverage_by_template["cover_fraction_all_8"] = (
        coverage_by_template["n_cover"] / coverage_by_template["n_descriptors"]
    )
    coverage_by_template["cover_fraction_available"] = (
        coverage_by_template["n_cover"] / coverage_by_template["n_available"]
    )
    coverage_by_template["range_template_pass"] = (
        coverage_by_template["cover_fraction_all_8"]
        >= float(config["gates"]["range_descriptor_fraction_min"])
    )

    coverage_by_descriptor = (
        coverage.groupby(["rugged_end_template", "descriptor"], as_index=False)
        .agg(
            n_templates=("template", "size"),
            n_available=("descriptor_available", "sum"),
            n_cover=("cover", "sum"),
            n_contains_median=("contains_empirical_median", "sum"),
            n_reaches_rugged_q90=("reaches_empirical_rugged_q90", "sum"),
            maximum_rugged_shortfall=("rugged_shortfall", "max"),
        )
    )
    coverage_by_descriptor["cover_fraction"] = (
        coverage_by_descriptor["n_cover"]
        / coverage_by_descriptor["n_templates"]
    )

    minimum_joint = int(config["gates"]["joint_descriptors_min"])
    accessibility_floor = float(config["gates"]["accessibility_floor"])
    reach_rows: list[dict[str, object]] = []
    for template in rugged_templates:
        template_bands = bands[bands["template"].eq(template)].set_index("descriptor")
        candidates: list[dict[str, object]] = []
        for row in primary[primary["template"].eq(template)].to_dict(orient="records"):
            names = [
                descriptor
                for descriptor in descriptors
                if finite(row[descriptor])
                and finite(template_bands.loc[descriptor, "empirical_q05"])
                and float(row[descriptor])
                >= float(template_bands.loc[descriptor, "empirical_q05"]) - 1e-12
                and float(row[descriptor])
                <= float(template_bands.loc[descriptor, "empirical_q95"]) + 1e-12
            ]
            candidates.append(
                {
                    **row,
                    "n_descriptors_in_band": len(names),
                    "descriptors_in_band": ";".join(sorted(names)),
                    "joint": len(names) >= minimum_joint,
                    "navigable": float(row["accessible_path_fraction"])
                    >= accessibility_floor,
                }
            )
        frame = pd.DataFrame(candidates)
        best_count = int(frame["n_descriptors_in_band"].max())
        reach_rows.append(
            {
                "template": template,
                "n_empirical_descriptors_available": int(
                    template_bands["descriptor_available"].sum()
                ),
                "maximum_descriptors_in_band": best_count,
                "any_joint_reach_pass": bool(frame["joint"].any()),
                "any_navigable_joint_reach_pass": bool(
                    (frame["joint"] & frame["navigable"]).any()
                ),
                "maximum_accessibility_among_joint_rows": float(
                    frame.loc[frame["joint"], "accessible_path_fraction"].max()
                )
                if bool(frame["joint"].any())
                else math.nan,
            }
        )
    reach = pd.DataFrame(reach_rows)
    reach_compare = reach.merge(recorded_reach, on="template")
    reach_match = bool(
        (
            reach_compare["maximum_descriptors_in_band"]
            == reach_compare["n_descriptors_in_band"]
        ).all()
        and (
            reach_compare["any_joint_reach_pass_x"]
            == reach_compare["any_joint_reach_pass_y"]
        ).all()
        and (
            reach_compare["any_navigable_joint_reach_pass_x"]
            == reach_compare["any_navigable_joint_reach_pass_y"]
        ).all()
    )

    increasing = [
        descriptor
        for descriptor, direction in config["descriptors"]["rugged_direction"].items()
        if direction == "increase"
    ]
    specificity_rows: list[dict[str, object]] = []
    for template in rugged_templates:
        template_rows = generated[generated["template"].eq(template)]
        masks = {
            "composite": (
                template_rows["representation"].eq("binary_r8")
                & template_rows["support_mode"].eq("structured")
                & template_rows["amplitude"].gt(0)
            ),
            "additive": (
                template_rows["representation"].eq("binary_r8")
                & template_rows["support_mode"].eq("structured")
                & template_rows["amplitude"].eq(0)
            ),
            "shuffled": (
                template_rows["representation"].eq("binary_r8")
                & template_rows["support_mode"].eq("shuffled")
                & template_rows["amplitude"].gt(0)
            ),
            "random": (
                template_rows["representation"].eq("binary_r8")
                & template_rows["support_mode"].eq("random")
                & template_rows["amplitude"].gt(0)
            ),
        }
        for descriptor in increasing:
            maxima = {
                name: float(template_rows.loc[mask, descriptor].max())
                for name, mask in masks.items()
            }
            required = float(config["gates"]["specificity_deltas"][descriptor])
            deltas = {
                "composite_minus_additive": maxima["composite"] - maxima["additive"],
                "composite_minus_shuffled": maxima["composite"] - maxima["shuffled"],
                "composite_minus_random": maxima["composite"] - maxima["random"],
            }
            specificity_rows.append(
                {
                    "template": template,
                    "descriptor": descriptor,
                    "composite_max": maxima["composite"],
                    "additive_core_max": maxima["additive"],
                    "shuffled_support_max": maxima["shuffled"],
                    "random_support_max": maxima["random"],
                    "required_delta": required,
                    **deltas,
                    "specificity_pass": all(
                        value >= required - 1e-12 for value in deltas.values()
                    ),
                }
            )
    specificity = pd.DataFrame(specificity_rows)
    specificity_numeric_difference = maximum_numeric_difference(
        specificity,
        recorded_specificity,
        ["template", "descriptor"],
        [
            "composite_max",
            "additive_core_max",
            "shuffled_support_max",
            "random_support_max",
            "composite_minus_additive",
            "composite_minus_shuffled",
            "composite_minus_random",
        ],
    )
    specificity_boolean_match = bool(
        specificity.merge(
            recorded_specificity,
            on=["template", "descriptor"],
            suffixes=("_independent", "_recorded"),
        )
        .eval("specificity_pass_independent == specificity_pass_recorded")
        .all()
    )
    specificity_by_descriptor = (
        specificity.groupby("descriptor", as_index=False)
        .agg(
            n_template_passes=("specificity_pass", "sum"),
            n_templates=("template", "size"),
            minimum_composite_minus_additive=("composite_minus_additive", "min"),
            minimum_composite_minus_shuffled=("composite_minus_shuffled", "min"),
            minimum_composite_minus_random=("composite_minus_random", "min"),
        )
    )
    specificity_by_descriptor["template_pass_fraction"] = (
        specificity_by_descriptor["n_template_passes"]
        / specificity_by_descriptor["n_templates"]
    )
    specificity_by_descriptor["descriptor_pass"] = (
        specificity_by_descriptor["template_pass_fraction"]
        >= float(config["gates"]["specificity_landscape_fraction_min"])
    )

    selections: list[dict[str, object]] = []
    for row in primary[primary["template"].isin(rugged_templates)].to_dict(
        orient="records"
    ):
        template_bands = bands[bands["template"].eq(row["template"])].set_index(
            "descriptor"
        )
        count = sum(
            finite(row[descriptor])
            and finite(template_bands.loc[descriptor, "empirical_q05"])
            and float(row[descriptor])
            >= float(template_bands.loc[descriptor, "empirical_q05"]) - 1e-12
            and float(row[descriptor])
            <= float(template_bands.loc[descriptor, "empirical_q95"]) + 1e-12
            for descriptor in descriptors
        )
        selections.append(
            {
                "template": row["template"],
                "beta": row["beta"],
                "p0": row["p0"],
                "amplitude": row["amplitude"],
                "n_descriptors_in_band": count,
                "navigable": row["accessible_path_fraction"] >= accessibility_floor,
            }
        )
    operating = (
        pd.DataFrame(selections)
        .groupby(["beta", "p0", "amplitude"], as_index=False)
        .agg(
            mean_descriptors_in_band=("n_descriptors_in_band", "mean"),
            maximum_descriptors_in_band=("n_descriptors_in_band", "max"),
            rugged_template_fraction_joint=(
                "n_descriptors_in_band",
                lambda values: float(np.mean(np.asarray(values) >= minimum_joint)),
            ),
            navigable_fraction=("navigable", "mean"),
            n_rows=("template", "size"),
        )
        .sort_values(
            [
                "rugged_template_fraction_joint",
                "mean_descriptors_in_band",
                "navigable_fraction",
            ],
            ascending=False,
        )
    )
    operating_difference = maximum_numeric_difference(
        operating,
        recorded_operating,
        ["beta", "p0", "amplitude"],
        [
            "mean_descriptors_in_band",
            "maximum_descriptors_in_band",
            "rugged_template_fraction_joint",
            "navigable_fraction",
            "n_rows",
        ],
    )

    rugged_cover = coverage_by_template[
        coverage_by_template["rugged_end_template"].astype(bool)
    ]
    g1_observed_fraction = float(rugged_cover["range_template_pass"].mean())
    rugged_joint_fraction = float(reach["any_joint_reach_pass"].mean())
    rugged_navigable_fraction = float(
        reach["any_navigable_joint_reach_pass"].mean()
    )
    specificity_descriptor_fraction = float(
        specificity_by_descriptor["descriptor_pass"].mean()
    )
    landscape_floor = float(
        config["gates"]["rugged_end_landscape_fraction_min"]
    )
    independent_gates = {
        "execution_and_pairing": bool(technical["pass"].astype(bool).all()),
        "empirical_band_reproduction": bool(
            technical.loc[
                technical["check"].astype(str).str.startswith(
                    (
                        "sha256_",
                        "tier1_scope_catalog",
                        "published_pooled_band_",
                        "published_dataset_median_",
                    )
                ),
                "pass",
            ]
            .astype(bool)
            .all()
        ),
        "G1_range_coverage": bool(g1_observed_fraction >= landscape_floor),
        "G2_joint_rugged_end_reachability": bool(
            rugged_joint_fraction >= landscape_floor
        ),
        "G3_specificity": bool(
            specificity_descriptor_fraction
            > float(config["gates"]["specificity_descriptor_majority_fraction"])
        ),
        "G4_navigable_rugged_sanity": bool(
            rugged_navigable_fraction >= landscape_floor
        ),
    }
    if not independent_gates["execution_and_pairing"]:
        independent_decision = "EXECUTION_OR_PAIRING_INVALID"
    elif not independent_gates["empirical_band_reproduction"]:
        independent_decision = "EMPIRICAL_BANDS_NOT_REPRODUCED"
    elif not (
        independent_gates["G1_range_coverage"]
        and independent_gates["G2_joint_rugged_end_reachability"]
    ):
        independent_decision = "OBSERVABLE_CONSTRUCTION_LIMIT"
    elif not independent_gates["G3_specificity"]:
        independent_decision = "COVERAGE_NONSPECIFIC"
    elif not independent_gates["G4_navigable_rugged_sanity"]:
        independent_decision = "NAVIGABILITY_DEGENERATE"
    else:
        independent_decision = "PROCEED_TO_FITTING"

    recorded_decision = pd.read_csv(RUN_ROOT / "decision_table.csv")
    acceptance_rows: list[dict[str, object]] = []
    observed = {
        "G1_range_coverage": g1_observed_fraction,
        "G2_joint_rugged_end_reachability": rugged_joint_fraction,
        "G3_specificity": specificity_descriptor_fraction,
        "G4_navigable_rugged_sanity": rugged_navigable_fraction,
    }
    for gate, independent_pass in independent_gates.items():
        registered_pass = bool(
            recorded_decision.loc[recorded_decision["gate"].eq(gate), "pass"].iloc[0]
        )
        acceptance_rows.append(
            {
                "gate": gate,
                "registered_pass": registered_pass,
                "independent_pass": independent_pass,
                "match": registered_pass == independent_pass,
                "independent_observed_fraction": observed.get(gate, math.nan),
            }
        )
    acceptance_rows.append(
        {
            "gate": "terminal_decision",
            "registered_pass": True,
            "independent_pass": True,
            "match": independent_decision == summary["decision"],
            "independent_observed_fraction": math.nan,
            "registered_decision": summary["decision"],
            "independent_decision": independent_decision,
        }
    )
    acceptance = pd.DataFrame(acceptance_rows)

    frozen_metrics_path = (
        REPO_ROOT
        / config["empirical"]["frozen_graphlet_metrics"]["path"]
    )
    frozen = pd.read_csv(frozen_metrics_path)
    rugged_target_rows: list[dict[str, object]] = []
    for template in rugged_templates:
        subset = frozen[frozen["dataset_id"].eq(template)].copy()
        density = subset["n_local_maxima"] / np.power(2, subset["dimension"])
        rugged_target_rows.append(
            {
                "template": template,
                "n_graphlets": len(subset),
                "flat_graphlet_fraction": float(
                    np.mean(
                        pd.to_numeric(subset["fitness_sd"], errors="coerce")
                        .fillna(0)
                        .le(1e-12)
                    )
                ),
                "median_local_maxima_density": float(density.median()),
                "maximum_local_maxima_density": float(density.max()),
                "finite_random_walk_fraction": float(
                    subset["random_walk_lag1_autocorrelation"].notna().mean()
                ),
                "finite_dirichlet_localisation_fraction": float(
                    subset["edge_energy_top10_fraction"].notna().mean()
                ),
            }
        )
    rugged_target_diagnostic = pd.DataFrame(rugged_target_rows)

    expected_outputs = list(config["outputs"]["expected"])
    protocol_rows = [
        {
            "check": f"declared_output_{name}",
            "expected": True,
            "observed": (RUN_ROOT / name).is_file(),
            "pass": (RUN_ROOT / name).is_file(),
        }
        for name in expected_outputs
    ]
    checkpoint_dirs = [path for path in (RUN_ROOT / "checkpoints").iterdir() if path.is_dir()]
    done_rows = [strict_json(path / "done.json") for path in checkpoint_dirs]
    protocol_rows.extend(
        [
            {
                "check": "generated_descriptor_row_count",
                "expected": 259200,
                "observed": len(generated),
                "pass": len(generated) == 259200,
            },
            {
                "check": "checkpoint_directory_count",
                "expected": 300,
                "observed": len(checkpoint_dirs),
                "pass": len(checkpoint_dirs) == 300,
            },
            {
                "check": "checkpoint_metric_row_sum",
                "expected": 259200,
                "observed": sum(int(row["n_metric_rows"]) for row in done_rows),
                "pass": sum(int(row["n_metric_rows"]) for row in done_rows)
                == 259200,
            },
            {
                "check": "complete_generation_grid",
                "expected": 259200,
                "observed": generated[
                    [
                        "template",
                        "seed",
                        "representation",
                        "support_mode",
                        "beta",
                        "p0",
                        "amplitude",
                    ]
                ].drop_duplicates().shape[0],
                "pass": generated[
                    [
                        "template",
                        "seed",
                        "representation",
                        "support_mode",
                        "beta",
                        "p0",
                        "amplitude",
                    ]
                ].drop_duplicates().shape[0]
                == 259200,
            },
        ]
    )
    protocol_audit = pd.DataFrame(protocol_rows)

    source_rows: list[dict[str, object]] = []
    for key in (
        "published_key_bands",
        "published_dataset_medians",
        "frozen_graphlet_metrics",
        "frozen_scope_audit",
    ):
        record = config["empirical"][key]
        path = REPO_ROOT / record["path"]
        actual = sha256(path)
        source_rows.append(
            {
                "source": key,
                "path": record["path"],
                "expected_sha256": record["sha256"],
                "actual_sha256": actual,
                "pass": actual == record["sha256"],
            }
        )
    source_audit = pd.DataFrame(source_rows)

    files = [path for path in RUN_ROOT.rglob("*") if path.is_file()]
    json_files = [path for path in files if path.suffix == ".json"]
    strict_json_pass = True
    for path in json_files:
        try:
            strict_json(path)
        except (json.JSONDecodeError, ValueError):
            strict_json_pass = False
    collection_rows = [
        {
            "check": "dvc_pointer_file_count",
            "expected": int(pointer["nfiles"]),
            "observed": len(files),
            "pass": int(pointer["nfiles"]) == len(files),
        },
        {
            "check": "dvc_pointer_byte_count",
            "expected": int(pointer["size"]),
            "observed": sum(path.stat().st_size for path in files),
            "pass": int(pointer["size"])
            == sum(path.stat().st_size for path in files),
        },
        {
            "check": "strict_json_parse",
            "expected": len(json_files),
            "observed": len(json_files) if strict_json_pass else -1,
            "pass": strict_json_pass,
        },
        {
            "check": "scheduler_exit_code",
            "expected": 0,
            "observed": int(run_record["attempts"][-1]["scheduler_exit_code"]),
            "pass": int(run_record["attempts"][-1]["scheduler_exit_code"]) == 0,
        },
        {
            "check": "stderr_empty",
            "expected": 0,
            "observed": (RUN_ROOT / "logs/stderr.txt").stat().st_size,
            "pass": (RUN_ROOT / "logs/stderr.txt").stat().st_size == 0,
        },
    ]
    collection_audit = pd.DataFrame(collection_rows)

    metrics = {
        "run_id": "RUN-162",
        "experiment_id": "EXP-064",
        "scheduler_job_id": str(run_record["scheduler_job_id"]),
        "scheduler_exit_code": int(run_record["attempts"][-1]["scheduler_exit_code"]),
        "decision_registered": summary["decision"],
        "decision_independent": independent_decision,
        "n_generated_rows": len(generated),
        "n_checkpoint_directories": len(checkpoint_dirs),
        "n_technical_checks": len(technical),
        "n_technical_failures": int((~technical["pass"].astype(bool)).sum()),
        "g1_rugged_template_pass_fraction": g1_observed_fraction,
        "g2_rugged_joint_fraction": rugged_joint_fraction,
        "g3_specificity_descriptor_pass_fraction": specificity_descriptor_fraction,
        "g4_rugged_navigable_fraction": rugged_navigable_fraction,
        "coverage_numeric_maximum_difference": coverage_numeric_difference,
        "coverage_boolean_match": coverage_boolean_match,
        "reachability_match": reach_match,
        "specificity_numeric_maximum_difference": specificity_numeric_difference,
        "specificity_boolean_match": specificity_boolean_match,
        "operating_region_numeric_maximum_difference": operating_difference,
        "dvc_directory_hash": pointer["md5"],
        "dvc_file_count": int(pointer["nfiles"]),
        "dvc_byte_count": int(pointer["size"]),
        "protocol_pass": bool(protocol_audit["pass"].all()),
        "source_hash_pass": bool(source_audit["pass"].all()),
        "collection_pass": bool(collection_audit["pass"].all()),
        "public_model_changed": bool(summary["public_model_changed"]),
        "uses_exp050_transfer_function": bool(summary["uses_exp050_transfer_function"]),
        "stage_4c_4d_4e_hard_stop": summary["stage_4c_4d_4e_hard_stop"],
    }

    write_csv(acceptance, "run-162_acceptance_audit.csv")
    write_csv(coverage_by_template, "run-162_coverage_by_template.csv")
    write_csv(coverage_by_descriptor, "run-162_coverage_by_descriptor.csv")
    write_csv(reach, "run-162_rugged_end_audit.csv")
    write_csv(
        rugged_target_diagnostic,
        "run-162_rugged_target_diagnostic.csv",
    )
    write_csv(specificity, "run-162_specificity_audit.csv")
    write_csv(
        specificity_by_descriptor,
        "run-162_specificity_by_descriptor.csv",
    )
    write_csv(operating, "run-162_operating_region_audit.csv")
    write_csv(protocol_audit, "run-162_protocol_output_audit.csv")
    write_csv(source_audit, "run-162_source_hash_audit.csv")
    write_csv(collection_audit, "run-162_collection_audit.csv")
    (TABLE_DIR / "metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (ANALYSIS_DIR / "source_tables").mkdir(exist_ok=True)
    for source, destination in (
        (RUN_ROOT / "summary.json", "run-162_summary.json"),
        (RUN_ROOT / "preflight.json", "run-162_preflight.json"),
        (
            RUN_ROOT / "checkpoint_resume_summary.json",
            "run-162_checkpoint_resume_summary.json",
        ),
    ):
        (ANALYSIS_DIR / "source_tables" / destination).write_bytes(
            source.read_bytes()
        )

    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
