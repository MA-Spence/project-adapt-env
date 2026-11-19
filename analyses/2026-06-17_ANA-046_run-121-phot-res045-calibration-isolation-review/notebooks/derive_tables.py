from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


RUN_ID = "RUN-121"
RUNTIME_OUTPUT_DIR = "/tank/data/labproj/runtime/project-adapt-env/runs/RUN-121/outputs"
STAGE_ORDER = ["A", "A+F", "A+F+B", "A+F+B+C", "A+F+B+C+D", "A+F+B+C+D+E"]
STAGE_RANK = {stage: idx for idx, stage in enumerate(STAGE_ORDER)}


def find_analysis_dir() -> Path:
    here = Path(__file__).resolve()
    if here.parent.name == "notebooks":
        return here.parents[1]
    return here.parent


def repo_root() -> Path:
    return find_analysis_dir().parents[1]


def short_id(value: object) -> object:
    prefix = "bounded_trait_space__"
    if isinstance(value, str) and value.startswith(prefix):
        return value[len(prefix) :]
    return value


def fmt(value: object, digits: int = 4) -> object:
    if pd.isna(value):
        return ""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return f"{value:.{digits}g}"
    return value


def write_table(path: Path, frame: pd.DataFrame, title: str, note: str | None = None) -> None:
    display = frame.copy()
    for column in display.columns:
        display[column] = display[column].map(fmt)
    parts = [f"# {title}", ""]
    if note:
        parts.extend([f"_{note}_", ""])
    parts.extend([display.to_markdown(index=False), ""])
    path.write_text("\n".join(parts))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


analysis_dir = find_analysis_dir()
source_dir = analysis_dir / "source_tables"
tables_dir = analysis_dir / "tables"
tables_dir.mkdir(exist_ok=True)

summary = read_json(source_dir / "run-121_summary.json")
manifest = pd.read_csv(source_dir / "run-121_remote_output_manifest.tsv", sep="\t")
equivalence = pd.read_csv(source_dir / "run-121_bounded_gate_equivalence.csv")
distances = pd.read_csv(source_dir / "run-121_bounded_gate_distances.csv")
graphlet_summary = pd.read_csv(source_dir / "run-121_graphlet_summary.csv")
selected = pd.read_csv(source_dir / "run-121_selected_bounded_gate_candidates.csv")
comparator = pd.read_csv(source_dir / "run-121_abc_comparator_selected_family_distances.csv")
probe = pd.read_csv(source_dir / "run-121_softmin_capability_probe.csv")
graphlet_catalog = pd.read_csv(source_dir / "run-121_graphlet_catalog.csv")
module_catalog = pd.read_csv(source_dir / "run-121_module_catalog.csv")
prior_117 = read_json(repo_root() / "analyses/2026-06-17_ANA-044_run-117-phot-cumulative-epistasis-calibration-fix-review/tables/metrics.json")
prior_119 = read_json(repo_root() / "analyses/2026-06-17_ANA-045_run-119-phot-raw-latent-softmin-pathology-fix-review/tables/metrics.json")

stage_objective = pd.DataFrame(summary["fix_stage_summaries"])
stage_equivalence = (
    equivalence.groupby("fix_stage", as_index=False)
    .agg(
        n_rows=("parameter_set_id", "size"),
        n_operationally_equivalent=("operationally_equivalent", "sum"),
        n_component_guardrail_pass=("no_component_distance_failure", "sum"),
        n_absolute_margin_pass=("within_absolute_margin", "sum"),
        n_relative_margin_pass=("within_relative_margin", "sum"),
        best_test_distance=("candidate_test_distance", "min"),
        median_test_distance=("candidate_test_distance", "median"),
        best_max_component_distance=("max_component_distance", "min"),
    )
)
stage_summary = stage_objective.merge(stage_equivalence, on="fix_stage", how="left")
stage_summary["sort_order"] = stage_summary["fix_stage"].map(STAGE_RANK)
stage_summary = stage_summary.sort_values("sort_order").drop(columns="sort_order")
stage_summary = stage_summary[
    [
        "fix_stage",
        "fixes_applied",
        "n_candidate_evaluations",
        "n_selected_candidates",
        "best_objective_distance",
        "best_test_distance",
        "median_test_distance",
        "n_operationally_equivalent",
        "n_component_guardrail_pass",
        "n_absolute_margin_pass",
        "n_relative_margin_pass",
        "best_max_component_distance",
    ]
]
write_table(
    tables_dir / "run-121_stage_equivalence_decision.md",
    stage_summary,
    "RUN-121 stage equivalence decision",
    "Operational equivalence is the same RUN-098 best-any-full-ABC decision used by the EXP-047 runner. Stage A+F+B is a probe-only stage and has no selected held-out candidate rows.",
)

best_test_idx = equivalence.groupby("fix_stage")["candidate_test_distance"].idxmin()
best_test = equivalence.loc[best_test_idx].copy()
best_test["sort_order"] = best_test["fix_stage"].map(STAGE_RANK)
best_test = best_test.sort_values("sort_order")
best_test["parameter_set_id"] = best_test["parameter_set_id"].map(short_id)
write_table(
    tables_dir / "run-121_stage_best_test_distances.md",
    best_test[
        [
            "fix_stage",
            "parameter_set_id",
            "candidate_test_distance",
            "reference_test_distance",
            "delta_test_distance",
            "relative_delta_test_distance",
            "max_component_distance",
            "within_absolute_margin",
            "within_relative_margin",
            "no_component_distance_failure",
            "operationally_equivalent",
        ]
    ],
    "RUN-121 best held-out candidate by stage",
)

all_test_rows = equivalence.copy()
all_test_rows["sort_order"] = all_test_rows["fix_stage"].map(STAGE_RANK)
all_test_rows["parameter_set_id"] = all_test_rows["parameter_set_id"].map(short_id)
all_test_rows = all_test_rows.sort_values(["sort_order", "candidate_test_distance"])
write_table(
    tables_dir / "run-121_bounded_gate_test_distances.md",
    all_test_rows[
        [
            "fix_stage",
            "parameter_set_id",
            "candidate_test_distance",
            "delta_test_distance",
            "relative_delta_test_distance",
            "max_component_distance",
            "no_component_distance_failure",
            "operationally_equivalent",
        ]
    ],
    "RUN-121 selected candidate held-out distances",
)

comparator_test = comparator[comparator["graphlet_split"] == "test"].sort_values(
    "distance_to_empirical_primary"
)
write_table(
    tables_dir / "run-121_abc_comparator_test_distances.md",
    comparator_test[
        [
            "candidate_family",
            "candidate_model",
            "parameter_set_id",
            "distance_to_empirical_primary",
            "selection_split",
        ]
    ],
    "RUN-098 comparator test distances copied into RUN-121",
)

best_ids = equivalence.loc[best_test_idx, ["fix_stage", "parameter_set_id"]]
best_distances = distances[distances["graphlet_split"] == "test"].merge(
    best_ids, on=["fix_stage", "parameter_set_id"], how="inner"
)
component_columns = [c for c in best_distances.columns if c.startswith("distance_component__")]
component_rows = []
for _, row in best_distances.iterrows():
    components = [
        (column.replace("distance_component__", ""), row[column])
        for column in component_columns
        if pd.notna(row[column])
    ]
    for component, value in sorted(components, key=lambda item: item[1], reverse=True)[:6]:
        component_rows.append(
            {
                "fix_stage": row["fix_stage"],
                "parameter_set_id": short_id(row["parameter_set_id"]),
                "component": component,
                "component_distance": value,
            }
        )
top_components = pd.DataFrame(component_rows)
top_components["sort_order"] = top_components["fix_stage"].map(STAGE_RANK)
top_components = top_components.sort_values(["sort_order", "component_distance"], ascending=[True, False])
write_table(
    tables_dir / "run-121_top_component_distances.md",
    top_components.drop(columns="sort_order"),
    "RUN-121 top component distances for each stage-best held-out candidate",
)

diagnostic_columns = [
    "accessible_path_fraction_mean",
    "raw_accessible_path_fraction_mean",
    "n_local_maxima_mean",
    "raw_n_local_maxima_mean",
    "magnitude_epistasis_mean_mean",
    "raw_magnitude_epistasis_mean_mean",
    "single_mutant_effect_sd_mean",
    "raw_single_mutant_effect_sd_mean",
    "magnitude_epistasis_to_single_mutant_sd_ratio_mean",
    "raw_magnitude_epistasis_to_single_mutant_sd_ratio_mean",
    "signed_pairwise_epistasis_sd_mean",
    "raw_signed_pairwise_epistasis_sd_mean",
    "sign_epistasis_fraction_local_mean",
    "reciprocal_sign_epistasis_fraction_local_mean",
    "diminishing_returns_epistasis_fraction_mean",
    "edge_energy_top10_fraction_mean",
    "site_energy_top1_fraction_mean",
    "walsh_high_order_fraction_mean",
]
empirical = graphlet_summary[
    (graphlet_summary["graphlet_split"] == "test")
    & (graphlet_summary["candidate_family"] == "empirical_target")
].iloc[0]
diagnostic_rows = [
    {
        "fix_stage": "empirical",
        "role": "empirical_PHOT",
        "parameter_set_id": "latent_phi",
        "test_distance": 0.0,
        **{column: empirical[column] for column in diagnostic_columns},
    }
]
best_lookup = equivalence.loc[
    best_test_idx, ["fix_stage", "parameter_set_id", "candidate_test_distance"]
]
best_diagnostics = best_lookup.merge(
    graphlet_summary[graphlet_summary["graphlet_split"] == "test"][
        ["fix_stage", "parameter_set_id", *diagnostic_columns]
    ],
    on=["fix_stage", "parameter_set_id"],
    how="left",
)
for _, row in best_diagnostics.iterrows():
    diagnostic_rows.append(
        {
            "fix_stage": row["fix_stage"],
            "role": "stage_best_test",
            "parameter_set_id": short_id(row["parameter_set_id"]),
            "test_distance": row["candidate_test_distance"],
            **{column: row[column] for column in diagnostic_columns},
        }
    )
diagnostics = pd.DataFrame(diagnostic_rows)
diagnostics["sort_order"] = diagnostics["fix_stage"].map(
    lambda value: -1 if value == "empirical" else STAGE_RANK[value]
)
write_table(
    tables_dir / "run-121_observable_vs_raw_epistasis_and_landscape.md",
    diagnostics.sort_values("sort_order").drop(columns="sort_order"),
    "RUN-121 observable, raw, and landscape diagnostics",
    "Rows use test graphlets; empirical PHOT is shown once as the target scale.",
)

empirical_targets = {
    "raw_accessible_path_fraction_mean": empirical["raw_accessible_path_fraction_mean"],
    "raw_n_local_maxima_mean": empirical["raw_n_local_maxima_mean"],
    "raw_single_mutant_effect_sd_mean": empirical["raw_single_mutant_effect_sd_mean"],
    "raw_magnitude_epistasis_mean_mean": empirical["raw_magnitude_epistasis_mean_mean"],
    "single_mutant_effect_sd_mean": empirical["single_mutant_effect_sd_mean"],
    "magnitude_epistasis_mean_mean": empirical["magnitude_epistasis_mean_mean"],
}
pathology = diagnostics[diagnostics["role"] == "stage_best_test"].copy()
pathology["raw_accessibility_fraction_of_empirical"] = (
    pathology["raw_accessible_path_fraction_mean"]
    / empirical_targets["raw_accessible_path_fraction_mean"]
)
pathology["raw_local_maxima_fraction_of_empirical"] = (
    pathology["raw_n_local_maxima_mean"] / empirical_targets["raw_n_local_maxima_mean"]
)
pathology["raw_single_mutant_sd_ratio_to_empirical"] = (
    pathology["raw_single_mutant_effect_sd_mean"]
    / empirical_targets["raw_single_mutant_effect_sd_mean"]
)
pathology["raw_magnitude_epistasis_ratio_to_empirical"] = (
    pathology["raw_magnitude_epistasis_mean_mean"]
    / empirical_targets["raw_magnitude_epistasis_mean_mean"]
)
pathology["observable_single_mutant_sd_ratio_to_empirical"] = (
    pathology["single_mutant_effect_sd_mean"]
    / empirical_targets["single_mutant_effect_sd_mean"]
)
pathology["observable_magnitude_epistasis_ratio_to_empirical"] = (
    pathology["magnitude_epistasis_mean_mean"]
    / empirical_targets["magnitude_epistasis_mean_mean"]
)
pathology["sort_order"] = pathology["fix_stage"].map(STAGE_RANK)
write_table(
    tables_dir / "run-121_raw_latent_pathology_diagnostics.md",
    pathology.sort_values("sort_order")[
        [
            "fix_stage",
            "parameter_set_id",
            "test_distance",
            "raw_accessible_path_fraction_mean",
            "raw_accessibility_fraction_of_empirical",
            "raw_n_local_maxima_mean",
            "raw_local_maxima_fraction_of_empirical",
            "raw_single_mutant_effect_sd_mean",
            "raw_single_mutant_sd_ratio_to_empirical",
            "raw_magnitude_epistasis_mean_mean",
            "raw_magnitude_epistasis_ratio_to_empirical",
            "observable_single_mutant_sd_ratio_to_empirical",
            "observable_magnitude_epistasis_ratio_to_empirical",
        ]
    ],
    "RUN-121 raw-latent pathology diagnostics",
    "Ratios compare each stage-best held-out row with the empirical PHOT test target.",
)

best_objective = selected.loc[selected.groupby("fix_stage")["objective_distance"].idxmin()].copy()
best_objective["selection_role"] = "best_objective"
heldout_ids = equivalence.loc[best_test_idx, ["fix_stage", "parameter_set_id"]]
best_heldout = selected.merge(heldout_ids, on=["fix_stage", "parameter_set_id"], how="inner")
best_heldout["selection_role"] = "best_heldout"
parameter_rows = pd.concat([best_objective, best_heldout], ignore_index=True)
parameter_rows = parameter_rows.merge(
    equivalence[["fix_stage", "parameter_set_id", "candidate_test_distance"]],
    on=["fix_stage", "parameter_set_id"],
    how="left",
)
parameter_rows["sort_order"] = parameter_rows["fix_stage"].map(STAGE_RANK)
parameter_rows["parameter_set_id"] = parameter_rows["parameter_set_id"].map(short_id)
parameter_columns = [
    "fix_stage",
    "selection_role",
    "parameter_set_id",
    "objective_distance",
    "calibration_graphlet_distance",
    "candidate_test_distance",
    "guardrail_penalty",
    "n_guardrail_violations",
    "guardrail_violations",
    "raw_latent_objective_distance",
    "raw_latent_objective_penalty",
    "raw_latent_objective_weight",
    "latent_generator",
    "raw__softmin_n_traits",
    "raw__softmin_auxiliary_order",
    "raw__softmin_structured_residual_scale",
    "raw__stability_margin",
    "raw__stability_scale",
    "raw__stability_conservation_power",
    "raw__threshold_steepness",
    "raw__blosum_blend",
    "raw__softmin_tau",
    "raw__softmin_auxiliary_scale",
    "raw__softmin_auxiliary_bias",
    "raw__softmin_auxiliary_correlation",
    "raw__softmin_output_scale",
    "readout_method",
    "readout_midpoint",
    "readout_sigmoid_slope",
    "readout_y_min",
    "readout_y_max",
    "readout_fit_rmse",
    "readout_fit_correlation",
]
write_table(
    tables_dir / "run-121_best_candidate_parameters.md",
    parameter_rows.sort_values(["sort_order", "selection_role"])[parameter_columns],
    "RUN-121 best objective and held-out candidate parameters",
)

probe_summary = (
    probe.groupby(["graphlet_split", "probe_k_traits"], as_index=False)
    .agg(
        n_graphlets=("graphlet_id", "nunique"),
        median_fit_rmse=("probe_fit_rmse", "median"),
        mean_fit_rmse=("probe_fit_rmse", "mean"),
        median_fit_correlation=("probe_fit_correlation", "median"),
        mean_fit_correlation=("probe_fit_correlation", "mean"),
        mean_accessible_path_fraction=("accessible_path_fraction", "mean"),
        mean_target_accessible_path_fraction=("target_accessible_path_fraction", "mean"),
        mean_n_local_maxima=("n_local_maxima", "mean"),
        mean_target_n_local_maxima=("target_n_local_maxima", "mean"),
        mean_single_mutant_sd=("single_mutant_effect_sd", "mean"),
        mean_target_single_mutant_sd=("target_single_mutant_effect_sd", "mean"),
        mean_magnitude_epistasis=("magnitude_epistasis_mean", "mean"),
        mean_target_magnitude_epistasis=("target_magnitude_epistasis_mean", "mean"),
        mean_local_sign_epistasis=("sign_epistasis_fraction_local", "mean"),
        mean_target_local_sign_epistasis=("target_sign_epistasis_fraction_local", "mean"),
    )
)
probe_summary["accessibility_fraction_of_target"] = (
    probe_summary["mean_accessible_path_fraction"] / probe_summary["mean_target_accessible_path_fraction"]
)
probe_summary["local_maxima_fraction_of_target"] = (
    probe_summary["mean_n_local_maxima"] / probe_summary["mean_target_n_local_maxima"]
)
probe_summary["single_mutant_sd_fraction_of_target"] = (
    probe_summary["mean_single_mutant_sd"] / probe_summary["mean_target_single_mutant_sd"]
)
probe_summary["magnitude_epistasis_fraction_of_target"] = (
    probe_summary["mean_magnitude_epistasis"] / probe_summary["mean_target_magnitude_epistasis"]
)
write_table(
    tables_dir / "run-121_softmin_capability_probe_summary.md",
    probe_summary,
    "RUN-121 model-free k-trait soft-min capability probe",
    "Stage A+F+B was intentionally probe-only; these rows summarize direct per-graphlet k-trait additive soft-min fits before sequence-conditioned fitting was re-enabled.",
)

failure_rows = []
best_by_stage = best_diagnostics.merge(
    equivalence[
        [
            "fix_stage",
            "parameter_set_id",
            "candidate_test_distance",
            "no_component_distance_failure",
            "max_component_distance",
        ]
    ],
    on=["fix_stage", "parameter_set_id", "candidate_test_distance"],
    how="left",
)
for _, row in best_by_stage.iterrows():
    stage = row["fix_stage"]
    if stage in {"A", "A+F"}:
        failure_mode = (
            "anchoring and decoupled objective did not repair held-out graphlet distance; "
            "component guardrails fail and raw topology remains weak"
        )
    elif stage == "A+F+B+C":
        failure_mode = (
            "sequence-conditioned k-trait/auxiliary sweep reduces aggregate distance relative "
            "to A/A+F, but remains worse than the full-ABC comparator and component guardrails fail"
        )
    elif stage == "A+F+B+C+D":
        failure_mode = (
            "structured residual and sigmoid readout give the best non-expanded search fit, "
            "but still fail absolute/relative margins and component equivalence"
        )
    else:
        failure_mode = (
            "expanded search slightly improves the structured-residual fit, but does not pass "
            "operational equivalence and remains above the best nested full-ABC context"
        )
    failure_rows.append(
        {
            "fix_stage": stage,
            "best_test_distance": row["candidate_test_distance"],
            "component_guardrails_pass": row["no_component_distance_failure"],
            "max_component_distance": row["max_component_distance"],
            "observable_magnitude_epistasis": row["magnitude_epistasis_mean_mean"],
            "raw_magnitude_epistasis": row["raw_magnitude_epistasis_mean_mean"],
            "raw_single_mutant_effect_sd": row["raw_single_mutant_effect_sd_mean"],
            "accessible_path_fraction": row["accessible_path_fraction_mean"],
            "raw_accessible_path_fraction": row["raw_accessible_path_fraction_mean"],
            "raw_n_local_maxima": row["raw_n_local_maxima_mean"],
            "failure_mode": failure_mode,
        }
    )
failure_modes = pd.DataFrame(failure_rows)
failure_modes["sort_order"] = failure_modes["fix_stage"].map(STAGE_RANK)
write_table(
    tables_dir / "run-121_stage_failure_modes.md",
    failure_modes.sort_values("sort_order").drop(columns="sort_order"),
    "RUN-121 stage failure modes",
)

copied_names = {path.name.replace("run-121_", "") for path in source_dir.iterdir() if path.is_file()}
collection = manifest.copy()
collection["copied_to_analysis"] = collection["file"].map(lambda name: name in copied_names)
collection["retention_note"] = collection["file"].map(
    lambda name: "left on remote runtime volume; aggregate summaries copied"
    if name in {"walsh_terms.csv", "gft_coefficients.csv", "PHOT_CHLRE_Chen_2023.alignment.fasta"}
    else "copied"
)
write_table(
    tables_dir / "run-121_output_collection_audit.md",
    collection[["file", "size_bytes", "sha256", "copied_to_analysis", "retention_note"]],
    "RUN-121 output collection audit",
    f"Remote runtime output directory: {RUNTIME_OUTPUT_DIR}",
)

catalog_summary = pd.DataFrame(
    [
        {
            "catalog": "graphlets",
            "rows": len(graphlet_catalog),
            "columns": len(graphlet_catalog.columns),
            "splits": ",".join(sorted(map(str, graphlet_catalog["graphlet_split"].unique()))),
        },
        {"catalog": "modules", "rows": len(module_catalog), "columns": len(module_catalog.columns), "splits": ""},
    ]
)
write_table(
    tables_dir / "run-121_graphlet_and_module_catalog_summary.md",
    catalog_summary,
    "RUN-121 graphlet and module catalog summary",
)

best_any = comparator_test.iloc[0]
best_nested = comparator_test[comparator_test["candidate_family"] != "empirical_readout_residual"].iloc[0]
best_run121 = equivalence.loc[equivalence["candidate_test_distance"].idxmin()]
run_comparison = pd.DataFrame(
    [
        {
            "run_id": "RUN-117",
            "experiment_id": prior_117["experiment_id"],
            "design": "observable cumulative A-D",
            "best_stage": prior_117["best_run117_stage"],
            "best_test_distance": prior_117["best_run117_test_distance"],
            "n_operationally_equivalent": prior_117["n_operationally_equivalent"],
        },
        {
            "run_id": "RUN-119",
            "experiment_id": prior_119["experiment_id"],
            "design": "raw-latent topology fixes",
            "best_stage": prior_119["best_run119_stage"],
            "best_test_distance": prior_119["best_run119_test_distance"],
            "n_operationally_equivalent": prior_119["n_operationally_equivalent"],
        },
        {
            "run_id": "RUN-121",
            "experiment_id": summary["experiment_id"],
            "design": "RES-045 calibration isolation",
            "best_stage": best_run121["fix_stage"],
            "best_test_distance": best_run121["candidate_test_distance"],
            "n_operationally_equivalent": summary["n_operationally_equivalent"],
        },
        {
            "run_id": "RUN-098",
            "experiment_id": "EXP-039",
            "design": "best nested full-ABC context",
            "best_stage": best_nested["candidate_family"],
            "best_test_distance": best_nested["distance_to_empirical_primary"],
            "n_operationally_equivalent": "",
        },
        {
            "run_id": "RUN-098",
            "experiment_id": "EXP-039",
            "design": "formal best-any full-ABC comparator",
            "best_stage": best_any["candidate_family"],
            "best_test_distance": best_any["distance_to_empirical_primary"],
            "n_operationally_equivalent": "",
        },
    ]
)
write_table(
    tables_dir / "run-121_prior_run_comparison.md",
    run_comparison,
    "RUN-121 comparison with prior PHOT bounded-gate runs",
)

metrics = {
    "run_id": RUN_ID,
    "experiment_id": summary["experiment_id"],
    "n_candidate_evaluations": int(summary["n_candidate_evaluations"]),
    "n_selected_candidates": int(summary["n_selected_candidates"]),
    "n_equivalence_rows": int(summary["n_equivalence_rows"]),
    "n_operationally_equivalent": int(summary["n_operationally_equivalent"]),
    "n_softmin_capability_probe_rows": int(summary["n_softmin_capability_probe_rows"]),
    "best_any_run098_test_distance": float(best_any["distance_to_empirical_primary"]),
    "best_any_run098_candidate_family": best_any["candidate_family"],
    "best_nested_run098_test_distance": float(best_nested["distance_to_empirical_primary"]),
    "best_nested_run098_candidate_family": best_nested["candidate_family"],
    "best_run117_test_distance": float(prior_117["best_run117_test_distance"]),
    "best_run119_test_distance": float(prior_119["best_run119_test_distance"]),
    "best_run121_stage": best_run121["fix_stage"],
    "best_run121_test_distance": float(best_run121["candidate_test_distance"]),
    "best_run121_parameter_set_id": best_run121["parameter_set_id"],
    "run121_improvement_vs_run119_fraction": float(
        (prior_119["best_run119_test_distance"] - best_run121["candidate_test_distance"])
        / prior_119["best_run119_test_distance"]
    ),
    "run121_distance_ratio_vs_best_run117": float(
        best_run121["candidate_test_distance"] / prior_117["best_run117_test_distance"]
    ),
    "run121_distance_ratio_vs_best_nested_run098": float(
        best_run121["candidate_test_distance"] / best_nested["distance_to_empirical_primary"]
    ),
    "stage_summary": json.loads(stage_summary.to_json(orient="records")),
    "raw_latent_pathology_summary": json.loads(
        pathology.sort_values("sort_order")[
            [
                "fix_stage",
                "raw_accessibility_fraction_of_empirical",
                "raw_local_maxima_fraction_of_empirical",
                "raw_single_mutant_sd_ratio_to_empirical",
                "raw_magnitude_epistasis_ratio_to_empirical",
                "observable_single_mutant_sd_ratio_to_empirical",
                "observable_magnitude_epistasis_ratio_to_empirical",
            ]
        ].to_json(orient="records")
    ),
    "softmin_capability_probe_summary": json.loads(probe_summary.to_json(orient="records")),
    "copied_source_payload_bytes": int(sum(path.stat().st_size for path in source_dir.iterdir() if path.is_file())),
}
(tables_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

print(f"Wrote ANA-046 derived tables for {RUN_ID} to {tables_dir}")
