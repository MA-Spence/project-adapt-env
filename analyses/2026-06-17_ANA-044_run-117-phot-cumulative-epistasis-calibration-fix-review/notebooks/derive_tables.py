from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


RUN_ID = "RUN-117"
RUNTIME_OUTPUT_DIR = "/tank/data/labproj/runtime/project-adapt-env/runs/RUN-117/outputs"
STAGE_ORDER = ["A", "A+B", "A+B+C", "A+B+C+D"]
STAGE_RANK = {stage: idx for idx, stage in enumerate(STAGE_ORDER)}


def find_analysis_dir() -> Path:
    here = Path(__file__).resolve()
    if here.parent.name == "notebooks":
        return here.parents[1]
    return here.parent


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


analysis_dir = find_analysis_dir()
source_dir = analysis_dir / "source_tables"
tables_dir = analysis_dir / "tables"
tables_dir.mkdir(exist_ok=True)

summary = json.loads((source_dir / "run-117_summary.json").read_text())
manifest = pd.read_csv(source_dir / "run-117_remote_output_manifest.tsv", sep="\t")
equivalence = pd.read_csv(source_dir / "run-117_bounded_gate_equivalence.csv")
distances = pd.read_csv(source_dir / "run-117_bounded_gate_distances.csv")
graphlet_summary = pd.read_csv(source_dir / "run-117_graphlet_summary.csv")
selected = pd.read_csv(source_dir / "run-117_selected_bounded_gate_candidates.csv")
comparator = pd.read_csv(source_dir / "run-117_abc_comparator_selected_family_distances.csv")
graphlet_catalog = pd.read_csv(source_dir / "run-117_graphlet_catalog.csv")
module_catalog = pd.read_csv(source_dir / "run-117_module_catalog.csv")

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
    tables_dir / "run-117_stage_equivalence_decision.md",
    stage_summary,
    "RUN-117 stage equivalence decision",
    "Operational equivalence is the same RUN-098 best-any-full-ABC decision used by the EXP-045 runner.",
)

best_test_idx = equivalence.groupby("fix_stage")["candidate_test_distance"].idxmin()
best_test = equivalence.loc[best_test_idx].copy()
best_test["sort_order"] = best_test["fix_stage"].map(STAGE_RANK)
best_test = best_test.sort_values("sort_order")
best_test["parameter_set_id"] = best_test["parameter_set_id"].map(short_id)
write_table(
    tables_dir / "run-117_stage_best_test_distances.md",
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
    "RUN-117 best held-out candidate by stage",
)

all_test_rows = equivalence.copy()
all_test_rows["sort_order"] = all_test_rows["fix_stage"].map(STAGE_RANK)
all_test_rows["parameter_set_id"] = all_test_rows["parameter_set_id"].map(short_id)
all_test_rows = all_test_rows.sort_values(["sort_order", "candidate_test_distance"])
write_table(
    tables_dir / "run-117_bounded_gate_test_distances.md",
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
    "RUN-117 selected candidate held-out distances",
)

comparator_test = comparator[comparator["graphlet_split"] == "test"].sort_values(
    "distance_to_empirical_primary"
)
write_table(
    tables_dir / "run-117_abc_comparator_test_distances.md",
    comparator_test[
        [
            "candidate_family",
            "candidate_model",
            "parameter_set_id",
            "distance_to_empirical_primary",
            "selection_split",
        ]
    ],
    "RUN-098 comparator test distances copied into RUN-117",
)

best_ids = set(equivalence.loc[best_test_idx, "parameter_set_id"])
best_distances = distances[
    (distances["graphlet_split"] == "test") & (distances["parameter_set_id"].isin(best_ids))
]
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
    tables_dir / "run-117_top_component_distances.md",
    top_components.drop(columns="sort_order"),
    "RUN-117 top component distances for each stage-best held-out candidate",
)

diagnostic_columns = [
    "accessible_path_fraction_mean",
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
    tables_dir / "run-117_observable_vs_raw_epistasis_and_landscape.md",
    diagnostics.sort_values("sort_order").drop(columns="sort_order"),
    "RUN-117 observable, raw, and landscape diagnostics",
    "Rows use test graphlets; empirical PHOT is shown once as the target scale.",
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
write_table(
    tables_dir / "run-117_best_candidate_parameters.md",
    parameter_rows.sort_values(["sort_order", "selection_role"])[
        [
            "fix_stage",
            "selection_role",
            "parameter_set_id",
            "objective_distance",
            "calibration_graphlet_distance",
            "candidate_test_distance",
            "guardrail_penalty",
            "n_guardrail_violations",
            "guardrail_violations",
            "raw__stability_margin",
            "raw__threshold_steepness",
            "raw__epistasis_strength",
            "raw__empirical_pairwise_strength",
            "raw__trait_space_epistasis_bound",
            "readout_method",
            "raw__global_readout_midpoint",
            "raw__global_readout_slope",
            "raw__global_readout_weight",
            "raw__specific_epistasis_variance_ratio",
        ]
    ],
    "RUN-117 best objective and held-out candidate parameters",
)

failure_rows = []
for _, row in best_diagnostics.iterrows():
    stage = row["fix_stage"]
    eq_row = equivalence[equivalence["parameter_set_id"] == row["parameter_set_id"]].iloc[0]
    if stage in {"A+B", "A+B+C"}:
        failure_mode = "collapsed zero-variation landscape; component guardrails fail"
    elif stage == "A":
        failure_mode = (
            "component guardrails pass but aggregate distance fails; "
            "accessibility low and raw latent epistasis remains large"
        )
    else:
        failure_mode = (
            "component guardrails pass but aggregate distance fails; "
            "accessibility low and raw latent variation is near zero"
        )
    failure_rows.append(
        {
            "fix_stage": stage,
            "best_test_distance": eq_row["candidate_test_distance"],
            "component_guardrails_pass": eq_row["no_component_distance_failure"],
            "observable_magnitude_epistasis": row["magnitude_epistasis_mean_mean"],
            "raw_magnitude_epistasis": row["raw_magnitude_epistasis_mean_mean"],
            "accessible_path_fraction": row["accessible_path_fraction_mean"],
            "failure_mode": failure_mode,
        }
    )
failure_modes = pd.DataFrame(failure_rows)
failure_modes["sort_order"] = failure_modes["fix_stage"].map(STAGE_RANK)
write_table(
    tables_dir / "run-117_stage_failure_modes.md",
    failure_modes.sort_values("sort_order").drop(columns="sort_order"),
    "RUN-117 stage failure modes",
)

copied_names = {path.name.replace("run-117_", "") for path in source_dir.iterdir() if path.is_file()}
collection = manifest.copy()
collection["copied_to_analysis"] = collection["file"].map(lambda name: name in copied_names)
collection["retention_note"] = collection["file"].map(
    lambda name: "left on remote runtime volume; aggregate summaries copied"
    if name in {"walsh_terms.csv", "gft_coefficients.csv", "PHOT_CHLRE_Chen_2023.alignment.fasta"}
    else "copied"
)
write_table(
    tables_dir / "run-117_output_collection_audit.md",
    collection[["file", "size_bytes", "sha256", "copied_to_analysis", "retention_note"]],
    "RUN-117 output collection audit",
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
    tables_dir / "run-117_graphlet_and_module_catalog_summary.md",
    catalog_summary,
    "RUN-117 graphlet and module catalog summary",
)

best_any = comparator_test.iloc[0]
best_nested = comparator_test[comparator_test["candidate_family"] != "empirical_readout_residual"].iloc[0]
best_run117 = equivalence.loc[equivalence["candidate_test_distance"].idxmin()]
metrics = {
    "run_id": RUN_ID,
    "experiment_id": summary["experiment_id"],
    "n_candidate_evaluations": int(summary["n_candidate_evaluations"]),
    "n_selected_candidates": int(summary["n_selected_candidates"]),
    "n_equivalence_rows": int(summary["n_equivalence_rows"]),
    "n_operationally_equivalent": int(summary["n_operationally_equivalent"]),
    "best_any_run098_test_distance": float(best_any["distance_to_empirical_primary"]),
    "best_any_run098_candidate_family": best_any["candidate_family"],
    "best_nested_run098_test_distance": float(best_nested["distance_to_empirical_primary"]),
    "best_nested_run098_candidate_family": best_nested["candidate_family"],
    "best_run117_stage": best_run117["fix_stage"],
    "best_run117_test_distance": float(best_run117["candidate_test_distance"]),
    "best_run117_parameter_set_id": best_run117["parameter_set_id"],
    "stage_summary": json.loads(stage_summary.to_json(orient="records")),
    "copied_source_payload_bytes": int(sum(path.stat().st_size for path in source_dir.iterdir() if path.is_file())),
}
(tables_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

print(f"Wrote ANA-044 derived tables for {RUN_ID} to {tables_dir}")
