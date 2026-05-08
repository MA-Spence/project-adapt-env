#!/usr/bin/env python3
"""Run a branch-comparison calibration panel for HYP-001."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import proteingym_latent_observation_calibration as base  # noqa: E402

ADAPT_ENV_ROOT = PROJECT_ROOT / "external" / "Adapt-Env"
if str(ADAPT_ENV_ROOT) not in sys.path:
    sys.path.insert(0, str(ADAPT_ENV_ROOT))

from adaptenv import LandscapeConfig  # noqa: E402
from adaptenv.calibration import calibrate_synthetic_landscape  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def flatten_mapping(prefix: str, payload: dict[str, Any]) -> dict[str, Any]:
    return {f"{prefix}__{key}": value for key, value in payload.items()}


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    config = base.load_config(Path(args.config))
    if args.quick and "quick" in config:
        config = base.merge_dict(config, config["quick"])

    proteingym_cfg = dict(config["proteingym"])
    panel_cfg = dict(config["panel"])
    mmseqs_cfg = dict(config["mmseqs"])
    mavenn_cfg = dict(config["mavenn"])
    base_config = dict(config["adapt_env"]["base_config"])
    branches_cfg = dict(config["branches"])

    reference_csv_path = project_root / proteingym_cfg["reference_csv_path"]
    substitutions_parquet_path = project_root / proteingym_cfg["substitutions_parquet_path"]
    alignment_dir = project_root / mmseqs_cfg["alignment_dir"]
    mmseqs_cfg["cache_dir"] = str(project_root / mmseqs_cfg["cache_dir"])

    base.download_file(proteingym_cfg["reference_csv_url"], reference_csv_path)
    base.download_file(proteingym_cfg["substitutions_parquet_url"], substitutions_parquet_path)

    reference_df = pd.read_csv(reference_csv_path)
    panel_df = base.select_assay_panel(reference_df, panel_cfg)
    panel_ids = panel_df["DMS_id"].tolist()

    assay_df = pd.read_parquet(
        substitutions_parquet_path,
        columns=["DMS_id", "mutant", "mutated_sequence", "DMS_score"],
        filters=[("DMS_id", "in", panel_ids)],
    )
    assay_df["mutated_sequence"] = assay_df["mutated_sequence"].map(base.canonical_sequence)

    raw_assay_dir = project_root / proteingym_cfg["assay_output_dir"]
    raw_assay_dir.mkdir(parents=True, exist_ok=True)

    shared_raw_landscapes = []
    shared_latent_landscapes = []
    alignment_paths: list[Path] = []
    wildtypes: list[str] = []
    per_assay_payload: list[dict[str, Any]] = []

    calibration_max_mutation_count = int(config["calibration_max_mutation_count"])

    for index, row in enumerate(panel_df.itertuples(index=False)):
        dms_id = str(row.DMS_id)
        wildtype_sequence = base.canonical_sequence(str(row.target_seq))
        subset = assay_df[assay_df["DMS_id"] == dms_id].copy()
        if subset.empty:
            raise ValueError(f"No ProteinGym rows found for selected assay {dms_id}")
        subset = subset[subset["mutated_sequence"].map(len) == len(wildtype_sequence)].copy()
        subset["mutation_count"] = base.mutation_counts(
            wildtype_sequence,
            subset["mutated_sequence"],
        )
        subset = subset[subset["mutation_count"] >= 1].copy()
        if subset.empty:
            raise ValueError(f"No valid ProteinGym rows remained for selected assay {dms_id}")
        singles = subset["mutation_count"] == 1
        subset.loc[singles, "mutation_position"] = base.mutation_positions(
            wildtype_sequence,
            subset.loc[singles, "mutated_sequence"],
        )

        assay_output_path = raw_assay_dir / f"{dms_id}.csv"
        subset.to_csv(assay_output_path, index=False)

        alignment_path = alignment_dir / f"{dms_id}.fasta"
        base.build_alignment_fasta(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            fasta_path=alignment_path,
            mmseqs_cfg=mmseqs_cfg,
        )

        mavenn_model, mavenn_metrics = base.fit_mavenn_model(
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

        calibration_subset = subset[
            subset["mutation_count"] <= calibration_max_mutation_count
        ].copy()
        if calibration_subset.empty:
            raise ValueError(f"No calibration subset remained for assay {dms_id}")

        raw_landscape = base.build_empirical_landscape(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            wildtype_score=float(mavenn_metrics["wildtype_yhat"]),
            assay_frame=calibration_subset,
            score_column="DMS_score",
        )
        latent_landscape = base.build_empirical_landscape(
            dms_id=dms_id,
            wildtype_sequence=wildtype_sequence,
            wildtype_score=float(mavenn_metrics["wildtype_phi"]),
            assay_frame=calibration_subset,
            score_column="latent_phi",
        )

        shared_raw_landscapes.append(raw_landscape)
        shared_latent_landscapes.append(latent_landscape)
        alignment_paths.append(alignment_path)
        wildtypes.append(wildtype_sequence)
        per_assay_payload.append(
            {
                "dms_id": dms_id,
                "uniprot_id": str(row.UniProt_ID),
                "taxon": str(row.taxon),
                "selection_type": str(row.selection_type),
                "sequence_length": int(row.seq_len),
                "n_single_mutants_reference": int(row.DMS_number_single_mutants),
                "n_multiple_mutants_reference": int(row.DMS_number_multiple_mutants),
                "n_variants_loaded": int(len(subset)),
                "n_variants_for_calibration": int(len(calibration_subset)),
                "alignment_path": str(alignment_path),
                "assay_csv_path": str(assay_output_path),
                "mavenn": base.to_builtin(mavenn_metrics),
                "raw_score_summary": base.assay_score_summary(
                    subset["DMS_score"].to_numpy(dtype=np.float64)
                ),
                "latent_phi_summary": base.assay_score_summary(
                    subset["latent_phi"].to_numpy(dtype=np.float64)
                ),
            }
        )

    landscapes_by_source = {
        "raw": shared_raw_landscapes,
        "latent": shared_latent_landscapes,
    }

    shared_branch_rows: list[dict[str, Any]] = []
    per_assay_branch_rows: list[dict[str, Any]] = []
    branch_summary: dict[str, Any] = {}

    for branch_name, branch_cfg in branches_cfg.items():
        source = str(branch_cfg["source"])
        fit_mode = str(branch_cfg["fit_mode"])
        if source not in landscapes_by_source:
            raise ValueError(f"Unsupported branch source: {source}")
        options = base.calibration_options_from_config(branch_cfg["calibration"])

        if fit_mode == "shared":
            result = calibrate_synthetic_landscape(
                functional_landscapes=landscapes_by_source[source],
                functional_alignment_profiles=alignment_paths,
                functional_wildtypes=wildtypes,
                base_config=LandscapeConfig(**base_config),
                options=options,
            )
            branch_summary[branch_name] = {
                "fit_mode": fit_mode,
                "source": source,
                **base.summarize_branch(branch_name=branch_name, calibration_result=result),
            }
            shared_branch_rows.append(
                {
                    "branch": branch_name,
                    "fit_mode": fit_mode,
                    "source": source,
                    **flatten_mapping("validation", base.to_builtin(result.validation)),
                    **flatten_mapping("fitted", base.to_builtin(result.fitted_parameters)),
                    **flatten_mapping("objective", base.to_builtin(result.objective_terms)),
                }
            )
            continue

        if fit_mode != "per_assay":
            raise ValueError(f"Unsupported fit_mode: {fit_mode}")

        branch_rows = []
        for payload, landscape, alignment_path, wildtype_sequence in zip(
            per_assay_payload,
            landscapes_by_source[source],
            alignment_paths,
            wildtypes,
        ):
            result = calibrate_synthetic_landscape(
                functional_landscapes=[landscape],
                functional_alignment_profiles=[alignment_path],
                functional_wildtypes=[wildtype_sequence],
                base_config=LandscapeConfig(**base_config),
                options=options,
            )
            row = {
                "branch": branch_name,
                "fit_mode": fit_mode,
                "source": source,
                "dms_id": payload["dms_id"],
                **flatten_mapping("validation", base.to_builtin(result.validation)),
                **flatten_mapping("fitted", base.to_builtin(result.fitted_parameters)),
                **flatten_mapping("objective", base.to_builtin(result.objective_terms)),
            }
            branch_rows.append(row)
            per_assay_branch_rows.append(row)
        branch_summary[branch_name] = {
            "fit_mode": fit_mode,
            "source": source,
            "rows": branch_rows,
        }

    panel_csv_path = output_dir / "selected_panel.csv"
    panel_df.to_csv(panel_csv_path, index=False)

    mavenn_metrics_path = output_dir / "mavenn_assay_metrics.csv"
    pd.DataFrame(
        [
            {
                "dms_id": payload["dms_id"],
                "taxon": payload["taxon"],
                "sequence_length": payload["sequence_length"],
                "fit_variant_count": payload["mavenn"]["fit_variant_count"],
                "train_variant_count": payload["mavenn"]["train_variant_count"],
                "test_variant_count": payload["mavenn"]["test_variant_count"],
                "wildtype_phi": payload["mavenn"]["wildtype_phi"],
                "wildtype_yhat": payload["mavenn"]["wildtype_yhat"],
                "test_nrmse": payload["mavenn"]["test_metrics"]["nrmse"],
                "test_spearman": payload["mavenn"]["test_metrics"]["spearman"],
                "I_var_test_bits": payload["mavenn"]["I_var_test_bits"],
                "I_pred_test_bits": payload["mavenn"]["I_pred_test_bits"],
            }
            for payload in per_assay_payload
        ]
    ).to_csv(mavenn_metrics_path, index=False)

    shared_branch_path = output_dir / "branch_validations.csv"
    pd.DataFrame(shared_branch_rows).to_csv(shared_branch_path, index=False)

    per_assay_branch_path = output_dir / "per_assay_branch_fits.csv"
    pd.DataFrame(per_assay_branch_rows).to_csv(per_assay_branch_path, index=False)

    summary_payload = {
        "experiment_id": config["experiment"]["id"],
        "proteingym_version": proteingym_cfg["version"],
        "proteingym_package_version": package_version("proteingym"),
        "mavenn_package_version": package_version("mavenn"),
        "reference_csv_path": str(reference_csv_path),
        "substitutions_parquet_path": str(substitutions_parquet_path),
        "alignment_dir": str(alignment_dir),
        "panel_csv_path": str(panel_csv_path),
        "branch_validation_csv_path": str(shared_branch_path),
        "per_assay_branch_fits_csv_path": str(per_assay_branch_path),
        "panel": base.to_builtin(per_assay_payload),
        "branches": base.to_builtin(branch_summary),
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")

    print(json.dumps({"summary_json": str(summary_path)}, indent=2))


if __name__ == "__main__":
    main()
