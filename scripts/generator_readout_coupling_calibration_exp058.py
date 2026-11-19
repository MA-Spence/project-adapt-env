#!/usr/bin/env python3
"""EXP-058 generator readout-coupling calibration.

This runner implements the single operator-free calibration experiment
registered from scratch/RES-056_assessment_readout_coupling_and_next_steps.md.
It reuses the EXP-057 support-aware graphlet machinery, but evaluates the three
registered readout-coupling fixes in one grid:

- A: non-saturated, re-centred single stability gate;
- B: mixed stability/score routing for GNK;
- C: two-channel product readout with GNK routed through function capacity.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
SCRIPT_DIR = PROJECT_ROOT / "scripts"
EXTERNAL_ADAPTENV = PROJECT_ROOT / "external" / "Adapt-Env"
for path in (SRC_ROOT, SCRIPT_DIR, EXTERNAL_ADAPTENV):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from phot_coupled_landscape_grid_search import load_config, write_csv, write_json  # noqa: E402
from design_matched_support_aware_probe_exp057 import (  # noqa: E402
    allele_choices,
    choose_support_sites,
    append_checkpoint_payload,
    checkpoint_root,
    copy_resume_checkpoints,
    evaluate_components,
    existing_parameter_checkpoints,
    factorial_sequences,
    finite,
    load_templates,
    robust_metrics,
    sample_anchor_rows,
    stable_uint64,
    support_catalog,
    write_parameter_checkpoint,
)
from hashed_gnk_benchmark_knob_validation import spearman  # noqa: E402

from adaptenv import FitnessLandscapeEnv, LandscapeConfig, PUBLIC_MODEL_SYSTEM  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args()


def parameter_grid(config: dict[str, Any]) -> list[dict[str, Any]]:
    gen = config["generator"]
    blocks = [
        ("A", "fix_a_non_saturated_single_gate", "single_stability_gate"),
        ("B", "fix_b_mixed_score_stability_routing", "mixed_score_stability"),
        ("C", "fix_c_two_channel_product", "two_channel_product"),
    ]
    rows: list[dict[str, Any]] = []
    for fix_label, block_name, architecture in blocks:
        block = config[block_name]
        for support_mode in gen["support_modes"]:
            for seed in gen["seeds"]:
                for steepness in block["threshold_steepness_grid"]:
                    for margin in block["stability_margin_grid"]:
                        for weight in block["gnk_epistasis_weight_grid"]:
                            for lambda_score in block.get("lambda_score_grid", [0.0]):
                                parameter_set_id = (
                                    f"fix{fix_label}__{support_mode}"
                                    f"__ts{float(steepness):.3f}"
                                    f"__sm{float(margin):.3f}"
                                    f"__gw{float(weight):.3f}"
                                    f"__ls{float(lambda_score):.3f}"
                                    f"__seed{int(seed)}"
                                )
                                rows.append(
                                    {
                                        "parameter_set_id": parameter_set_id,
                                        "fix_label": fix_label,
                                        "fix_block": block_name,
                                        "architecture": architecture,
                                        "support_mode": str(support_mode),
                                        "seed": int(seed),
                                        "sequence_length": int(gen["sequence_length"]),
                                        "threshold_steepness": float(steepness),
                                        "stability_margin": float(margin),
                                        "gnk_epistasis_weight": float(weight),
                                        "lambda_score": float(lambda_score),
                                        "lambda_stability": float(1.0 - float(lambda_score)),
                                        "gnk_terms": int(gen["gnk_terms"]),
                                        "observed_fitness_combine_mode": str(
                                            block.get("observed_fitness_combine_mode", "stability_gate")
                                        ),
                                        "generic_epistasis_target": str(
                                            block.get("gnk_target", "stability")
                                        ),
                                    }
                                )
    return rows


def build_config(
    config: dict[str, Any],
    parameter: dict[str, Any],
    *,
    gnk_weight: float | None = None,
    generic_target: str | None = None,
    combine_mode: str | None = None,
) -> LandscapeConfig:
    gen = config["generator"]
    order_weights = tuple(float(x) for x in gen.get("gnk_order_weights", [0.35, 0.45, 0.20]))
    mode = str(combine_mode or parameter["observed_fitness_combine_mode"])
    target = str(generic_target or parameter["generic_epistasis_target"])
    kwargs: dict[str, Any] = {
        "model_system": PUBLIC_MODEL_SYSTEM,
        "observed_fitness_combine_mode": mode,
        "L": int(parameter["sequence_length"]),
        "alphabet_size": 20,
        "seed": int(parameter["seed"]),
        "stability_margin": float(parameter["stability_margin"]),
        "stability_scale": float(gen.get("stability_scale", 1.0)),
        "threshold_steepness": float(parameter["threshold_steepness"]),
        "fitness_output_scale": float(gen.get("fitness_output_scale", 10.0)),
        "gnk_epistasis_weight": float(
            parameter["gnk_epistasis_weight"] if gnk_weight is None else gnk_weight
        ),
        "gnk_terms": int(parameter["gnk_terms"]),
        "gnk_min_order": int(gen.get("gnk_min_order", 2)),
        "gnk_max_order": int(gen.get("gnk_max_order", 4)),
        "gnk_order_weights": order_weights,
        "gnk_support_mode": str(parameter["support_mode"]),
        "gnk_seed": int(parameter["seed"]) + int(gen.get("gnk_seed_offset", 10000)),
        "n_epistatic_terms": 0,
        "noise_amplitude": float(gen.get("noise_scale", 0.0)),
        "generic_epistasis_target": target,
    }
    if mode == "product":
        kwargs.update(
            {
                "observed_fitness_terms": (
                    {"source": "stability_gate"},
                    {"source": "function_capacity"},
                ),
                # Keep the function channel finite and broad so GNK can shift
                # capacity without introducing a narrow independent bottleneck.
                "n_functional_dims": 1,
                "functional_sigma_base": 20.0,
                "functional_sigma_anisotropy": 0.0,
                "functional_weight": 1.0,
            }
        )
    return LandscapeConfig(**kwargs)


def build_evaluation_envs(
    config: dict[str, Any],
    parameter: dict[str, Any],
) -> dict[str, FitnessLandscapeEnv]:
    envs = {
        "core": FitnessLandscapeEnv(build_config(config, parameter, gnk_weight=0.0)),
        "support": FitnessLandscapeEnv(build_config(config, parameter, generic_target="stability")),
    }
    weight = float(parameter["gnk_epistasis_weight"])
    if str(parameter["architecture"]) == "mixed_score_stability":
        envs["stability"] = FitnessLandscapeEnv(
            build_config(
                config,
                parameter,
                gnk_weight=weight * float(parameter["lambda_stability"]),
                generic_target="stability",
            )
        )
        envs["gnk"] = FitnessLandscapeEnv(
            build_config(config, parameter, gnk_weight=weight, generic_target="stability")
        )
    else:
        envs["full"] = FitnessLandscapeEnv(build_config(config, parameter))
    return envs


def evaluate_layers(
    config: dict[str, Any],
    parameter: dict[str, Any],
    seqs: np.ndarray,
    envs: dict[str, FitnessLandscapeEnv],
) -> dict[str, np.ndarray]:
    output_scale = float(config["generator"].get("fitness_output_scale", 10.0))
    lambda_score = float(parameter.get("lambda_score", 0.0))

    core_components = evaluate_components(
        envs["core"],
        seqs,
        chunk_size=int(config["probe"].get("eval_chunk_size", 65536)),
    )

    if str(parameter["architecture"]) == "mixed_score_stability":
        stability_components = evaluate_components(
            envs["stability"],
            seqs,
            chunk_size=int(config["probe"].get("eval_chunk_size", 65536)),
        )
        gnk_components = evaluate_components(
            envs["gnk"],
            seqs,
            chunk_size=int(config["probe"].get("eval_chunk_size", 65536)),
        )
        full_values = (
            np.asarray(stability_components["fitness"], dtype=np.float64)
            + output_scale * lambda_score * np.asarray(gnk_components["gnk_epistatic_score"], dtype=np.float64)
        )
        return {
            "full_generator_observable": full_values,
            "additive_core_only_observable": np.asarray(core_components["fitness"], dtype=np.float64),
            "bare_gnk_texture": np.asarray(gnk_components["gnk_epistatic_score"], dtype=np.float64),
            "effective_stability_margin": np.asarray(
                stability_components["effective_stability_margin_used"],
                dtype=np.float64,
            ),
        }

    full_components = evaluate_components(
        envs["full"],
        seqs,
        chunk_size=int(config["probe"].get("eval_chunk_size", 65536)),
    )
    return {
        "full_generator_observable": np.asarray(full_components["fitness"], dtype=np.float64),
        "additive_core_only_observable": np.asarray(core_components["fitness"], dtype=np.float64),
        "bare_gnk_texture": np.asarray(full_components["gnk_epistatic_score"], dtype=np.float64),
        "effective_stability_margin": np.asarray(full_components["effective_stability_margin_used"], dtype=np.float64),
    }


def parameter_tables(
    config: dict[str, Any],
    parameter: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    templates = config["_templates"]
    probe_cfg = config["probe"]
    envs = build_evaluation_envs(config, parameter)
    support_env = envs["support"]
    metric_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    graphlet_rows: list[dict[str, Any]] = []
    anchor_rows_out: list[dict[str, Any]] = []

    terms = support_catalog(support_env, str(parameter["parameter_set_id"]))
    support_rows.extend(terms)
    anchor_rows = sample_anchor_rows(support_env, config, parameter)
    for anchor_row in anchor_rows:
        anchor_rows_out.append(
            {
                key: value
                for key, value in anchor_row.items()
                if key != "anchor_array"
            }
            | {
                "parameter_set_id": parameter["parameter_set_id"],
                "fix_label": parameter["fix_label"],
                "architecture": parameter["architecture"],
                "support_mode": parameter["support_mode"],
                "stability_margin": parameter["stability_margin"],
                "threshold_steepness": parameter["threshold_steepness"],
                "gnk_epistasis_weight": parameter["gnk_epistasis_weight"],
                "lambda_score": parameter["lambda_score"],
            }
        )

    for template in templates:
        dim = int(template["dimension"])
        graphlet_count = int(probe_cfg.get("support_graphlets_per_template", 4))
        if int(template["max_alphabet_cardinality"]) >= 10:
            graphlet_count = int(probe_cfg.get("potts_graphlets_per_template", graphlet_count))
        for anchor_row in anchor_rows:
            anchor = np.asarray(anchor_row["anchor_array"], dtype=np.int64)
            for graphlet_index in range(graphlet_count):
                sites = choose_support_sites(
                    terms,
                    dimension=dim,
                    graphlet_index=graphlet_index,
                    sequence_length=int(parameter["sequence_length"]),
                    seed=int(parameter["seed"]) + int(probe_cfg.get("site_seed_offset", 51000)),
                )
                choices = allele_choices(
                    anchor,
                    sites,
                    [int(x) for x in template["alphabet_cardinalities"]],
                    alphabet_size=int(support_env.cfg.alphabet_size),
                    seed=int(parameter["seed"]) + graphlet_index,
                )
                seqs, shape = factorial_sequences(anchor, sites, choices)
                layers = evaluate_layers(config, parameter, seqs, envs)
                graphlet_id = (
                    f"{parameter['parameter_set_id']}__{template['dataset_id']}"
                    f"__{anchor_row['anchor_regime']}__g{graphlet_index:03d}"
                )
                graphlet_rows.append(
                    {
                        "graphlet_id": graphlet_id,
                        "parameter_set_id": parameter["parameter_set_id"],
                        "dataset_id": template["dataset_id"],
                        "anchor_regime": anchor_row["anchor_regime"],
                        "graphlet_index": int(graphlet_index),
                        "dimension": dim,
                        "alphabet_cardinalities": ";".join(str(k) for k in shape),
                        "n_nodes": int(seqs.shape[0]),
                        "varied_sites": ";".join(str(pos) for pos in sites),
                        "support_aware": True,
                    }
                )
                for layer, values in layers.items():
                    row = {
                        "source": "generator",
                        "model_family": layer,
                        "graphlet_id": graphlet_id,
                        "dataset_id": template["dataset_id"],
                        "anchor_regime": anchor_row["anchor_regime"],
                        "graphlet_index": int(graphlet_index),
                        "template_dimension": dim,
                        "template_max_alphabet_cardinality": int(template["max_alphabet_cardinality"]),
                        "node_count": int(seqs.shape[0]),
                        "varied_sites": ";".join(str(pos) for pos in sites),
                        **parameter,
                        **robust_metrics(
                            values,
                            shape,
                            config["descriptors"],
                            seed=stable_uint64(graphlet_id, layer, seed=int(parameter["seed"])) % (2**32),
                        ),
                    }
                    metric_rows.append(row)
                    if layer != "full_generator_observable":
                        control_rows.append(row)
    return {
        "metric_rows": metric_rows,
        "control_rows": control_rows,
        "support_rows": support_rows,
        "graphlet_rows": graphlet_rows,
        "anchor_audit_rows": anchor_rows_out,
    }


def worker_parameter_tables(args: tuple[dict[str, Any], dict[str, Any]]) -> tuple[str, dict[str, list[dict[str, Any]]]]:
    config, parameter = args
    return str(parameter["parameter_set_id"]), parameter_tables(config, parameter)


def generator_tables(
    config: dict[str, Any],
    *,
    output_dir: Path,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    parameters = parameter_grid(config)
    copied_checkpoints = copy_resume_checkpoints(output_dir)
    completed_checkpoints = existing_parameter_checkpoints(output_dir)
    metric_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    graphlet_rows: list[dict[str, Any]] = []
    anchor_rows_out: list[dict[str, Any]] = []
    resumed_parameter_sets: list[str] = []
    computed_parameter_sets: list[str] = []
    pending_parameters: list[tuple[int, dict[str, Any]]] = []

    for index, parameter in enumerate(parameters, start=1):
        parameter_set_id = str(parameter["parameter_set_id"])
        if parameter_set_id in completed_checkpoints:
            append_checkpoint_payload(
                completed_checkpoints[parameter_set_id],
                metric_rows=metric_rows,
                control_rows=control_rows,
                support_rows=support_rows,
                graphlet_rows=graphlet_rows,
                anchor_audit_rows=anchor_rows_out,
            )
            resumed_parameter_sets.append(parameter_set_id)
            print(f"[EXP-058] resume {index}/{len(parameters)} {parameter_set_id}", flush=True)
        else:
            pending_parameters.append((index, parameter))

    def finish_parameter(parameter: dict[str, Any], payload: dict[str, list[dict[str, Any]]]) -> None:
        write_parameter_checkpoint(
            output_dir,
            parameter,
            metric_rows=payload["metric_rows"],
            control_rows=payload["control_rows"],
            support_rows=payload["support_rows"],
            graphlet_rows=payload["graphlet_rows"],
            anchor_audit_rows=payload["anchor_audit_rows"],
        )
        append_checkpoint_payload(
            payload,
            metric_rows=metric_rows,
            control_rows=control_rows,
            support_rows=support_rows,
            graphlet_rows=graphlet_rows,
            anchor_audit_rows=anchor_rows_out,
        )
        computed_parameter_sets.append(str(parameter["parameter_set_id"]))

    requested_workers = int(os.environ.get("LABPROJ_WORKERS") or os.environ.get("SLURM_CPUS_PER_TASK") or "1")
    workers = max(1, min(requested_workers, len(pending_parameters) or 1))
    if workers == 1:
        for index, parameter in pending_parameters:
            parameter_set_id = str(parameter["parameter_set_id"])
            print(f"[EXP-058] compute {index}/{len(parameters)} {parameter_set_id}", flush=True)
            finish_parameter(parameter, parameter_tables(config, parameter))
    else:
        print(
            f"[EXP-058] compute {len(pending_parameters)} pending parameter sets with {workers} workers",
            flush=True,
        )
        parameter_by_id = {str(parameter["parameter_set_id"]): parameter for _, parameter in pending_parameters}
        index_by_id = {str(parameter["parameter_set_id"]): index for index, parameter in pending_parameters}
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(worker_parameter_tables, (config, parameter))
                for _, parameter in pending_parameters
            ]
            for future in as_completed(futures):
                parameter_set_id, payload = future.result()
                print(
                    f"[EXP-058] complete {index_by_id[parameter_set_id]}/{len(parameters)} {parameter_set_id}",
                    flush=True,
                )
                finish_parameter(parameter_by_id[parameter_set_id], payload)

    resume_summary = {
        "checkpoint_root": str(checkpoint_root(output_dir)),
        "copied_checkpoint_dirs": copied_checkpoints,
        "n_resumed_parameter_sets": len(resumed_parameter_sets),
        "n_computed_parameter_sets": len(computed_parameter_sets),
        "resumed_parameter_sets": resumed_parameter_sets,
        "computed_parameter_sets": computed_parameter_sets,
        "workers": workers,
    }
    return parameters, metric_rows, support_rows, graphlet_rows, anchor_rows_out, resume_summary


def core_guardrail_checks(config: dict[str, Any], metric_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(metric_rows)
    if frame.empty:
        return []
    core = frame[frame["model_family"] == "additive_core_only_observable"].copy()
    out: list[dict[str, Any]] = []
    scale = float(config["generator"].get("fitness_output_scale", 10.0))
    low_max = float(config["gates"]["core_guardrail"].get("low_plateau_fitness_max", 0.1 * scale))
    high_min = float(config["gates"]["core_guardrail"].get("high_plateau_fitness_min", 0.9 * scale))
    rho_max = float(config["gates"]["core_guardrail"]["beneficial_fraction_peak_decline_rho_max"])
    group_cols = [
        "parameter_set_id",
        "fix_label",
        "architecture",
        "support_mode",
        "threshold_steepness",
        "stability_margin",
        "gnk_epistasis_weight",
        "lambda_score",
    ]
    for keys, group in core.groupby(group_cols, dropna=False):
        record = dict(zip(group_cols, keys))
        low = finite(group.loc[group["anchor_regime"] == "low_plateau", "fitness_mean"].median())
        high = finite(group.loc[group["anchor_regime"] == "high_plateau", "fitness_mean"].median())
        collapsed = group.groupby("anchor_regime", as_index=False).mean(numeric_only=True)
        order = {"low_plateau": -2.0, "fold_transition": 0.0, "high_plateau": 2.0}
        x = [order.get(str(v), math.nan) for v in collapsed["anchor_regime"].tolist()]
        y = [finite(v) for v in collapsed.get("beneficial_fraction", pd.Series(dtype=float)).tolist()]
        rho = spearman(x, y)
        low_pass = bool(np.isfinite(low) and low <= low_max)
        high_pass = bool(np.isfinite(high) and high >= high_min)
        decline_pass = bool(np.isfinite(rho) and rho < rho_max)
        out.append(
            record
            | {
                "low_plateau_fitness_median": low,
                "high_plateau_fitness_median": high,
                "low_plateau_fitness_max": low_max,
                "high_plateau_fitness_min": high_min,
                "beneficial_fraction_peak_decline_rho": rho,
                "beneficial_fraction_peak_decline_rho_max": rho_max,
                "low_plateau_pass": low_pass,
                "high_plateau_pass": high_pass,
                "beneficial_decline_pass": decline_pass,
                "core_guardrail_pass": bool(low_pass and high_pass and decline_pass),
            }
        )
    return out


def gnk_separation_checks(
    config: dict[str, Any],
    metric_rows: list[dict[str, Any]],
    core_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    frame = pd.DataFrame(metric_rows)
    if frame.empty:
        return []
    core_pass = {
        str(row["parameter_set_id"]): bool(row["core_guardrail_pass"])
        for row in core_rows
    }
    gate_cfg = config["gates"]["gnk_separation"]
    metrics = list(gate_cfg["metrics"])
    thresholds = {str(k): float(v) for k, v in dict(gate_cfg["min_abs_delta"]).items()}
    rows: list[dict[str, Any]] = []
    key_cols = [
        "parameter_set_id",
        "fix_label",
        "architecture",
        "support_mode",
        "threshold_steepness",
        "stability_margin",
        "gnk_epistasis_weight",
        "lambda_score",
    ]
    for keys, group in frame.groupby(key_cols, dropna=False):
        record = dict(zip(key_cols, keys))
        metric_passes: list[bool] = []
        for metric in metrics:
            full = finite(group.loc[group["model_family"] == "full_generator_observable", metric].mean())
            core = finite(group.loc[group["model_family"] == "additive_core_only_observable", metric].mean())
            delta = full - core if np.isfinite(full) and np.isfinite(core) else math.nan
            threshold = float(thresholds[metric])
            record[f"{metric}_full_minus_core"] = delta
            record[f"{metric}_min_abs_delta"] = threshold
            record[f"{metric}_separation_pass"] = bool(np.isfinite(delta) and abs(delta) >= threshold)
            metric_passes.append(bool(record[f"{metric}_separation_pass"]))
        guardrail_pass = bool(core_pass.get(str(record["parameter_set_id"]), False))
        record["matched_core_guardrail_pass"] = guardrail_pass
        record["gnk_separation_pass"] = bool(
            guardrail_pass
            and finite(record.get("gnk_epistasis_weight"), 0.0) > 0.0
            and all(metric_passes)
        )
        rows.append(record)
    return rows


def knob_monotonicity_checks(config: dict[str, Any], metric_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    frame = pd.DataFrame(metric_rows)
    if frame.empty:
        return []
    sub = frame[frame["model_family"] == "full_generator_observable"].copy()
    expected = dict(config["gates"]["knob_monotonicity"]["expected_directions"])
    rho_floor = float(config["gates"]["knob_monotonicity"]["abs_rho_floor"])
    rows: list[dict[str, Any]] = []
    group_cols = [
        "fix_label",
        "architecture",
        "support_mode",
        "threshold_steepness",
        "stability_margin",
        "lambda_score",
        "dataset_id",
        "anchor_regime",
    ]
    for keys, group in sub.groupby(group_cols, dropna=False):
        record_base = dict(zip(group_cols, keys))
        collapsed = group.groupby("gnk_epistasis_weight", as_index=False).mean(numeric_only=True)
        for metric, direction in expected.items():
            if metric not in collapsed:
                continue
            x = [finite(v) for v in collapsed["gnk_epistasis_weight"].tolist()]
            y = [finite(v) for v in collapsed[metric].tolist()]
            rho = spearman(x, y)
            sign = 1.0 if str(direction) == "increase" else -1.0
            direction_pass = bool(np.isfinite(rho) and rho * sign >= rho_floor)
            wrong_direction = bool(np.isfinite(rho) and rho * sign < 0.0)
            rows.append(
                record_base
                | {
                    "knob": "gnk_epistasis_weight",
                    "metric": metric,
                    "expected_direction": direction,
                    "spearman_rho": rho,
                    "abs_rho_floor": rho_floor,
                    "direction_pass": direction_pass,
                    "wrong_direction": wrong_direction,
                    "n_weight_levels": int(len(collapsed)),
                }
            )
    return rows


def decision_table(
    config: dict[str, Any],
    core_rows: list[dict[str, Any]],
    separation_rows: list[dict[str, Any]],
    monotonicity_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    decisions = dict(config["decision_rule"])
    stages = [
        ("A", decisions["fix_a_alone_pass"]),
        ("B", decisions["fix_a_needs_b"]),
        ("C", decisions["fix_a_plus_b_fail"]),
    ]
    core = pd.DataFrame(core_rows)
    sep = pd.DataFrame(separation_rows)
    mono = pd.DataFrame(monotonicity_rows)
    rows: list[dict[str, Any]] = []
    for fix_label, success_decision in stages:
        csub = core[core["fix_label"] == fix_label] if not core.empty else pd.DataFrame()
        ssub = sep[sep["fix_label"] == fix_label] if not sep.empty else pd.DataFrame()
        msub = mono[mono["fix_label"] == fix_label] if not mono.empty else pd.DataFrame()
        n_mono = int(len(msub))
        n_mono_pass = int(msub["direction_pass"].sum()) if n_mono and "direction_pass" in msub else 0
        any_wrong = bool(msub["wrong_direction"].any()) if n_mono and "wrong_direction" in msub else False
        mono_pass = bool(n_mono > 0 and n_mono_pass / n_mono > 0.5 and not any_wrong)
        core_pass_ids = set(csub.loc[csub["core_guardrail_pass"], "parameter_set_id"].astype(str)) if not csub.empty else set()
        sep_pass_ids = set(ssub.loc[ssub["gnk_separation_pass"], "parameter_set_id"].astype(str)) if not ssub.empty else set()
        joint_ids = sorted(core_pass_ids & sep_pass_ids)
        all_gates_pass = bool(joint_ids and mono_pass)
        rows.append(
            {
                "fix_label": fix_label,
                "registered_success_decision": success_decision,
                "n_core_guardrail_pass": int(len(core_pass_ids)),
                "n_gnk_separation_pass": int(len(sep_pass_ids)),
                "n_joint_core_and_separation_pass": int(len(joint_ids)),
                "n_monotonicity_checks": n_mono,
                "n_monotonicity_pass": n_mono_pass,
                "monotonicity_majority_no_wrong_direction_pass": mono_pass,
                "all_registered_gates_pass": all_gates_pass,
                "example_passing_parameter_set": joint_ids[0] if joint_ids else "",
            }
        )
    if rows[0]["all_registered_gates_pass"]:
        final = decisions["fix_a_alone_pass"]
    elif rows[1]["all_registered_gates_pass"]:
        final = decisions["fix_a_needs_b"]
    elif rows[2]["all_registered_gates_pass"]:
        final = decisions["fix_a_plus_b_fail"]
    else:
        final = decisions["fix_c_fail"]
    for row in rows:
        row["final_preregistered_decision_if_run_completed"] = final
        row["hard_stop"] = decisions["hard_stop"]
    return rows


def summarize(
    config: dict[str, Any],
    parameters: list[dict[str, Any]],
    metric_rows: list[dict[str, Any]],
    core_rows: list[dict[str, Any]],
    separation_rows: list[dict[str, Any]],
    monotonicity_rows: list[dict[str, Any]],
    decisions: list[dict[str, Any]],
    resume_summary: dict[str, Any],
) -> dict[str, Any]:
    return {
        "experiment_id": config["experiment"]["id"],
        "mode": config["experiment"]["mode"],
        "public_model_system": PUBLIC_MODEL_SYSTEM,
        "n_templates": len(config["_templates"]),
        "n_parameter_sets": len(parameters),
        "n_generator_metric_rows": len(metric_rows),
        "n_core_guardrail_checks": len(core_rows),
        "n_core_guardrail_pass": int(sum(bool(row.get("core_guardrail_pass")) for row in core_rows)),
        "n_gnk_separation_checks": len(separation_rows),
        "n_gnk_separation_pass": int(sum(bool(row.get("gnk_separation_pass")) for row in separation_rows)),
        "n_monotonicity_checks": len(monotonicity_rows),
        "n_monotonicity_pass": int(sum(bool(row.get("direction_pass")) for row in monotonicity_rows)),
        "final_preregistered_decision_if_run_completed": (
            decisions[0]["final_preregistered_decision_if_run_completed"] if decisions else ""
        ),
        "hard_stop": config["decision_rule"]["hard_stop"],
        "resume": resume_summary,
    }


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = load_config(Path(args.config), quick=bool(args.quick))
    config["_templates"] = load_templates(config, project_root)
    write_json(
        output_dir / "preflight.json",
        {
            "experiment_id": config["experiment"]["id"],
            "mode": config["experiment"]["mode"],
            "quick": bool(args.quick),
            "public_model_system": PUBLIC_MODEL_SYSTEM,
            "project_root": str(project_root),
            "scratch_design_source": config["experiment"]["advisory_source"],
            "operator_free": True,
            "uses_exp050_transfer_function": False,
        },
    )
    parameters, metric_rows, support_rows, graphlet_rows, anchor_rows, resume_summary = generator_tables(
        config,
        output_dir=output_dir,
    )
    core_rows = core_guardrail_checks(config, metric_rows)
    separation_rows = gnk_separation_checks(config, metric_rows, core_rows)
    monotonicity_rows = knob_monotonicity_checks(config, metric_rows)
    decisions = decision_table(config, core_rows, separation_rows, monotonicity_rows)

    write_csv(output_dir / "generator_calibration_grid.csv", parameters)
    write_csv(output_dir / "support_catalog.csv", support_rows)
    write_csv(output_dir / "support_aware_graphlet_audit.csv", graphlet_rows)
    write_csv(output_dir / "anchor_regime_audit.csv", anchor_rows)
    write_csv(output_dir / "generator_descriptor_metrics.csv", metric_rows)
    write_csv(output_dir / "core_guardrail_checks.csv", core_rows)
    write_csv(output_dir / "gnk_separation_checks.csv", separation_rows)
    write_csv(output_dir / "knob_monotonicity_checks.csv", monotonicity_rows)
    write_csv(output_dir / "calibration_decision_table.csv", decisions)
    write_json(
        output_dir / "summary.json",
        summarize(
            config,
            parameters,
            metric_rows,
            core_rows,
            separation_rows,
            monotonicity_rows,
            decisions,
            resume_summary,
        ),
    )


if __name__ == "__main__":
    main()
