#!/usr/bin/env python3
"""Independent engineering checks for the matched full-pairwise TrpB repeat.

This reads saved artifacts without importing the runner or model. Passing means
reproducibility/design checks passed; it is not statistical model validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import numbers
from pathlib import Path
import time

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
RTOL, ATOL, ADDITIVE_ATOL = 1e-10, 1e-12, 1e-14
RATIOS = (0., .5, 1., 2., 4., 8.)
MODES = ("none", "independent_difference", "shared_stop_sensitivity")
RESOLUTIONS = (0., .001, .01)
KEYS = ["protein", "seed", "cohort", "ratio", "observation", "layer",
        "relative_resolution"]
MODEL = dict(interaction_backend="full_pairwise", additive_prior="profile",
             representation_rank=None, orders=[2], order_weights=[1.0],
             target_degrees=[6.0], reference_rank=.5)
CONFIG_CHANGES = {"protocol", "interpretation", "model", "profile_control_run",
                  "profile_control_bundle", "profile_control_manifest_sha256"}


class ValidationError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise ValidationError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def check_files(directory, files):
    require(bool(files), f"Empty files manifest: {directory}")
    directory = Path(directory).resolve()
    for relative, entry in files.items():
        path = (directory / relative).resolve()
        require(path.is_relative_to(directory), f"Manifest path escapes root: {relative}")
        expected = entry["sha256"] if isinstance(entry, dict) else entry
        require(path.is_file(), f"Missing manifest member: {path}")
        require(sha(path) == expected, f"Checksum mismatch: {path}")
    return dict(root=str(directory), files=len(files))


def manifested_root(manifest, candidates):
    """Prefer a complete frozen snapshot; never hide its corruption by fallback."""
    for candidate in candidates:
        candidate = Path(candidate)
        if candidate.is_dir() and all((candidate / name).is_file()
                                      for name in manifest["files"]):
            return check_files(candidate, manifest["files"])
    raise ValidationError(f"No complete payload for manifest among {list(map(str, candidates))}")


def optional_archival_sources(manifest, candidates):
    for candidate in candidates:
        if Path(candidate).is_dir() and all((Path(candidate)/name).is_file()
                                            for name in manifest["files"]):
            return dict(available=True, **check_files(candidate, manifest["files"]))
    return dict(available=False, required=False,
                reason="Original source payload is not included here; the copied manifest and all control outputs are pinned. Restore the original RUN bundle for source replay.")


def check_array(actual, expected, label, *, rtol=RTOL, atol=ATOL):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape, f"Shape mismatch: {label}")
    if actual.dtype.kind in "buifc" and expected.dtype.kind in "buifc":
        require(np.allclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True),
                f"Numeric mismatch: {label} (rtol={rtol}, atol={atol})")
        finite = np.isfinite(actual) & np.isfinite(expected)
        difference = np.abs(actual[finite].astype(float) - expected[finite].astype(float))
        maximum = float(difference.max()) if difference.size else 0.
        exact = bool(np.array_equal(actual, expected, equal_nan=True))
    else:
        require(np.array_equal(actual, expected), f"Categorical mismatch: {label}")
        maximum, exact = 0., True
    return dict(maximum_absolute_difference=maximum, array_equal=exact,
                shape=list(actual.shape))


def compare_frames(actual, expected, label, keys=None, *, rtol=RTOL, atol=ATOL):
    require(set(actual.columns) == set(expected.columns), f"Column mismatch: {label}")
    if keys:
        require(not actual.duplicated(keys).any() and not expected.duplicated(keys).any(),
                f"Duplicate keys: {label}")
        actual = actual.sort_values(keys).reset_index(drop=True)
        expected = expected.sort_values(keys).reset_index(drop=True)
    require(len(actual) == len(expected), f"Row-count mismatch: {label}")
    columns = {}
    for column in expected:
        left, right = actual[column], expected[column]
        if pd.api.types.is_numeric_dtype(left) and pd.api.types.is_numeric_dtype(right):
            columns[column] = check_array(left, right, f"{label}/{column}",
                                          rtol=rtol, atol=atol)
        else:
            require(left.fillna("<NA>").tolist() == right.fillna("<NA>").tolist(),
                    f"Categorical mismatch: {label}/{column}")
    return dict(rows=len(actual), numeric_columns=columns)


def numeric_leaves(value, prefix=""):
    """Keep all numerical fit outputs, excluding only timing/provenance strings."""
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in {"elapsed_seconds", "wall_seconds", "runtime_seconds"}:
                continue
            result.update(numeric_leaves(item, f"{prefix}/{key}"))
        return result
    if isinstance(value, list):
        result = {}
        for index, item in enumerate(value):
            result.update(numeric_leaves(item, f"{prefix}/{index}"))
        return result
    return {prefix: value} if isinstance(value, numbers.Real) else {}


def compare_numeric_json(actual, expected, label):
    left, right = numeric_leaves(actual), numeric_leaves(expected)
    require(set(left) == set(right), f"Numerical JSON schema mismatch: {label}")
    names = sorted(left)
    return dict(fields=names, **check_array([left[k] for k in names],
                [right[k] for k in names], label))


def validate_config(config, control):
    current = {k: v for k, v in config.items() if k not in CONFIG_CHANGES}
    historical = {k: v for k, v in control.items() if k not in CONFIG_CHANGES}
    changed = [k for k in set(current) | set(historical)
               if current.get(k) != historical.get(k)]
    require(not changed, f"Undeclared configuration change: {sorted(changed)}")
    require(config["model"] == MODEL, "Unexpected full-pairwise model configuration")
    require(config["profile_control_run"] == "RUN-178", "Control must be RUN-178")
    require(config["proteins"] == ["TrpB"], "Expected TrpB-only repeat")
    require(tuple(config["interaction_ratios"]) == RATIOS, "Interaction-ratio grid changed")
    require(tuple(config["observation_modes"]["TrpB"]) == MODES, "Observation grid changed")
    require(tuple(config["measurement"]["relative_resolutions"]) == RESOLUTIONS,
            "Resolution grid changed")
    return dict(unchanged_keys=sorted(current), declared_changed_keys=sorted(CONFIG_CHANGES))


def coverage(frame, seeds, *, latent=False):
    require(set(KEYS).issubset(frame.columns), "Missing summary identity columns")
    require(not frame.duplicated(KEYS).any(), "Duplicate summary design cell")
    expected = set()
    for seed, cohort in seeds:
        for ratio in RATIOS:
            if latent:
                expected.add(("TrpB", seed, cohort, ratio, "none", "latent", 0.))
            else:
                for mode in MODES:
                    for layer in (["intrinsic"] if mode == "none" else ["intrinsic", "observed"]):
                        for epsilon in RESOLUTIONS:
                            expected.add(("TrpB", seed, cohort, ratio, mode, layer, epsilon))
    if not latent:
        expected.update(("TrpB", -1, "empirical", -1., "measured", "empirical", r)
                        for r in RESOLUTIONS)
    observed = set(frame[KEYS].itertuples(index=False, name=None))
    require(observed == expected,
            f"Summary design coverage differs: missing={len(expected-observed)}, extra={len(observed-expected)}")
    return dict(rows=len(frame), expected_rows=len(expected), duplicate_cells=0)


def check_settings(settings, control_settings, config):
    mapping = {(s["observation"], float(s["ratio"])): s for s in settings}
    expected = {(mode, ratio) for mode in MODES for ratio in RATIOS}
    require(len(settings) == 18 and set(mapping) == expected, "Selected-setting grid incomplete/duplicated")
    bounds = np.asarray([config["fit"]["offset_bounds"],
                         np.log(config["fit"]["scale_bounds"]),
                         config["fit"]["position_bounds"],
                         np.log(config["fit"]["steepness_bounds"])])
    for key, setting in mapping.items():
        parameters = np.asarray(setting["parameters"])
        require(parameters.shape == (4,) and np.isfinite(parameters).all(), f"Invalid fit parameters: {key}")
        check_array(setting["bounds"], bounds, f"fit bounds {key}", rtol=0, atol=1e-14)
        require(np.all(parameters >= bounds[:, 0]-ATOL)
                and np.all(parameters <= bounds[:, 1]+ATOL), f"Fit outside declared bounds: {key}")
        objective = setting["objective"]
        require(set(objective) == {"quantile", "mean_sd", "floor_cdf", "total"}, "Objective schema changed")
        require(all(np.isfinite(v) and v >= 0 for v in objective.values()), f"Invalid objective: {key}")
        check_array(objective["total"], sum(objective[k] for k in ("quantile", "mean_sd", "floor_cdf")),
                    f"Objective decomposition {key}")
    controls = {(s["observation"], float(s["ratio"])): s for s in control_settings}
    zero = {mode: compare_numeric_json(mapping[mode, 0.], controls[mode, 0.],
                                       f"r0 selected parameters/{mode}") for mode in MODES}
    return dict(settings=18, all_parameters_within_bounds=True, r0=zero)


def check_npz(actual_path, expected_path, label, *, keys=None, rtol=0, atol=0):
    with np.load(actual_path, allow_pickle=False) as actual, np.load(expected_path, allow_pickle=False) as expected:
        if keys is None:
            require(set(actual.files) == set(expected.files), f"NPZ key mismatch: {label}")
            keys = expected.files
        return {key: check_array(actual[key], expected[key], f"{label}/{key}", rtol=rtol, atol=atol)
                for key in keys}


def _compare_scientific_payloads(output, control, config, report):
    """Payload-only helper also permits RUN-178 self-checks; production mode is separate."""
    audit, old_audit = (pd.read_csv(path / "seed_audit.csv") for path in (output, control))
    require(len(audit) == 24 and not audit.seed.duplicated().any(), "Expected 24 distinct confirmation seeds")
    require(audit.cohort.value_counts().to_dict() == {"reused": 16, "fresh": 8}, "Seed cohorts differ from 16+8")
    columns = ["protein", "seed", "cohort", "noise_seed"]
    report["seed_design"] = compare_frames(audit[columns], old_audit[columns], "seed identity", ["seed"], rtol=0, atol=0)
    seeds = [(int(row.seed), row.cohort) for row in old_audit.itertuples()]
    for name in ("latent_fields", "cube_values"):
        expected = {f"TrpB_{seed}.npz" for seed, _ in seeds}
        require({p.name for p in (output / name).glob("*.npz")} == expected,
                f"Seed payload coverage differs: {name}")
    old_refs = {int(row.seed): json.loads(row.reference) for row in old_audit.itertuples()}
    for row in audit.itertuples():
        require(json.loads(row.reference) == old_refs[int(row.seed)], f"Changed reference: {row.seed}")
    report["confirmation_references_equal"] = True
    for name in ("TrpB_empirical.npz", "TrpB_empirical_cubes.npz"):
        report[name] = check_npz(output/name, control/name, name)
    report["additive_matrix"] = check_array(np.load(output/"additive_matrix.npy", allow_pickle=False),
        np.load(control/"additive_matrix.npy", allow_pickle=False), "frozen additive matrix", rtol=0, atol=0)
    report["coverage"] = compare_frames(pd.read_csv(output/"coverage.csv"),
        pd.read_csv(control/"coverage.csv"), "coverage metadata", ["protein"], rtol=0, atol=0)
    report["summaries"] = {}
    for name in ("graph_summary.csv", "cube_summary.csv", "latent_graph_summary.csv"):
        actual, old = (pd.read_csv(path/name) for path in (output, control))
        result = coverage(actual, seeds, latent=name.startswith("latent"))
        result["r0"] = compare_frames(actual[actual.ratio.eq(0)], old[old.ratio.eq(0)], f"r0/{name}", KEYS)
        if not name.startswith("latent"):
            result["empirical"] = compare_frames(actual[actual.cohort.eq("empirical")],
                old[old.cohort.eq("empirical")], f"empirical/{name}", KEYS)
        report["summaries"][name] = result
    report["calibration_arrays"] = check_npz(output/"latent_calibration.npz", control/"latent_calibration.npz",
        "calibration", keys=["a", "genotypes"], rtol=0, atol=ADDITIVE_ATOL)
    indices = [pd.read_csv(path/"latent_calibration_index.csv").drop(columns="fingerprint")
               for path in (output, control)]
    report["calibration_index"] = compare_frames(*indices, "calibration index", ["seed", "cube"], rtol=0, atol=0)
    cal_audit, old_cal_audit = [read_json(path/"calibration_prior_audit.json") for path in (output, control)]
    current = {item["seed"]: item for item in cal_audit}
    previous = {item["seed"]: item for item in old_cal_audit}
    require(len(current) == len(cal_audit) == len(previous), "Calibration audit seed coverage differs")
    require(set(current) == set(previous), "Calibration seed identities changed")
    require(not set(current) & {seed for seed, _ in seeds}, "Calibration/confirmation seeds overlap")
    for seed in current:
        require(current[seed]["reference"] == previous[seed]["reference"], f"Calibration reference changed: {seed}")
        require(current[seed]["cubes"] == previous[seed]["cubes"], f"Calibration cube count changed: {seed}")
        check_array(current[seed]["sigma_a"], previous[seed]["sigma_a"], "calibration sigma_a", rtol=0, atol=ADDITIVE_ATOL)
    report["calibration_references_equal"] = True
    report["confirmation_fields"] = []
    for seed, cohort in seeds:
        name = f"TrpB_{seed}.npz"
        with np.load(output/"latent_fields"/name) as actual, np.load(control/"latent_fields"/name) as old:
            require(set(actual.files) == set(old.files) == {"a", "g"}, f"Field schema changed: {seed}")
            require(actual["a"].shape == actual["g"].shape == (160000,), f"Full field coverage missing: {seed}")
            additive = check_array(actual["a"], old["a"], f"additive/{seed}", rtol=0, atol=ADDITIVE_ATOL)
            require(np.isfinite(actual["g"]).all(), f"Nonfinite interaction field: {seed}")
            difference = float(np.max(np.abs(actual["g"] - old["g"])))
        with np.load(output/"cube_values"/name) as actual, np.load(control/"cube_values"/name) as old:
            require(set(actual.files) == set(old.files), f"Saved cube key mismatch: {seed}")
            check_array(actual["genotypes"], old["genotypes"], f"cube geometry/{seed}", rtol=0, atol=0)
            zero = {key: check_array(actual[key], old[key], f"r0 cube values/{seed}/{key}")
                    for key in old.files if key.endswith("_r0")}
            require(len(zero) == 5, f"Missing r0 observation/layer cube values: {seed}")
        report["confirmation_fields"].append(dict(seed=seed, cohort=cohort, additive=additive,
            maximum_interaction_difference=difference, r0_cube_values=zero))
    report["interaction_change_is_descriptive_only"] = True
    report["settings"] = check_settings(read_json(output/"selected_parameters.json"),
        read_json(control/"selected_parameters.json"), config)
    report["noise"] = {}
    names = {p.name for p in (control/"noise").iterdir() if p.is_file()}
    require({p.name for p in (output/"noise").iterdir() if p.is_file()} == names, "Noise artifact set changed")
    for name in sorted(names):
        actual, old = output/"noise"/name, control/"noise"/name
        result = dict(byte_identical=sha(actual) == sha(old))
        if name.endswith(".json"):
            require(read_json(actual) == read_json(old), f"Noise model changed: {name}")
        elif ".csv" in name:
            result["semantic_comparison"] = compare_frames(pd.read_csv(actual), pd.read_csv(old),
                f"noise/{name}", rtol=0, atol=ADDITIVE_ATOL)
        else:
            require(result["byte_identical"], f"Unknown changed noise artifact: {name}")
        report["noise"][name] = result


def validate_run(output, control, report):
    output, control = Path(output).resolve(), Path(control).resolve()
    manifests = {"candidate": read_json(output/"manifest.json"), "control": read_json(control/"manifest.json")}
    candidate, old = manifests["candidate"], manifests["control"]
    require(candidate.get("mode") == "MATCHED_PAIRWISE_REPEAT", "Candidate mode must be MATCHED_PAIRWISE_REPEAT")
    require(candidate.get("public_model") == "thermodynamic_pairwise_v3", "Candidate model ID mismatch")
    require(candidate.get("interaction_backend") == "full_pairwise", "Candidate backend mismatch")
    require(candidate.get("source_run") == "RUN-178", "Candidate source_run must be RUN-178")
    require(old.get("mode") == "MATCHED_PRIOR_REPEAT" and old.get("public_model") == "thermodynamic_gnk_v2",
            "Control is not the historical profile repeat")
    for key, expected in dict(settings=18, seed_audits=24, graph_rows=2163,
                              cube_summary_rows=2163, latent_graph_rows=144).items():
        require(candidate.get(key) == old.get(key) == expected, f"Manifest count mismatch: {key}")
    config, old_config = [yaml.safe_load((path/"config.yaml").read_text()) for path in (output, control)]
    report["configuration"] = validate_config(config, old_config)
    require(candidate.get("model") == config["model"], "Manifest/config model disagreement")
    require(candidate["protocol"] == config["protocol"], "Manifest/config protocol disagreement")
    require(candidate["additive_prior"] == config["additive_prior"], "Manifest/config prior disagreement")
    require(candidate["reference_policy"] == old["reference_policy"], "Reference policy changed")
    control_hash = sha(control/"manifest.json")
    require(candidate["profile_control_manifest_sha256"] == config["profile_control_manifest_sha256"] == control_hash,
            "Wrong profile-control manifest digest")
    require(sha(output/"profile_control_manifest.json") == control_hash, "Copied profile control manifest differs")
    report["manifests"] = {}
    for label, path, cfg in (("candidate", output, config), ("control", control, old_config)):
        manifest = manifests[label]
        entry = {"outputs": check_files(path, manifest["files"])}
        for field, filename in (("input_manifest_sha256", "input_manifest.json"),
                                ("source_manifest_sha256", "source_manifest.json"),
                                ("config_sha256", "config.yaml")):
            require(manifest[field] == sha(path/filename), f"Manifest digest mismatch: {label}/{filename}")
        source_manifest = read_json(path/"source_manifest.json")
        entry["source"] = (manifested_root(source_manifest, [path.parent/"source", ROOT])
                           if label == "candidate" else
                           optional_archival_sources(source_manifest, [path.parent/"source"]))
        entry["inputs"] = manifested_root(read_json(path/"input_manifest.json"),
            [path.parent/"inputs"/Path(cfg["input_bundle"]).name,
             output.parent/"inputs"/Path(cfg["input_bundle"]).name, ROOT/cfg["input_bundle"]])
        prior = cfg["additive_prior"]
        entry["prior"] = manifested_root(read_json(path/"prior_manifest.json"),
            [path.parent/"inputs"/Path(prior["bundle"]).name,
             output.parent/"inputs"/Path(prior["bundle"]).name, ROOT/prior["bundle"]])
        historical = ROOT/"data/processed/aggregate-thermodynamic-gnk-calibration/RUN-175/source"
        entry["historical_source"] = optional_archival_sources(read_json(path/"historical_source_manifest.json"), [historical])
        for key, filename in (("matrix_sha256", "additive_matrix.npy"), ("metadata_sha256", "prior_manifest.json")):
            require(prior[key] == sha(path/filename), f"Prior digest mismatch: {label}/{filename}")
        environment = read_json(path/"worker_environment.json")
        require(environment == cfg["required_worker_environment"], f"Worker environment mismatch: {label}")
        report["manifests"][label] = entry
    for filename in ("input_manifest.json", "prior_manifest.json", "worker_environment.json"):
        require(read_json(output/filename) == read_json(control/filename), f"Frozen metadata changed: {filename}")
    report["worker_environment"] = read_json(output/"worker_environment.json")
    _compare_scientific_payloads(output, control, config, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--control-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    report = dict(schema="exp065_pairwise_engineering_validation_v1", passed=False,
                  scope="Engineering reproducibility only; no statistical adequacy or superiority claim",
                  output_dir=str(args.output_dir.resolve()), control_dir=str(args.control_dir.resolve()),
                  tolerances=dict(rtol=RTOL, atol=ATOL, additive_atol=ADDITIVE_ATOL),
                  validator_sha256=sha(Path(__file__)))
    try:
        validate_run(args.output_dir, args.control_dir, report)
        report["passed"] = True
    except (ValidationError, KeyError, ValueError, OSError, TypeError) as error:
        report["error"] = f"{type(error).__name__}: {error}"
    report["elapsed_seconds"] = time.perf_counter() - started
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(passed=report["passed"], report=str(args.report),
                         elapsed_seconds=report["elapsed_seconds"], error=report.get("error"))))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
