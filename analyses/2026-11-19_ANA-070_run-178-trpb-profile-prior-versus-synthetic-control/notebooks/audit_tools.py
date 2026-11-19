"""Reproducible numerical/provenance checks; no scientific acceptance gates.

Only archived RUN bundles are read. No optimization, synthetic control run or
whole-graph recomputation is performed. Walsh powers and selected objective
components are recomputed independently; observation draws use the hash-checked
archived noise kernel to preserve the declared measurement model.
"""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import tempfile

import numpy as np
import pandas as pd
import yaml
from scipy.special import expit


ANA = Path(__file__).resolve().parents[1]
ROOT = ANA.parents[1]
TABLES = ANA / "tables"
BASE = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration"
KEYS = ["protein", "seed", "cohort", "ratio", "observation", "layer", "relative_resolution"]
MODES = ["none", "independent_difference", "shared_stop_sensitivity"]
RATIOS = [0, 0.5, 1, 2, 4, 8]
PARAMETERS = ["p0", "p1", "p2", "p3"]
DATASET = "trpb_johnston2024_four_site_graphfla"
AA = "ACDEFGHIKLMNPQRSTVWY"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def verify_hashes(root, files):
    for relative, expected in files.items():
        if isinstance(expected, dict):
            expected = expected["sha256"]
        assert sha(root / relative) == expected, str(root / relative)
    return len(files)


def independent_wht(values):
    """Butterfly transform, independent of the runner's Hadamard matrix product."""
    coefficients = np.asarray(values, dtype=float).copy()
    assert coefficients.ndim == 2 and coefficients.shape[1] == 16
    assert np.isfinite(coefficients).all()
    for stride in (1, 2, 4, 8):
        blocks = coefficients.reshape(len(values), -1, 2 * stride)
        left, right = blocks[:, :, :stride].copy(), blocks[:, :, stride:].copy()
        blocks[:, :, :stride], blocks[:, :, stride:] = left + right, left - right
    coefficients /= 16
    return np.stack([(coefficients[:, [i for i in range(16) if i.bit_count() == k]] ** 2).sum(axis=1)
                     for k in range(1, 5)], axis=1)


def _wht_audit(run, output, seeds, settings, cubes):
    lookup = cubes.set_index(KEYS)
    scale = settings[0]["calibration_scale"]
    records = []

    def check(values, base):
        power = independent_wht(values)
        total = power.sum(axis=1)
        # Parseval also checks normalization without the runner's WHT code.
        np.testing.assert_allclose(total, np.var(values, axis=1), rtol=1e-10, atol=1e-18)
        for resolution in (0.0, 0.001, 0.01):
            threshold = max((resolution * scale) ** 2, 1e-24)
            eligible = total > threshold
            metrics = {
                "epistatic_fraction": float(np.mean(power[eligible, 1:].sum(axis=1) / total[eligible])) if eligible.any() else np.nan,
                "ho_fraction": float(np.mean(power[eligible, 2:].sum(axis=1) / total[eligible])) if eligible.any() else np.nan,
                "power_availability": float(eligible.mean()),
            }
            for k in range(4):
                metrics[f"pooled_power_{k + 1}"] = float(power[:, k].sum() / total.sum()) if total.sum() > len(values) * threshold else np.nan
            metrics["pooled_epistatic_fraction"] = sum(metrics[f"pooled_power_{k}"] for k in (2, 3, 4))
            metrics["pooled_ho_fraction"] = sum(metrics[f"pooled_power_{k}"] for k in (3, 4))
            key = tuple(base[k] for k in KEYS[:-1]) + (resolution,)
            saved = lookup.loc[key]
            errors = []
            for name, value in metrics.items():
                np.testing.assert_allclose(value, saved[name], rtol=1e-7, atol=1e-8, equal_nan=True,
                                           err_msg=f"{run} {key} {name}")
                if np.isfinite(value):
                    errors.append(abs(value - float(saved[name])))
            records.append(dict(run=run, **base, relative_resolution=resolution, n_cubes=len(values),
                                n_power_available=int(eligible.sum()), maximum_absolute_difference=max(errors),
                                **{f"recomputed_{name}": value for name, value in metrics.items()}))

    with np.load(output / "TrpB_empirical_cubes.npz") as empirical:
        check(empirical["fitness"], dict(protein="TrpB", seed=-1, cohort="empirical", ratio=-1.0,
                                        observation="measured", layer="empirical"))
    for row in seeds.itertuples():
        with np.load(output / "cube_values" / f"TrpB_{row.seed}.npz") as values:
            for label in values.files:
                if label == "genotypes":
                    continue
                head, ratio = label.rsplit("_r", 1)
                mode, layer = head.rsplit("_", 1)
                check(values[label], dict(protein="TrpB", seed=int(row.seed), cohort=row.cohort,
                                         ratio=float(ratio), observation=mode, layer=layer))
    assert len(records) == len(cubes) == 2163
    return records


def _noise_banks(source, inputs, old_config):
    path = source / "scripts/replicate_observation_exp065.py"
    spec = importlib.util.spec_from_file_location("ana070_frozen_noise", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory(prefix="ana070-noise-") as temporary:
        noise, assignments = module.calibrate_noise(
            pd.read_csv(inputs / "TrpB/donor_assignments.csv.gz"),
            set(json.loads((inputs / "TrpB/excluded_noise_donors.json").read_text())),
            old_config, Path(temporary))
        frozen = json.loads((inputs / "TrpB/noise_model.json").read_text())
        regenerated = json.loads((Path(temporary) / "noise_model.json").read_text())
        for key in ("depth_breaks", "fitness_breaks", "shared_covariance_raw"):
            np.testing.assert_allclose(regenerated[key], frozen[key], rtol=0, atol=1e-14)
    original_assignments = pd.read_csv(inputs / "TrpB/donor_assignments.csv.gz")
    pd.testing.assert_frame_equal(assignments, original_assignments, check_dtype=False)
    uniforms = np.load(inputs / "TrpB/fit_observation_uniforms.npz")["uniforms"]
    breaks = noise["fitness_breaks"]
    representatives = (breaks[0] - 1, (breaks[0] + breaks[1]) / 2,
                       (breaks[1] + breaks[2]) / 2, breaks[2] + 1)
    banks = {"none": None}
    for mode in MODES[1:]:
        values = []
        for value in representatives:
            first, second, _, _ = module.observe(np.full(uniforms.shape[1:-1], value), noise,
                                                 uniforms, mode == "shared_stop_sensitivity")
            values.append(np.stack([first - value, second - value]))
        banks[mode] = (np.stack(values), breaks)
    return banks


def _objective(setting, additive, interaction, calibration, bank, quantiles):
    parameters = np.asarray(setting["parameters"])
    origin, scale = setting["calibration_origin"], setting["calibration_scale"]
    target = (calibration - origin) / scale
    raw = origin + scale * (parameters[0] + np.exp(parameters[1]) * expit(
        parameters[2] + np.exp(parameters[3]) * (additive + setting["ratio"] * interaction)))
    if bank is None:
        values = (raw - origin) / scale
    else:
        errors, breaks = bank
        bins = np.searchsorted(breaks, raw, side="right")
        selected = np.take_along_axis(errors, np.broadcast_to(bins[None, None, None],
                                      (1, errors.shape[1], errors.shape[2], *bins.shape)), axis=0)[0]
        values = (((raw + selected[0]) + (raw + selected[1])) / 2 - origin) / scale
    probabilities = np.linspace(0, 1, quantiles)
    floor = (setting["near_floor_raw"] - origin) / scale
    parts = dict(quantile=float(np.mean((np.quantile(values, probabilities) - np.quantile(target, probabilities)) ** 2) / target.var()),
                 mean_sd=float((values.std(axis=-1).mean() / target.std(axis=-1).mean() - 1) ** 2),
                 floor_cdf=float(((values <= floor).mean() - (target <= floor).mean()) ** 2))
    parts["total"] = sum(parts.values())
    return parts


def _trace_record(run, path, setting, expected_count, phase):
    trace = pd.read_csv(path)
    assert len(trace) == expected_count == setting["evaluations"], str(path)
    np.testing.assert_array_equal(trace.evaluation, np.arange(len(trace)))
    np.testing.assert_allclose(trace.total, trace[["quantile", "mean_sd", "floor_cdf"]].sum(axis=1),
                               rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(trace.total.min(), setting["objective"]["total"], rtol=1e-9, atol=1e-10)
    parameters = np.asarray(setting["parameters"])
    matches = np.isclose(trace[PARAMETERS].to_numpy(), parameters, rtol=1e-12, atol=1e-12).all(axis=1)
    assert matches.any(), str(path)
    np.testing.assert_allclose(trace.loc[matches, "total"], setting["objective"]["total"], rtol=1e-9, atol=1e-10)
    bounds = np.asarray(setting["bounds"])
    assert np.all(trace[PARAMETERS].to_numpy() >= bounds[:, 0] - 1e-12)
    assert np.all(trace[PARAMETERS].to_numpy() <= bounds[:, 1] + 1e-12)
    return trace, dict(run=run, phase=phase, trace=path.name, evaluations=len(trace),
                       expected_evaluations=expected_count, selected_objective=setting["objective"]["total"],
                       minimum_trace_objective=float(trace.total.min()), optimizer_converged=setting["optimizer_converged"])


def _fit_audit(run, output, inputs, config, old_config, settings, banks):
    calibration_file = output / "latent_calibration.npz" if run == "RUN-178" else inputs / "TrpB/latent_calibration.npz"
    with np.load(calibration_file) as latent:
        additive, interaction = latent["a"], latent["g"]
    with np.load(inputs / "TrpB/empirical_panel.npz") as panel:
        index = pd.read_csv(inputs / "TrpB/empirical_panel_index.csv")
        calibration = panel["fitness"][index.split.eq("calibration")]
    old_settings = json.loads((inputs / "selected_parameters.json").read_text())
    full_budget = (config["fit"]["maxiter"] + 1) * config["fit"]["popsize"] * 4
    assert full_budget == 1952
    dataset_seed = int.from_bytes(hashlib.sha256(DATASET.encode()).digest()[:4], "little") % 1000000000
    fit_records, trace_records, intermediate = [], [], {}
    for setting in settings:
        mode, ratio = setting["observation"], setting["ratio"]
        expected_bounds = [config["fit"]["offset_bounds"], np.log(config["fit"]["scale_bounds"]).tolist(),
                           config["fit"]["position_bounds"], np.log(config["fit"]["steepness_bounds"]).tolist()]
        np.testing.assert_array_equal(setting["bounds"], expected_bounds)
        parts = _objective(setting, additive, interaction, calibration, banks[mode], config["fit"]["quantiles"])
        for name, value in parts.items():
            np.testing.assert_allclose(value, setting["objective"][name], rtol=1e-9, atol=1e-10)
        fit_records.append(dict(run=run, observation=mode, ratio=ratio, objective_recomputed=parts["total"],
                                maximum_component_difference=max(abs(parts[k] - setting["objective"][k]) for k in parts),
                                optimizer_converged=setting["optimizer_converged"], selected_source=setting["selected_source"],
                                **{f"objective_{key}": value for key, value in parts.items()}))
        trace, record = _trace_record(run, output / "fits" / f"TrpB_{mode}_r{ratio:g}.csv.gz", setting,
                                     1 if ratio in (0, 4) else full_budget + 2, "selected")
        trace_records.append(record)
        if ratio in (0, 4):
            if run == "RUN-175":
                old = next(x for x in old_settings if x["dataset_id"] == DATASET and
                           x["observation"] == mode and x["arm"] == f"r{int(ratio)}_sigmoid")
                np.testing.assert_array_equal(setting["parameters"], old["parameters"])
            else:
                endpoint = json.loads((output / "fits" / f"endpoint_parameters_{mode}.json").read_text())[str(int(ratio))]
                np.testing.assert_array_equal(setting["parameters"], endpoint["parameters"])
                expected_seed = dataset_seed + old_config["fit"]["seed"] + MODES.index(mode) * 10 + old_config["model"]["arms"].index(f"r{int(ratio)}_sigmoid")
                assert setting["optimizer_seed"] == endpoint["optimizer_seed"] == expected_seed
                _, record = _trace_record(run, output / "fits" / f"endpoint_{mode}_r{int(ratio)}.csv.gz",
                                          endpoint, full_budget, "endpoint_refit")
                assert setting["endpoint_optimizer_evaluations"] == full_budget
                trace_records.append(record)
        else:
            assert (trace.candidate.iloc[:full_budget] == "differential_evolution").all()
            for offset, endpoint_ratio in enumerate((0, 4)):
                endpoint = next(x for x in settings if x["observation"] == mode and x["ratio"] == endpoint_ratio)
                np.testing.assert_allclose(trace[PARAMETERS].iloc[full_budget + offset], endpoint["parameters"], rtol=1e-12, atol=1e-12)
            if run == "RUN-178":
                assert setting["optimizer_seed"] == dataset_seed + config["fit"]["seed"] + 1000 + MODES.index(mode) * 100 + RATIOS.index(ratio)
            intermediate[(mode, ratio)] = trace[PARAMETERS].iloc[:32].to_numpy()
    return fit_records, trace_records, intermediate


def verify():
    """Write audit tables and return a JSON-compatible engineering summary."""
    TABLES.mkdir(parents=True, exist_ok=True)
    bundles = {name: BASE / name for name in ("RUN-175", "RUN-178")}
    outputs = {name: path / "outputs" for name, path in bundles.items()}
    inputs = {"RUN-175": bundles["RUN-175"] / "inputs",
              "RUN-178": bundles["RUN-178"] / "inputs/interaction-strength-exp065-inputs"}
    configs, settings, seeds, cubes, graphs, provenance = {}, {}, {}, {}, {}, []
    for run in bundles:
        output, bundle = outputs[run], bundles[run]
        manifest = json.loads((output / "manifest.json").read_text())
        source_manifest = json.loads((output / "source_manifest.json").read_text())
        input_manifest = json.loads((inputs[run] / "manifest.json").read_text())
        provenance.append(dict(run=run, output_hashes=verify_hashes(output, manifest["files"]),
                               source_hashes=verify_hashes(bundle / "source", source_manifest["files"]),
                               input_hashes=verify_hashes(inputs[run], input_manifest["files"])))
        assert manifest["mode"] == ("MATCHED_PRIOR_REPEAT" if run == "RUN-178" else "PAIRED_SWEEP")
        for key, path in (("source_manifest_sha256", output / "source_manifest.json"),
                          ("input_manifest_sha256", inputs[run] / "manifest.json"),
                          ("config_sha256", output / "config.yaml")):
            assert manifest[key] == sha(path)
        configs[run] = yaml.safe_load((output / "config.yaml").read_text())
        settings[run] = [x for x in json.loads((output / "selected_parameters.json").read_text()) if x["protein"] == "TrpB"]
        assert len(settings[run]) == 18
        assert {(x["observation"], x["ratio"]) for x in settings[run]} == {(m, r) for m in MODES for r in RATIOS}
        seeds[run] = pd.read_csv(output / "seed_audit.csv").query("protein == 'TrpB'")
        assert seeds[run].groupby("cohort").size().to_dict() == {"fresh": 8, "reused": 16}
        for name, collection in (("cube_summary.csv", cubes), ("graph_summary.csv", graphs)):
            collection[run] = pd.read_csv(output / name).query("protein == 'TrpB'")
            assert len(collection[run]) == 2163 and not collection[run].duplicated(KEYS).any()
        assert json.loads((output / "kernel_validation.json").read_text())["small_graph_cases"] == 20
    assert sha(inputs["RUN-175"] / "manifest.json") == sha(inputs["RUN-178"] / "manifest.json")
    control = json.loads((bundles["RUN-178"] / "source/experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/alignment_prior_control.json").read_text())
    control_hashes = verify_hashes(ROOT, control["control_output_hashes"])
    profile_input = bundles["RUN-178"] / "inputs/trpb-profile-additive-inputs"
    prior = json.loads((profile_input / "prior_manifest.json").read_text())
    profile_hashes = verify_hashes(profile_input, prior["files"])
    assert prior["alphabet"] == AA and prior["assay_positions_1based"] == [183, 184, 227, 228]
    assert prior["biological_parent_at_assay_positions"] == "VFVS"
    for key in ("interaction_ratios", "seed_design", "model", "fit", "sampling", "measurement"):
        assert configs["RUN-175"][key] == configs["RUN-178"][key], key
    assert configs["RUN-175"]["observation_modes"]["TrpB"] == configs["RUN-178"]["observation_modes"]["TrpB"] == MODES
    old_config = yaml.safe_load((inputs["RUN-178"] / "config.yaml").read_text())
    assert old_config["fit"] == configs["RUN-178"]["fit"]
    pd.testing.assert_frame_equal(seeds["RUN-175"][["seed", "cohort", "noise_seed"]].reset_index(drop=True),
                                  seeds["RUN-178"][["seed", "cohort", "noise_seed"]].reset_index(drop=True))
    paired_records = []
    for filename in ("TrpB_empirical.npz", "TrpB_empirical_cubes.npz"):
        with np.load(outputs["RUN-175"] / filename) as before, np.load(outputs["RUN-178"] / filename) as after:
            assert set(before.files) == set(after.files)
            for key in before.files:
                np.testing.assert_array_equal(before[key], after[key])
    empirical = np.load(outputs["RUN-178"] / "TrpB_empirical.npz")
    assert empirical["valid"].sum() == 159129
    variants = pd.read_parquet(inputs["RUN-178"] / "TrpB/variants.parquet")
    codes = np.column_stack([pd.Categorical(variants[f"pos{j + 1}"].astype(str), categories=list(AA)).codes for j in range(4)])
    assert (codes >= 0).all()
    reconstructed_mask = np.zeros(20 ** 4, dtype=bool)
    reconstructed_mask[np.ravel_multi_index(codes.T, (20,) * 4)] = True
    np.testing.assert_array_equal(reconstructed_mask, empirical["valid"])
    for collection in (graphs, cubes):
        a, b = [collection[run].query("cohort == 'empirical'").reset_index(drop=True) for run in bundles]
        pd.testing.assert_frame_equal(a, b, check_exact=True)
        assert set(map(tuple, collection["RUN-175"][KEYS].to_numpy())) == set(map(tuple, collection["RUN-178"][KEYS].to_numpy()))
    for row in seeds["RUN-178"].itertuples():
        filename = f"TrpB_{row.seed}.npz"
        with np.load(outputs["RUN-175"] / "latent_fields" / filename) as a, np.load(outputs["RUN-178"] / "latent_fields" / filename) as b:
            np.testing.assert_allclose(a["g"], b["g"], rtol=1e-12, atol=1e-12)
            paired_records.append(dict(seed=int(row.seed), cohort=row.cohort, nodes=len(a["g"]),
                                       maximum_g_difference=float(abs(a["g"] - b["g"]).max()),
                                       additive_rms_difference=float(np.sqrt(np.mean((a["a"] - b["a"]) ** 2)))))
        with np.load(outputs["RUN-175"] / "cube_values" / filename) as a, np.load(outputs["RUN-178"] / "cube_values" / filename) as b:
            assert set(a.files) == set(b.files)
            np.testing.assert_array_equal(a["genotypes"], b["genotypes"])
    unchanged_methods = []
    def methods(path):
        model = next(node for node in ast.parse(path.read_text()).body if isinstance(node, ast.ClassDef) and node.name == "ThermodynamicLandscape")
        return {node.name: ast.dump(node, include_attributes=False) for node in model.body if isinstance(node, ast.FunctionDef)}
    before, after = [methods(bundles[run] / "source/external/Adapt-Env/adaptenv/thermodynamic.py") for run in bundles]
    for name in before:
        if name != "summary":
            assert before[name] == after[name], name
            unchanged_methods.append(name)
    banks = _noise_banks(bundles["RUN-178"] / "source", inputs["RUN-178"], old_config)
    fit_records, trace_records, wht_records, populations = [], [], [], {}
    for run in bundles:
        fit, trace, populations[run] = _fit_audit(run, outputs[run], inputs[run], configs[run], old_config, settings[run], banks)
        fit_records.extend(fit)
        trace_records.extend(trace)
        wht_records.extend(_wht_audit(run, outputs[run], seeds[run], settings[run], cubes[run]))
    for key in populations["RUN-175"]:
        np.testing.assert_array_equal(populations["RUN-175"][key], populations["RUN-178"][key])
    for filename, records in (("fit_audit.csv", fit_records), ("fit_trace_audit.csv", trace_records),
                               ("wht_audit.csv", wht_records), ("paired_field_audit.csv", paired_records),
                               ("provenance_audit.csv", provenance)):
        pd.DataFrame(records).to_csv(TABLES / filename, index=False)
    result = dict(
        scope="TrpB RUN-178 profile versus archived RUN-175 synthetic control",
        provenance=provenance, frozen_control_hashes_verified=control_hashes,
        profile_input_hashes_verified=profile_hashes, measured_genotypes=159129,
        paired_seeds=24, paired_full_field_nodes=24 * 160000,
        maximum_normalized_g_difference=max(x["maximum_g_difference"] for x in paired_records),
        unchanged_computational_methods=unchanged_methods,
        empirical_graph_and_cube_rows_exact=True, paired_cube_geometries_exact=True,
        selected_objectives_recomputed=len(fit_records), traces_audited=len(trace_records),
        intermediate_initial_populations_exact=len(populations["RUN-175"]),
        independent_wht_summary_rows=len(wht_records),
        maximum_wht_difference=max(x["maximum_absolute_difference"] for x in wht_records),
        maximum_objective_difference=max(x["maximum_component_difference"] for x in fit_records),
        finite_budget_evaluations=1952, scientific_gates=False,
        limitations=[
            "Numerical consistency is not empirical adequacy or global optimizer convergence.",
            "Objective formulas are independent; observation draws reuse the verified archived noise kernel.",
            "Trace minima, components, parameters and budgets are checked; the optimizer is not rerun.",
            "All three WHT resolutions and all saved TrpB conditions are audited; whole-graph dynamic programming is not rerun.",
            "Shared inputs and fresh generator seeds do not create an independent biological test.",
            "RUN-178 uses macOS/Clang and RUN-175 used Linux/Zig; exact numerical package versions do not imply identical platforms.",
        ])
    dump(TABLES / "audit.json", result)
    return result


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
