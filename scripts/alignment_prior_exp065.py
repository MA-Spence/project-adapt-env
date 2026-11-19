#!/usr/bin/env python3
"""Matched TrpB RUN-175 repeat with only an MSA-derived additive prior changed.

The recorded synthetic controls are reused, never recalibrated here. RUN-172's
endpoint fitting and RUN-175's intermediate fitting, observations, graph kernel,
and cube summaries are imported unchanged. Endpoint readouts must be refitted
for the new prior before they can anchor the intermediate-strength searches.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys
import time

import numpy as np
import pandas as pd
from scipy.special import expit
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "external/Adapt-Env"))
import interaction_strength_exp065 as sweep
import production_landscapes_exp066 as production
from adaptenv.amino_acids import AA_ORDER

DATASET = sweep.IDS["TrpB"]
ASSAY_POSITIONS = [183, 184, 227, 228]
REFERENCE_POLICY = (
    "Preserve each seed's original automatically generated reference. It is a "
    "computational centering anchor, not the biological wild type or consensus."
)
REUSED_SOURCE_FILES = (
    "scripts/interaction_strength_exp065.py",
    "scripts/production_landscapes_exp066.py",
    "scripts/observation_calibration_exp065.py",
    "scripts/replicate_observation_exp065.py",
    "scripts/landscape_graph.py",
    "scripts/landscape_graph_kernel.cpp",
    "external/Adapt-Env/adaptenv/extensive_gnk.py",
    "external/Adapt-Env/adaptenv/ranked_categorical_gnk.py",
    "external/Adapt-Env/adaptenv/amino_acids.py",
)


def verify_sources(source_manifest: Path, historical_manifest: Path) -> dict:
    current = json.loads(source_manifest.read_text())
    historical = json.loads(historical_manifest.read_text())
    for relative, expected in current["files"].items():
        if sweep.sha(ROOT / relative) != expected:
            raise ValueError(f"Current source checksum mismatch: {relative}")
    audited = {}
    for relative in REUSED_SOURCE_FILES:
        digest = sweep.sha(ROOT / relative)
        if digest != historical["files"][relative]:
            raise ValueError(f"Frozen scientific kernel changed: {relative}")
        audited[relative] = digest
    return audited


def load_prior(matrix_path: Path, metadata_path: Path, config: dict) -> tuple:
    expected = config["additive_prior"]
    for path, key in [(matrix_path, "matrix_sha256"),
                      (metadata_path, "metadata_sha256")]:
        if sweep.sha(path) != expected[key]:
            raise ValueError(f"Additive-prior checksum mismatch: {path}")
    matrix = np.load(matrix_path, allow_pickle=False)
    if matrix.shape != (4, 20) or not np.isfinite(matrix).all():
        raise ValueError("The matched TrpB prior must be a finite (4, 20) matrix")
    if not np.any(np.ptp(matrix, axis=1) > 0):
        raise ValueError("The additive prior has zero mutation variance")
    metadata = json.loads(metadata_path.read_text())
    if list(metadata["alphabet"]) != list(AA_ORDER):
        raise ValueError("Prior columns must use the public model's AA_ORDER")
    for key in ["assay_positions_1based", "matrix_rows_correspond_to_assay_positions_1based"]:
        if metadata.get(key) != ASSAY_POSITIONS:
            raise ValueError(f"Prior rows must map to ordered TrpB positions {ASSAY_POSITIONS}: {key}")
    matrix = np.asarray(matrix, dtype=np.float64)
    matrix.setflags(write=False)
    return matrix, metadata


def build_model(seed: int, model_config: dict, matrix: np.ndarray, **readout):
    """Keep the original seed/reference stream; replace only the additive matrix."""
    config = sweep.ThermodynamicConfig(
        length=4, seed=int(seed), **model_config, **readout
    )
    return sweep.ThermodynamicLandscape(config, additive_matrix=matrix)


def regenerate_panel(index: pd.DataFrame, frozen: dict, model_config: dict,
                     matrix: np.ndarray) -> tuple:
    """Use every saved seed/cube in its saved row order, auditing unchanged GNK."""
    genotypes = np.asarray(frozen["genotypes"])
    if genotypes.shape != (len(index), 16, 4):
        raise ValueError("Frozen latent panel and its index disagree")
    if index.duplicated(["seed", "cube"]).any():
        raise ValueError("Duplicate seed/cube in frozen latent panel")
    additive = np.empty(genotypes.shape[:2])
    interactions = np.empty_like(additive)
    updated_index = index.copy()
    audits = []
    for seed in index.seed.unique():
        selected = np.flatnonzero(index.seed.eq(seed))
        model = build_model(int(seed), model_config, matrix)
        a, g = model.latent_components(genotypes[selected].reshape(-1, 4))
        additive[selected] = a.reshape(len(selected), 16)
        interactions[selected] = g.reshape(len(selected), 16)
        original_g = np.asarray(frozen["g"])[selected]
        np.testing.assert_allclose(interactions[selected], original_g,
                                   rtol=1e-10, atol=1e-10)
        updated_index.loc[selected, "fingerprint"] = model.fingerprint
        audits.append(dict(
            seed=int(seed), cubes=len(selected), reference=model.reference.tolist(),
            maximum_normalized_g_difference=float(
                np.max(np.abs(interactions[selected] - original_g))
            ), fingerprint=model.fingerprint, sigma_a=float(model.sigma_a),
        ))
    return updated_index, additive, interactions, audits


def fit_settings(ca: np.ndarray, cg: np.ndarray, target: np.ndarray, banks: dict,
                 origin: float, scale: float, floor: float, config: dict,
                 original_config: dict, output: Path, quick: bool = False) -> list:
    """Repeat the endpoint-then-intermediate fitting lineage on calibration only."""
    output.mkdir(parents=True, exist_ok=True)
    modes = config["observation_modes"]["TrpB"]
    ratios = config["interaction_ratios"]
    if config["fit"] != original_config["fit"]:
        raise ValueError("Matched repeat cannot change the fitting configuration")
    original_arms = original_config["model"]["arms"]
    selected = []
    for mode_index, mode in enumerate(modes):
        anchors = {}
        # RUN-175 inherited these endpoint fits from RUN-172. They must be
        # recalibrated under the new prior, with the original RUN-172 seeds.
        for ratio in [0, 4]:
            arm = f"r{ratio}_sigmoid"
            optimizer_seed = (production.seed_for(DATASET)
                              + original_config["fit"]["seed"]
                              + mode_index * 10 + original_arms.index(arm))
            setting, trace = production.fit_map(
                arm, ca, cg, target, banks[mode], origin, scale,
                original_config["fit"], optimizer_seed, quick, floor,
            )
            if setting is None:
                raise ValueError("Unidentifiable endpoint calibration target")
            setting.update(optimizer_seed=optimizer_seed,
                           calibration_lineage="RUN172_endpoint_prior_only_refit")
            anchors[ratio] = setting
            trace.to_csv(output / f"endpoint_{mode}_r{ratio}.csv.gz", index=False)
        sweep.dump(output / f"endpoint_parameters_{mode}.json", anchors)
        for ratio_index, ratio in enumerate(ratios):
            optimizer_seed = (production.seed_for(DATASET) + config["fit"]["seed"]
                              + 1000 + mode_index * 100 + ratio_index)
            setting, trace = sweep.fit_readout(
                ratio, mode, ca, cg, target, banks[mode], origin, scale, floor,
                config["fit"], optimizer_seed, anchors, quick,
            )
            # Rename provenance labels only; the original selection algorithm
            # receives the newly fitted endpoints and is otherwise unchanged.
            labels = {
                "frozen_RUN172_anchor": "matched_alignment_endpoint_refit",
                "RUN172_r0_candidate": "matched_alignment_r0_endpoint_candidate",
                "RUN172_r4_candidate": "matched_alignment_r4_endpoint_candidate",
            }
            trace["candidate"] = trace.candidate.replace(labels)
            setting["selected_source"] = labels.get(
                setting["selected_source"], setting["selected_source"]
            )
            setting.update(
                protein="TrpB", dataset_id=DATASET,
                additive_prior="empirical_msa_log_profile",
                optimizer_seed=(anchors[int(ratio)]["optimizer_seed"]
                                if ratio in [0, 4] else optimizer_seed),
                calibration_lineage=("RUN172_endpoint_prior_only_refit"
                                     if ratio in [0, 4] else
                                     "RUN175_intermediate_with_refitted_endpoints"),
            )
            if ratio in [0, 4]:
                setting["endpoint_optimizer_evaluations"] = anchors[int(ratio)]["evaluations"]
            trace.to_csv(output / f"TrpB_{mode}_r{ratio:g}.csv.gz", index=False)
            selected.append(setting)
            print("FIT", mode, ratio, setting["objective"]["total"], flush=True)
    return selected


def validate_config(config: dict, historical_config: dict):
    if config["proteins"] != ["TrpB"]:
        raise ValueError("This matched intervention is restricted to TrpB")
    for key in ["interaction_ratios", "seed_design", "model", "fit", "sampling", "measurement"]:
        if config[key] != historical_config[key]:
            raise ValueError(f"Matched RUN-175 design changed: {key}")
    if config["observation_modes"]["TrpB"] != historical_config["observation_modes"]["TrpB"]:
        raise ValueError("Matched RUN-175 observation modes changed")


def run(config_path: Path, output: Path, data: Path, matrix_path: Path,
        metadata_path: Path, source_manifest: Path, historical_manifest: Path,
        workers: int = 4, quick: bool = False):
    start = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("Use an empty output directory; prior results must not be overwritten")
    config = yaml.safe_load(config_path.read_text())
    historical_config = yaml.safe_load(
        (historical_manifest.parent / "interaction_strength.yaml").read_text()
    )
    validate_config(config, historical_config)
    kernels = verify_sources(source_manifest, historical_manifest)
    matrix, metadata = load_prior(matrix_path, metadata_path, config)
    inputs = json.loads((data / "manifest.json").read_text())
    for relative, expected in inputs["files"].items():
        if sweep.sha(data / relative) != expected:
            raise ValueError(f"Frozen input checksum mismatch: {relative}")
    for source, name in [(config_path, "config.yaml"),
                         (source_manifest, "source_manifest.json"),
                         (historical_manifest, "historical_source_manifest.json"),
                         (data / "manifest.json", "input_manifest.json"),
                         (matrix_path, "additive_matrix.npy"),
                         (metadata_path, "prior_manifest.json")]:
        shutil.copyfile(source, output / name)
    environment = dict(
        python=platform.python_version(),
        packages={name: importlib.metadata.version(name)
                  for name in ["numpy", "pandas", "scipy", "pyarrow", "PyYAML"]},
    )
    sweep.dump(output / "worker_environment.json", environment)
    if not quick and environment != config["required_worker_environment"]:
        raise ValueError("Worker Python/numerical dependencies differ from recorded RUN-175; "
                         "restore the recorded environment before the matched scientific run")
    sweep.dump(output / "intervention.json", dict(
        prior="empirical_msa_log_profile", reference_policy=REFERENCE_POLICY,
        unchanged_scientific_sources=kernels, prior_metadata=metadata,
        synthetic_control="Reuse recorded RUN-175 outputs; no synthetic refit",
        endpoint_policy="Repeat original RUN-172 fitting under changed prior, then use those endpoints in RUN-175 fitting",
    ))
    sweep.graph.validate_kernel(output)
    sweep.graph.KERNEL = sweep.graph.kernel()
    original_config = yaml.safe_load((data / "config.yaml").read_text())
    original_settings = json.loads((data / "selected_parameters.json").read_text())
    anchor = next(s for s in original_settings if s["dataset_id"] == DATASET
                  and s["arm"] == "r4_sigmoid" and s["observation"] == "none")
    origin, scale, floor = (anchor[k] for k in
                            ["calibration_origin", "calibration_scale", "near_floor_raw"])
    sub = data / "TrpB"
    geometry = json.loads((sub / "geometry.json").read_text())
    if geometry["alphabets"] != [list(AA_ORDER)] * 4:
        raise ValueError("Frozen assay codes do not match the prior's amino-acid identities")
    table = pd.read_parquet(sub / "variants.parquet")
    measured_codes = np.column_stack([
        pd.Categorical(table[f"pos{j+1}"].astype(str), categories=list(AA_ORDER)).codes
        for j in range(4)
    ])
    if (measured_codes < 0).any():
        raise ValueError("Unknown empirical amino acid")
    measured_indices = np.ravel_multi_index(measured_codes.T, (20,) * 4)
    if len(np.unique(measured_indices)) != len(measured_indices):
        raise ValueError("Duplicate empirical genotype")
    codes = np.array(np.unravel_index(np.arange(20**4), (20,) * 4)).T.astype(np.uint8)
    empirical = np.full(len(codes), np.nan)
    empirical[measured_indices] = table.fitness.to_numpy(float)
    valid = np.isfinite(empirical)
    if valid.sum() != geometry["unique_genotypes"]:
        raise ValueError("Empirical mask differs from the frozen geometry")
    np.savez_compressed(output / "TrpB_empirical.npz", fitness=empirical, valid=valid, genotypes=codes)
    with np.load(sub / "empirical_panel.npz") as panel:
        index = pd.read_csv(sub / "empirical_panel_index.csv")
        calibration_mask = index.split.eq("calibration").to_numpy()
        cal = panel["fitness"][calibration_mask]
        conf = panel["fitness"][~calibration_mask]
        conf_genotypes = panel["genotypes"][~calibration_mask]
    np.savez_compressed(output / "TrpB_empirical_cubes.npz", fitness=conf, genotypes=conf_genotypes)
    calibration_index = pd.read_csv(sub / "latent_calibration_index.csv")
    with np.load(sub / "latent_calibration.npz") as latent:
        cindex, ca, cg, calibration_audit = regenerate_panel(
            calibration_index, latent, config["model"], matrix
        )
        np.savez_compressed(output / "latent_calibration.npz", a=ca, g=cg, genotypes=latent["genotypes"])
    cindex.to_csv(output / "latent_calibration_index.csv", index=False)
    sweep.dump(output / "calibration_prior_audit.json", calibration_audit)
    noise_dir = output / "noise"
    noise_dir.mkdir()
    noise, _ = sweep.calibrate_noise(
        pd.read_csv(sub / "donor_assignments.csv.gz"),
        set(json.loads((sub / "excluded_noise_donors.json").read_text())),
        original_config, noise_dir,
    )
    old_noise = json.loads((sub / "noise_model.json").read_text())
    new_noise = json.loads((noise_dir / "noise_model.json").read_text())
    for key in ["depth_breaks", "fitness_breaks", "shared_covariance_raw"]:
        np.testing.assert_allclose(new_noise[key], old_noise[key], rtol=0, atol=1e-14)
    with np.load(sub / "fit_observation_uniforms.npz") as uniform_file:
        uniforms = uniform_file["uniforms"]
    banks = {"none": None}
    for mode in config["observation_modes"]["TrpB"][1:]:
        banks[mode] = sweep.observation_bank(noise, uniforms, mode == "shared_stop_sensitivity")
    fitted = fit_settings(ca, cg, (cal-origin)/scale, banks, origin, scale, floor,
                          config, original_config, output / "fits", quick)
    sweep.dump(output / "selected_parameters.json", fitted)
    relative = config["measurement"]["relative_resolutions"]
    ratios = config["interaction_ratios"]
    base = dict(protein="TrpB", seed=-1, cohort="empirical", ratio=-1.,
                observation="measured", layer="empirical")
    graphs = [sweep.graph_record(
        np.nan_to_num(empirical), valid, rel*scale,
        float(np.quantile(empirical[valid], .99)), floor,
        dict(**base, relative_resolution=rel),
    ) for rel in relative]
    spectra = sweep.cube_records(conf, pd.DataFrame({"cube": np.arange(len(conf))}),
                                 relative, scale, floor, base)
    latent_graphs, seed_audits = [], []
    old_index = pd.read_csv(sub / "latent_confirmation_index.csv")
    with np.load(sub / "latent_confirmation.npz") as values:
        old_latent = {key: values[key] for key in values.files}
    old_seeds = sorted(old_index.seed.unique().tolist())
    fresh = [production.seed_for(DATASET) + n for n in config["seed_design"]["fresh_seed_offsets"]]
    if set(old_seeds) & set(fresh) or set(calibration_index.seed) & set(old_seeds + fresh):
        raise ValueError("Calibration and confirmation seed cohorts overlap")
    if quick:
        old_seeds, fresh = old_seeds[:1], fresh[:1]
    tasks = ([(seed, "reused", i) for i, seed in enumerate(old_seeds)]
             + [(seed, "fresh", i) for i, seed in enumerate(fresh)])
    fields, cubes_dir = output / "latent_fields", output / "cube_values"
    fields.mkdir(); cubes_dir.mkdir()

    def one_seed(task):
        seed, cohort, seed_index = task
        model = build_model(seed, config["model"], matrix)
        a, g = model.latent_components(codes)
        np.savez_compressed(fields / f"TrpB_{seed}.npz", a=a, g=g)
        maximum_g_difference = None
        if cohort == "reused":
            selected = np.flatnonzero(old_index.seed.eq(seed))
            genotype_cubes = old_latent["genotypes"][selected]
            cube_indices = np.ravel_multi_index(genotype_cubes.reshape(-1, 4).T, (20,)*4)
            expected_g = old_latent["g"][selected].ravel()
            np.testing.assert_allclose(g[cube_indices], expected_g, rtol=1e-10, atol=1e-10)
            maximum_g_difference = float(np.max(np.abs(g[cube_indices] - expected_g)))
        else:
            count = config["sampling"]["synthetic_cubes_per_seed"]
            selected = (np.arange(count) + seed_index*count) % len(conf_genotypes)
            genotype_cubes = conf_genotypes[selected]
            cube_indices = np.ravel_multi_index(genotype_cubes.reshape(-1, 4).T, (20,)*4)
        noise_seed = config["seed_design"]["noise_seed_offset"] + int(seed)
        uniforms = np.random.default_rng(noise_seed).random((len(codes), 6))
        graph_rows, cube_rows, latent_rows = [], [], []
        saved = {"genotypes": genotype_cubes}
        audit = dict(protein="TrpB", seed=int(seed), cohort=cohort, noise_seed=noise_seed,
                     reference=model.reference.tolist(), fingerprint=model.fingerprint,
                     normalized_g_checked_nodes=len(cube_indices) if cohort == "reused" else 0,
                     maximum_normalized_g_difference=maximum_g_difference)
        for ratio in ratios:
            latent = a + ratio*g
            label = dict(protein="TrpB", seed=int(seed), cohort=cohort, ratio=float(ratio),
                         observation="none", layer="latent", relative_resolution=0.)
            latent_rows.append(sweep.graph_record(
                latent, valid, 0., float(np.quantile(latent[valid], .99)),
                float(latent[valid].min()-1), label,
            ))
        for setting in fitted:
            ratio, mode = setting["ratio"], setting["observation"]
            raw = origin + scale*sweep.predict(setting["parameters"], a, g, ratio)
            replay = build_model(seed, config["model"], matrix, interaction_ratio=ratio,
                                 reference_occupancy=float(expit(setting["position"])),
                                 steepness=setting["steepness"], fitness_min=setting["offset"],
                                 fitness_max=setting["offset"]+setting["scale"])
            np.testing.assert_allclose(replay.evaluate_batch(codes[::211]), raw[::211],
                                       rtol=1e-10, atol=1e-10*scale)
            maps = [("intrinsic", raw)]
            if mode != "none":
                r1, r2, _, _ = sweep.observe(raw, noise, uniforms, shared=mode == "shared_stop_sensitivity")
                maps.append(("observed", (r1+r2)/2))
            for layer, values in maps:
                label = dict(protein="TrpB", seed=int(seed), cohort=cohort,
                             ratio=ratio, observation=mode, layer=layer)
                for rel in relative:
                    row = sweep.graph_record(values, valid, rel*scale,
                                             float(np.quantile(values[valid], .99)), floor,
                                             dict(**label, relative_resolution=rel))
                    row["unique_fitness_fraction"] = len(np.unique(values[valid]))/valid.sum()
                    graph_rows.append(row)
                cube_values = values[cube_indices].reshape(genotype_cubes.shape[:2])
                saved[f"{mode}_{layer}_r{ratio:g}"] = cube_values
                cube_rows += sweep.cube_records(cube_values,
                    pd.DataFrame({"cube": np.arange(len(cube_values))}), relative, scale, floor, label)
        np.savez_compressed(cubes_dir / f"TrpB_{seed}.npz", **saved)
        return graph_rows, cube_rows, latent_rows, audit

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for grows, crows, lrows, audit in pool.map(one_seed, tasks):
            graphs += grows; spectra += crows; latent_graphs += lrows; seed_audits.append(audit)
            pd.DataFrame(graphs).to_csv(output / "graph_summary.csv", index=False)
            print("GRAPHS", audit["cohort"], audit["seed"], "complete", flush=True)
    pd.DataFrame(spectra).to_csv(output / "cube_summary.csv", index=False)
    pd.DataFrame(latent_graphs).to_csv(output / "latent_graph_summary.csv", index=False)
    pd.DataFrame(seed_audits).to_csv(output / "seed_audit.csv", index=False)
    pd.DataFrame([dict(protein="TrpB", nodes=int(valid.sum()), possible_nodes=len(codes),
                       origin=origin, scale=scale, floor=floor,
                       calibration_cubes=len(cal), confirmation_cubes=len(conf))]).to_csv(
        output / "coverage.csv", index=False)
    files = {str(p.relative_to(output)): sweep.sha(p) for p in sorted(output.rglob("*"))
             if p.is_file() and p.name != "manifest.json"}
    sweep.dump(output / "manifest.json", dict(
        protocol=config["protocol"], mode="SMOKE_ONLY" if quick else "MATCHED_PRIOR_REPEAT",
        public_model=sweep.PUBLIC_MODEL_SYSTEM, elapsed_seconds=time.perf_counter()-start,
        source_run="RUN-175", precursor_endpoint_run="RUN-172",
        input_manifest_sha256=sweep.sha(data / "manifest.json"),
        config_sha256=sweep.sha(config_path), source_manifest_sha256=sweep.sha(source_manifest),
        additive_prior=config["additive_prior"], reference_policy=REFERENCE_POLICY,
        settings=len(fitted), graph_rows=len(graphs), latent_graph_rows=len(latent_graphs),
        cube_summary_rows=len(spectra), seed_audits=len(seed_audits), files=files,
        interpretation=config["interpretation"],
    ))
    print("COMPLETE", len(fitted), "settings", len(seed_audits), "seeds", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--prior-matrix", type=Path, required=True)
    parser.add_argument("--prior-metadata", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--historical-source-manifest", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    historical = (args.historical_source_manifest or
                  args.config.parent / "interaction_strength_source.json")
    run(args.config, args.output_dir, args.data_root, args.prior_matrix,
        args.prior_metadata, args.source_manifest, historical, args.workers, args.quick)
