#!/usr/bin/env python3
"""Replay intrinsic release numerics against frozen RUN-175/RUN-178 evidence.

This is a bounded engineering check, not calibration or empirical validation.
Independent subprocesses prevent current and archived package imports mixing.
The shared-stop observation layer is not replayed: packaged presets are intrinsic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import scipy
import yaml


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration"
PRESETS = {
    "RUN-175": "trpb_run175_r05_shared_stop_intrinsic",
    "RUN-178": "trpb_run178_r05_profile_shared_stop_intrinsic",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_hashes(root, files):
    for relative, expected in files.items():
        if isinstance(expected, dict):
            expected = expected["sha256"]
        assert sha(root / relative) == expected, str(root / relative)


def worker(source: Path, run: str, destination: Path, current: bool):
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(source))
    from adaptenv import ThermodynamicConfig, ThermodynamicLandscape, __version__

    output = BASE / run / "outputs"
    config = yaml.safe_load((output / "config.yaml").read_text())
    settings = json.loads((output / "selected_parameters.json").read_text())
    selected = next(row for row in settings if row["protein"] == "TrpB"
                    and row["ratio"] == 0.5 and row["observation"] == "shared_stop_sensitivity")
    seed_table = pd.read_csv(output / "seed_audit.csv").query("protein == 'TrpB'").sort_values("seed")
    assert len(seed_table) == 24
    with np.load(output / "TrpB_empirical.npz") as empirical:
        sampled_codes = empirical["genotypes"][::157]
    matrix = np.load(output / "additive_matrix.npy") if run == "RUN-178" else None
    parameters = selected["parameters"]
    occupancy = float(1.0 / (1.0 + np.exp(-parameters[2])))
    latent_a, latent_g, cube_fitness, references, matrices, fingerprints = [], [], [], [], [], []
    for row in seed_table.itertuples():
        fields = dict(length=4, seed=int(row.seed), **config["model"])
        if current and run == "RUN-175":
            fields["additive_prior"] = "synthetic"
        neutral = ThermodynamicLandscape(ThermodynamicConfig(**fields), additive_matrix=matrix)
        a, g = neutral.latent_components(sampled_codes)
        latent_a.append(a)
        latent_g.append(g)
        if current:
            from adaptenv import build_landscape_from_preset
            model = build_landscape_from_preset(PRESETS[run], seed=int(row.seed))
            assert model.config.additive_prior == ("synthetic" if run == "RUN-175" else "profile")
            if run == "RUN-178":
                assert "profile_prior" in model.summary()
        else:
            fields.update(reference_occupancy=occupancy, steepness=selected["steepness"],
                          interaction_ratio=selected["ratio"], fitness_min=selected["offset"],
                          fitness_max=selected["offset"] + selected["scale"])
            model = ThermodynamicLandscape(ThermodynamicConfig(**fields), additive_matrix=matrix)
        with np.load(output / "cube_values" / f"TrpB_{row.seed}.npz") as cubes:
            scores = model.evaluate_batch(cubes["genotypes"].reshape(-1, 4)).reshape(8, 16)
        cube_fitness.append(scores)
        references.append(model.reference)
        matrices.append(model.additive_matrix)
        fingerprints.append(model.fingerprint)
    np.savez_compressed(destination, seeds=seed_table.seed.to_numpy(),
                        a=np.asarray(latent_a), g=np.asarray(latent_g),
                        cube_fitness=np.asarray(cube_fitness), references=np.asarray(references),
                        matrices=np.asarray(matrices), sampled_codes=sampled_codes)
    destination.with_suffix(".json").write_text(json.dumps(dict(
        run=run, package_version=__version__, package_source=str(source),
        factory="build_landscape_from_preset" if current else "archived constructor and saved readout",
        preset=PRESETS[run], fingerprints=fingerprints), indent=2) + "\n")


def replay():
    records = []
    current_package = ROOT / "external/Adapt-Env"
    current_inputs = sorted((current_package / "adaptenv").glob("*.py")) + [
        current_package / "adaptenv/presets.yaml",
        current_package / "adaptenv/trpb_profile_run178.json",
        current_package / "pyproject.toml",
    ]
    current_hashes = {str(path.relative_to(current_package)): sha(path) for path in current_inputs}
    with tempfile.TemporaryDirectory(prefix="adaptenv-release-replay-") as temporary:
        temp = Path(temporary)
        for run in PRESETS:
            output = BASE / run / "outputs"
            manifest = json.loads((output / "manifest.json").read_text())
            check_hashes(output, manifest["files"])
            source_manifest = json.loads((output / "source_manifest.json").read_text())
            check_hashes(BASE / run / "source", source_manifest["files"])
            profile_payload_hashes = 0
            if run == "RUN-178":
                inputs = BASE / run / "inputs/trpb-profile-additive-inputs"
                profile_files = json.loads((inputs / "prior_manifest.json").read_text())["files"]
                check_hashes(inputs, profile_files)
                profile_payload_hashes = len(profile_files)
            products = {}
            for label, source in (("archived", BASE / run / "source/external/Adapt-Env"),
                                  ("current", ROOT / "external/Adapt-Env")):
                path = temp / f"{run}-{label}.npz"
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--run", run,
                           "--source", str(source), "--destination", str(path)]
                if label == "current":
                    command.append("--current")
                result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=180)
                if result.returncode:
                    raise RuntimeError(f"{label} replay failed for {run}:\n{result.stdout}\n{result.stderr}")
                with np.load(path) as arrays:
                    products[label] = {key: arrays[key] for key in arrays.files}
                products[label]["metadata"] = json.loads(path.with_suffix(".json").read_text())
            before, after = products["archived"], products["current"]
            for key in ("seeds", "references", "sampled_codes"):
                np.testing.assert_array_equal(before[key], after[key])
            maxima = {}
            for key in ("a", "g", "cube_fitness", "matrices"):
                np.testing.assert_allclose(before[key], after[key], rtol=1e-12, atol=1e-12)
                maxima[f"current_vs_archived_{key}"] = float(np.max(np.abs(before[key] - after[key])))
            for i, seed in enumerate(after["seeds"]):
                with np.load(output / "latent_fields" / f"TrpB_{seed}.npz") as field:
                    for key in ("a", "g"):
                        saved = field[key][::157]
                        np.testing.assert_allclose(after[key][i], saved, rtol=1e-12, atol=1e-12)
                        name = f"current_vs_saved_{key}"
                        maxima[name] = max(maxima.get(name, 0.0), float(np.max(abs(after[key][i] - saved))))
                with np.load(output / "cube_values" / f"TrpB_{seed}.npz") as cube:
                    saved = cube["shared_stop_sensitivity_intrinsic_r0.5"]
                    np.testing.assert_allclose(after["cube_fitness"][i], saved, rtol=1e-12, atol=1e-12)
                    name = "current_preset_vs_saved_intrinsic_fitness"
                    maxima[name] = max(maxima.get(name, 0.0), float(np.max(abs(after["cube_fitness"][i] - saved))))
            records.append(dict(run=run, preset=PRESETS[run], seeds=24,
                                latent_genotypes_per_seed=len(after["sampled_codes"]),
                                cube_genotypes_per_seed=128, output_hashes=len(manifest["files"]),
                                source_hashes=len(source_manifest["files"]), maximum_absolute_differences=maxima,
                                profile_payload_hashes=profile_payload_hashes,
                                archived_version=before["metadata"]["package_version"],
                                current_version=after["metadata"]["package_version"],
                                changed_fingerprints=sum(a != b for a, b in zip(
                                    before["metadata"]["fingerprints"], after["metadata"]["fingerprints"]))))
    check_hashes(current_package, current_hashes)
    artifact = dict(
        purpose="Engineering replay of release numerics and prior-matched intrinsic presets",
        passed=True, runs=records, replay_script_sha256=sha(Path(__file__)),
        current_package_input_sha256=current_hashes,
        environment=dict(python=sys.version, executable=sys.executable,
                         platform=platform.platform(), numpy=np.__version__,
                         pandas=pd.__version__, scipy=scipy.__version__, pyyaml=yaml.__version__),
        limitations=[
            "No optimization, scientific control campaign or observed-map reconstruction was performed.",
            "The packaged readouts are intrinsic; shared-stop calibration provenance does not enable observation noise.",
            "Numeric agreement does not establish profile superiority, biological realism or ML performance.",
            "A new prior-policy configuration field may change fingerprints even when the numerical map is unchanged.",
        ])
    destination = Path(__file__).with_name("replay.json")
    destination.write_text(json.dumps(artifact, indent=2, allow_nan=False) + "\n")
    print(json.dumps(artifact, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--run", choices=tuple(PRESETS))
    parser.add_argument("--source", type=Path)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--current", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.source, args.run, args.destination, args.current)
    else:
        replay()
