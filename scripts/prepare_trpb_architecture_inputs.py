#!/usr/bin/env python3
"""Freeze bounded, genotype-disjoint TrpB panels for architecture calibration.

Selection uses only the historical measurement mask, genotype order and fixed
random seed. Actual published replicate columns are retained without averaging.
This creates inputs, not a scientific RUN or a calibration result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from itertools import combinations
import json
from pathlib import Path
import shutil

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
REPLICATES = ROOT / "data/processed/trpb-replicate-observation-exp065"
HISTORICAL = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs"
PINNED = {
    "source/replicates/replicates.csv": "d3c19fa231fbbb20150c304d32f93d757898921651da0eaff8d6e6711461515a",
    "source/replicates/source_manifest.json": "d9af3a8cd15da4769d1c2f076f5ff07030eaf9a3d2785b79bf1859ff628410bc",
    "source/historical/TrpB_empirical.npz": "5eb0efdb3ca65589b561bbb8a68a7e990d6fd66263f94e249b26cce01cb7e5ac",
    "source/profile/additive_matrix.npy": "d5ac4d8724fc753709887afbe5a82a172c47c5a5bf18cc8a13dd1977875f7be3",
    "source/profile/prior_manifest.json": "408ffe4a5052fe0d2525ed5c52113381753ba6cfb9fce3d775885648a804abd9",
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sources(source_bundle=None):
    paths = {
        "source/replicates/replicates.csv": REPLICATES / "replicates.csv",
        "source/replicates/input_manifest.json": REPLICATES / "input_manifest.json",
        "source/replicates/source_manifest.json": REPLICATES / "source_manifest.json",
        "source/historical/TrpB_empirical.npz": HISTORICAL / "TrpB_empirical.npz",
        "source/historical/run180_manifest.json": HISTORICAL / "manifest.json",
        "source/profile/additive_matrix.npy": HISTORICAL / "additive_matrix.npy",
        "source/profile/prior_manifest.json": HISTORICAL / "prior_manifest.json",
    }
    if source_bundle is not None:
        paths = {name: Path(source_bundle) / name for name in paths}
    for relative, expected in PINNED.items():
        if sha(paths[relative]) != expected:
            raise ValueError(f"Frozen input checksum mismatch: {relative}")
    source_manifest = json.loads(paths["source/replicates/source_manifest.json"].read_text())
    for record in source_manifest["files"]:
        relative = "source/replicates/" + record["file"]
        path = Path(source_bundle) / relative if source_bundle is not None else REPLICATES / record["file"]
        if sha(path) != record["sha256"]:
            raise ValueError(f"Author source checksum mismatch: {record['file']}")
        paths[relative] = path
    return paths


def load_inputs(paths):
    """Check canonical library identity and historical mean agreement."""
    with np.load(paths["source/historical/TrpB_empirical.npz"], allow_pickle=False) as archive:
        all_sequences = archive["genotypes"].copy()
        historical_fitness = archive["fitness"].copy()
        valid = archive["valid"].copy()
    indices = np.flatnonzero(valid)
    if (all_sequences.shape != (160000, 4) or valid.shape != (160000,)
            or valid.dtype != np.bool_ or historical_fitness.shape != (160000,)
            or len(indices) != 159129):
        raise ValueError("Expected the unchanged RUN-180 four-site TrpB measurement mask")
    code = all_sequences.astype(np.int64) @ np.array([8000, 400, 20, 1])
    if not np.array_equal(code, np.arange(160000)):
        raise ValueError("Historical sequence table is not canonical q20 lexicographic order")
    sequences = all_sequences[indices]
    lookup = {}
    rows = 0
    with paths["source/replicates/replicates.csv"].open(newline="") as stream:
        for row in csv.DictReader(stream):
            rows += 1
            genotype = row["genotype"]
            if len(genotype) == 4 and all(aa in AA_ORDER for aa in genotype):
                if genotype in lookup:
                    raise ValueError("Duplicate genotype in published replicate table")
                lookup[genotype] = (float(row["replicate_1"]), float(row["replicate_2"]))
    if rows != 193170 or len(lookup) != 159129:
        raise ValueError("Unexpected source replicate population")
    names = ["".join(AA_ORDER[int(allele)] for allele in row) for row in sequences]
    if set(names) != set(lookup):
        raise ValueError("Replicate and historical measured genotype masks differ")
    fitness = np.asarray([lookup[name] for name in names], dtype=float)
    if not np.isfinite(fitness).all():
        raise ValueError("Expected two finite actual measurements per selected genotype")
    mean_error = float(np.max(np.abs(fitness.mean(axis=1) - historical_fitness[indices])))
    if mean_error > 5.1e-8:
        raise ValueError("Published replicate means disagree with historical fitness beyond rounding")
    matrix = np.load(paths["source/profile/additive_matrix.npy"], allow_pickle=False)
    prior = json.loads(paths["source/profile/prior_manifest.json"].read_text())
    if (matrix.shape != (4, 20) or not np.isfinite(matrix).all()
            or prior["alphabet"] != AA_ORDER
            or prior["assay_positions_1based"] != [183, 184, 227, 228]):
        raise ValueError("Prior must preserve the archived q20 alphabet and TrpB site mapping")
    return sequences, fitness, matrix, indices, mean_error


def fixed_folds(sequences, n_folds, seed):
    """Same lexicographic genotype-group assignment as adaptenv.calibration."""
    unique, inverse = np.unique(sequences, axis=0, return_inverse=True)
    if n_folds < 2 or n_folds > len(unique):
        raise ValueError("n_folds must be between two and the number of unique genotypes")
    labels = np.empty(len(unique), dtype=np.int64)
    labels[np.random.default_rng(seed).permutation(len(unique))] = np.arange(len(unique)) % n_folds
    return labels[inverse]


def choose_panels(sequences, *, n_calibration, n_confirmation, selection_seed, n_folds, inference_seed):
    if (min(n_calibration, n_confirmation) < 4
            or n_calibration + n_confirmation > len(sequences)):
        raise ValueError("Panel sizes must be >=4 and fit within the measured genotype population")
    order = np.random.default_rng(selection_seed).permutation(len(sequences))
    selected = {"calibration": order[:n_calibration],
                "confirmation": order[n_calibration:n_calibration + n_confirmation]}
    return {name: (rows, fixed_folds(sequences[rows], n_folds, inference_seed))
            for name, rows in selected.items()}


def coverage(sequences, folds):
    result = []
    for i, j in combinations(range(4), 2):
        codes = sequences[:, i].astype(np.int64) * 20 + sequences[:, j]
        counts = np.bincount(codes, minlength=400)
        result.append(dict(sites_zero_based=[i, j], observed_cells=int(np.count_nonzero(counts)),
                           minimum_cell_count=int(counts.min()), maximum_cell_count=int(counts.max()),
                           minimum_training_cell_count=min(int(np.bincount(codes[folds != f], minlength=400).min())
                                                           for f in np.unique(folds))))
    return result


def prepare(output, *, n_calibration=8192, n_confirmation=8192,
            selection_seed=20261121, n_folds=10, inference_seed=20261122, source_bundle=None):
    output = Path(output)
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be new or empty; existing scientific inputs are not overwritten")
    paths = _sources(source_bundle)
    sequences, fitness, matrix, library_indices, mean_error = load_inputs(paths)
    panels = choose_panels(sequences, n_calibration=n_calibration, n_confirmation=n_confirmation,
                           selection_seed=selection_seed, n_folds=n_folds, inference_seed=inference_seed)
    output.mkdir(parents=True, exist_ok=True)
    for name, source in paths.items():
        destination = output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    script = output / "source/preparation_script.py"
    shutil.copyfile(Path(__file__), script)
    panel_records = {}
    with (output / "selection.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["panel", "panel_row", "genotype", "canonical_library_index", "measured_population_index", "fold"])
        for name, (rows, folds) in panels.items():
            np.savez_compressed(output / f"{name}.npz", sequences=sequences[rows], fitness=fitness[rows],
                                additive_matrix=matrix, folds=folds)
            for row, (index, fold) in enumerate(zip(rows, folds)):
                writer.writerow([name, row, "".join(AA_ORDER[int(a)] for a in sequences[index]),
                                 int(library_indices[index]), int(index), int(fold)])
            panel_records[name] = dict(unique_genotypes=len(rows), actual_replicates_per_genotype=2,
                                       observed_measurements=int(fitness[rows].size),
                                       fold_sizes=np.bincount(folds, minlength=n_folds).tolist(),
                                       pair_allele_coverage=coverage(sequences[rows], folds))
    manifest = dict(
        schema="trpb_architecture_calibration_inputs_v1", dataset="Johnston 2024 four-site TrpB",
        model_alphabet=AA_ORDER, assay_positions_1based=[183, 184, 227, 228],
        population_unique_genotypes=len(sequences), disjoint_calibration_confirmation=True,
        selection=dict(method="single PCG64 uniform permutation of lexicographically ordered measured canonical genotypes",
                       selection_seed=selection_seed, uses_fitness_values=False,
                       calibration_count=n_calibration, confirmation_count=n_confirmation,
                       genotype_overlap=0, inference_seed=inference_seed, n_folds=n_folds,
                       folds_rule="same genotype-grouped permutation rule as adaptenv.calibration._folds"),
        measurements=dict(source="published mu_1-bg/max and mu_2-bg/max columns, preserved as two observations",
                          transformed=False, averaged_for_fitting=False,
                          maximum_replicate_mean_vs_historical_fitness_error=mean_error,
                          historical_mean_rounding_tolerance=5.1e-8,
                          shared_normalization="Both published replicates share average input-count normalization; independence of assay errors is not inferred."),
        latent_coefficient_degrees_of_freedom=1 + 4 * 19 + 6 * 19 ** 2,
        panels=panel_records,
        interpretation=["New fitting protocol; not an exact repeat of RUN-180's objective or sampling design.",
                        "Confirmation is disjoint for this fit; the full public TrpB landscape was inspected in historical work, so it is not previously unseen biological evidence.",
                        "Selection is fixed without conditioning on fitness, activity labels or confirmation outcomes.",
                        "Pair-cell coverage and sample counts do not prove design rank, coefficient identifiability or numerical convergence.",
                        "This bundle prepares inputs only; it contains no calibration result or new biological adequacy claim."],
        source_origins={name: str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
                        for name, path in paths.items()},
        preparation_script_sha256=sha(script),
        files={str(path.relative_to(output)): sha(path) for path in sorted(output.rglob("*")) if path.is_file()},
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/processed/trpb-architecture-calibration-inputs")
    parser.add_argument("--calibration-count", type=int, default=8192)
    parser.add_argument("--confirmation-count", type=int, default=8192)
    parser.add_argument("--selection-seed", type=int, default=20261121)
    parser.add_argument("--inference-seed", type=int, default=20261122)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument("--source-bundle", type=Path, help="Replay from a previously prepared self-contained input bundle")
    args = parser.parse_args(argv)
    manifest = prepare(args.output_dir, n_calibration=args.calibration_count, n_confirmation=args.confirmation_count,
                       selection_seed=args.selection_seed, n_folds=args.folds, inference_seed=args.inference_seed,
                       source_bundle=args.source_bundle)
    print(json.dumps({"output": str(args.output_dir.resolve()), "panels": manifest["panels"],
                      "files": len(manifest["files"])}, indent=2))


if __name__ == "__main__":
    main()
