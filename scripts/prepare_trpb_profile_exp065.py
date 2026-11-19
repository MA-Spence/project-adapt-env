#!/usr/bin/env python3
"""Freeze the author's TrpB frequency profile; no calibration or queue activity.

The deposited raw A2M is NUL-filled. The complete deposited alignment-frequency
table is usable, but rounded and has no recoverable weighting/N_eff metadata.
It is therefore treated as a frequency profile, never as reconstructed counts.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "external/Adapt-Env"))

import numpy as np
from adaptenv.alignment import additive_prior_from_frequencies
from adaptenv.amino_acids import AA_ORDER

DEFAULT_OUTPUT = ROOT / "data/processed/trpb-profile-additive-inputs"
POSITIONS = [183, 184, 227, 228]
UNIFORM_MIXTURE = 0.01
SOURCE_ALPHABET = "-ACDEFGHIKLMNPQRSTVWY"
PRECISION = "Author alignment CSV decimal frequencies rounded to 0.001; sequence weighting and N_eff unavailable."
BASE_URL = "https://data.caltech.edu/records/h5rah-5z170/files/"
SOURCES = {
    "TARGET_b0.3.a2m": ("data", "data/EVmutation/TARGET_b0.3/align/TARGET_b0.3.a2m", "fea9258294fa60028e8ad93df807937f18157056e9b64311577ee27d5a2d0455"),
    "TARGET_b0.3.fa": ("data", "data/EVmutation/TARGET_b0.3/align/TARGET_b0.3.fa", "351a567fc1e709242d19f1c95af2cf2a8d5c27e9892886c19df96083920ceaa4"),
    "TARGET_b0.3_alignment_statistics.csv": ("data", "data/EVmutation/TARGET_b0.3/align/TARGET_b0.3_alignment_statistics.csv", "6eb9f5ac3426528173f09cd34aa92930f2a226e4f8a2823b44ac5519ecbe98a0"),
    "TARGET_b0.3_frequencies.csv": ("data", "data/EVmutation/TARGET_b0.3/align/TARGET_b0.3_frequencies.csv", "5093d7d9b7658a94af5264702386e773ece1358a73cbc963eb771d0f2114b6a5"),
    "TARGET_b0.3_single_mutant_matrix.csv": ("data", "data/EVmutation/TARGET_b0.3/mutate/TARGET_b0.3_single_mutant_matrix.csv", "4386083de6e627257e97c52084c9c64f872fc2edb7fbf479161552059edd0e68"),
    "globals.py": ("code", "code/scripts/growth_analysis/growth_analysis/globals.py", "1566f17e53ed5d4a90fe07ec6f8026ad400562030f56a83d30a92cf66e27e7fd"),
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def fetch_sources(source: Path) -> None:
    """Read published ZIP members by HTTP ranges; never execute external code."""
    from fetch_trpb_replicates_exp065 import RemoteZip

    source.mkdir(parents=True, exist_ok=True)
    for archive in ["data", "code"]:
        selected = [(name, member, expected) for name, (kind, member, expected)
                    in SOURCES.items() if kind == archive
                    and (not (source / name).exists() or sha(source / name) != expected)]
        if not selected:
            continue
        with zipfile.ZipFile(RemoteZip(BASE_URL + archive + ".zip?download=1")) as remote:
            for name, member, expected in selected:
                content = remote.read(member)  # zipfile verifies the member CRC.
                if hashlib.sha256(content).hexdigest() != expected:
                    raise ValueError(f"Published source checksum changed: {member}")
                (source / name).write_bytes(content)


def read_profile(source: Path):
    """Validate author sequence and assay coordinates independently of outcomes."""
    for name, (_, _, expected) in SOURCES.items():
        if sha(source / name) != expected:
            raise ValueError(f"Source checksum mismatch: {name}")
    assignments = {}
    for node in ast.parse((source / "globals.py").read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {"TM9D8S_AA", "LIBRARY_INFO_DICT"}:
                    assignments[target.id] = ast.literal_eval(node.value)
    sequence = "".join(line.strip() for line in (source / "TARGET_b0.3.fa").read_text().splitlines()
                       if not line.startswith(">"))
    if sequence != assignments["TM9D8S_AA"].removesuffix("*") or len(sequence) != 397:
        raise ValueError("EVmutation target does not match the author's 397-residue TrpB construct")
    library = assignments["LIBRARY_INFO_DICT"]["4-site"]
    positions = [library["positions"][i] for i in range(1, 5)]
    parent = "".join(library["AAs"][i] for i in range(1, 5))
    if positions != POSITIONS or parent != "VFVS" or "".join(sequence[p-1] for p in positions) != parent:
        raise ValueError("Physical TrpB assay position/parent mapping mismatch")
    with (source / "TARGET_b0.3_frequencies.csv").open() as handle:
        records = list(csv.DictReader(handle))
    if [int(row["i"]) for row in records] != list(range(1, 398)):
        raise ValueError("Frequency rows must cover all author query positions in order")
    if "".join(row["A_i"].upper() for row in records) != sequence:
        raise ValueError("Frequency-table reference residues disagree with the author target")
    frequencies = np.asarray([[float(row[a]) for a in SOURCE_ALPHABET] for row in records])
    return sequence, frequencies


def prepare(output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    source = output / "source"
    sequence, frequencies = read_profile(source)
    rows = [p - 1 for p in POSITIONS]
    matrix, conversion = additive_prior_from_frequencies(
        frequencies, alphabet=SOURCE_ALPHABET, rows=rows,
        uniform_mixture=UNIFORM_MIXTURE, source_precision=PRECISION,
    )
    np.save(output / "additive_matrix.npy", matrix, allow_pickle=False)
    np.save(output / "source_frequencies.npy", frequencies, allow_pickle=False)
    selected = frequencies[rows, 1:]
    probabilities = (1 - UNIFORM_MIXTURE) * selected / selected.sum(axis=1, keepdims=True) + UNIFORM_MIXTURE / 20
    np.save(output / "selected_probabilities.npy", probabilities, allow_pickle=False)
    source_manifest = {
        "source_doi": "10.22002/h5rah-5z170", "paper_doi": "10.1073/pnas.2400439121",
        "source_selection": "Published Methods specifies Tm9D8* EVmutation bitscore 0.3; original author profile selected before calibration outcomes.",
        "retrieval": "Original ZIP members by HTTP ranges; ZIP member CRC checked. Complete archives not downloaded or MD5-verified.",
        "files": [{"url": BASE_URL + archive + ".zip?download=1", "member": member,
                   "path": "source/" + name, "bytes": (source / name).stat().st_size, "sha256": expected}
                  for name, (archive, member, expected) in SOURCES.items()],
        "raw_alignment_audit": {
            "status": "unusable_source_member_all_NUL_bytes", "bytes": 8101144, "crc32": 2853580136,
            "independent_verification": "Alternate /api/records/h5rah-5z170/files/data.zip/content endpoint, no-cache byte range, manual raw-DEFLATE decompression reproduced identical SHA/CRC. Author ftmlde_data duplicate has identical NUL bytes.",
            "used_to_estimate_frequencies": False,
        },
        "other_sources_rejected": [
            "Original full-precision single-mutant CSV gives only 19 non-reference amino acids per site. Missing WT frequency and uncertain gap conditioning prevent safe reconstruction.",
            "Zenodo13910506 SSMuLA TrpB.model header/full frequencies inspected: later 22746 valid +169 invalid sequences, not original18857 alignment; values differ. Not used.",
        ],
    }
    dump(output / "source_manifest.json", source_manifest)
    with (source / "TARGET_b0.3_alignment_statistics.csv").open() as handle:
        statistics = next(csv.DictReader(handle))
    manifest = {
        "schema": "trpb_frequency_profile_additive_prior_v1", "dataset": "TrpB",
        "input_kind": "author_alignment_frequency_table_not_raw_MSA",
        "alphabet": AA_ORDER, "assay_positions_1based": POSITIONS,
        "matrix_rows_correspond_to_assay_positions_1based": POSITIONS,
        "selected_query_columns_zero_based": rows, "profile_rows_zero_based": rows,
        "biological_parent_at_assay_positions": "VFVS",
        "profile_consensus_at_assay_positions": "".join(AA_ORDER[i] for i in selected.argmax(axis=1)),
        "author_target_sequence": sequence,
        "author_target_length": len(sequence), "author_target_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "terminal_construct_tag": "LEHHHHHH at positions390-397; retained in author query; none of the four selected rows overlap it",
        "author_alignment_statistics": statistics,
        "sequence_weighting": "unknown for deposited alignment frequency table; not inferred from EVmutation default settings",
        "n_effective_sequences": None, "sequence_weights_available": False,
        "raw_MSA_available": False, "source_frequency_precision": PRECISION,
        "additive_score_formula": "ln(20 * ((1 - 0.01) * f_ia / sum_b_in_20AA(f_ib) + 0.01 / 20))",
        "uniform_probability_mixture": UNIFORM_MIXTURE, "count_pseudocount": None,
        "calibration_outcomes_used_to_construct_prior": False,
        "synthetic_reference_policy": "unchanged seed-specific historical reference; biological VFVS is a coordinate check, not a replacement reference",
        "interaction_prior": "unchanged procedural GNK; no covariance or Potts couplings imported",
        "interpretation": "Family conservation prior; evolutionary frequency is not the measured assay fitness or a validated stability energy.",
        "conversion": conversion,
        "preparation_script": "scripts/prepare_trpb_profile_exp065.py",
        "preparation_script_sha256": sha(Path(__file__)),
    }
    manifest["files"] = {str(path.relative_to(output)): sha(path) for path in sorted(output.rglob("*"))
                         if path.is_file() and path.name != "prior_manifest.json"}
    dump(output / "prior_manifest.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fetch-source", action="store_true", help="Retrieve/check the exact public ZIP members first")
    args = parser.parse_args()
    if args.fetch_source:
        fetch_sources(args.output / "source")
    manifest = prepare(args.output)
    print(json.dumps({"output": str(args.output), "matrix_sha256": manifest["files"]["additive_matrix.npy"],
                      "prior_manifest_sha256": sha(args.output / "prior_manifest.json"),
                      "positions": POSITIONS, "profile_consensus": manifest["profile_consensus_at_assay_positions"]}, indent=2))


if __name__ == "__main__":
    main()
