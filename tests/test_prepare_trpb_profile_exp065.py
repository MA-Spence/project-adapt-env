"""Scientific input controls for the frozen author-profile intervention."""
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import prepare_trpb_profile_exp065 as prior


@pytest.fixture
def prepared(tmp_path):
    shutil.copytree(prior.DEFAULT_OUTPUT / "source", tmp_path / "source")
    metadata = prior.prepare(tmp_path)
    return tmp_path, metadata


def test_source_coordinates_and_hand_checked_author_frequencies(prepared):
    path, metadata = prepared
    matrix = np.load(path / "additive_matrix.npy", allow_pickle=False)
    assert matrix.shape == (4, 20)
    assert metadata["alphabet"] == "ACDEFGHIKLMNPQRSTVWY"
    assert metadata["assay_positions_1based"] == [183, 184, 227, 228]
    assert metadata["matrix_rows_correspond_to_assay_positions_1based"] == [183, 184, 227, 228]
    assert metadata["biological_parent_at_assay_positions"] == "VFVS"
    assert metadata["profile_consensus_at_assay_positions"] == "ILVG"
    assert "".join(metadata["author_target_sequence"][p-1] for p in prior.POSITIONS) == "VFVS"
    # Author CSV has 183I=.286 with non-gap sum=.999, and 228G=.989
    # with non-gap sum=.993. These independently check physical row/alphabet mapping.
    assert matrix[0, prior.AA_ORDER.index("I")] == pytest.approx(np.log(20*(.99*.286/.999+.0005)))
    assert matrix[3, prior.AA_ORDER.index("G")] == pytest.approx(np.log(20*(.99*.989/.993+.0005)))
    # Published rounded 227K is zero. The declared uniform mixture gives finite mass.
    assert matrix[2, prior.AA_ORDER.index("K")] == pytest.approx(np.log(.01))
    np.testing.assert_allclose(np.exp(matrix).sum(axis=1), 20)


def test_manifest_does_not_invent_alignment_counts_and_is_reproducible(prepared):
    path, metadata = prepared
    assert metadata["raw_MSA_available"] is False
    assert metadata["sequence_weights_available"] is False
    assert metadata["n_effective_sequences"] is None
    assert metadata["count_pseudocount"] is None
    assert metadata["uniform_probability_mixture"] == .01
    assert metadata["conversion"]["input_kind"] == "frequency_table_not_raw_alignment"
    assert metadata["author_alignment_statistics"]["num_seqs"] == "18857"
    for name, digest in metadata["files"].items():
        assert prior.sha(path / name) == digest
    before = prior.sha(path / "prior_manifest.json")
    prior.prepare(path)
    assert prior.sha(path / "prior_manifest.json") == before


def test_published_source_tampering_fails_before_prior_build(prepared):
    path, _ = prepared
    source = path / "source/TARGET_b0.3_frequencies.csv"
    source.write_text(source.read_text().replace("0.286", "0.999", 1))
    with pytest.raises(ValueError, match="Source checksum mismatch"):
        prior.prepare(path)
