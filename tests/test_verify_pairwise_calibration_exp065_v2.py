"""Keep gzip-provenance relaxation narrower than scientific-data equality."""
import gzip
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pairwise_validator_v2", ROOT/"scripts/verify_pairwise_calibration_exp065_v2.py")
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


def donor_fixture(root, payload=b"genotype,split\nAAAA,calibration\n", mtime=1):
    directory = root/"noise"
    directory.mkdir(parents=True)
    path = directory/"donor_assignments.csv.gz"
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="donor_assignments.csv", mode="wb", fileobj=raw, mtime=mtime) as handle:
            handle.write(payload)
    model = {"membership_sha256": validator.sha(path), "depth_breaks": [1., 2.],
             "conditioning": {"method": "frozen"}}
    (directory/"noise_model.json").write_text(json.dumps(model))
    return root


def test_different_gzip_timestamp_accepts_identical_decompressed_bytes(tmp_path):
    left = donor_fixture(tmp_path/"left", mtime=1)
    right = donor_fixture(tmp_path/"right", mtime=2)
    result = validator.verify_membership_amendment(left, right)
    assert result["candidate"]["compressed_sha256"] != result["control"]["compressed_sha256"]
    assert result["decompressed_bytes_equal"]
    assert result["other_noise_model_fields_equal"]


@pytest.mark.parametrize("field", ["missing", "incorrect"])
def test_each_compressed_digest_is_mandatory_and_verified(tmp_path, field):
    left, right = donor_fixture(tmp_path/"left"), donor_fixture(tmp_path/"right")
    path = left/"noise/noise_model.json"
    model = json.loads(path.read_text())
    if field == "missing":
        del model["membership_sha256"]
    else:
        model["membership_sha256"] = "0"*64
    path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match="own compressed"):
        validator.verify_membership_amendment(left, right)


def test_changed_donor_data_is_rejected_even_with_correct_own_digest(tmp_path):
    left = donor_fixture(tmp_path/"left")
    right = donor_fixture(tmp_path/"right", payload=b"genotype,split\nAAAA,confirmation\n")
    with pytest.raises(ValueError, match="Decompressed donor membership"):
        validator.verify_membership_amendment(left, right)


def test_changed_non_digest_noise_metadata_is_rejected(tmp_path):
    left, right = donor_fixture(tmp_path/"left"), donor_fixture(tmp_path/"right")
    path = right/"noise/noise_model.json"
    model = json.loads(path.read_text())
    model["depth_breaks"][0] = 3.
    path.write_text(json.dumps(model))
    with pytest.raises(ValueError, match="other than membership_sha256"):
        validator.verify_membership_amendment(left, right)


def test_frozen_v1_hash_is_pinned_and_tampering_rejected(tmp_path, monkeypatch):
    assert validator.load_v1().RTOL == 1e-10
    changed = tmp_path/"changed.py"
    changed.write_text(validator.V1_PATH.read_text()+"\n# changed\n")
    monkeypatch.setattr(validator, "V1_PATH", changed)
    with pytest.raises(ValueError, match="frozen v1 validator changed"):
        validator.load_v1()
