"""Rollback configuration and dispatch safeguards; no scientific acceptance gates."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import optuna_pairwise_calibration_exp065 as protocol

EXP = protocol.EXP
BASELINE = ROOT / "data/processed/aggregate-thermodynamic-gnk-calibration/RUN-180/outputs"


def configs():
    return (yaml.safe_load((EXP / "optuna_readout.yaml").read_text()),
            yaml.safe_load((BASELINE / "config.yaml").read_text()))


@pytest.mark.parametrize("key", ["fit", "model", "interaction_ratios", "seed_design",
                                     "sampling", "measurement", "observation_modes",
                                     "additive_prior", "required_worker_environment"])
def test_rollback_refuses_changed_scientific_design(key):
    current, old = configs()
    protocol.validate_config(current, old)
    changed = copy.deepcopy(current)
    changed[key] = None
    with pytest.raises(ValueError, match="RUN-180 design changed"):
        protocol.validate_config(changed, old)


def test_baseline_authenticates_saved_parameters_and_config():
    current, _ = configs()
    settings, manifest = protocol.load_baseline(BASELINE, current)
    assert len(settings) == 18
    assert manifest["protocol"] == "trpb_full_pairwise_profile_calibration_v1"
    current["warm_start"]["manifest_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="manifest checksum"):
        protocol.load_baseline(BASELINE, current)


@pytest.mark.parametrize("explicit,runner,config,source", [
    (None, "optuna_joint_pairwise_calibration_exp065.py", "joint_readout.yaml", "joint_readout_source.json"),
    ("optuna-readout", "optuna_pairwise_calibration_exp065.py", "optuna_readout.yaml", "optuna_readout_source.json"),
    ("full-pairwise", "pairwise_calibration_exp065.py", "full_pairwise.yaml", "full_pairwise_source.json"),
])
def test_launcher_default_and_explicit_historical_route(tmp_path, explicit, runner, config, source):
    recorder = tmp_path / "record-argv"
    output = tmp_path / "argv.json"
    recorder.write_text(f"#!{sys.executable}\nimport json, sys\nfrom pathlib import Path\n"
                        f"Path({str(output)!r}).write_text(json.dumps(sys.argv[1:]))\n")
    recorder.chmod(0o755)
    env = dict(os.environ, PYTHON=str(recorder))
    env.pop("ADAPTENV_PROTOCOL", None)
    for key in ["ADAPTENV_DATA_ROOT", "ADAPTENV_PROFILE_ROOT", "ADAPTENV_PROFILE_CONTROL_ROOT",
                "ADAPTENV_BASELINE_ROOT", "LABPROJ_QUICK"]:
        env.pop(key, None)
    if explicit:
        env["ADAPTENV_PROTOCOL"] = explicit
    subprocess.run(["bash", str(EXP / "run.sh")], env=env, check=True)
    args = json.loads(output.read_text())
    assert Path(args[0]).name == runner
    assert Path(args[args.index("--config") + 1]).name == config
    assert Path(args[args.index("--source-manifest") + 1]).name == source
    assert ("--baseline-root" in args) == (explicit in (None, "optuna-readout"))
