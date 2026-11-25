"""Lock/check the unpublished Gamma and indirect-inference implementation."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFEST = HERE / "manifest.json"


def paths():
    package = ROOT / "external/Adapt-Env"
    result = [package / p for p in ("pyproject.toml", "README.md", "LICENSE")]
    for name, suffixes in (("adaptenv", {".py", ".yaml", ".json", ".md"}),
                           ("tests", {".py", ".json"}), ("docs", {".md"})):
        result.extend(p for p in (package / name).rglob("*")
                      if p.is_file() and p.suffix in suffixes and "__pycache__" not in p.parts)
    result.extend((ROOT / p) for p in ("scripts/calibrate_architecture.py", "tests/test_calibrate_architecture.py",
                  "environment/calibration-py311.lock.txt", "environment/calibration-py313.lock.txt"))
    result.extend(p for p in (ROOT / "docs/calibration").glob("*") if p.is_file())
    result.extend(p for p in HERE.iterdir() if p.is_file() and p.name != "manifest.json")
    return sorted(set(result))


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    files = {str(p.relative_to(ROOT)): {"sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                                      "bytes": p.stat().st_size} for p in paths()}
    if args.check:
        if json.loads(MANIFEST.read_text())["files"] != files:
            raise SystemExit("Candidate source or verification record differs from the lock")
        print(f"Verified {len(files)} files; candidate remains unpublished and not empirically calibrated")
    else:
        if MANIFEST.exists():
            raise SystemExit("Candidate already locked; do not overwrite historical manifests")
        MANIFEST.write_text(json.dumps(dict(version="0.7.0rc1", status="unpublished_implementation_candidate",
             calibration_status="new_architecture_objective_not_empirically_calibrated",
             model="thermodynamic_pairwise_v3", default_additive_prior="profile",
             default_interaction_concentration=None, pair_orders=[2],
             method="thornton_auxiliary_optuna_indirect_inference_v1", files=files), indent=2) + "\n")
        print(f"Locked {len(files)} files")


if __name__ == "__main__":
    main()
