"""Freeze/check the unpublished full-pairwise candidate and its engineering record."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFEST = HERE / "manifest.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paths():
    package = ROOT / "external/Adapt-Env"
    result = [package / name for name in ["pyproject.toml", "README.md", "LICENSE"]]
    for subdir, suffixes in [("adaptenv", {".py", ".yaml", ".json", ".md"}),
                            ("tests", {".py", ".json"}), ("docs", {".md"})]:
        result.extend(p for p in (package / subdir).rglob("*")
                      if p.is_file() and p.suffix in suffixes and "__pycache__" not in p.parts)
    result.append(ROOT / "src/project_adapt_env/__init__.py")
    result.extend(HERE / name for name in ["decision.md", "tractability.py", "tractability.json",
                                         "verification.json", "lock.py"])
    return sorted(set(result))


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        manifest = json.loads(MANIFEST.read_text())
        expected = manifest["files"]
        if {str(p.relative_to(ROOT)) for p in paths()} != set(expected):
            raise SystemExit("Candidate file set changed")
        for relative, detail in expected.items():
            if sha(ROOT / relative) != detail["sha256"]:
                raise SystemExit(f"Candidate changed: {relative}")
        print(f"Verified {len(expected)} candidate files; 0.6.0rc1 remains unpublished and uncalibrated")
        return
    if MANIFEST.exists():
        raise SystemExit("Candidate already locked; review changes before deliberately replacing the lock")
    manifest = {
        "version": "0.6.0rc1", "status": "unpublished_release_candidate",
        "model_family": "thermodynamic_pairwise_v3", "default_additive_prior": "profile",
        "interaction_backend": "full_pairwise", "orders": [2], "target_degrees": [6.0],
        "calibration_status": "uncalibrated", "supersedes_default": "0.5.0rc1",
        "historical_backend": "ranked_gnk",
        "historical_source_commit": "745c00b3a369da33e7d7c7b9850ccbc212fac714",
        "no_release_tag_deployment_or_publication": True,
        "files": {str(p.relative_to(ROOT)): {"sha256": sha(p), "bytes": p.stat().st_size}
                  for p in paths()},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Locked {len(manifest['files'])} candidate files")


if __name__ == "__main__":
    main()
