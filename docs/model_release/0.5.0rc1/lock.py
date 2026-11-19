"""Freeze or verify the unpublished candidate's source, data and evidence hashes."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MANIFEST = HERE / 'manifest.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paths():
    package = ROOT / 'external/Adapt-Env'
    result = [package / name for name in ['pyproject.toml', 'README.md', 'LICENSE']]
    for subdir, suffixes in [('adaptenv', {'.py', '.yaml', '.json'}),
                            ('tests', {'.py'}), ('docs', {'.md'})]:
        result.extend(p for p in (package / subdir).rglob('*')
                      if p.is_file() and p.suffix in suffixes and '__pycache__' not in p.parts)
    result.append(ROOT / 'src/project_adapt_env/__init__.py')
    result.extend(HERE / name for name in ['decision.md', 'replay.py', 'replay.json',
                                         'verification.json', 'lock.py'])
    return sorted(set(result))


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        manifest = json.loads(MANIFEST.read_text())
        expected = manifest['files']
        assert {str(p.relative_to(ROOT)) for p in paths()} == set(expected), 'Candidate file set changed'
        for relative, detail in expected.items():
            path = ROOT / relative
            assert sha(path) == detail['sha256'], f'Candidate changed: {relative}'
        print(f"Verified {len(expected)} locked candidate files; version {manifest['version']} remains unpublished")
        return
    if MANIFEST.exists():
        raise SystemExit('A candidate lock already exists. Review changes and intentionally remove/update it before freezing again.')
    manifest = {'version': '0.5.0rc1', 'status': 'unpublished_release_candidate',
        'model_family': 'thermodynamic_gnk_v2', 'default_additive_prior': 'profile',
        'default_requires_supplied_profile': True,
        'synthetic_control_requires_explicit_selection': True,
        'evidence': ['EXP-065', 'RUN-178', 'RUN-175', 'ANA-070'],
        'no_release_tag_deployment_or_publication': True,
        'files': {str(p.relative_to(ROOT)): {'sha256': sha(p), 'bytes': p.stat().st_size} for p in paths()}}
    MANIFEST.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"Locked {len(manifest['files'])} candidate files")


if __name__ == '__main__':
    main()
