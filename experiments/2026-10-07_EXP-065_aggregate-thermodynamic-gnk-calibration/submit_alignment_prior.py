"""Prepare/submit the self-contained TrpB profile repeat through labproj.

The default is a nonexecuting dry-run. --execute requires an explicit invocation;
the source/profile manifests are checked before either operation. Project-local
transport, overlay and staging settings are restored even if submission fails.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys

import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--dry-run', action='store_true')
    action.add_argument('--execute', action='store_true')
    parser.add_argument('--resubmit', metavar='RUN_ID')
    parser.add_argument('--doctor', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    exp = Path(__file__).resolve().parent.relative_to(root)
    source_path = root / exp / 'alignment_prior_source.json'
    source = json.loads(source_path.read_text())
    for relative, expected in source['files'].items():
        actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Frozen source changed: {relative}')
    for bundle in ['interaction-strength-exp065-inputs', 'trpb-profile-additive-inputs']:
        directory = root / 'data/processed' / bundle
        manifest = json.loads((directory / ('manifest.json' if bundle.startswith('interaction') else 'prior_manifest.json')).read_text())
        for relative, expected in manifest['files'].items():
            if isinstance(expected, dict):
                expected = expected['sha256']
            if hashlib.sha256((directory / relative).read_bytes()).hexdigest() != expected:
                raise ValueError(f'Frozen input changed: {bundle}/{relative}')

    files = set(source['files']) | {
        str(exp / 'alignment_prior_source.json'),
        str(exp / 'metadata.yaml'),
        str(exp / 'alignment_prior_validation.json'),
        'AGENTS.md', 'AIM.md', 'PROJECT_STATE.md', 'docs/hypotheses.md',
        'labproject.yaml', '.gitignore', '.gitmodules', '.labproj-stageignore',
        '.agents/skills/labproj/SKILL.md',
        'data/processed/interaction-strength-exp065-inputs.dvc',
        'data/processed/trpb-profile-additive-inputs.dvc',
    }
    trees = {'environment', 'registry', str(exp / 'runs'),
             'data/processed/interaction-strength-exp065-inputs',
             'data/processed/trpb-profile-additive-inputs'}
    excludes = []

    def inspect(folder):
        for entry in os.scandir(folder):
            path = Path(entry.path)
            relative = path.relative_to(root).as_posix()
            if relative in files | trees:
                continue
            if entry.is_dir(follow_symlinks=False) and any(
                    item.startswith(relative + '/') for item in files | trees):
                inspect(path)
            else:
                excludes.append(relative)

    inspect(root)
    filter_text = ('# TrpB alignment-profile repeat; frozen source and both input bundles.\n'
                   + '\n'.join(sorted(excludes)) + '\n')
    (root / exp / 'alignment_prior.staging-filter.txt').write_text(filter_text)
    changing = [root / '.labproj-stageignore', root / 'labproject.yaml', root / '.env']
    originals = {path: path.read_bytes() for path in changing}
    project = yaml.safe_load(originals[changing[1]])
    project['environment']['overlay']['enabled'] = False
    transport = re.sub(rb'(?m)^LAB_SSH_HOST=.*$', b'LAB_SSH_HOST=dutch', originals[changing[2]])
    if not re.search(rb'(?m)^LAB_SSH_HOST=dutch\s*$', transport):
        raise ValueError('Cannot select the recorded Dutch submission transport')
    try:
        changing[0].write_text(filter_text)
        changing[1].write_text(yaml.safe_dump(project, sort_keys=False))
        changing[2].write_bytes(transport)
        if args.doctor:
            command = ['labproj', 'target', 'doctor', 'lab-slurm', '--network']
        else:
            command = (['labproj', 'run', 'resubmit', args.resubmit] if args.resubmit
                       else ['labproj', 'submit', 'EXP-065', '--target', 'lab-slurm'])
            command.append('--execute' if args.execute else '--dry-run')
        print('Self-contained input/source staging; transport=dutch; action='
              + ('execute' if args.execute else 'dry-run'), flush=True)
        result = subprocess.run(command, cwd=root, timeout=300)
        return result.returncode
    finally:
        for path, content in originals.items():
            path.write_bytes(content)
        print('Restored project settings, transport and staging filter.', flush=True)


if __name__ == '__main__':
    sys.exit(main())
