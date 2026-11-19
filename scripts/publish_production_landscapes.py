#!/usr/bin/env python3
"""Publish a completed EXP-066 output bundle atomically into processed data."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


def publish(source, destination):
    manifest=json.loads((source/'manifest.json').read_text())
    if manifest['mode']!='PRODUCTION' or manifest['execution_errors']:
        raise ValueError('Only a complete production execution can be published')
    def verify(root):
        for name,expected in manifest['files'].items():
            if hashlib.sha256((root/name).read_bytes()).hexdigest()!=expected:
                raise ValueError(f'Output checksum mismatch: {name}')
    verify(source)
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists():
        if (destination/'manifest.json').read_bytes()==(source/'manifest.json').read_bytes():
            verify(destination);return
        raise FileExistsError(f'Preserving earlier production bundle: {destination}')
    temporary=Path(tempfile.mkdtemp(prefix=destination.name+'.publishing-',dir=destination.parent))
    try:
        shutil.copytree(source,temporary,dirs_exist_ok=True)
        verify(temporary)
        os.replace(temporary,destination)
    except BaseException:
        shutil.rmtree(temporary)
        raise
    print(f'PUBLISHED {destination}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--destination',type=Path,required=True)
    a=p.parse_args();publish(a.source,a.destination)
