#!/usr/bin/env python3
"""Prepare an isolated pinned runtime for the authorized joint-readout Slurm RUN."""
import base64
import json
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
TARGET = "/tank/data/labproj/runtime/project-adapt-env/envs/joint-readout-py311"
BASE = "/tank/data/labproj/stacks/STACK-000/bin/python"
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=12", "dutch"]


def main():
    lock = (ROOT / "environment/readout-py311.lock.txt").read_bytes()
    (HERE / "readout-py311.lock.txt").write_bytes(lock)
    script = f'''
import base64,datetime,hashlib,json,pathlib,platform,subprocess,sys
target=pathlib.Path({TARGET!r})
lock_bytes=base64.b64decode({base64.b64encode(lock).decode()!r})
digest=hashlib.sha256(lock_bytes).hexdigest()
evidence=target/'preparation'
marker=evidence/'ownership.json'
if platform.python_version()!='3.11.15': raise RuntimeError('Expected STACK-000 Python 3.11.15')
if target.exists():
    if not marker.is_file() or json.loads(marker.read_text()).get('lock_sha256')!=digest:
        raise RuntimeError('Refusing to alter an unowned or differently pinned environment')
else:
    evidence.mkdir(parents=True)
    marker.write_text(json.dumps(dict(scope='EXP065 joint continuous-r readout fitting',lock_sha256=digest,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())))
if not (target/'bin/python').is_file(): subprocess.run([sys.executable,'-m','venv',str(target)],check=True)
lock=evidence/'readout-py311.lock.txt';lock.write_bytes(lock_bytes)
with (evidence/'installation.log').open('a') as stream:
    p=subprocess.run([str(target/'bin/python'),'-m','pip','install','-r',str(lock)],stdout=stream,stderr=subprocess.STDOUT)
if p.returncode: raise RuntimeError('Installation failed; inspect '+str(evidence/'installation.log'))
print(json.dumps(dict(target=str(target),lock_sha256=digest,installation_returncode=p.returncode)))
'''
    result = subprocess.run([*SSH, shlex.quote(BASE) + " -"], input=script,
                            text=True, capture_output=True)
    (HERE / "provision.stdout.txt").write_text(result.stdout)
    (HERE / "provision.stderr.txt").write_text(result.stderr)
    log = subprocess.run([*SSH, "cat " + shlex.quote(TARGET + "/preparation/installation.log")],
                         capture_output=True, text=True)
    (HERE / "installation.log").write_text(log.stdout + log.stderr)
    if result.returncode:
        raise RuntimeError("Remote environment preparation failed; inspect recorded evidence")
    print(result.stdout, end="")


if __name__ == "__main__":
    main()
