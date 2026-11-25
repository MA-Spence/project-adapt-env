#!/usr/bin/env python3
"""Prepare only the isolated Dutch runtime for the authorized RUN-182 migration.

No job submission, empirical fitting, service restart or node-state mutation.
The source protocol, data and exact dependency lock come from RUN-182's snapshot.
"""
from __future__ import annotations

import base64
import hashlib
import io
import json
from pathlib import Path
import shlex
import subprocess
import tarfile


ROOT = Path(__file__).resolve().parents[3]
EVIDENCE = Path(__file__).resolve().parent
SNAPSHOT = ROOT / ".labproj/local-preparations/architecture-20261121T013859Z-fxagfbej"
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", "matthew-spence@dutch"]
BASE = "/tank/data/labproj/stacks/STACK-000/bin/python"
TARGET = "/tank/data/labproj/runtime/project-adapt-env/envs/architecture-0.7.0rc1-py311"
EXP_REL = "experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration"


def remote_python(script, *, interpreter=BASE):
    compile(script, "<remote-environment-script>", "exec")
    result = subprocess.run([*SSH, f"{shlex.quote(interpreter)} -"], input=script,
                            text=True, capture_output=True, check=False)
    return result


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def provision():
    lock_bytes = (SNAPSHOT / "source/environment/calibration-py311.lock.txt").read_bytes()
    (EVIDENCE / "calibration-py311.lock.txt").write_bytes(lock_bytes)
    encoded = base64.b64encode(lock_bytes).decode()
    script = f'''
import base64, datetime, hashlib, json, pathlib, platform, subprocess, sys
target=pathlib.Path({TARGET!r})
lock_bytes=base64.b64decode({encoded!r})
lock_hash=hashlib.sha256(lock_bytes).hexdigest()
evidence=target/'architecture_environment_evidence'
marker=evidence/'ownership.json'
if platform.python_version()!='3.11.15':
    raise RuntimeError('Expected existing STACK-000 Python 3.11.15')
if target.exists():
    if not marker.is_file() or json.loads(marker.read_text()).get('lock_sha256')!=lock_hash:
        raise RuntimeError('Refusing to alter an existing unowned or differently pinned environment')
else:
    evidence.mkdir(parents=True)
    marker.write_text(json.dumps(dict(scope='RUN-182 unchanged protocol migration to Slurm',
        lock_sha256=lock_hash, created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2)+chr(10))
if not (target/'bin/python').is_file():
    subprocess.run([sys.executable,'-m','venv',str(target)],check=True)
lock=evidence/'calibration-py311.lock.txt'
lock.write_bytes(lock_bytes)
with (evidence/'installation.log').open('a') as stream:
    process=subprocess.run([str(target/'bin/python'),'-m','pip','install','--requirement',str(lock)],
                           stdout=stream,stderr=subprocess.STDOUT,check=False)
if process.returncode:
    raise RuntimeError('Pinned installation failed; see '+str(evidence/'installation.log'))
print(json.dumps(dict(target=str(target),lock_sha256=lock_hash,installation_returncode=process.returncode)))
'''
    result = remote_python(script)
    (EVIDENCE / "provision.stdout.txt").write_text(result.stdout)
    (EVIDENCE / "provision.stderr.txt").write_text(result.stderr)
    log = subprocess.run([*SSH, "cat " + shlex.quote(TARGET + "/architecture_environment_evidence/installation.log")],
                         capture_output=True, text=True, check=False)
    (EVIDENCE / "installation.log").write_text(log.stdout + log.stderr)
    if result.returncode:
        raise RuntimeError(f"Remote provisioning failed ({result.returncode}); inspect local evidence logs")
    return lock_bytes


def stage_preflight():
    package = SNAPSHOT / "source/external/Adapt-Env"
    paths = [path for path in package.rglob("*") if path.is_file()]
    paths += [SNAPSHOT / "source/scripts/calibrate_architecture.py",
              SNAPSHOT / f"source/{EXP_REL}/architecture_calibration.json",
              SNAPSHOT / "inputs/trpb-architecture-calibration-inputs/calibration.npz",
              SNAPSHOT / "inputs/trpb-architecture-calibration-inputs/confirmation.npz"]
    files = {str(path.relative_to(SNAPSHOT)): digest(path.read_bytes()) for path in sorted(paths)}
    manifest = dict(source_snapshot=str(SNAPSHOT.relative_to(ROOT)), files=files,
                    scope="Small structural preflight copy only; complete worker snapshot is staged separately")
    (EVIDENCE / "preflight_sources.json").write_text(json.dumps(manifest, indent=2) + "\n")
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as archive:
        for path in sorted(paths):
            archive.add(path, arcname=str(path.relative_to(SNAPSHOT)), recursive=False)
        record = json.dumps(manifest, indent=2).encode()
        info = tarfile.TarInfo("preflight_sources.json")
        info.size = len(record)
        archive.addfile(info, io.BytesIO(record))
    extract = f'''
import json,pathlib,shutil,sys,tarfile
root=pathlib.Path({(TARGET + '/architecture_environment_evidence/preflight_snapshot')!r})
root.mkdir(parents=True,exist_ok=True)
with tarfile.open(fileobj=sys.stdin.buffer,mode='r|gz') as archive:
    for member in archive:
        destination=(root/member.name).resolve()
        if not member.isfile() or not destination.is_relative_to(root.resolve()):
            raise RuntimeError('Invalid structural preflight archive member')
        destination.parent.mkdir(parents=True,exist_ok=True)
        with archive.extractfile(member) as source,destination.open('wb') as target:
            shutil.copyfileobj(source,target)
print(json.dumps(dict(staged=str(root))))
'''
    command = [*SSH, shlex.quote(TARGET + "/bin/python") + " -c " + shlex.quote(extract)]
    result = subprocess.run(command, input=payload.getvalue(), capture_output=True, check=False)
    (EVIDENCE / "staging.stdout.txt").write_bytes(result.stdout)
    (EVIDENCE / "staging.stderr.txt").write_bytes(result.stderr)
    if result.returncode:
        raise RuntimeError("Remote structural preflight staging failed")


def capture(lock_bytes):
    expected = dict(line.split("==", 1) for line in lock_bytes.decode().splitlines()
                    if line.strip() and not line.startswith("#"))
    script = f'''
import contextlib,hashlib,importlib.metadata as metadata,io,json,os,pathlib,platform,re,socket,subprocess,sys
expected={expected!r}
target=pathlib.Path({TARGET!r})
evidence=target/'architecture_environment_evidence'
snapshot=evidence/'preflight_snapshot'
manifest=json.loads((snapshot/'preflight_sources.json').read_text())
for relative,expected_hash in manifest['files'].items():
    actual_hash=hashlib.sha256((snapshot/relative).read_bytes()).hexdigest()
    if actual_hash!=expected_hash: raise RuntimeError('Preflight source checksum mismatch: '+relative)
actual={{d.metadata['Name']:d.version for d in metadata.distributions()}}
normalized={{re.sub(r'[-_.]+','-',key).lower():value for key,value in actual.items()}}
mismatch={{name:dict(expected=value,actual=normalized.get(re.sub(r'[-_.]+','-',name).lower()))
          for name,value in expected.items() if normalized.get(re.sub(r'[-_.]+','-',name).lower())!=value}}
if platform.python_version()!='3.11.15' or mismatch:
    raise RuntimeError('Pinned runtime mismatch: '+str(mismatch))
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
sys.dont_write_bytecode=True
sys.path.insert(0,str(snapshot/'source/external/Adapt-Env'))
import numpy,scipy,optuna
from adaptenv import __version__
from adaptenv.inference import RFAConfig,RFADesign,fit_reference_free
from adaptenv.calibration import CalibrationConfig,CalibrationData,calibrate,build_calibrated_landscape
configuration=RFAConfig()
design=RFADesign.from_sequences(numpy.array([[0,0],[0,1],[1,0],[1,1]],dtype=numpy.int64))
blas=io.StringIO()
with contextlib.redirect_stdout(blas): numpy.show_config()
check=subprocess.run([sys.executable,'-m','pip','check'],capture_output=True,text=True,check=False)
if check.returncode: raise RuntimeError('pip check failed: '+check.stdout+check.stderr)
command=[sys.executable,str(snapshot/'source/scripts/calibrate_architecture.py'),
         '--input',str(snapshot/'inputs/trpb-architecture-calibration-inputs/calibration.npz'),
         '--confirmation',str(snapshot/'inputs/trpb-architecture-calibration-inputs/confirmation.npz'),
         '--config',str(snapshot/'source/{EXP_REL}/architecture_calibration.json'),'--validate-only']
validation=subprocess.run(command,env=os.environ.copy(),capture_output=True,text=True,check=False)
if validation.returncode: raise RuntimeError('Structural preflight failed: '+validation.stderr)
record=dict(status='prepared_and_structurally_valid',python=platform.python_version(),
    executable=sys.executable,base_prefix=sys.base_prefix,prefix=sys.prefix,platform=platform.platform(),
    machine=platform.machine(),hostname=socket.gethostname(),python_build=platform.python_build(),
    packages=actual,pinned_packages=expected,pinned_mismatches=mismatch,
    environment_lock_sha256=hashlib.sha256((evidence/'calibration-py311.lock.txt').read_bytes()).hexdigest(),
    numpy_configuration=blas.getvalue(),pip_check=dict(returncode=check.returncode,stdout=check.stdout,stderr=check.stderr),
    package_version=__version__,package_source_root=str(snapshot/'source/external/Adapt-Env'),
    structural_design=dict(shape=list(design.matrix.shape),nonzeros=int(design.matrix.nnz)),
    structural_input_validation=json.loads(validation.stdout),source_files_verified=len(manifest['files']),
    empirical_fitting_performed=False,job_submitted=False,
    interpretation='Same Python and pinned dependency versions as RUN-182; Linux architecture/BLAS differ from the MacBook and cross-platform bitwise equivalence is not asserted.')
(evidence/'environment.json').write_text(json.dumps(record,indent=2)+chr(10))
print(json.dumps(record,indent=2))
'''
    result = remote_python(script, interpreter=TARGET + "/bin/python")
    (EVIDENCE / "capture.stderr.txt").write_text(result.stderr)
    if result.returncode:
        (EVIDENCE / "capture.stdout.txt").write_text(result.stdout)
        raise RuntimeError("Remote environment validation failed; inspect capture logs")
    record = json.loads(result.stdout)
    (EVIDENCE / "environment.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    lock_bytes = provision()
    stage_preflight()
    record = capture(lock_bytes)
    files = {path.name: digest(path.read_bytes()) for path in sorted(EVIDENCE.iterdir())
             if path.is_file() and path.name != "manifest.json"}
    (EVIDENCE / "manifest.json").write_text(json.dumps(dict(files=files, remote_environment=TARGET,
        source_run="RUN-182", source_snapshot=str(SNAPSHOT.relative_to(ROOT))), indent=2) + "\n")
    print(json.dumps({key: record[key] for key in ("status", "python", "platform", "package_version",
                     "source_files_verified", "empirical_fitting_performed", "job_submitted")}, indent=2))


if __name__ == "__main__":
    main()
