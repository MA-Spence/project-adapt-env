"""Narrow, reversible source staging for the EXP-065 interaction sweep."""
from pathlib import Path
import json,os,re,subprocess,sys,yaml
root=Path.cwd();exp=Path('experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration')
manifest=json.loads((root/exp/'interaction_strength_source.json').read_text())
files=set(manifest['files'])|{'AGENTS.md','AIM.md','PROJECT_STATE.md','labproject.yaml','.gitignore','.gitmodules','.labproj-stageignore','.agents/skills/labproj/SKILL.md','external/Adapt-Env/README.md','data/processed/interaction-strength-exp065-inputs.dvc'}
files|={str(exp/p) for p in ['metadata.yaml','interaction_strength_source.json','interaction_strength_validation.json','submit_interaction_strength.py','interaction_strength.staging-filter.txt']}
trees={'environment','registry',str(exp/'runs')};excludes=[]
def inspect(folder):
 for entry in os.scandir(folder):
  p=Path(entry.path);rel=p.relative_to(root).as_posix()
  if rel in files|trees:continue
  if entry.is_dir(follow_symlinks=False) and any(s.startswith(rel+'/') for s in files|trees):inspect(p)
  else:excludes.append(rel)
inspect(root)
filtertext='# EXP-065 interaction-strength source snapshot; frozen input bundle staged separately.\n'+'\n'.join(sorted(excludes))+'\n'
filterpath=root/'.labproj-stageignore';original_filter=filterpath.read_bytes()
projectpath=root/'labproject.yaml';original_project=projectpath.read_bytes();project=yaml.safe_load(original_project);project['environment']['overlay']['enabled']=False
# labproj reads .env directly, overriding inherited environment variables.
envpath=root/'.env';original_env=envpath.read_bytes();changed_env=re.sub(rb'(?m)^LAB_SSH_HOST=.*$',b'LAB_SSH_HOST=dutch',original_env)
assert changed_env!=original_env or b'LAB_SSH_HOST=dutch' in original_env
(root/exp/'interaction_strength.staging-filter.txt').write_text(filtertext)
try:
 filterpath.write_text(filtertext);projectpath.write_text(yaml.safe_dump(project,sort_keys=False));envpath.write_bytes(changed_env)
 print('Submitting through dutch; source snapshot excludes',len(excludes),'unrelated prefixes.',flush=True)
 command=(['labproj','run','resubmit',sys.argv[2],'--execute'] if sys.argv[1]=='--resubmit' else ['labproj','submit','EXP-065','--target','lab-slurm',sys.argv[1]])
 result=subprocess.run(command,timeout=240)
 sys.exit(result.returncode)
finally:
 filterpath.write_bytes(original_filter);projectpath.write_bytes(original_project);envpath.write_bytes(original_env)
 print('Restored project config, transport setting and staging filter.',flush=True)
