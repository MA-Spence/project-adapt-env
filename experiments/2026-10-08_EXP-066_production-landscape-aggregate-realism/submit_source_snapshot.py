from pathlib import Path
import os,subprocess,sys,yaml
root=Path.cwd();exp=Path('experiments/2026-10-08_EXP-066_production-landscape-aggregate-realism')
meta=yaml.safe_load((exp/'metadata.yaml').read_text())
files=set(meta['inputs']['declared_files'])|{'AGENTS.md','AIM.md','PROJECT_STATE.md','labproject.yaml','.gitignore','.gitmodules','.labproj-stageignore','.agents/skills/labproj/SKILL.md','external/Adapt-Env/README.md'}
files|={str(p) for p in exp.iterdir() if p.is_file()};files.add(str(exp/'production.staging-filter.txt'))
trees={'environment','registry',str(exp/'runs')};excludes=[]
def inspect(folder):
 for entry in os.scandir(folder):
  path=Path(entry.path);rel=path.relative_to(root).as_posix()
  if rel in trees or rel in files:continue
  if entry.is_dir(follow_symlinks=False) and any(s.startswith(rel+'/') for s in files|trees):inspect(path)
  else:excludes.append(rel)
inspect(root)
text='# Source-only snapshot for EXP-066; exact literal prefixes.\n# Empirical bundles are separately SHA-256 verified.\n'+'\n'.join(sorted(excludes))+'\n'
filterpath=root/'.labproj-stageignore';original=filterpath.read_bytes()
projectpath=root/'labproject.yaml';project_original=projectpath.read_bytes()
project_config=yaml.safe_load(project_original);project_config['environment']['overlay']['enabled']=False
(root/exp/'production.staging-filter.txt').write_text(text)
print('Temporary staging filter:',len(excludes),'excluded prefixes',flush=True)
try:
 filterpath.write_text(text)
 projectpath.write_text(yaml.safe_dump(project_config,sort_keys=False))
 command=(['labproj','run','resubmit','RUN-172','--execute'] if sys.argv[1]=='--resubmit' else ['labproj','submit','EXP-066','--target','lab-slurm',sys.argv[1]])
 result=subprocess.run(command,timeout=240)
 sys.exit(result.returncode)
finally:
 filterpath.write_bytes(original)
 projectpath.write_bytes(project_original)
 print('Original staging filter restored.',flush=True)
