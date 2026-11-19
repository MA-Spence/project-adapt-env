#!/usr/bin/env python3
"""Collect the durable EXP-066 bundle into its local processed-data namespace."""
import argparse
from pathlib import Path
import subprocess
from publish_production_landscapes import publish
p=argparse.ArgumentParser();p.add_argument('run');a=p.parse_args()
if not a.run.startswith('RUN-') or not a.run[4:].isdigit():raise ValueError('Expected RUN-NNN')
root=Path(__file__).resolve().parents[1]/'data/processed/production_landscapes'
staging=root/(a.run+'.collecting');staging.mkdir(parents=True,exist_ok=True)
remote=f'matthew-spence@raspi-1:/tank/data/labproj/projects/project-adapt-env/data/processed/production_landscapes/{a.run}/'
subprocess.run(['rsync','-a',remote,str(staging)+'/'],check=True)
publish(staging,root/a.run)
print(f'Verified local bundle: {root/a.run}')
