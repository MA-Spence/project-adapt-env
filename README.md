# Project Adapt-Env

This repository is the scientific project workspace for this study.

It is intended to hold:

- the scientific record and lineage registries
- project code, notebooks, and small configuration files
- DVC pointer files for durable data and promoted artifacts

It is not intended to hold:

- machine-specific credentials in `.env`
- the only copy of important data in scratch or runtime paths
- lab-wide execution policy that belongs in `lab-infra-config`

## Current Model

The current `labproj` model is:

- the project repository lives in Git
- infrastructure policy is loaded from a separate `lab-infra-config`
- durable data is backed up with DVC
- remote jobs can be submitted to configured targets such as Slurm
- `data/processed` is the default durable data root

## Key Configuration

- Project ID: `project-adapt-env`
- Default target: `lab-slurm`
- Default stack: `STACK-000`
- Infra config source: `../lab-infra-config-homelab`

## Getting Started

Initialize or refresh the local machine configuration from the selected lab infra config:

```bash
labproj env apply-infra
```

If this project needs to be repointed to a different infra config directory or target YAML:

```bash
labproj env apply-infra /path/to/lab-infra-config
labproj env apply-infra /path/to/lab-infra-config/targets/lab-slurm.yaml
```

Configure durable data backup for this project:

```bash
labproj data configure-backup
labproj data status
```

Check the project state:

```bash
labproj status
labproj refresh --write --strict
labproj check --strict
```

## Common Workflow

Create the first scientific records:

```bash
labproj new aim --yes
labproj new hypothesis --yes
labproj new experiment --yes
```

Track durable project data and back it up with DVC:

```bash
labproj data track data/processed/example/output.tsv
labproj data push
```

Prepare and submit runs:

```bash
labproj target doctor lab-slurm --network
labproj submit EXP-001 --target lab-slurm --dry-run
labproj submit EXP-001 --target lab-slurm --execute
labproj run reconcile
```

Update statuses when a record has reached a real milestone:

```bash
labproj set-status RUN-001 completed
labproj set-status ANA-001 completed
```

## Repository Layout

- `registry/`: machine-readable lineage and record state
- `experiments/`: experiment definitions and RUN records
- `analyses/`: notebook-backed analysis records
- `results/`: interpreted result records
- `claims/`: manuscript-facing claims
- `data/`: durable project data, typically backed up with DVC
- `scratch_space/`: exploratory, non-record work only
- `environment/`: lightweight project overlay dependencies

## Working Rules

- Keep `.env` local and untracked.
- Treat `data/processed` as durable by default.
- Track promoted data with DVC rather than Git.
- Promote important outputs out of scratch/runtime into tracked project locations.
- Keep infrastructure changes in `lab-infra-config`, not hardcoded into the project.
