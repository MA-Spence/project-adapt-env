---
name: labproj
description: Use the labproj CLI to inspect the scientific record, create records, track durable data with DVC, and submit or reconcile runs for this repository.
---

# labproj

Use this skill when working inside a labproj scientific project repository.

## First Steps

Before doing anything substantive:

1. Read `AGENTS.md`.
2. Read `AIM.md`, `docs/hypotheses.md`, relevant `results/*/result.md`,
   relevant `claims/*`, and `PROJECT_STATE.md`.
3. Read the `README.md` files under `external/` for any package relevant to
   the task.
4. Check the current project state with:

```bash
labproj status
labproj refresh --write --strict
labproj check --strict
```

Do not suggest new directions unless explicitly asked.

## Project Setup and Configuration

Create or refresh a project bootstrap:

```bash
labproj setup . --infra-config /path/to/lab-infra-config
labproj env bootstrap
labproj env apply-infra
labproj doctor local --network
```

## Scientific Record Commands

Create records only when explicitly asked:

```bash
labproj new aim --yes
labproj new hypothesis --yes
labproj new experiment --yes
labproj new analysis --yes
labproj new result --yes
labproj new claim --yes
```

Update status only when the milestone is real:

```bash
labproj set-status RUN-001 completed
labproj set-status ANA-001 completed
labproj set-status RES-001 completed
```

## Data and DVC

Durable data should usually live under `data/processed`.

```bash
labproj data configure-backup
labproj data status --detailed
labproj data track data/processed/example-output
labproj data push
labproj data restore data/processed/example-output --pull
```

The `track`, `push`, and `restore` commands accept either a file or a
directory.

## External Packages

If the project uses an independent source package, keep it independent:

```bash
labproj external add https://github.com/owner/package.git --editable --sync
```

Treat code in `external/` as its own package history. Do not collapse it into
the scientific record.

## Run Submission

When asked to execute an experiment:

```bash
labproj target doctor lab-slurm --network
labproj submit EXP-001 --target lab-slurm --dry-run
labproj submit EXP-001 --target lab-slurm --execute
labproj run reconcile
```

Use dry-runs first when changing execution config.

## Analysis and Reproducibility

- Keep `PROJECT_STATE.md` current if you materially change the project.
- Keep claims, results, and analyses aligned with actual evidence.
- Flag scope drift when requests move beyond the registered aims or hypotheses.
