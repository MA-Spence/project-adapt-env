# Example chemical biology computational project

This repository is a project-facing example of the `labproj` scientific record
model for a computational chemical biology / protein engineering workflow.

It is intended to show what a student-facing scientific project repo should
look like when:

- the scientific record lives in the project repository
- infrastructure policy lives in a separate `lab-infra-config`
- durable data is backed up with DVC
- external source packages remain independently installable

## Scientific Record Model

This project follows the lineage chain:

`AIM -> HYP -> EXP -> RUN -> ANA -> RES -> CLAIM`

Each layer should correspond to a real scientific object, not a convenience
label.

## What Belongs Here

- aims, hypotheses, experiments, analyses, results, and claims
- project code, notebooks, and small configuration files
- DVC pointer files for durable data and promoted artifacts

## What Does Not Belong Here

- machine-specific credentials in `.env`
- the only copy of important evidence in scratch or runtime paths
- lab-wide execution policy that belongs in `lab-infra-config`

## Project Configuration

- Project ID: `example-chemical-biology-project`
- Default target: `local`
- Default stack: `STACK-000`
- Infra config source: `../lab-infra-config`

## Scientific and Agent Guidance

Agent-facing operating rules live in [AGENTS.md](/Users/matthewspence/Documents/Documents/labproj_collection/scientific-project-repo/AGENTS.md).
Project-local agent tooling guidance lives in [.agents/skills/labproj/SKILL.md](/Users/matthewspence/Documents/Documents/labproj_collection/scientific-project-repo/.agents/skills/labproj/SKILL.md).

## Repository Layout

- `registry/`: machine-readable lineage and record state
- `experiments/`: experiment definitions and RUN records
- `analyses/`: notebook-backed analysis records
- `results/`: interpreted result records
- `claims/`: manuscript-facing claims
- `data/`: durable project data, typically backed up with DVC
- `external/`: external source packages that remain independently installable
- `scratch_space/`: exploratory, non-record work only
- `environment/`: lightweight project overlay dependencies

## Data and Reproducibility

- Keep `.env` local and untracked.
- Treat `data/processed` as durable by default.
- Track promoted data with DVC rather than Git.
- Promote important outputs out of scratch/runtime into tracked project locations.
- Keep infrastructure changes in `lab-infra-config`, not hardcoded into the project.
