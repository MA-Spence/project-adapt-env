#!/usr/bin/env bash
set -euo pipefail
export LABPROJ_RUN_ATTEMPT=1
export LABPROJ_RESUME_RUN=0

export LABPROJ_TARGET=local
export LABPROJ_EXPERIMENT_ID=EXP-065
export LABPROJ_RUN_ID=RUN-179
export LABPROJ_OUTPUT_DIR="/Users/matthewspence/Documents/Documents/project-adapt-env/experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/outputs/RUN-179"
export LABPROJ_WORKDIR="/Users/matthewspence/Documents/Documents/project-adapt-env/.labproj/work/RUN-179"

mkdir -p "${LABPROJ_OUTPUT_DIR}" "${LABPROJ_WORKDIR}"

bash "/Users/matthewspence/Documents/Documents/project-adapt-env/experiments/2026-10-07_EXP-065_aggregate-thermodynamic-gnk-calibration/run.sh"
