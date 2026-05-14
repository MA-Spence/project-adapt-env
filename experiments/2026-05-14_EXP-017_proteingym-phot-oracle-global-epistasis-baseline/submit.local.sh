#!/usr/bin/env bash
set -euo pipefail

export LABPROJ_TARGET=local
export LABPROJ_EXPERIMENT_ID=EXP-017
export LABPROJ_RUN_ID=RUN-TEMPLATE
export LABPROJ_OUTPUT_DIR="/Users/matthewspence/Documents/Documents/project-adapt-env/experiments/2026-05-14_EXP-017_proteingym-phot-oracle-global-epistasis-baseline/outputs/RUN-TEMPLATE"
export LABPROJ_WORKDIR="/Users/matthewspence/Documents/Documents/project-adapt-env/.labproj/work/RUN-TEMPLATE"

mkdir -p "${LABPROJ_OUTPUT_DIR}" "${LABPROJ_WORKDIR}"

bash "/Users/matthewspence/Documents/Documents/project-adapt-env/experiments/2026-05-14_EXP-017_proteingym-phot-oracle-global-epistasis-baseline/run.sh"
