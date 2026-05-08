#!/usr/bin/env bash
#SBATCH --job-name=EXP-001RUN-002
#SBATCH --partition=all
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/tank/data/labproj/runtime/project-adapt-env/runs/RUN-002/logs/stdout.txt
#SBATCH --error=/tank/data/labproj/runtime/project-adapt-env/runs/RUN-002/logs/stderr.txt
set -euo pipefail
if ! command -v module >/dev/null 2>&1; then
  set +u
  source /etc/profile
  set -u
fi
module purge
export LABPROJ_TARGET=lab-slurm
export LABPROJ_CLUSTER_NAME=homelab
export LABPROJ_EXPERIMENT_ID=EXP-001
export LABPROJ_RUN_ID=RUN-002
export LABPROJ_OUTPUT_DIR="/tank/data/labproj/runtime/project-adapt-env/runs/RUN-002/outputs"
export LABPROJ_LOG_DIR="/tank/data/labproj/runtime/project-adapt-env/runs/RUN-002/logs"
export LABPROJ_WORKDIR="${TMPDIR:-/tmp/labproj/project-adapt-env/RUN-002}"

mkdir -p "${LABPROJ_OUTPUT_DIR}" "${LABPROJ_LOG_DIR}" "${LABPROJ_WORKDIR}"

bash "/tank/data/labproj/staging/project-adapt-env/RUN-002/project/experiments/2026-05-07_EXP-001_proteingym-dms-distributional-realism-panel/run.sh"
