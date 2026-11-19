#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT_DIR="${LABPROJ_OUTPUT_DIR:-${SCRIPT_DIR}/outputs/local}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/adaptenv-mpl-${UID}}"
case "${ADAPTENV_PROTOCOL:-alignment-prior}" in
  alignment-prior)
    RUNNER="alignment_prior_exp065.py"
    CONFIG="alignment_prior.yaml"
    BUNDLE="interaction-strength-exp065-inputs"
    REQUIRED="manifest.json"
    if [[ -x /tank/data/labproj/runtime/project-adapt-env/tools/graph-compiler/ziglang/zig ]]; then
      export ADAPTENV_GRAPH_COMPILER=/tank/data/labproj/runtime/project-adapt-env/tools/graph-compiler/ziglang/zig
    fi
    ;;
  interaction-strength)
    RUNNER="interaction_strength_exp065.py"
    CONFIG="interaction_strength.yaml"
    BUNDLE="interaction-strength-exp065-inputs"
    REQUIRED="manifest.json"
    if [[ -x /tank/data/labproj/runtime/project-adapt-env/tools/graph-compiler/ziglang/zig ]]; then
      export ADAPTENV_GRAPH_COMPILER=/tank/data/labproj/runtime/project-adapt-env/tools/graph-compiler/ziglang/zig
    fi
    ;;
  initial)
    RUNNER="aggregate_thermodynamic_calibration.py"
    CONFIG="config.yaml"
    BUNDLE="ensemble-landscape-validation-datasets"
    REQUIRED="datasets"
    ;;
  paired-readout)
    RUNNER="paired_readout_exp065.py"
    CONFIG="paired_readout.yaml"
    BUNDLE="ensemble-landscape-validation-datasets"
    REQUIRED="datasets"
    ;;
  observation-calibration)
    RUNNER="observation_calibration_exp065.py"
    CONFIG="observation_calibration.yaml"
    BUNDLE="trpb-replicate-observation-exp065"
    REQUIRED="input_manifest.json"
    ;;
  replicate-observation)
    RUNNER="replicate_observation_exp065.py"
    CONFIG="replicate_observation.yaml"
    BUNDLE="trpb-replicate-observation-exp065"
    REQUIRED="input_manifest.json"
    ;;
  *) echo "Unknown EXP-065 protocol: ${ADAPTENV_PROTOCOL}" >&2; exit 1 ;;
esac
DATA_ROOT="${ADAPTENV_DATA_ROOT:-${PROJECT_ROOT}/data/processed/${BUNDLE}}"
SHARED_DATA_ROOT="/tank/data/labproj/runtime/project-adapt-env/shared/${BUNDLE}"
if [[ ! -e "${DATA_ROOT}/${REQUIRED}" && -z "${ADAPTENV_DATA_ROOT:-}" && -e "${SHARED_DATA_ROOT}/${REQUIRED}" ]]; then
  DATA_ROOT="${SHARED_DATA_ROOT}"
fi
if [[ ! -e "${DATA_ROOT}/${REQUIRED}" ]]; then
  echo "Empirical input bundle not found: ${DATA_ROOT}" >&2
  exit 1
fi
set --
if [[ "${ADAPTENV_PROTOCOL:-alignment-prior}" == "alignment-prior" ]]; then
  PRIOR_ROOT="${ADAPTENV_PROFILE_ROOT:-${PROJECT_ROOT}/data/processed/trpb-profile-additive-inputs}"
  if [[ ! -f "${PRIOR_ROOT}/additive_matrix.npy" || ! -f "${PRIOR_ROOT}/prior_manifest.json" ]]; then
    echo "Frozen alignment prior not found: ${PRIOR_ROOT}" >&2
    exit 1
  fi
  set -- --prior-matrix "${PRIOR_ROOT}/additive_matrix.npy" \
    --prior-metadata "${PRIOR_ROOT}/prior_manifest.json" \
    --source-manifest "${SCRIPT_DIR}/alignment_prior_source.json"
fi
if [[ "${LABPROJ_QUICK:-0}" == "1" ]]; then set -- "$@" --quick; fi
exec "${PYTHON:-python3}" "${PROJECT_ROOT}/scripts/${RUNNER}" \
  --config "${SCRIPT_DIR}/${CONFIG}" --output-dir "${OUTPUT_DIR}" --data-root "${DATA_ROOT}" "$@"
