#!/usr/bin/env bash
set -euo pipefail
# EXP-067 is a prospective protocol, not an executable benchmark yet.
# Never create placeholder results or report success from this planning scaffold.
printf '%s\n' 'EXP-067 is planned only. Implement and freeze the benchmark harness, counter contract, isolated environment and correctness gates described in experiment.md before creating a RUN.' >&2
exit 64
