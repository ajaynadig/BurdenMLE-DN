#!/usr/bin/env bash
set -euo pipefail
ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# No work is launched without --run or an explicit legacy RUN_* selection.
exec "${PYTHON_BIN:-python3}" "${ANALYSIS_DIR}/workflow/run.py" "$@"
