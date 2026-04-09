#!/usr/bin/env bash
# Run the mosdepth coverage summary R script.
# Execute from the project root directory.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

Rscript "${SCRIPT_DIR}/summarize_mosdepth.R" "$@"
