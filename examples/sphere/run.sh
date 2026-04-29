#!/usr/bin/env bash
# Run the cosmo sphere example.
#
# Execute from any directory; the script resolves the repo root automatically.
# Prerequisites: Julia installed and on PATH.
#
# One-time setup (install Julia dependencies):
#   julia --project=cosmo -e 'using Pkg; Pkg.instantiate()'

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
INPUT="${SCRIPT_DIR}/input.json"

cd "${REPO_ROOT}"

julia --project=cosmo cosmo/run.jl "${INPUT}"
