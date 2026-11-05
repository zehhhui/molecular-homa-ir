#!/usr/bin/env bash
# Reduced synthetic workflow only; this does not reproduce manuscript results.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"
PYTHON_BIN="${PYTHON:-${PYTHON_BIN:-python3}}"
RSCRIPT_BIN="${RSCRIPT:-${RSCRIPT_BIN:-Rscript}}"
"$PYTHON_BIN" examples/generate_synthetic_data.py --output-dir data/synthetic
"$RSCRIPT_BIN" scripts/train_models.R --config config/smoke.json
"$RSCRIPT_BIN" scripts/extract_contributions.R --model-run results/smoke --output results/smoke/extracted
"$PYTHON_BIN" scripts/treatment_contributions.py --input results/smoke/extracted --output results/smoke/contributions --bootstrap 100 --min-folds 1 --top-n 10
"$RSCRIPT_BIN" tests/test_model_invariants.R results/smoke
"$PYTHON_BIN" -m unittest discover -s tests -p 'test_*.py' -v
