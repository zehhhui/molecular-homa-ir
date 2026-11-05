#!/usr/bin/env python3
"""Run the complete synthetic demonstration with Code Ocean filesystem rules.

No inputs are downloaded, no packages are installed and no study data are read.
The source tree is treated as read-only. Set CAPSULE_RESULTS_DIR and
CAPSULE_SCRATCH_DIR to ordinary writable folders to run outside Code Ocean.
"""
from __future__ import annotations

import argparse
from importlib import metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap", type=int, default=100,
                        help="Synthetic-demo bootstrap count; manuscript settings use 2000.")
    args = parser.parse_args()
    if args.bootstrap < 2:
        parser.error("Use at least two bootstrap replicates.")
    repo = Path(__file__).resolve().parents[1]
    results = Path(os.environ.get("CAPSULE_RESULTS_DIR", "/results" if Path("/results").is_dir()
                                  else str(repo.parent / "results"))).resolve()
    scratch = Path(os.environ.get("CAPSULE_SCRATCH_DIR", "/scratch" if Path("/scratch").is_dir()
                                  else tempfile.gettempdir())).resolve()
    if results == repo or repo in results.parents:
        raise ValueError("CAPSULE_RESULTS_DIR must be outside the read-only code tree")
    if scratch == repo or repo in scratch.parents:
        raise ValueError("CAPSULE_SCRATCH_DIR must be outside the read-only code tree")
    results.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    if (results / "synthetic_model").exists() or (results / "run_summary.json").exists():
        raise ValueError("Choose an empty results location for a fresh run; existing checkpoints are not overwritten.")
    python = os.environ.get("PYTHON", os.environ.get("PYTHON_BIN", sys.executable))
    rscript = os.environ.get("RSCRIPT", os.environ.get("RSCRIPT_BIN", "Rscript"))
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[key] = "1"
    summaries = {"status": "running", "synthetic_only": True, "original_study_data_included": False,
                 "reproduces_manuscript_numerical_results": False, "bootstrap_per_trial": args.bootstrap,
                 "runtime_network_access_required": False, "steps": []}

    def save_summary() -> None:
        (results / "run_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")

    with (results / "run.log").open("w") as log:
        def execute(label: str, command: list[str]) -> None:
            message = f"\n=== {label} ===\n"
            print(message, end="", flush=True)
            log.write(message)
            log.flush()
            process = subprocess.Popen(command, cwd=repo, env=env, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True)
            assert process.stdout is not None
            for line in process.stdout:
                print(line, end="", flush=True)
                log.write(line)
                log.flush()
            code = process.wait()
            summaries["steps"].append({"step": label, "exit_code": code})
            save_summary()
            if code:
                raise subprocess.CalledProcessError(code, command)

        try:
            with tempfile.TemporaryDirectory(prefix="molecular-homa-ir-", dir=scratch) as workdir:
                work = Path(workdir)
                generated = work / "synthetic_inputs"
                execute("Generate entirely synthetic paired inputs", [python, str(repo / "examples/generate_synthetic_data.py"),
                        "--output-dir", str(generated)])
                config = json.loads((repo / "config/smoke.json").read_text())
                for key in ("phenotype", "metabolites", "proteins", "features"):
                    config["inputs"][key] = str(generated / f"{key}.csv")
                config["output_dir"] = str(results / "synthetic_model")
                config_file = work / "capsule_config.json"
                config_file.write_text(json.dumps(config, indent=2) + "\n")
                # Preserve the exact demonstration settings and generation description.
                (results / "demonstration_config.json").write_text(json.dumps(config, indent=2) + "\n")
                shutil.copyfile(generated / "synthetic_generation.json", results / "synthetic_generation.json")
                execute("Train nested cross-validated models and predict week 24", [rscript,
                        "scripts/train_models.R", "--config", str(config_file)])
                execute("Extract frozen molecular contributions", [rscript,
                        "scripts/extract_contributions.R", "--model-run", str(results / "synthetic_model"),
                        "--output", str(results / "decomposition")])
                execute("Analyze trial-specific contributions", [python,
                        str(repo / "scripts/treatment_contributions.py"), "--input", str(results / "decomposition"),
                        "--output", str(results / "contributions"), "--bootstrap", str(args.bootstrap),
                        "--min-folds", "1", "--top-n", "10"])
                execute("Verify held-out and frozen-model invariants", [rscript,
                        "tests/test_model_invariants.R", str(results / "synthetic_model")])
                execute("Verify contribution mathematics and selection", [python, "-m", "unittest", "discover",
                        "-s", str(repo / "tests"), "-p", "test_*.py", "-v"])
                execute("Record R runtime", [rscript, "-e", "for(p in c('glmnet','matrixStats','dplyr','jsonlite')) cat(p,as.character(packageVersion(p)),'\\n'); sessionInfo()"])
            summaries["status"] = "passed"
            summaries["python_version"] = sys.version.split()[0]
            summaries["python_packages"] = {name: metadata.version(name) for name in ("numpy", "pandas", "scipy")}
            if Path("/opt/molecular-homa-ir-explicit.txt").exists():
                shutil.copyfile("/opt/molecular-homa-ir-explicit.txt", results / "environment-conda-explicit.txt")
            save_summary()
            (results / "README.txt").write_text(
                "Synthetic workflow demonstration only. No original participant data are included.\n"
                "The generated observations, scores and analyte names do not reproduce clinical findings.\n"
                "run_summary.json reports completion; run.log records the full workflow and tests.\n"
                "synthetic_model/: fitted models and predictions; decomposition/: exact additive inputs;\n"
                "contributions/: adjusted contrasts, percentages, selection and conservation audits.\n")
        except Exception as exc:
            summaries["status"] = "failed"
            summaries["error"] = str(exc)
            save_summary()
            raise


if __name__ == "__main__":
    main()
