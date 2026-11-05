# Code Ocean capsule

This repository includes a headless, network-free **synthetic workflow demonstration** for Code Ocean. It runs the same scientific model and contribution code as the repository. It does not reproduce the paper's numerical results because the clinical participant data and original fold assignments are not publicly distributed.

## Build and upload a capsule bundle

From the repository root, choose a new output folder outside this repository:

```bash
python codeocean/build_bundle.py --output ../molecular-homa-ir-codeocean
```

The builder creates an upload directory and ZIP containing:

```text
code/           Source code, documentation, tests and run entry point
environment/    Build-only postInstall and environment specification
data/           A note explaining that inputs are simulated at runtime
REPRODUCING.md  Execution instructions
```

The source allowlist excludes raw data, generated results, fitted models, caches and Git internals. The scientific files are copied unchanged. The derived bundle also records their SHA-256 checksums.

Create a new capsule or copy the public GitHub repository. Imported files may require rearrangement into the Code folder; use the bundle's `code/` contents and set `code/run` as the file to run. Code Ocean's documented Git import does not automatically guarantee the intended folder placement for files that do not match its default folders. See [importing Git repositories](https://docs.codeocean.com/osl-guide/version-control/moving-repositories-in-and-out-of-code-ocean/importing-git-repositories).

## Environment setup

Select a Linux Python starter with Conda available. Paste `codeocean/postInstall` into the Environment editor's **postInstall** field. It creates `/opt/molecular-homa-ir` with R 4.5, glmnet 5.0 and the Python/R dependencies from `environment.yml`. The entry point uses this environment when present. Record the actual solved package versions from the build and run logs; compatible version ranges are not a complete immutable lockfile.

All installation happens during the environment build. The script is self-contained because Code Ocean's post-install stage must not depend on reading the Code or Data folders. A fresh reproducible run installs nothing and downloads no data. See [environment configuration](https://docs.codeocean.com/osl-guide/getting-started/the-computational-environment/configuring-your-computational-environment-an-overview) and [the post-install stage](https://docs.codeocean.com/user-guide/v3.8.0/setting-up-the-environment/the-post-install-script).

The local package preparation validates the run script using an existing tested R/Python environment. It does **not** claim that the Linux Conda solve, Code Ocean environment build or remote reproducible run has already succeeded. Complete those checks in the target account before publishing; do not substitute package versions silently if the build fails.

## Run and inspect

The master script simulates 120 paired participants with 14 synthetic analytes, uses the reduced smoke configuration, extracts frozen contributions, calculates trial-specific contribution percentages with 100 bootstrap replicates, and runs the R and Python scientific invariant tests.

The source code is read-only. Synthetic inputs and temporary configuration go to scratch storage. Models, predictions, contribution tables, logs and checks are written to `/results`, where Code Ocean preserves output files. See [saving results](https://docs.codeocean.com/osl-guide/getting-started/uploading-code-and-data/saving-files) and [writing a run script](https://docs.codeocean.com/osl-guide/user-manual/reproducible-code-execution/how-to-write-and-set-a-run-script).

Expected output files include:

- `run_summary.json`: overall status and exit codes.
- `run.log`: complete step-by-step console output and R session information.
- `synthetic_generation.json` and `demonstration_config.json`: simulated population and exact reduced settings.
- `synthetic_model/`: synthetic fitted models, hold-out predictions and performance.
- `decomposition/`: synthetic participant-level contributions and reconstruction audits.
- `contributions/`: adjusted effects, bootstrap uncertainty, ranked features and numerical conservation checks.

Results with synthetic participant IDs are part of the demonstration and contain no real patient observations. These should not be described as the manuscript's source data. The manuscript 705/593-person, 85/80/52-feature values remain aggregate reference checks for authorized original inputs.

## Local capsule-entry test

Use empty output folders and your installed interpreters:

```bash
CAPSULE_RESULTS_DIR=/path/to/empty/results \
CAPSULE_SCRATCH_DIR=/path/to/scratch \
PYTHON=/path/to/python RSCRIPT=/path/to/Rscript \
bash codeocean/run
```

This checks capsule-style paths and persisted outputs without requiring root-level `/code`, `/data`, `/scratch` or `/results` folders on a local machine. Existing model output is not overwritten. A remote Code Ocean reproducible run and publication are separate steps and must be confirmed by the platform itself.
