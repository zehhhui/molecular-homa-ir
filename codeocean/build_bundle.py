#!/usr/bin/env python3
"""Build a Code Ocean upload directory from an explicit source allowlist.

This packages code only, never data/, results/, model objects or local caches.
The output directory must not exist and must be outside the source repository.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    if output.exists() or output == repo or repo in output.parents:
        parser.error("Use a new output directory outside the repository")
    archive = output.with_suffix(".zip")
    if archive.exists():
        parser.error("Archive already exists; choose a new output location")
    output.mkdir(parents=True)
    code = output / "code"
    code.mkdir()
    roots = ("R", "scripts", "config", "examples", "tests", "docs", "reference", "codeocean")
    extensions = {".R", ".py", ".json", ".md", ".sh", ".yml", ".txt"}
    root_files = ("README.md", "README_zh.md", "VALIDATION.md", "LICENSE", "CITATION.cff", "VERSION", "CHANGELOG.md",
                  ".zenodo.json", ".gitignore", "environment.yml", "requirements.txt", "RELEASE_MANIFEST.json")
    paths = [repo / name for name in root_files if (repo / name).is_file()]
    for root in roots:
        for path in (repo / root).rglob("*"):
            if path.is_file() and "__pycache__" not in path.parts and (
                    path.suffix in extensions or path.name in {"run", "postInstall"}):
                paths.append(path)
    for path in sorted(set(paths)):
        if path.is_symlink():
            raise ValueError(f"Symlink not allowed in release source: {path.name}")
        dest = code / path.relative_to(repo)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dest)
    # The capsule's actual master entry point lives at /code/run.
    (code / "run").write_text('#!/usr/bin/env bash\nset -euo pipefail\nSCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"\nexec bash "$SCRIPT_DIR/codeocean/run" "$@"\n')
    for path in (code / "run", code / "codeocean/run", code / "codeocean/postInstall"):
        path.chmod(0o755)
    environment = output / "environment"
    environment.mkdir()
    shutil.copyfile(repo / "codeocean/postInstall", environment / "postInstall")
    (environment / "postInstall").chmod(0o755)
    shutil.copyfile(repo / "environment.yml", environment / "environment.yml")
    (environment / "README.md").write_text(
        "Select a Conda-based Python/Linux starter environment in Code Ocean.\n"
        "Paste postInstall into the Environment editor and build it before running.\n"
        "The recipe has not been built on Code Ocean by the local bundle builder.\n"
        "No Dockerfile or published capsule environment digest is fabricated.\n")
    data = output / "data"
    data.mkdir()
    (data / "README.txt").write_text(
        "No original individual data are distributed. Synthetic inputs are generated\n"
        "deterministically in scratch storage during the reproducible run.\n")
    (output / "REPRODUCING.md").write_text(
        "# Molecular HOMA-IR synthetic demonstration\n\n"
        "Upload code/ to the Code Ocean code folder and select code/run as the entry point.\n"
        "Configure the environment using environment/postInstall, then run once and inspect\n"
        "results/run_summary.json. This demonstrates the workflow, not the clinical estimates.\n"
        "See code/docs/code_ocean.md for configuration and local execution.\n")
    records = {str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(output.rglob("*")) if p.is_file()}
    (output / "bundle_manifest.json").write_text(json.dumps({
        "kind": "Code Ocean upload bundle", "synthetic_only": True,
        "individual_study_data_included": False, "code_ocean_remote_run_status": "not_run",
        "source_files_sha256": records}, indent=2) + "\n")
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
        for path in sorted(output.rglob("*")):
            if path.is_file():
                z.write(path, path.relative_to(output))
    print(json.dumps({"bundle_directory": str(output), "archive": str(archive),
                      "files": len(records) + 1, "individual_study_data_included": False}, indent=2))


if __name__ == "__main__":
    main()
