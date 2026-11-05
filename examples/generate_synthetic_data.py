#!/usr/bin/env python3
"""Generate entirely artificial paired inputs for the end-to-end smoke run.

These observations are simulated from a fixed random seed. They contain no
ChiHOPE data and are not intended to reproduce the manuscript's estimates.
Molecular values are already on an analysis scale; do not log them again.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def generate(output_dir: Path, n_per_arm: int = 20, seed: int = 20261105) -> dict:
    if n_per_arm < 12:
        raise ValueError("Use at least 12 participants per trial/treatment arm.")
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    treatment_groups = [("CMAP", a) for a in ("CG32", "CG48", "PB")]
    treatment_groups += [("CMAS", a) for a in ("CG32", "CG48", "SP")]
    n = len(treatment_groups) * n_per_arm
    trial = np.repeat([g[0] for g in treatment_groups], n_per_arm)
    arm = np.repeat([g[1] for g in treatment_groups], n_per_arm)
    sample_id = [f"SYNTHETIC_{i + 1:04d}" for i in range(n)]
    sex = np.tile(np.resize(np.array([1, 2]), n_per_arm), len(treatment_groups))
    age = rng.integers(35, 70, size=n)
    bmi0 = np.clip(rng.normal(27, 3.2, n), 19, 38)
    latent = rng.normal(size=n)
    met0 = latent[:, None] * np.array([0.9, 0.65, -0.75, 0.4, -0.3, 0.2, 0.1, -0.1])
    met0 += rng.normal(0, 0.65, (n, 8))
    prot0 = latent[:, None] * np.array([0.8, -0.7, 0.5, -0.25, 0.15, -0.1])
    prot0 += rng.normal(0, 0.75, (n, 6))
    log_homa0 = 1.3 + 0.50 * latent + 0.11 * met0[:, 0] - 0.10 * prot0[:, 1]
    log_homa0 += rng.normal(0, 0.12, n)
    drug = np.where(arm == "CG32", 0.55, np.where(arm == "CG48", 0.75, 0.0))
    comparator = np.where(arm == "SP", 0.06, 0.0)
    met_delta = -(drug + comparator)[:, None] * np.array([0.8, 0.6, -0.65, 0.3, -0.2, 0.0, -0.12, 0.0])
    prot_delta = -(drug + comparator)[:, None] * np.array([0.75, -0.6, 0.45, -0.2, 0.0, 0.1])
    met24 = met0 + met_delta + rng.normal(0, 0.20, met0.shape)
    prot24 = prot0 + prot_delta + rng.normal(0, 0.22, prot0.shape)
    log_homa24 = log_homa0 - 0.52 * drug - comparator - 0.04 + rng.normal(0, 0.10, n)
    bmi24 = bmi0 + 0.85 * drug + rng.normal(0, 0.30, n)
    met_features = [f"SYN_MET_{i:02d}" for i in range(1, 9)]
    prot_features = [f"SYN_PROT_{i:02d}" for i in range(1, 7)]
    pheno, metabolites, proteins = [], [], []
    for week, bmi, homa, met, prot in [
        (0, bmi0, np.exp(log_homa0), met0, prot0),
        (24, bmi24, np.exp(log_homa24), met24, prot24),
    ]:
        pheno.append(pd.DataFrame({"sample_id": sample_id, "week": week,
            "trial": trial, "arm": arm, "sex": sex, "age": age,
            "homa_ir": homa, "bmi": bmi}))
        for dest, values, names in [(metabolites, met, met_features), (proteins, prot, prot_features)]:
            frame = pd.DataFrame(values, columns=names)
            frame.insert(0, "week", week)
            frame.insert(0, "sample_id", sample_id)
            dest.append(frame)
    pd.concat(pheno, ignore_index=True).to_csv(output_dir / "phenotype.csv", index=False)
    pd.concat(metabolites, ignore_index=True).to_csv(output_dir / "metabolites.csv", index=False)
    pd.concat(proteins, ignore_index=True).to_csv(output_dir / "proteins.csv", index=False)
    features = pd.DataFrame({"feature_id": met_features + prot_features,
        "omics": ["Metabolomics"] * 8 + ["Proteomics"] * 6,
        "label": [f"Synthetic metabolite {i}" for i in range(1, 9)] + [f"Synthetic protein {i}" for i in range(1, 7)]})
    features.to_csv(output_dir / "features.csv", index=False)
    description = {"synthetic_only": True, "seed": seed, "baseline_participants": n,
        "paired_week24_participants": n, "n_per_trial_arm": n_per_arm,
        "metabolites": 8, "proteins": 6, "sex_coding": {"1": "male", "2": "female"},
        "molecular_scale": "already-transformed continuous assay values",
        "note": "Simulation for workflow testing; not a reproduction of clinical findings."}
    (output_dir / "synthetic_generation.json").write_text(json.dumps(description, indent=2) + "\n")
    return description


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("data/synthetic"))
    parser.add_argument("--n-per-arm", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()
    print(json.dumps(generate(args.output_dir, args.n_per_arm, args.seed), indent=2))


if __name__ == "__main__":
    main()
