#!/usr/bin/env python3
"""
mahalanobis.py --folder SCZ_corr \
               --prefix scz_ \
               --multphen full_embeddings_IID.csv \
               --out results
"""

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2

def extract_z_score(log_path):
    """Extract Z-score from an LDSC .log file."""
    with open(log_path, 'r') as f:
        text = f.read()
    match = re.search(r"Z-score:\s*([-+]?[0-9]*\.?[0-9]+)", text)
    if match:
        return float(match.group(1))
    else:
        print(f"Warning: Z-score not found in {log_path}")
        return np.nan

def load_z_scores(folder, prefix):
    """Load Z-scores from all .log files sharing prefix."""
    z_scores = []
    pheno_names = []
    for file in sorted(os.listdir(folder)):
        if file.startswith(prefix) and file.endswith(".log"):
            z = extract_z_score(os.path.join(folder, file))
            pheno_name = re.sub(r'\.log$', '', file.replace(prefix, ''))
            z_scores.append(z)
            pheno_names.append(pheno_name)
    return np.array(z_scores), pheno_names

def compute_correlation(csv_path, pheno_names):
    """Compute correlation matrix among phenotypes from CSV."""
    df = pd.read_csv(csv_path)
    df = df[pheno_names]  # select only relevant columns
    corr_matrix = df.corr().values
    return corr_matrix

def mahalanobis_chi2(z_scores, corr_matrix):
    """Compute Mahalanobis chi-square statistic."""
    inv_corr = np.linalg.inv(corr_matrix)
    chi2_val = float(z_scores.T @ inv_corr @ z_scores)
    df = len(z_scores)
    p_val = 1 - chi2.cdf(chi2_val, df)
    return chi2_val, df, p_val

def main():
    parser = argparse.ArgumentParser(description="Meta LDSC Mahalanobis test across phenotypes")
    parser.add_argument("--folder", required=True, help="Folder containing LDSC .log files")
    parser.add_argument("--prefix", required=True, help="Common prefix for LDSC files (e.g., 'scz_')")
    parser.add_argument("--multphen", required=True, help="CSV file with phenotype values")
    parser.add_argument("--out", required=False, help="Output folder to write the results", default=None)
    args = parser.parse_args()

    out = args.out 
    if not out:
        out = args.folder

    print("Loading Z-scores...")
    z_scores, pheno_names = load_z_scores(args.folder, args.prefix)
    print(f"Found {len(z_scores)} phenotypes:\n {pheno_names}")
    print(f"List of the z-scores:\n", z_scores)

    print("Computing phenotype correlation...")
    corr_matrix = compute_correlation(args.multphen, pheno_names)

    print("Performing Mahalanobis test...")
    chi2_val, df, p_val = mahalanobis_chi2(z_scores, corr_matrix)

    print("\n=== Meta LDSC Mahalanobis Test ===")
    print(f"Chi² statistic = {chi2_val:.4f}")
    print(f"Degrees of freedom = {df}")
    print(f"P-value = {p_val:.4g}")

    results = {
        "chi2_statistic": chi2_val,
        "degrees_of_freedom": df,
        "p_value": p_val,
        "phenotypes": pheno_names
    }

    json_path = os.path.join(out, "meta_ldsc_results.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=4)

    print(f"Results saved at: {json_path}")

if __name__ == "__main__":
    main()

