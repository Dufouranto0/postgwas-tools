#!/usr/bin/env python3
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
# Inspired from Vincent Frouin's code
##########################################################################

import os
import json
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2
from postgwas_tools.annot.utils import find_files


def _get_h2(path):
    h2_target = "Total Observed scale h2:"
    h2_values = []

    with open(path, "rt") as f:
        lines = f.readlines()

    for line in lines:
        if h2_target in line:
            try:
                # Extract the float before the parentheses
                h2_str = line.replace(h2_target, "").strip().split(" ", 1)[0]
                h2_values.append(float(h2_str))
            except ValueError:
                print(f"Warning: Could not extract h2 from line: {line.strip()}")

    # Ensure we got exactly two phenotypes
    if len(h2_values) != 2:
        print(f"Warning: Expected 2 h2 values, got {len(h2_values)}")

    return tuple(h2_values)

def _get_gencorr(path):
    gencov_target = "Total Observed scale gencov:"
    gencov = None
    gencorr_target = "Genetic Correlation:"
    gencorr = None
    se = None
    zscore_target = "Z-score:"
    zscore = None
    P_target = "P:"
    P = None
    with open(path, "rt") as of:
        lines = of.readlines()
    for line in lines:
        if gencov_target in line:
            try:
                gencov = float(line.replace(gencov_target, "").strip().split(" ", 1)[0])
            except ValueError:
                print(f"Warning: Could not extract gencov from {path}")
        if gencorr_target in line:
            try:
                gencorr = float(line.replace(gencorr_target, "").strip().split(" ", 1)[0])
                floatse = line.replace(gencorr_target, "").strip().split(" ", 1)[1]
                se = float(floatse[1:-1])
            except ValueError:
                print(f"Warning: Could not extract gencorr from {path} \n (likely h2 out of bounds)")
        if zscore_target in line:
            try:
                zscore = float(line.replace(zscore_target, "").strip().split(" ", 1)[0])
            except ValueError:
                print(f"Warning: Could not extract zscore from {path}")
        if P_target in line:
            try:
                P = float(line.replace(P_target, "").strip().split(" ", 1)[0])
            except ValueError:
                print(f"Warning: Could not extract P from {path}")
            break
    return gencov,gencorr,se,zscore,P

def create_df(file_paths, prefix, h2threshold):
    gencorr_dic = {"pheno":[], "gencov":[], "gencorr":[],"se":[], "zscore":[], "P":[]}
    for file_path in file_paths:
        print(f"Working with file: {file_path}")
        base_name = os.path.basename(file_path)
        pheno = base_name.replace(".log", "")
        pheno = pheno.replace(prefix, "")
        h2_ph1, h2_ph2 = _get_h2(file_path)
        if h2_ph1 >h2threshold and h2_ph2 > h2threshold:
            gencov,gencorr,se,zscore,P = _get_gencorr(file_path)
            gencorr_dic["pheno"].append(pheno)
            gencorr_dic["gencov"].append(gencov)
            gencorr_dic["gencorr"].append(gencorr)
            gencorr_dic["se"].append(se)
            gencorr_dic["zscore"].append(zscore)
            gencorr_dic["P"].append(P)
     
    gencorr_df = pd.DataFrame(gencorr_dic)
    return gencorr_df

def compute_correlation(multphen_path, pheno_names):
    """Compute correlation matrix among phenotypes from a delimited file (auto-detects delimiter)."""
    if multphen_path and os.path.exists(multphen_path):
        df = pd.read_csv(multphen_path, sep=None, engine='python')  # auto-detects delimiter
        df = df[pheno_names]  # select only relevant columns
        corr_matrix = df.corr().values
    else:
        corr_matrix = np.identity(len(pheno_names))
    return corr_matrix

def omnibus(z_scores, corr_matrix):
    """Compute Mahalanobis chi-square statistic."""
    inv_corr = np.linalg.inv(corr_matrix)
    chi2_val = float(z_scores.T @ inv_corr @ z_scores)
    df = len(z_scores)
    p_val = 1 - chi2.cdf(chi2_val, df)
    return chi2_val, df, p_val

def main():
    parser = argparse.ArgumentParser(description="Extract gencorr estimation from ldsc logs.")
    parser.add_argument('-p', '--paths', nargs='+', type=str, 
                        help="List of ldsc logs paths.")
    parser.add_argument('-o', '--out', type=str, default="./", 
                        help="Directory of the output")
    parser.add_argument("--prefix", default="", help="Common prefix for LDSC files (e.g., 'scz_')")
    parser.add_argument("--omnibus", action="store_true",
                        help="Optional flag to do omnibus test as well")
    parser.add_argument("--multphen", default=None, help="Path to the file with the phenotypes")
    parser.add_argument("--h2threshold", type=float, default=0, help="h2 threshold to decide which phenotype to keep based on the h2 estimation in the .log")
    args = parser.parse_args()

    # Find all relevant files
    file_paths = find_files(args.paths)
    out = args.out
    out = out.replace('/', '') if out.endswith('/') else out
     
    gencorr_df = create_df(file_paths, args.prefix, args.h2threshold)

    gencorr_df.to_csv(f"{out}/gencorr_summary.tsv", sep='\t', index=False)
    print("Summary of gencorr saved at:","\n", f"{out}/gencorr_summary.tsv")

    if args.omnibus:
        gencorr_df = gencorr_df.dropna()
        pheno_names = gencorr_df["pheno"].to_list()
        z_scores = np.array(gencorr_df["zscore"].to_list())
        if not args.multphen:
            print(f"Specify the path to the file with the phenotypes just before the gwas ot the mostest with --multphen")
        corr_matrix = compute_correlation(args.multphen, pheno_names)
        chi2_val, df, p_val = omnibus(z_scores, corr_matrix)  
        results = {
        "chi2_statistic": chi2_val,
        "degrees_of_freedom": df,
        "p_value": p_val,
        "phenotypes": pheno_names
        }  
        os.makedirs(out, exist_ok=True)
        json_path = os.path.join(out, f"meta_{args.prefix}ldsc.json")
        with open(json_path, "w") as f:
            json.dump(results, f, indent=4)

        print(f"Results saved at: {json_path}")


if __name__ == "__main__":
    main()
