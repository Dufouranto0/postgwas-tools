# coding: utf-8
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
# Inspired from Vincent Frouin's code
##########################################################################
"""
Example of usage:

ROOTDIR=/ccc/workflash/cont003/n4h00001/n4h00001
RESULTSDIR=25irene_AD_allRegionMOSTEST/results/Champollion_V1_32/*/*/32PCs/white.British.ancestry

     python3 $ROOTDIR/LDSC/extract_gencorr.py \
            -p $ROOTDIR/$RESULTSDIR/SCZ_corr/scz_dim1.log \
            -o $ROOTDIR/25irene_AD_allRegionMOSTEST/results/Champollion_V1_32/SCZ_corr/.
"""
##########################################################################

import os
import argparse
import pandas as pd
from postgwas_tools.annot.utils import find_files

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
                print(f"Warning: Could not extract gencorr from {path}")
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

def main():
    parser = argparse.ArgumentParser(description="Extract gencorr estimation from ldsc logs.")
    parser.add_argument('-p', '--paths', nargs='+', type=str, 
                        help="List of ldsc logs paths.")
    parser.add_argument('-o', '--out', type=str, default="./", 
                        help="Directory of the output")
    args = parser.parse_args()

    # Find all relevant files
    file_paths = find_files(args.paths)
    out = args.out
    out = os.path.dirname(out)

    # Call the plotting function with the found file paths
    gencorr_dic = {"pheno":[], "dim":[], "gencov":[], "gencorr":[],"se":[], "zscore":[], "P":[]}
    for file_path in file_paths:
        print(f"Working with file: {file_path}")
        base_name = os.path.basename(file_path)
        pheno = file_path.split("/")[-6]
        dim = base_name.replace(".log", "")
        gencov,gencorr,se,zscore,P = _get_gencorr(file_path)
        gencorr_dic["pheno"].append(pheno)
        gencorr_dic["dim"].append(dim)
        gencorr_dic["gencov"].append(gencov)
        gencorr_dic["gencorr"].append(gencorr)
        gencorr_dic["se"].append(se)
        gencorr_dic["zscore"].append(zscore)
        gencorr_dic["P"].append(P)
     
    gencorr_dic = pd.DataFrame(gencorr_dic)
    gencorr_dic.to_csv(f"{out}/gencorr_summary.tsv", sep='\t', index=False)
    print("Summary of gencorr saved at:","\n", f"{out}/gencorr_summary.tsv")

if __name__ == "__main__":
    main()
