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

pcocc-rs run --env UKB_HOME n4h00001rs:quant-genetics-0.3 \
     python3 $ROOTDIR/LDSC/extract_gencorr.py --\
            -p $ROOTDIR/$RESULTSDIR/SCZ_corr/ \
            -f scz_dim1.log \
            -o $ROOTDIR/25irene_AD_allRegionMOSTEST/results/Champollion_V1_32/SCZ_corr/.
"""
##########################################################################

import os
import glob
import argparse
import pandas as pd

def find_files(paths, file_name='scz_dim1.log'):
    """Find all files matching the pattern 'file_name' under the given base paths."""
    all_files = []
    for path in paths:
        if path.endswith(file_name):
            # If the path directly points to a file
            matched_files = glob.glob(path)
        else:
            # If the path is a directory or contains wildcards, search for matching files
            search_pattern = f"{path.rstrip('/')}/**/*{file_name}"
            matched_files = glob.glob(search_pattern, recursive=True)
        
        all_files.extend(matched_files)
    return all_files

def _get_gencorr(path, file_name='scz_dim1.log'):
    print("Working with file:", path)
    gencov_target = "Total Observed scale gencov:"
    gencov = None
    gencorr_target = "Genetic Correlation:"
    gencorr = None
    se = None
    zscore_target = "Z-score:"
    zscore = None
    P_target = "P:"
    P = None
    region_name = path.split("/")[9]
    model_name = path.split("/")[10]
    dim = file_name.replace("scz_", "")
    dim = dim.replace(".log", "")
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
    return region_name,dim,gencov,gencorr,se,zscore,P

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract gencorr estimation from ldsc logs.")
    parser.add_argument('-p', '--paths', nargs='+', type=str, 
                        help="List of base folder paths or file paths to ldsc logs files.")
    parser.add_argument('-f', '--filename', type=str, default='scz_dim1.log', 
                        help="Name of ldsc log files.")
    parser.add_argument('-o', '--out', type=str, default="./", 
                        help="Directory of the output")
    args = parser.parse_args()

    # Find all relevant files
    file_name = args.filename
    file_paths = find_files(args.paths, file_name)
    out = args.out
    out = os.path.dirname(out)

    if not file_paths:
        print(f"No files found matching the pattern '{file_name}'. Please check your paths.")
        return

    # Call the plotting function with the found file paths
    gencorr_dic = {"pheno":[], "dim":[], "gencov":[], "gencorr":[],"se":[], "zscore":[], "P":[]}
    for file_path in file_paths:
        pheno,dim,gencov,gencorr,se,zscore,P = _get_gencorr(file_path, file_name)
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
