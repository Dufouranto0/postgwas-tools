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
RESULTSDIR=25irene_AD_CING_rIG/results/ChampollionV0/CINGULATE_left/09-35-58_201/non.white.British.ancestry

pcocc-rs run --env UKB_HOME n4h00001rs:quant-genetics-0.3 \
     python3 $ROOTDIR/LDSC/extract_h2.py --\
            -p $ROOTDIR/$RESULTSDIR/h2/ \
            -f orig_munge_h2.log \
            -o $ROOTDIR/$RESULTSDIR/h2/.
"""
##########################################################################

import os
import glob
import argparse
import pandas as pd

def find_files(paths, file_name='orig_munge_h2.log'):
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

def _get_h2(path, file_name='orig_munge_h2.log'):
    print("Working with file", path)
    h2_target = "Total Observed scale h2:"
    h2 = None
    se = None
    lambdaGC_target = "Lambda GC:"
    lambdaGC = None
    Mean_Chi2_target = "Mean Chi^2:"
    Mean_Chi2 = None
    Intercept_target = "Intercept:"
    Intercept = None
    pheno = (os.path.basename(path)).replace(f'.{file_name}', "")
    with open(path, "rt") as of:
        lines = of.readlines()
    for line in lines:
        if h2_target in line:
            try:
                h2 = float(line.replace(h2_target, "").strip().split(" ", 1)[0])
                floatse = line.replace(h2_target, "").strip().split(" ", 1)[1]
                se = float(floatse[1:-1])
            except ValueError:
                print(f"Warning: Could not extract h2 from {path}")
        if lambdaGC_target in line:
            try:
                lambdaGC = float(line.replace(lambdaGC_target, "").strip().split(" ", 1)[0])
            except ValueError:
                print(f"Warning: Could not extract lambdaGC from {path}")
        if Mean_Chi2_target in line:
            try:
                Mean_Chi2 = float(line.replace(Mean_Chi2_target, "").strip().split(" ", 1)[0])
            except ValueError:
                print(f"Warning: Could not extract Mean_Chi2 from {path}")
        if Intercept_target in line:
            try:
                Intercept = float(line.replace(Intercept_target, "").strip().split(" ", 1)[0])
            except ValueError:
                print(f"Warning: Could not extract Intercept from {path}")
            break
    return pheno,h2,se,lambdaGC,Mean_Chi2,Intercept

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract h2 estimation from ldsc logs for univariate results.")
    parser.add_argument('-p', '--paths', nargs='+', type=str, 
                        help="List of base folder paths or file paths to ldsc logs files.")
    parser.add_argument('-f', '--filename', type=str, default='orig_munge_h2.log', 
                        help="End of the name of the ldsc log files, which does not contain the phenotype name.")
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
    h2_dic = {"pheno":[], "h2":[], "se":[], "lambdaGC":[], "Mean_Chi2":[], "Intercept":[]}
    for file_path in file_paths:
        pheno,h2,se,lambdaGC,Mean_Chi2,Intercept = _get_h2(file_path, file_name)
        h2_dic["pheno"].append(pheno)
        h2_dic["h2"].append(h2)
        h2_dic["se"].append(se)
        h2_dic["lambdaGC"].append(lambdaGC)
        h2_dic["Mean_Chi2"].append(Mean_Chi2)
        h2_dic["Intercept"].append(Intercept)
     
    h2_df = pd.DataFrame(h2_dic)
    h2_df.to_csv(f"{out}/h2_summary.tsv", sep='\t', index=False)
    print("Summary of h2 saved at:","\n", f"{out}/h2_summary.tsv")

if __name__ == "__main__":
    main()
