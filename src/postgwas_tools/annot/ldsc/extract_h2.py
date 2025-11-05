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

##########################################################################

import os
import argparse
import pandas as pd
from postgwas_tools.annot.utils import find_files

def _get_h2(path):
    h2_target = "Total Observed scale h2:"
    h2 = None
    se = None
    lambdaGC_target = "Lambda GC:"
    lambdaGC = None
    Mean_Chi2_target = "Mean Chi^2:"
    Mean_Chi2 = None
    Intercept_target = "Intercept:"
    Intercept = None
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
    return h2,se,lambdaGC,Mean_Chi2,Intercept

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract h2 estimation from ldsc logs for univariate results.")
    parser.add_argument('-p', '--paths', nargs='+', type=str, 
                        help="List of ldsc logs files.")
    parser.add_argument('-o', '--out', type=str, default="./", 
                        help="Directory of the output")
    args = parser.parse_args()

    # Find all relevant files
    file_paths = find_files(args.paths)
    out = args.out
    out = out.replace('/', '') if out.endswith('/') else out

    # Call the plotting function with the found file paths
    h2_dic = {"pheno":[], "h2":[], "se":[], "lambdaGC":[], "Mean_Chi2":[], "Intercept":[]}
    for file_path in file_paths:
        print(f"Working with file: {file_path}")
        base_name = os.path.basename(file_path)
        pheno = base_name.split('.')[0]
        h2,se,lambdaGC,Mean_Chi2,Intercept = _get_h2(file_path)
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
