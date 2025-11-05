#!/usr/bin/env python3
"""
Inspired from:

https://github.com/ShujiaHuang/qmplot

"""

import argparse
import matplotlib.pyplot as plt
from qmplot import qqplot
import os
from postgwas_tools.annot.utils import read_sumstats

def main():
    parser = argparse.ArgumentParser(description="Generate Manhattan plot for a given summary statistic file.")
    parser.add_argument('-p', '--path', type=str, 
                        help="File path to the summary statistic file.")
    parser.add_argument('-o', '--out', type=str, default=None,
                        help="Output folder for saving the plot.")
    

    args = parser.parse_args()
    file_path = args.path
    output_folder= args.out

    print(f"Working with file: {file_path}")

    df = read_sumstats(file_path)
    df = df.dropna(how="any", axis=0)  # clean data

    # Create a Q-Q plot
    f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qqplot(data=df["P"],
           marker="o",
           title=f"QQ plot",
           xlabel=r"Expected $-log_{10}{(p)}$",
           ylabel=r"Observed $-log_{10}{(p)}$",
           ax=ax)

    if not output_folder:
       path_to_save = f"{os.path.dirname(file_path)}/QQplot.png"
    else:
       path_to_save = f"{output_folder}/QQplot.png"
    plt.savefig(path_to_save)
    print(f"File saved at: {path_to_save}")

if __name__ == "__main__":
    main()