"""
Inspired from:

https://github.com/ShujiaHuang/qmplot


Example of use:

python3 QQ_plot.py \
            -f maf-0.01.most_orig.sumstats \
            -p PVAL
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from qmplot import qqplot
import os


if __name__ == "__main__":

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate Manhattan plot for a given summary statistic file.")
    parser.add_argument('-f', '--filepath', type=str, 
                        help="File path to the summary statistic file.")
    parser.add_argument('-p', '--pvalue', type=str, 
                        help="Name of the column with the P-values.")

    args = parser.parse_args()
    file_path = args.filepath
    Pvalue_column = args.pvalue

    print(f"Working with file: {file_path}")

    model = file_path.split('/')[-3]  
    region = file_path.split('/')[-4]

    df = pd.read_table(file_path, sep="\t")
    df = df.dropna(how="any", axis=0)  # clean data

    # Create a Q-Q plot
    f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qqplot(data=df[Pvalue_column],
           marker="o",
           title=f"QQ plot for {region} {model}",
           xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$",
           ax=ax)

    path_to_save = f"{os.path.dirname(file_path)}/QQplot.png"
    plt.savefig(path_to_save)
    print(f"File saved at: {path_to_save}")