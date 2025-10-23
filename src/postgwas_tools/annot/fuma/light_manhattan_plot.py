# coding: utf-8
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
##########################################################################

from postgwas_tools.annot.utils import (
    adjust_color_brightness,
    read_sumstats,
)

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_manhattan(file_path, endoftitle=None):
    base_colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

    plt.figure(figsize=(19.20,10.80))

    # Initialize chrom_offsets for alignment
    chrom_offsets = None
    
    print(f"Working with file: {file_path}")
    df = read_sumstats(file_path)
        
    dic_start_end_chr = {}
    chromosomes = sorted(df['CHR'].unique())


    # Calculate cumulative position for x-axis alignment
    df.sort_values(by=['CHR', 'BP'], inplace=True)
    chrom_offsets = df.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)

    padding_values = [30_000_000]*11 + [3_000_000]*3 + [30_000_000]*5 + [3_000_000]*3
    padding_series = pd.Series(padding_values, index=chrom_offsets.index)
    padding_cumsum = padding_series.cumsum().shift(fill_value=0)
    # Final offsets with padding between chromosomes
    chrom_offsets += padding_cumsum

    for chrom in chromosomes:
        dic_start_end_chr[chrom] = [df[df["CHR"]==chrom]["BP"].min(),df[df["CHR"]==chrom]["BP"].max()]

    # Filter significant SNPs
    df = df[df["P"] < 1e-2] #1e-3
    df['neg_log10_pval'] = -np.log10(df['P'])
    df['x_val'] = df.apply(lambda row: row['BP'] + chrom_offsets.loc[row['CHR']], axis=1)

    # Loop through chromosomes to alternate colors
    for chrom_idx, chrom in enumerate(chromosomes):
        chrom_data = df[df['CHR'] == chrom]

        # Adjust color brightness based on chromosome index
        base_color = base_colors[0]  # choose the color
        brightness_factor = 1.5 if chrom_idx % 2 == 0 else 0.5  # Stronger light or dark shades
        adjusted_color = adjust_color_brightness(base_color, brightness_factor)

        # Plot data for the chromosome
        plt.scatter(chrom_data['x_val'], chrom_data['neg_log10_pval'], 
                    color=adjusted_color, s=7) 

        start = df[df["CHR"]==chrom]["BP"].min() + chrom_offsets[chrom]
        end = df[df["CHR"]==chrom]["BP"].max() + chrom_offsets[chrom] + 1_000_000
        mid = (start + end) / 2
        width = end - start + 6_000_000
        # Height for the bar below -log10(p)=3
        plt.bar(mid, 2, width=width, color=adjusted_color, align='center')

    # Add a horizontal significance threshold line
    plt.axhline(y=-np.log10(0.05/1e6), color='r', linestyle='--')
    plt.axhline(y=-np.log10(0.05/1e4), color='g', linestyle='--')

    # Customize plot labels and title
    plt.xlabel('Chromosome', fontsize=22)
    plt.ylabel('-log10(p-value)', fontsize=22)
    plt.title(f'Manhattan Plot {endoftitle}', fontsize=24)

    # Add chromosome labels at their midpoints
    chromosome_ticks = [chrom_offsets[chrom] + df[df['CHR'] == chrom]['BP'].max() / 2 for chrom in sorted(df['CHR'].unique())]
    chromosome_labels = [f"{chrom}" for chrom in sorted(df['CHR'].unique())]
    plt.xticks(chromosome_ticks, chromosome_labels, fontsize=16) #rotation=45,
    plt.yticks(fontsize=16)
    plt.xlim([dic_start_end_chr[1][0] + chrom_offsets[1] - 13_000_000, dic_start_end_chr[22][1] + chrom_offsets[22]+ 13_000_000])
    plt.ylim(bottom=0)
    #plt.legend(loc='upper right')
    plt.tight_layout()
    
    output_path = f'{os.path.dirname(file_path)}/manhattan_plot.png'
    plt.savefig(output_path, format='png', transparent=True)
    print(f"Plot saved to: {output_path}")
    plt.show()

def main():

    parser = argparse.ArgumentParser(description="Generate Manhattan plot for a given summary statistic file.")
    parser.add_argument('-p', '--path', type=str, 
                        help="File path to the summary statistic file.")
    parser.add_argument('-t', '--endoftitle', type=str, default=None,
                        help="End of the title you want to plot.")
    
    args = parser.parse_args()
    file_path = args.path
    endoftitle = args.endoftitle

    if not endoftitle:
        model = file_path.split('/')[-3]  
        region = file_path.split('/')[-4] 
        endoftitle = f"of the {region} region with model {model}"

    # Call the plotting function with the found file paths
    plot_manhattan(file_path, endoftitle)

if __name__ == "__main__":
    main()
