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
    find_files,
    adjust_color_brightness,
    read_sumstats,
)

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_manhattan(file_paths, plot_type, output_folder):
    # List of distinct base colors for each model
    base_colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

    plt.figure(figsize=(21, 10.5))

    # Initialize chrom_offsets for alignment
    chrom_offsets = None
    dic_start_end_chr = {}

    # Loop through each file and plot data
    for i, file in enumerate(file_paths):
        print(f"Working with file: {file}")
        model = file.split('/')[-2]
        
        df = read_sumstats(file)
        
        if i == 0:
            # Calculate cumulative position for x-axis alignment only once
            df.sort_values(by=['CHR', 'BP'], inplace=True)
            chromosomes = sorted(df['CHR'].unique())
            chrom_offsets = df.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)
            padding_values = [30_000_000]*11 + [3_000_000]*3 + [30_000_000]*5 + [3_000_000]*3
            padding_series = pd.Series(padding_values, index=chrom_offsets.index)
            padding_cumsum = padding_series.cumsum().shift(fill_value=0)
            # Final offsets with padding between chromosomes
            chrom_offsets += padding_cumsum

            for chrom in chromosomes:
                dic_start_end_chr[chrom] = [df[df["CHR"]==chrom]["BP"].min(),df[df["CHR"]==chrom]["BP"].max()]

        # Filter significant SNPs
        df = df[df["P"] < 1] #1e-3

        if plot_type == 'manhattan':
            df['neg_log10_pval'] = -np.log10(df['P'])

        elif plot_type == 'miami':
            if "non.white.British.ancestry" in file:#i%2==1:
                df['neg_log10_pval'] = np.log10(df['P'])
            else:
                df['neg_log10_pval'] = -np.log10(df['P'])

        df['x_val'] = df.apply(lambda row: row['BP'] + chrom_offsets.loc[row['CHR']], axis=1)

        # Loop through chromosomes to alternate colors
        chromosomes = sorted(df['CHR'].unique())
        for chrom_idx, chrom in enumerate(chromosomes):
            chrom_data = df[df['CHR'] == chrom]

            # Adjust color brightness based on chromosome index
            base_color = base_colors[i % len(base_colors)]  # Base color based on the model
            brightness_factor = 1.5 if chrom_idx % 2 == 0 else 0.5  # Stronger light or dark shades
            adjusted_color = adjust_color_brightness(base_color, brightness_factor)

            # Plot data for the chromosome
            plt.scatter(chrom_data['x_val'], chrom_data['neg_log10_pval'], 
                        color=adjusted_color, s=7,
                        label=model if chrom_idx == 0 else "")  # Add label only once per model

    # Add a horizontal significance threshold line
    plt.axhline(y=-np.log10(0.05/1000000), color='r', linestyle='--')
    plt.axhline(y=-np.log10(0.05/10000), color='g', linestyle='--')

    # Customize plot labels and title
    plt.xlabel('Chromosome', fontsize=22)
    plt.ylabel(r"$-log_{10}{(p)}$", fontsize=22)
    plt.title('Manhattan Plot', fontsize=24)

    # Add chromosome labels at their midpoints
    chromosome_ticks = [chrom_offsets[chrom] + df[df['CHR'] == chrom]['BP'].max() / 2 for chrom in sorted(df['CHR'].unique())]
    chromosome_labels = [f"{chrom}" for chrom in sorted(df['CHR'].unique())]
    plt.xticks(chromosome_ticks, chromosome_labels, fontsize=16) #rotation=45
    plt.yticks(fontsize=16)
    plt.xlim([dic_start_end_chr[1][0] + chrom_offsets[1] - 13_000_000, dic_start_end_chr[22][1] + chrom_offsets[22]+ 13_000_000])
    if plot_type == 'manhattan':
        plt.ylim(bottom=0)
    #plt.legend(loc='upper right',
    #           fontsize=16,
    #           frameon=False     )
    plt.tight_layout()

    output_path = f'{plot_type}_plot.png'
    plt.savefig(f"{output_folder}/{output_path}", format='png', transparent=True)
    print(f"Plot saved to: {output_folder}/{output_path}")
    #plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Generate Manhattan plots for one or multiple GWAS summary statistics files."
    )
    parser.add_argument(
        '-p', '--paths', nargs='+', type=str, required=True,
        help="List of paths, file names, or wildcard patterns (e.g. /data/*.sumstats)"
    )
    parser.add_argument(
        '-k', '--kind', type=str, default='manhattan',
        choices=['manhattan', 'miami'],
        help="Type of plot to generate: 'manhattan' (default) or 'miami'"
    )
    parser.add_argument(
        '-o', '--out', type=str, default='.',
        help="Output folder for saving the plot."
    )
    args = parser.parse_args()

    # Find files matching any path or pattern
    file_paths = find_files(args.paths)
    plot_type = args.kind
    output_folder = args.out

    if not file_paths:
        print("No files found. Please check your paths or patterns.")
        return

    print(f"Found {len(file_paths)} file(s):")
    for f in file_paths:
        print(f"   - {f}")

    plot_manhattan(file_paths, plot_type, output_folder)

if __name__ == "__main__":
    main()
