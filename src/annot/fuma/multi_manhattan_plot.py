# coding: utf-8
##########################################################################
# NSAp - Copyright (C) CEA, 2023
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
# Antoine Dufournet
##########################################################################

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import glob

def find_files(paths, file_name='most_orig.sumstats'):
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
    
    return sorted(all_files)

def plot_manhattan(file_paths, plot_type):
    # List of distinct base colors for each model
    base_colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

    plt.figure(figsize=(21, 10.5))

    # Initialize chrom_offsets for alignment
    chrom_offsets = None

    # Function to modify color brightness
    def adjust_color_brightness(color, factor):
        """Adjusts brightness of a color by blending with white (factor > 1) or black (factor < 1)."""
        color = mcolors.to_rgb(color)  # Convert color to RGB
        return tuple(min(1, max(0, c * factor)) for c in color)

    # Loop through each file and plot data
    for i, file in enumerate(file_paths):
        print(f"Working with file: {file}")
        model = file.split('/')[-2]
        
        df = pd.read_csv(file, sep='\t')
	
        if file.endswith('glm.linear'):
            df = df.rename(columns={"#CHROM":"CHR",
                                    "POS":"BP", 
                                    "ID":"rsID", 
                                    "REF":"A2", 
                                    "ALT":"ALT", 
                                    "A1":"A1", 
                                    "TEST":"TEST", 
                                    "OBS_CT":"OR", 
                                    "P":"PVAL"})
            df = df.drop(["rsID","A2","ALT","A1","TEST","OR"], axis=1)

        elif file.endswith('.txt'):
            df = df.rename(columns={"POS":"BP",
                                    "SNP":"rsID",
                                    "P":"PVAL"})

        elif file.endswith('.mqfam.total'):
            df = df.rename(columns={"SNP":"rsID", "P":"PVAL"})

        elif file.endswith('.sumstats'):
            df = df.drop(["SNP","A1","A2","Z_FAKE","N"], axis=1)
        
        if i == 0:
            # Calculate cumulative position for x-axis alignment only once
            df.sort_values(by=['CHR', 'BP'], inplace=True)
            chrom_offsets = df.groupby('CHR')['BP'].max().cumsum().shift(fill_value=0)

        # Filter significant SNPs
        df = df[df["PVAL"] < 1] #1e-3

        if plot_type == 'manhattan':
            df['neg_log10_pval'] = -np.log10(df['PVAL'])
        elif plot_type == 'miami':
            if "non.white.British.ancestry" in file:#i%2==1:
                df['neg_log10_pval'] = np.log10(df['PVAL'])
            else:
                df['neg_log10_pval'] = -np.log10(df['PVAL'])

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
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot of Multiple Models')

    # Add chromosome labels at their midpoints
    chromosome_ticks = [chrom_offsets[chrom] + df[df['CHR'] == chrom]['BP'].max() / 2 for chrom in sorted(df['CHR'].unique())]
    chromosome_labels = [f"Chr {chrom}" for chrom in sorted(df['CHR'].unique())]
    plt.xticks(chromosome_ticks, chromosome_labels, rotation=45)

    # Add legend
    plt.legend(loc='upper right')
    plt.tight_layout()

    # Save the plot to the current directory
    output_path = f'{plot_type}_plot.eps'
    plt.savefig(output_path, format='eps')
    print(f"Plot saved to: {output_path}")
    #plt.show()

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate multi-model Manhattan plots for given summary statistic files.")
    parser.add_argument('-p', '--paths', nargs='+', type=str, 
                        help="List of base folder paths or file paths to the summary statistic files.")
    parser.add_argument('-f', '--filename', type=str, default='most_orig.sumstats', 
                        help="End of the name of the summary statistic files.")
    parser.add_argument('-t', '--plottype', type=str, default='manhattan', 
                        help="manhattan or miami")
    args = parser.parse_args()

    # Find all relevant files
    file_name = args.filename
    file_paths = find_files(args.paths, file_name)
    plot_type = args.plottype

    if not file_paths:
        print(f"No files found matching the pattern '{file_name}'. Please check your paths.")
        return

    # Call the plotting function with the found file paths
    plot_manhattan(file_paths,plot_type)

if __name__ == "__main__":
    main()
