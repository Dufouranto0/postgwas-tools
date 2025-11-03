import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy.stats as stats
import argparse
import h5py
import sys
import os

def correlation_heatmap(z_score, path_to_save, name):
    '''
    This function calculates and visualizes the correlation matrix of SNPs from a given region,
    creating a heatmap without using seaborn.
    
    Parameters:
    z_score (DataFrame): The input DataFrame containing SNP data. 
    It should include 'CHR', 'SNP', 'PVAL', 'N', 'FREQ' and additional dimensions (dim1, dim2, ...).

    Returns:
    correlation_matrix (DataFrame): The correlation matrix of SNPs based on their dimensions.
    
    Visualization:
    The function will plot a heatmap of the correlation matrix.
    '''


    df_numerical = z_score.drop(['CHR', 'SNP', 'PVAL', 'N', 'FREQ'], axis=1) #['CHR', 'SNP', 'PVAL']
    df_transposed = df_numerical.T

    correlation_matrix = df_transposed.corr()
    
    labels = z_score['CHR'].astype(str) + " - " + z_score['SNP']
    correlation_matrix.columns = labels
    correlation_matrix.index = labels
    

    fig, ax = plt.subplots(figsize=(12, 12))
    cax = ax.matshow(correlation_matrix, cmap='coolwarm')
    cbar = fig.colorbar(cax, ax=ax, shrink=0.7, fraction=0.046, pad=0.04)
    #ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))

    # Display every nth label to avoid clutter (adjust `step` based on the number of labels)
    step = max(1, len(labels) // 40)
    #ax.set_xticks(np.arange(0, len(labels), step))
    ax.set_yticks(np.arange(0, len(labels), step))

    # Update tick labels and alignment
    #ax.set_xticklabels(labels[::step], rotation=45, ha='center', fontsize=10)  # Center alignment
    ax.set_yticklabels(labels[::step], fontsize=10)

    # Adjust the spacing of the labels
    #ax.tick_params(axis='x', which='major', pad=10)  # Increase space between labels and axis

    # Title and axis labels
    ax.set_title("Correlation Matrix Between SNPs", pad=40, fontsize=14)
    ax.set_xlabel("SNP (CHR - SNP ID)", fontsize=12)
    ax.set_ylabel("SNP (CHR - SNP ID)", fontsize=12)

    # Adjust layout for better fit
    plt.tight_layout()
    clean_path = os.path.dirname(path_to_save)
    plt.savefig(f'{clean_path}/Correlation_Matrix_SNPs_{name}.eps', format='eps')
    #plt.show()
    
    return correlation_matrix

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract mostest results.\n"
            "Usage: process_results_ext.py --bim <bim> --fname <fname> --out <out> --uniGWAS\n"
            "where:"
        )
    )
    parser.add_argument("--bim", type=str,
                        help="Path to .bim file (reference set of SNPs)", required=True)
    parser.add_argument("--fname", type=str,
                        help="Prefix of .mat files output by mostest.m", required=True)
    parser.add_argument("--out", type=str,
                        help="Optional suffix for output files (defaults to fname)",
                        default=None)
    parser.add_argument("--uniGWAS", action="store_true",
                        help="Optional flag to save the univariate GWAS as well")

    args = parser.parse_args()
    bim_file = args.bim
    fname = args.fname
    out = args.out if args.out else fname
    uniGWAS = args.uniGWAS

    """if len(sys.argv) <= 2:
        print('Usage: process_results_ext.py <bim> <fname> [<out>], where')
        print(' bim   - path to bim file (reference set of SNPs')
        print(' fname - prefix of .mat files output by mostest.m, ie. fname should be the same as "out" argument of the mostest.m')
        print(' out   - optional suffix for output files, by defautl fname will be used')
        print(' uniGWAS   - save univariate GWAS')
        sys.exit()

    bim_file = sys.argv[1] #'UKB26502_QCed_230519_maf0p005_chr21.bim'
    fname = sys.argv[2]    # 'all_chr21'
    out = sys.argv[3] if (len(sys.argv) > 3) else sys.argv[2]
    uniGWAS = sys.argv[4] if (len(sys.argv) > 4) else False"""

    # read .bim file (reference set of SNPs)
    print('Load {}...'.format(bim_file))
    bim = pd.read_csv(bim_file, sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
    del bim['GP']

    mat = sio.loadmat(fname + '.mat')
    bim['N'] = mat['nvec']

    # save matrix of z scores for SNPs passing 5e-08 threshold
    SNP_threshold = 5e-08
    print('Generate {}_***.zmat.tsv files...'.format(out))
    with h5py.File(fname + '_zmat.mat', 'r') as h5file:
        # Properly dereference the object references
        measures = []
        for i in range(h5file['measures'].shape[0]):
            ref = h5file['measures'][i, 0]
            measure = ''.join(chr(c[0]) for c in h5file[ref][:])
            measures.append(measure)
        #print(measures)
        freqvec = np.array(h5file['freqvec'])[0]

        zmat_orig = np.array(h5file['zmat_orig']) #if you want to work with minp ad it to the list below
        for test in ['most_log10pval_orig', 'minp_log10pval_orig']: 
            # save matrix of z scores for SNPs passing 5e-08 threshold
            # for both mostest and minp methods
            pval = np.power(10, -mat[test].flatten())
            df_zmat_orig = pd.DataFrame(np.transpose(zmat_orig[:, pval < SNP_threshold]), columns=measures)
            df_zmat_orig.insert(0, 'FREQ', freqvec[pval < SNP_threshold])
            df_zmat_orig.insert(0, 'N', bim.N.values[pval < SNP_threshold])
            df_zmat_orig.insert(0, 'PVAL', pval[pval < SNP_threshold])
            df_zmat_orig.insert(0, 'SNP', bim.SNP.values[pval < SNP_threshold])
            df_zmat_orig.insert(0, 'CHR', bim.CHR.values[pval < SNP_threshold])
            #print(out + "_" + test.replace('_log10pval', '') + '.zmat.tsv')
            #print(df_zmat_orig.head())
            df_zmat_orig.to_csv(out + "_" + test.replace('_log10pval', '') + '.zmat.tsv', index=False, sep='\t')
            correlation_heatmap(df_zmat_orig, out, test.replace('_log10pval_orig', ''))

	# save the z scores from the permuted genotype 
        zmat_perm = np.array(h5file['zmat_perm'])

        for test in ['most_log10pval_perm', 'minp_log10pval_perm']:
            pval = np.power(10, -mat[test].flatten())
            df_zmat_perm = pd.DataFrame(np.transpose(zmat_perm[:, pval < SNP_threshold]), columns=measures)
            df_zmat_perm.insert(0, 'FREQ', freqvec[pval < SNP_threshold])
            df_zmat_perm.insert(0, 'N', bim.N.values[pval < SNP_threshold])
            df_zmat_perm.insert(0, 'PVAL', pval[pval < SNP_threshold])
            df_zmat_perm.insert(0, 'SNP', bim.SNP.values[pval < SNP_threshold])
            df_zmat_perm.insert(0, 'CHR', bim.CHR.values[pval < SNP_threshold])
            df_zmat_perm.to_csv(out + "_" + test.replace('_log10pval', '') + '.zmat.tsv', index=False, sep='\t')

        del zmat_perm

        # save individual GWAS results ('freqvec' is an indicator that we've saved individual GWAS beta's)
        if uniGWAS:
            if 'freqvec' in h5file:
                bim['FRQ'] = np.transpose(np.array(h5file['freqvec']))

                # Create the folder for the individual GWAS results 
                folder_uni = "GWAS_uni"
                os.makedirs(f"{os.path.dirname(out)}/{folder_uni}", exist_ok=True)
                print(f"Folder '{folder_uni}' created successfully.")

                beta_orig = np.array(h5file['beta_orig'])
                # division by 0 occured when zmat_orig was to close to 0
                # therefore the default value for the standard error will
                # be set as 999
                se_orig = 999*np.ones(shape=beta_orig.shape, dtype=float)
                se_orig = np.divide(beta_orig, zmat_orig, out=se_orig, where=(np.abs(zmat_orig)>0.01))
                pval_orig = stats.norm.sf(np.abs(zmat_orig)) * 2.0

                for measure_index, measure in enumerate(measures):
                    fname = '{}/{}/{}.orig.sumstats.gz'.format(os.path.dirname(out),folder_uni, measure)
                    print('Generate {}...'.format(fname))
                    bim['PVAL'] = np.transpose(pval_orig[measure_index, :])
                    bim['Z'] = np.transpose(zmat_orig[measure_index, :])
                    bim['BETA'] = np.transpose(beta_orig[measure_index, :])
                    bim['SE'] = np.transpose(se_orig[measure_index, :])
                    bim.to_csv(fname, compression='gzip', sep='\t', index=False)

                del beta_orig
                del se_orig
                del pval_orig
        del zmat_orig


if __name__ == '__main__':
    main()
