import csv
import pandas as pd
import numpy as np
import os
from argparse import ArgumentParser

# Get the process name, use it in the output directory
def get_argument_parser():
    p = ArgumentParser()
    p.add_argument("--process", "-p",
        help="The name of the process")
    return p

def adjust_pvalue(df, col):

    # get indices of null values
    na_idx = df[df[col].isnull()].index.tolist()

    # adjust non-null p values
    pvals = df[col].values.tolist()
    pvals = [x for x in pvals if str(x) != 'nan']
    adj_pval = p_adjust_bh(pvals).tolist()

    # insert null at original indices
    for idx in na_idx:
        adj_pval.insert(idx, None)

    # adjusted p values as new df
    df.loc[:, 'adj_pval'] = adj_pval
    return df


def p_adjust_bh(p):
    ## multiple hypothesis testing
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def process_results(result_dirs, type, proc):
    proteins = ['BRUNOL4(Hs/Mm)', 'BRUNOL5(Hs/Mm)', 'BRUNOL6(Hs/Mm)', 'DAZAP1(Hs/Mm)', 'ESRP2(Hs/Mm)', 'FMR1(Hs/Mm)', 'FUS(Hs/Mm)', 'FXR1(Hs/Mm)', 'FXR2(Hs/Mm)', 'HNRNPA1(Hs/Mm)', 'HNRNPA1L2(Hs/Mm)', 'HNRNPA2B1(Hs/Mm)', 'HNRNPC(Hs/Mm)', 'HNRNPF(Hs/Mm)', 'HNRNPH1(Hs/Mm)', 'HNRNPH2(Hs/Mm)', 'HNRNPK(Hs/Mm)', 'HNRNPL(Hs/Mm)', 'HNRNPM(Hs/Mm)', 'HNRNPU(Hs/Mm)', 'HuR(Hs/Mm)', 'KHDRBS1(Hs/Mm)', 'KHDRBS2(Hs/Mm)', 'KHDRBS3(Hs/Mm)', 'MBNL1(Hs/Mm)', 'PABPC1(Hs/Mm)', 'PABPN1(Hs/Mm)', 'PCBP1(Hs/Mm)', 'PCBP2(Hs/Mm)', 'PTB3(Hs/Mm)', 'QKI(Hs/Mm)', 'RALY(Hs/Mm)', 'RBFOX1(Hs/Mm)', 'RBM24(Hs/Mm)', 'RBM28(Hs/Mm)', 'RBM3(Hs/Mm)', 'RBM4(Hs/Mm)', 'RBM42(Hs/Mm)', 'RBM5(Hs/Mm)', 'RBM8A(Hs/Mm)', 'SART3(Hs/Mm)', 'SFPQ(Hs/Mm)', 'SNRNP70(Hs/Mm)', 'SNRPA(Hs/Mm)', 'SRSF1(Hs/Mm)', 'SRSF10(Hs/Mm)', 'SRSF2(Hs/Mm)', 'SRSF7(Hs/Mm)', 'SRSF9(Hs/Mm)', 'TARDBP(Hs/Mm)', 'TIA1(Hs/Mm)', 'U2AF2(Hs/Mm)', 'YBX1(Hs/Mm)', 'ZC3H10(Hs/Mm)', 'ZCRB1(Hs/Mm)', 'ZNF638(Hs/Mm)']
    dfs_collection = []

    zscore_collection = []  # holds 'zscore_list(s)' of all flanks
    pval_collection = []

    # STEP 1: Read RBPmap results into a dataframe
    for current_dir in result_dirs:

        file = os.path.join(current_dir, 'Predictions.txt')

        # temp files
        subdir = file.split('Predictions.txt')[0]
        new_file = os.path.join(subdir, 'Predictions_new.txt')

        # list of proteins with a putative binding site in the current flank
        proteins_file = os.path.join(subdir, 'proteins.txt')

        try:

            os.system('grep -vwE "(Protein:)" ' + file + ' > ' + new_file)  # removing lines with protein names
            os.system('sed -i "1,8d" ' + new_file)
            os.system(
                'grep -oP "(?<=Protein: ).*" ' + file + ' > ' + proteins_file)  # "Protein: BRUNOL4" -> BRUNOL4

            result_df = pd.read_fwf(new_file, skip_blank_lines=True)
            result_df = result_df[['Genomic Coordinate', 'Motif', 'Z-score', 'P-value']]

            # STEP 2: Mark each Z score with its corresponding protein name

            # indices of first and last occurrences of each protein's Zscores
            idx = result_df.index[result_df['Motif'] == 'Motif'].tolist()
            idx.insert(0, -1)
            idx.append(len(result_df))

            # list of proteins with binding sites in current flank
            with open(proteins_file) as f:
                binding_proteins = f.read().splitlines()

            # construct column with appropriate num repeats of proteins
            protein_col = []
            i = 1
            for p in binding_proteins:
                if i <= len(binding_proteins):
                    for j in range(idx[i] - idx[i - 1] - 1):
                        protein_col.append(p)

                    i += 1
                else:
                    break

            # remove rows with text
            result_df = result_df[result_df.Motif != 'Motif']
            result_df.reset_index(drop=True, inplace=True)

            # add protein column
            result_df['protein'] = protein_col
            dfs_collection.append(result_df)

            # STEP 3: get largest Z score of each protein for the current flank
            zscore_list = []
            pval_list = []
            for p in proteins:
                protein_df = result_df.loc[result_df['protein'] == p].copy()
                if protein_df.empty:  # protein has no binding site in flank
                    zscore_list.append(0)
                    pval_list.append(None)
                else:  # assign z score based on significant adj p value
                    protein_df.loc[:, 'P-value'] = pd.to_numeric(protein_df['P-value'])
                    protein_df = adjust_pvalue(protein_df, 'P-value')
                    protein_df = protein_df[protein_df['adj_pval'] <= 0.05]
                    minidx = protein_df['P-value'].idxmin()
                    significant_z = protein_df.loc[minidx]['Z-score']
                    significant_p = protein_df.loc[minidx]['P-value']
                    zscore_list.append(significant_z)
                    pval_list.append(significant_p)

            # no protein has binding site in flanks
            if zscore_list.count(0) == len(zscore_list):
                zscore_list = [None] * len(zscore_list)

            zscore_collection.append(zscore_list)
            pval_collection.append(pval_list)

            # delete temp files
            os.remove(proteins_file)
            os.remove(new_file)

        except:
            pass

    col_names = [p.split('(Hs/Mm)')[0] for p in proteins]
    zscore_collection.insert(0, col_names)
    pval_collection.insert(0,  col_names)

    with open(f"{proc}_0_Files/FilteredZscores_" + type + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(zscore_collection)

    with open(f"{proc}_0_Files/FilteredPvalues_" + type + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(pval_collection)


def main(args):
    proc = args.process
    path = '../RBPmap'
    results_dirs = [x[0] for x in os.walk(path)]
    epi_dirs = [p for p in results_dirs if 'resultsrbp_input_epi' in p and 'sequence' in p]
    nonepi_dirs = [p for p in results_dirs if 'resultsrbp_input_nonepi' in p and 'sequence' in p]

    process_results(epi_dirs, 'epi', proc)
    process_results(nonepi_dirs, 'nonepi', proc)

    print('processed')


if __name__ == "__main__":    
      p = get_argument_parser()
      args = p.parse_args()
      main(args)

