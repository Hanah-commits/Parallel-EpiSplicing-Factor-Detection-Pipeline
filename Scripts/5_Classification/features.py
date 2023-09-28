import pandas as pd
import sys
import matplotlib.pyplot as plt
from argparse import ArgumentParser

# Get the process name, use it in the output directory
def get_argument_parser():
    p = ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("weights")
    p.add_argument("--process", "-p",
        help="The name of the process")
    return p


def feature_matrix(filename1, filename2, filename3, output_dir, proc, i, weighted=False):

    prefix = f'{proc}_0_Files/'
    features = pd.read_csv(prefix+filename1, delimiter='\t')
    rbp = pd.read_csv(prefix+filename2, delimiter=',')

    flanks = pd.read_csv(prefix+filename3, delimiter='\t')
    rbp = pd.concat([flanks, rbp], axis=1)

    name = ''
    if i == 0:
        name = 'epi'
    else:
        name = 'nonepi'

    # # FILTER 1: If RBP motif has 1+ PSSMs, keep only one
    rbp = rbp.loc[:, ~rbp.columns.duplicated()]

    if weighted:
        rbp.to_csv(output_dir+'/flanks_rbp_' + name + '.csv', sep='\t', index=False)

    cols = features.columns.tolist()
    cols.remove('gene_id')
    features = pd.merge(features[cols], rbp, on='flanks', how='outer')
    features.fillna(0, inplace=True)  # non-epigene flanks with no annotated peaks
    features['gene_id'] = features['gene_id'].str.split('.').str[0]  # ENSG00000116691.11 -> ENSG00000116691

    names = pd.read_csv('HelperFunctions/GeneID_Name.csv', delimiter='\t')
    names.columns = ['gene_id', 'gene']
    features = pd.merge(features, names, on='gene_id')

    del features['gene_id']
    cols = features.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    features = features[cols]

    with open(f'{prefix}{name}genes.txt', 'w') as f:
        for item in set(features['gene'].values.tolist()):
            f.write("%s\n" % item)

    if weighted:
        features.to_csv(f'{prefix}features_scaled_' + name + '.csv', sep='\t', index=False)
    else:
        features.to_csv(f'{prefix}features_' + name + '.csv', sep='\t', index=False)



    # plt.hist(features['mean_dpsi_per_lsv_junction'])
    # plt.title('Distribution of dPSI values')
    # plt.show()


def main(args):
    proc = args.process
    output_dir = args.output_dir
    if args.weights == 'True':
        weights = True
    else:
        weights = False

    dPSI_Mval_files = ['dPSI_Mval_epi.csv', 'dPSI_Mval_nonepi.csv']
    Zscore_files = ['FilteredZscores_epi.csv', 'FilteredZscores_nonepi.csv']
    logFCZscore_files = ['WeightedZscores_epi.csv', 'WeightedZscores_nonepi.csv']
    query_files = ['query_flanks_epi.csv', 'query_flanks_nonepi.csv']

    for i in range(len(query_files)):
        feature_matrix(dPSI_Mval_files[i], Zscore_files[i], query_files[i], output_dir, proc, i, weighted=False)

    if weights:
        
        for i in range(len(query_files)):
            feature_matrix(dPSI_Mval_files[i], logFCZscore_files[i], query_files[i], output_dir, proc, i, weighted=True)


if __name__ == "__main__":    
      p = get_argument_parser()
      args = p.parse_args()
      main(args)
