import pandas as pd
import json
import natsort
from glob import glob
import matplotlib.pyplot as plt


def feature_matrix(filename1, filename2, filename3, weighted=False):

    prefix = '0_Files/'
    features = pd.read_csv(prefix+filename1, delimiter='\t')
    rbp = pd.read_csv(prefix+filename2, delimiter=',')

    if not weighted: # TO DO: NAME RBP columns in post-rbp
        rbp.columns = ['A1CF', 'ANKHD1', 'BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'CNOT4', 'CPEB2', 'CPEB4', 'DAZAP1', 'ENOX1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'G3BP2', 'HNRNPA1', 'HNRNPA1', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPA2B1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPCL1', 'HNRNPF', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRPLL', 'HuR', 'IGF2BP2', 'IGF2BP3', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'LIN28A', 'MATR3', 'MBNL1', 'MSI1', 'PABPC1', 'PABPC3', 'PABPC4', 'PABPC5', 'PABPN1', 'PCBP1', 'PCBP2', 'PPRC1', 'PTBP1', 'PUM2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM41', 'RBM42', 'RBM45', 'RBM46', 'RBM5', 'RBM6', 'RBM8A', 'RBMS1', 'RBMS3', 'SAMD4A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'SRSF9', 'TARDBP', 'TIA1', 'TUT1', 'U2AF2', 'YBX1', 'YBX2', 'ZC3H10', 'ZC3H14', 'ZCRB1', 'ZNF638']


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
        rbp.to_csv('0_Files/flanks_rbp_' + name + '.csv', sep='\t', index=False)

    cols = features.columns.tolist()
    cols.remove('gene_id')
    features = pd.merge(features[cols], rbp, on='flanks', how='outer')
    features.fillna(0, inplace=True)  # non-epigene flanks with no annotated peaks
    features['gene_id'] = features['gene_id'].str.split('.').str[0]  # ENSG00000116691.11 -> ENSG00000116691

    names = pd.read_csv('0_Files/GeneID_Name.csv', delimiter='\t')
    names.columns = ['gene_id', 'gene']
    features = pd.merge(features, names, on='gene_id')

    # if '_epi' in filename1:
    #     clusters = json.load(open("clusters.txt"))
    #     for label, cluster in clusters.items():
    #         file = 'clustermatrix' + str(label) + '.csv'
    #         cluster_features = features[features['gene_id'].isin(cluster)]
    #         del cluster_features['gene_id']
    #         cols = cluster_features.columns.tolist()
    #         cols = cols[-1:] + cols[:-1]
    #         cluster_features = cluster_features[cols]
    #         cluster_features.to_csv(file, sep='\t', index=False)

    del features['gene_id']
    cols = features.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    features = features[cols]

    with open('0_Files/' + name+'genes.txt', 'w') as f:
        for item in set(features['gene'].values.tolist()):
            f.write("%s\n" % item)

    if weighted:
        features.to_csv('0_Files/features_scaled_' + name + '.csv', sep='\t', index=False)
    else:
        features.to_csv('0_Files/features_' + name + '.csv', sep='\t', index=False)



    # plt.hist(features['mean_dpsi_per_lsv_junction'])
    # plt.title('Distribution of dPSI values')
    # plt.show()


if __name__ == "__main__":


    dPSI_Mval_files = ['dPSI_Mval_epi.csv', 'dPSI_Mval_nonepi.csv']
    Zscore_files = ['FilteredZscores_epi.csv', 'FilteredZscores_nonepi.csv']
    logFCZscore_files = ['WeightedZscores_epi.csv', 'WeightedZscores_nonepi.csv']
    query_files = ['query_flanks_epi.csv', 'query_flanks_nonepi.csv']

    for i in range(len(query_files)):
        feature_matrix(dPSI_Mval_files[i], Zscore_files[i], query_files[i], weighted=False)

    for i in range(len(query_files)):
        feature_matrix(dPSI_Mval_files[i], logFCZscore_files[i], query_files[i], weighted=True)





