import pandas as pd
import matplotlib.pyplot as plt


def feature_matrix(filename1, filename2, filename3, weighted=False):
    features = pd.read_csv(filename1, delimiter='\t')
    rbp = pd.read_csv(filename2, delimiter=',')

    if not weighted: # TO DO: NAME RBP columns in post-rbp
        rbp.columns = ['A1CF', 'ANKHD1', 'BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'CNOT4', 'CPEB2', 'CPEB4', 'DAZAP1', 'ENOX1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'G3BP2', 'HNRNPA1', 'HNRNPA1', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPA2B1', 'HNRNPA2B1', 'HNRNPC', 'HNRNPCL1', 'HNRNPF', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HNRPLL', 'HuR', 'IGF2BP2', 'IGF2BP3', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'LIN28A', 'MATR3', 'MBNL1', 'MSI1', 'PABPC1', 'PABPC3', 'PABPC4', 'PABPC5', 'PABPN1', 'PCBP1', 'PCBP2', 'PPRC1', 'PTBP1', 'PUM2', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM41', 'RBM42', 'RBM45', 'RBM46', 'RBM5', 'RBM6', 'RBM8A', 'RBMS1', 'RBMS3', 'SAMD4A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'SRSF9', 'TARDBP', 'TIA1', 'TUT1', 'U2AF2', 'YBX1', 'YBX2', 'ZC3H10', 'ZC3H14', 'ZCRB1', 'ZNF638']

    flanks = pd.read_csv(filename3, delimiter='\t')
    rbp = pd.concat([flanks, rbp], axis=1)

    name = ''
    if i == 0:
        name = 'epi'
    else:
        name = 'nonepi'

    # # FILTER 1: If RBP motif has 1+ PSSMs, keep only one
    rbp = rbp.loc[:, ~rbp.columns.duplicated()]


    cols = features.columns.tolist()
    cols.remove('gene_id')
    features = pd.merge(features[cols], rbp, on='flanks', how='outer')
    features.fillna(0, inplace=True)  # non-epigene flanks with no annotated peaks
    features['gene_id'] = features['gene_id'].str.split('.').str[0]  # ENSG00000116691.11 -> ENSG00000116691

    names = pd.read_csv('0_Files/GeneID_Name.csv', delimiter='\t')
    names.columns = ['gene_id', 'gene']
    features = pd.merge(features, names, on='gene_id')

    del features['gene_id']
    cols = features.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    features = features[cols]

    features.to_csv('pvals_rbp' + name + '.csv', sep='\t', index=False)


if __name__ == "__main__":

    dPSI_Mval_files = ['0_Files/dPSI_Mval_epi.csv', '0_Files/dPSI_Mval_nonepi.csv']
    Zscore_files = ['0_Files/FilteredPvalues_epi.csv', '0_Files/FilteredPvalues_nonepi.csv']
    query_files = ['0_Files/query_flanks_epi.csv', '0_Files/query_flanks_nonepi.csv']

    for i in range(len(query_files)):
        feature_matrix(dPSI_Mval_files[i], Zscore_files[i], query_files[i])





