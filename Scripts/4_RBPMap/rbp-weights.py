import pandas as pd

proteins = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTB3', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']
weights = pd.read_csv('0_Files/rbp_logFC.csv', delimiter='\t')
weights = weights.set_index('rbp').T.to_dict('list')  # {'TARDBP': [-15400.0],  'LIN28A': [-162471.666666667]}

types = ['epi', 'nonepi']
for type in types:
    file = '0_Files/FilteredZscores_'+type+'.csv'
    rbp_scores = pd.read_csv(file, delimiter=',')
    rbp_scores = rbp_scores.loc[:, ~rbp_scores.columns.duplicated()]
    
    # drop_proteins = ['HNRPLL']
    # rbp_scores = rbp_scores.drop(drop_proteins, axis=1)

    for protein in rbp_scores.columns:
        if protein in weights:
            wt = weights[protein][0]
            rbp_scores[protein] = rbp_scores[protein].apply(lambda x: x * abs(weights[protein][0]))  # abs log FC
        else:
            print(protein)

    rbp_scores.to_csv('0_Files/WeightedZscores_'+type+'.csv', sep=',', index=False)


