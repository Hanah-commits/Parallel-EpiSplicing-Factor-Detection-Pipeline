import pandas as pd

proteins = ['A1CF(Hs/Mm)', 'ANKHD1(Hs/Mm)', 'BRUNOL4(Hs/Mm)', 'BRUNOL5(Hs/Mm)', 'BRUNOL6(Hs/Mm)', 'CNOT4(Hs/Mm)',
            'CPEB2(Hs/Mm)', 'CPEB4(Hs/Mm)', 'DAZAP1(Hs/Mm)', 'ENOX1(Hs/Mm)', 'ESRP2(Hs/Mm)', 'FMR1(Hs/Mm)',
            'FUS(Hs/Mm)', 'FXR1(Hs/Mm)', 'FXR2(Hs/Mm)', 'G3BP2(Hs/Mm)', 'HNRNPA1(Hs/Mm)', 'HNRNPA1(Hs/Mm)',
            'HNRNPA1(Hs/Mm)', 'HNRNPA1L2(Hs/Mm)', 'HNRNPA2B1(Hs/Mm)', 'HNRNPA2B1(Hs/Mm)', 'HNRNPA2B1(Hs/Mm)',
            'HNRNPC(Hs/Mm)', 'HNRNPCL1(Hs/Mm)', 'HNRNPF(Hs/Mm)', 'HNRNPF(Hs/Mm)', 'HNRNPH1(Hs/Mm)', 'HNRNPH2(Hs/Mm)',
            'HNRNPK(Hs/Mm)', 'HNRNPL(Hs/Mm)', 'HNRNPL(Hs/Mm)', 'HNRNPM(Hs/Mm)', 'HNRNPU(Hs/Mm)', 'HNRPLL(Hs/Mm)',
            'HuR(Hs/Mm)', 'IGF2BP2(Hs/Mm)', 'IGF2BP3(Hs/Mm)', 'KHDRBS1(Hs/Mm)', 'KHDRBS2(Hs/Mm)', 'KHDRBS3(Hs/Mm)',
            'LIN28A(Hs/Mm)', 'MATR3(Hs/Mm)', 'MBNL1(Hs/Mm)', 'MSI1(Hs/Mm)', 'PABPC1(Hs/Mm)', 'PABPC3(Hs/Mm)',
            'PABPC4(Hs/Mm)', 'PABPC5(Hs/Mm)', 'PABPN1(Hs/Mm)', 'PCBP1(Hs/Mm)', 'PCBP2(Hs/Mm)', 'PPRC1(Hs/Mm)',
            'PTBP1(Hs/Mm)', 'PUM2(Hs/Mm)', 'QKI(Hs/Mm)', 'RALY(Hs/Mm)', 'RBFOX1(Hs/Mm)', 'RBM24(Hs/Mm)', 'RBM28(Hs/Mm)',
            'RBM3(Hs/Mm)', 'RBM4(Hs/Mm)', 'RBM41(Hs/Mm)', 'RBM42(Hs/Mm)', 'RBM45(Hs/Mm)', 'RBM46(Hs/Mm)', 'RBM5(Hs/Mm)',
            'RBM6(Hs/Mm)', 'RBM8A(Hs/Mm)', 'RBMS1(Hs/Mm)', 'RBMS3(Hs/Mm)', 'SAMD4A(Hs/Mm)', 'SART3(Hs/Mm)',
            'SFPQ(Hs/Mm)', 'SNRNP70(Hs/Mm)', 'SNRPA(Hs/Mm)', 'SRSF1(Hs/Mm)', 'SRSF10(Hs/Mm)', 'SRSF2(Hs/Mm)',
            'SRSF7(Hs/Mm)', 'SRSF9(Hs/Mm)', 'SRSF9(Hs/Mm)', 'TARDBP(Hs/Mm)', 'TIA1(Hs/Mm)', 'TUT1(Hs/Mm)',
            'U2AF2(Hs/Mm)', 'YBX1(Hs/Mm)', 'YBX2(Hs/Mm)', 'ZC3H10(Hs/Mm)', 'ZC3H14(Hs/Mm)', 'ZCRB1(Hs/Mm)',
            'ZNF638(Hs/Mm)']
proteins = [p.replace('(Hs/Mm)', '').strip() for p in proteins]
weights = pd.read_csv('0_Files/rbp_logFC.csv', delimiter='\t')
weights = weights.set_index('rbp').T.to_dict('list')  # {'TARDBP': [-15400.0],  'LIN28A': [-162471.666666667]}

types = ['epi', 'nonepi']
for type in types:
    file = '0_Files/FilteredZscores_'+type+'.csv'
    rbp_scores = pd.read_csv(file, delimiter=',', header=None)
    rbp_scores.columns = proteins
    rbp_scores = rbp_scores.loc[:, ~rbp_scores.columns.duplicated()]

    # ## extra vela
    drop_proteins = ['HNRPLL']
    rbp_scores = rbp_scores.drop(drop_proteins, axis=1)

    for protein in rbp_scores.columns:
        if protein in weights:
            wt = weights[protein][0]
            rbp_scores[protein] = rbp_scores[protein].apply(lambda x: x * abs(weights[protein][0]))  # abs log FC
        else:
            print(protein)

    rbp_scores.to_csv('0_Files/WeightedZscores_'+type+'.csv', sep=',', index=False)


