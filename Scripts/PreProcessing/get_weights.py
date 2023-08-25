import pandas as pd
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("--process", "-p",
    help="The name of the process")

args = p.parse_args()
proc = args.process
tmp_out_dir = proc + '_0_Files'

logFC_Scores = pd.read_csv(f'{tmp_out_dir}/logFC.csv', delimiter='\t',
                           header=None, names=['Gene stable ID', 'logFC'])
names = pd.read_csv('HelperFunctions/GeneID_Name.csv', delimiter='\t')

rbp = ['A1CF', 'ANKHD1', 'CELF4', 'CELF5', 'CELF6', 'CNOT4', 'CPEB2', 'CPEB4', 'DAZAP1', 'ENOX1', 'ESRP2', 'FMR1',
       'FUS', 'FXR1', 'FXR2', 'G3BP2', 'HNRNPA1', 'HNRNPA1', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPA2B1',
       'HNRNPA2B1', 'HNRNPC', 'HNRNPCL1', 'HNRNPF', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPL',
       'HNRNPM', 'HNRNPU', 'HNRPLL', 'ELAVL1', 'IGF2BP2', 'IGF2BP3', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'LIN28A', 'MATR3',
       'MBNL1', 'MSI1', 'PABPC1', 'PABPC3', 'PABPC4', 'PABPC5', 'PABPN1', 'PCBP1', 'PCBP2', 'PPRC1', 'PTB3', 'PUM2',
       'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM41', 'RBM42', 'RBM45', 'RBM46', 'RBM5', 'RBM6',
       'RBM8A', 'RBMS1', 'RBMS3', 'SAMD4A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7',
       'SRSF9', 'SRSF9', 'TARDBP', 'TIA1', 'TUT1', 'U2AF2', 'YBX1', 'YBX2', 'ZC3H10', 'ZC3H14', 'ZCRB1', 'ZNF638']
names = names[names['Gene name'].isin(rbp)]
logFC_Scores['Gene stable ID'] = logFC_Scores['Gene stable ID'].str.split(
    '.').str[0]
logFC_Scores = pd.merge(logFC_Scores, names, on='Gene stable ID')

logFC_Scores['Gene name'] = logFC_Scores['Gene name'].replace({'CELF4': 'BRUNOL4',
                                                               'CELF5': 'BRUNOL5',
                                                               'CELF6': 'BRUNOL6',
                                                               'ELAVL1': 'HuR'})
logFC_Scores[['Gene name', 'logFC']].to_csv(
    f'{tmp_out_dir}/rbp_logFC.csv', sep='\t', index=False, header=['rbp', 'logFC'])
