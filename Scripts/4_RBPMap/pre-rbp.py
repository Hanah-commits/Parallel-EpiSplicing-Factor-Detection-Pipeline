import pandas as pd
import os

flanks = pd.read_csv('0_Files/all_flanks.csv', delimiter='\t')
epigene_flanks = pd.read_csv('0_Files/dPSI_Mval_epi.csv', delimiter='\t')
nonepigene_flanks = pd.read_csv('0_Files/dPSI_Mval_nonepi.csv', delimiter='\t')

epi_flanks = flanks[flanks['flanks'].isin(list(set(epigene_flanks['flanks'].values)))]
nonepi_flanks = flanks[flanks['flanks'].isin(list(set(nonepigene_flanks['flanks'].values)))]

input_files = []

# Get input sequences for RBPmap
n = 0
i = 1
while n < len(epi_flanks):
        name = '0_Files/rbp_input_epi'+str(i)+'.csv'
        input_files.append(name)
        if n+5000 <= len(epi_flanks):
                epi_flanks[['seqid', 'flanks', 'strand']].iloc[n:n+5000].to_csv(name, index=False, sep=':', header=False)
                n += 5000
        else:
                epi_flanks[['seqid', 'flanks', 'strand']].iloc[n:len(epi_flanks)+1].to_csv(name, index=False, sep=':', header=False)
                break
        i += 1

# Get input sequences for RBPmap
n = 0
i = 1
while n < len(nonepi_flanks):
        name = '0_Files/rbp_input_nonepi'+str(i)+'.csv'
        input_files.append(name)
        if n+5000 <= len(nonepi_flanks):
                nonepi_flanks[['seqid', 'flanks', 'strand']].iloc[n:n+5000].to_csv(name, index=False, sep=':', header=False)
                n += 5000
        else:
                nonepi_flanks[['seqid', 'flanks', 'strand']].iloc[n:len(nonepi_flanks)+1].to_csv(name, index=False, sep=':', header=False)
                break
        i += 1

# # get flanks for feature matrix preparation step
epi_flanks[['gene_id', 'flanks']].to_csv('0_Files/query_flanks_epi.csv', sep='\t', index=False)
nonepi_flanks[['gene_id', 'flanks']].to_csv('0_Files/query_flanks_nonepi.csv', sep='\t', index=False)

input_files = [os.getcwd() + '/' + file for file in input_files]
with open('0_Files/input.txt', 'w') as f:
    for item in input_files:
        f.write("%s\n" % item)
