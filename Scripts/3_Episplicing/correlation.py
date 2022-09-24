import os
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import rankdata
import matplotlib.pyplot as plt


def pearsonr_coeff(x, y):
    return pearsonr(x, y)[0]


def pearsonr_pval(x, y):
    return pearsonr(x, y)[1]


def adjust_pvalue(df):
    pval_cols = df.columns.tolist()[1:]  # skipping gene-id column
    new_cols = []
    col_names = []
    for col in pval_cols:

        # get indices of null values
        na_idx = df[df[col].isnull()].index.tolist()

        # adjust non-null p values
        pvals = df[col].values.tolist()
        pvals = [x for x in pvals if str(x) != 'nan']
        adj_pval = p_adjust_bh(pvals).tolist()

        # insert null at original indices
        for idx in na_idx:
            adj_pval.insert(idx, None)

        new_cols.append(adj_pval)
        col_names.append(col + '_adj')

    # adjusted p values as new df
    df1 = pd.DataFrame(columns=col_names)
    df1['gene_id'] = df['gene_id'].values.tolist()
    for i in range(len(new_cols)):
        df1[col_names[i]] = new_cols[i]

    return df1


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def fdr(p_vals):
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1

    return fdr


def plot_histogram(df, columns, status=0):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
    col = 0
    for i in range(2):
        for j in range(2):
            ax = axes[i][j]
            ax.hist(df[columns[col]], bins=10, color='blue', alpha=0.5, label='{}'.format(columns[col]))
            ax.set_xlabel('P-values')
            ax.set_ylabel('count')
            ax.set_ylim([0, 250])
            leg = ax.legend(loc='upper left')
            leg.draw_frame(False)
            col += 1
    if status == 0:
        plt.suptitle('Non-adjusted P-value histogram')
    else:
        plt.suptitle('FDR adjusted P-value histogram')

    plt.show()


def find_epigenes(df):
    pval_cols = df.columns.tolist()[:-1]  # skipping gene-id column
    epi_genes = []
    for col in pval_cols:
        epi_genes.extend((df.loc[df[col] <= 0.05, 'gene_id']).values.tolist())

    epi_genes = sorted(list(set(epi_genes)))
    return epi_genes


with open('paths.json') as f:
    d = json.load(f)

hms = d["Histone modifications"]

dPSI = pd.read_csv('0_Files/Filtered_dPSI.csv', delimiter='\t')
peaks = pd.read_csv('0_Files/Filtered_MValues.csv', delimiter='\t')
dPSI.drop_duplicates(inplace=True)
peaks.drop_duplicates(inplace=True)
flanks = pd.merge(dPSI, peaks, how="outer")

# FILTER 1: drop genes with less than 10 flanks
flanks = flanks[flanks.groupby('gene_id').gene_id.transform(len) > 2]

# # FILTER 2: drop genes with dPSI values but no peak
cols = ['gene_id'] + hms
grouped = flanks[cols].groupby('gene_id')
non_epi = []
for gene, group in grouped:
    if group[hms].isnull().all().all():
        non_epi.append(gene)

filtered_flanks = flanks[~flanks['gene_id'].isin(non_epi)]
filtered_flanks.set_index('flanks', inplace=True)

# if one or more flanks of a gene have less than four peaks
filtered_flanks.fillna(0, inplace=True)

# num = flanks.groupby(['gene_id']).size().reset_index(name='no_flanks') # # number of flanks for each gene
# num = num[num['no_flanks'] > 9]
# print(num)

# Step 1: Calculate correlation coefficients and p values
pval = filtered_flanks.groupby('gene_id').corr(method=pearsonr_pval)

# Step 2: Keep only relevant correlations
pval.to_csv('0_Files/pvals.csv', sep='\t')
new_pvals = pd.read_csv('0_Files/pvals.csv', delimiter='\t')
new_pvals.drop(['Unnamed: 1', 'mean_dpsi_per_lsv_junction'], axis=1, inplace=True)

# dropping p-values of hm-hm correlations
new_pvals = new_pvals.iloc[::5, :]

# drop genes where no dPSI-HM correlations exist
new_pvals.dropna(subset=hms, how='all', inplace=True)

# cleanup
new_pvals.reset_index(inplace=True)
del new_pvals['index']

# Step 3: Visualise p-value distribution
# plot_histogram(new_pvals, columns=["H3K4me3", "H3K27me3", "H3K9me3", "H3K27ac"])

# Step 4: Adjust the p values using Benjamini-Hochberg method
adj_pvals = adjust_pvalue(new_pvals)

# plot_histogram(adj_pvals, columns=["H3K4me3_adj", "H3K27me3_adj", "H3K9me3_adj", "H3K27ac_adj"])
# print(adj_pvals[adj_pvals['gene_id'] == 'ENSG00000154134.15'])

# Step 5: Obtain genes where adjusted_pval < 0.05
epigenes = find_epigenes(adj_pvals)
print('Epigenes ', len(epigenes))
print('Non-Epigenes ', len(non_epi))


# get flanks of epispliced genes
flanks[flanks['gene_id'].isin(epigenes)].to_csv('0_Files/dPSI_Mval_epi.csv', sep='\t', index=False)

# # get flanks of non-epispliced genes
flanks[flanks['gene_id'].isin(non_epi)].to_csv('0_Files/dPSI_Mval_nonepi.csv', sep='\t', index=False)

#clean-up
os.remove('0_Files/pvals.csv')
os.remove('0_Files/Filtered_dPSI.csv')
os.remove('0_Files/Filtered_MValues.csv')
