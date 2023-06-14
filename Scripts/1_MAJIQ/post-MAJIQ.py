import pandas as pd
import numpy as np
from pandas import Series
import sys


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
    df['adj_pval'] = adj_pval
    return df


def p_adjust_bh(p):
    ## multiple hypothesis testing
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


# STEP 1: Extract required columns and split individual dpsi values, their probabilities and junction coords

# Keep relevant columns
file = sys.argv[1]+'MAJIQ/majiq_output'
voila = pd.read_csv(file, delimiter='\t', skiprows=10)
col_list = ['gene_id', 'lsv_id', 'seqid', 'mean_dpsi_per_lsv_junction', 'probability_changing', 'junctions_coords', 'num_exons', 'strand'] #, 'exons_coords']
voila = voila[col_list]

# FILTER 1: remove LSVs with 2 exons
voila["num_exons"] = voila["num_exons"].replace('na' ,'0')
voila["num_exons"] = pd.to_numeric(voila["num_exons"])
voila = voila[voila['num_exons'] > 2]

# split column values to multiple lines
coord = voila['junctions_coords'].str.split(';').apply(Series, 1).stack()
coord.index = coord.index.droplevel(-1)
s1 = voila.junctions_coords.str.split(';', expand=True).stack().str.strip().reset_index(level=1, drop=True)

dpsi = voila['mean_dpsi_per_lsv_junction'].str.split(';').apply(Series, 1).stack()
dpsi.index = dpsi.index.droplevel(-1)
s2 = voila.mean_dpsi_per_lsv_junction.str.split(';', expand=True).stack().str.strip().reset_index(level=1, drop=True)

prob = voila['probability_changing'].str.split(';').apply(Series, 1).stack()
prob.index = prob.index.droplevel(-1)
s3 = voila.probability_changing.str.split(';', expand=True).stack().str.strip().reset_index(level=1, drop=True)

new_cols = pd.concat([coord, dpsi, prob], axis =1, keys=['junctions_coords', 'mean_dpsi_per_lsv_junction', 'probability_changing'])
voila = voila.drop(['junctions_coords', 'mean_dpsi_per_lsv_junction', 'probability_changing'], axis=1).join(new_cols).reset_index(drop=True)

voila = voila[col_list]
# skipping nan -> 99464127-nan
voila = voila[~voila['junctions_coords'].str.contains("nan")]

# FILTER 2: Drop rows with changing probability > 0.05
voila['pval'] = 1- pd.to_numeric(voila['probability_changing'])
voila = adjust_pvalue(voila, col='pval')
voila = voila[pd.to_numeric(voila['adj_pval']) <= 0.05]


# STEP 3: Get the source and target indices of splice junctions
voila[['source', 'target']] = voila['junctions_coords'].str.split('-', 1, expand=True)

del(voila['junctions_coords'])
voila_temp = voila.copy()
del(voila_temp['target'])
del(voila['source'])

voila.rename(columns={'target': 'junction0'}, inplace = True)
voila_temp.rename(columns={'source': 'junction0'}, inplace=True)

voila = pd.concat([voila_temp, voila]).sort_index(kind='merge')

keep_cols = ['seqid', 'junction0']
majiq_bed = voila[keep_cols]
majiq_bed = majiq_bed.drop_duplicates()
majiq_bed['junction1'] = pd.to_numeric(majiq_bed['junction0']) + 1  # to fit bedtools input requirements
majiq_bed.to_csv('0_Files/majiq.bed', index=False, sep='\t', header=False)  # input for bedtools intersect
voila.to_csv('0_Files/majiq_junctions.csv', index=False, sep='\t', header=True)

