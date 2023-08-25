import pandas as pd
import numpy as np
import json
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("output_dir")
p.add_argument("--process", "-p",
    help="The name of the process")

args = p.parse_args()
proc = args.process

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


prefix = args.output_dir + 'MANorm/'

with open('paths.json') as f:
    data = json.load(f)
d = data[proc]

hms = d["Histone modifications"]
peaksfiles = [hm+'_flanks.bed' for hm in hms]
peak_dfs = []

flanks = pd.read_csv(f'{proc}_0_Files/filtered_flanks.bed', delimiter='\t', header=None)
flanks.columns = ['seqid', 'start', 'stop']
flanks['flanks'] = flanks[['start', 'stop']].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
flanks.drop_duplicates(inplace=True)

for file in peaksfiles:

    hm = file.split('_flanks')[0]
    peaks = pd.read_csv(prefix+file, delimiter='\t', header=None)
    peaks.drop([3, 8, 10, 11, 12], axis=1, inplace=True)

    # assign 0 to flanks that have no peaks
    peaks.replace([-1, '.'], [0, 0], inplace=True)
    peaks[3] = peaks[[1, 2]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
    peaks.set_axis(['seqid', 'flank0', 'flank1', 'peak0', 'peak1', 'summit', 'M_value', 'p_value', 'flanks'],
                   axis=1, inplace=True)
    peaks['M_value_abs'] = pd.to_numeric(peaks['M_value']).abs()
    peaks.loc[peaks['M_value'] == 0, 'p_value'] = np.NaN # null pvalues if no peak in flank

    # keep M-values of peaks with adj P-value <= 0.05
    peaks['p_value'] = pd.to_numeric(peaks['p_value'])
    peaks = adjust_pvalue(peaks, col='p_value')
    peaks.loc[pd.to_numeric(peaks['adj_pval']) > 0.05, 'M_value_abs'] = 0 

    # # get all peaks that belong to each flank
    flank_peaks_group = peaks.groupby(['seqid', 'flanks'])['M_value_abs'] \
        .apply(lambda val: ','.join(str(v) for v in val)).reset_index()

    # # FILTER 1: if flank has 1+ peaks, keep peak with highest abs M-value
    flank_peaks_group['M_value_abs'] = flank_peaks_group['M_value_abs'].str.split(',')  # string -> list of strings
    flank_peaks_group['max_' + hm] = flank_peaks_group['M_value_abs'].apply(lambda x: max(map(float, x)))  # max MValue
    # # flank_peaks_group['#peaks_'+hm] = flank_peaks_group['M_value'].apply(lambda x: len(x))  # no of peaks/ flank

    # # get the corresponding peak for each flank's max M-value
    flank_peaks_group = pd.merge(flank_peaks_group[['flanks', 'max_' + hm]], peaks, on=['flanks'],
                                 how='inner')
    # # flank_peaks_group = pd.merge(flank_peaks_group[['flanks', 'max_'+hm, '#peaks_'+hm]], peaks, on=['flanks'], how='inner')
    flank_peaks_group = flank_peaks_group[flank_peaks_group['M_value_abs'] == flank_peaks_group['max_' + hm]]

    # FILTER 4: If flank has 1+ peaks with same max |Mvalue|, keep one
    flank_peaks_group.drop_duplicates(subset='flanks', keep='first', inplace=True)

    flank_peaks_group.rename(columns={'M_value': hm}, inplace=True)
    peak_dfs.append(flank_peaks_group[['flanks', hm]])  # 'max_' + hm, '#peaks_'+hm]])

peak_dfs = [df.set_index('flanks') for df in peak_dfs]
peak_dfs = pd.concat(peak_dfs, axis=1)
peak_dfs = peak_dfs.loc[~(peak_dfs==0).all(axis=1)]

peak_dfs.to_csv(f'{proc}_0_Files/Filtered_MValues.csv', sep='\t')
