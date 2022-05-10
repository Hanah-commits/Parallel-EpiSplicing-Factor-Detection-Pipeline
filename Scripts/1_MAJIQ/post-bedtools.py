import pandas as pd

junctions = pd.read_csv('0_Files/majiq_junctions.csv', delimiter='\t')
flanks = pd.read_csv('0_Files/majiq_flanks.bed', delimiter='\t', header=None)

# drop flanks that have no junction
flanks = flanks[flanks[4] != -1]

# merge the flanks df with the jns df
flanks[6] = flanks[[1, 2]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
flanks.drop([3, 5],  axis=1, inplace=True)
flanks.set_axis(['seqid', 'start', 'stop', 'junction0', 'flanks'], axis=1, inplace=True)

junctions['index'] = junctions.index
flank_jns = pd.merge(junctions, flanks, on=['junction0', 'seqid'])
flank_jns.sort_values(['index'], inplace=True)
flank_jns = flank_jns.reset_index(drop=True)
del(flank_jns['index'])
flank_jns['mean_dpsi_per_lsv_junction'] = pd.to_numeric(flank_jns['mean_dpsi_per_lsv_junction'])

# # get all junctions that belong to each flank
flank_jns_group = flank_jns.groupby(['seqid', 'flanks'])['mean_dpsi_per_lsv_junction']\
        .apply(lambda val: ','.join(str(v) for v in val)).reset_index()

# FILTER 3: if flank has 1+ junctions, keep junction with highest dPSI value
flank_jns_group['max_dPSI'] = flank_jns_group['mean_dpsi_per_lsv_junction'].str.split(',')\
        .apply(lambda x: max(map(float, x)))  # string -> list of strings -> list of floats -> max float


# get the corresponding junction for each flank's max dPSI value
flank_jns_group = pd.merge(flank_jns_group[['flanks', 'max_dPSI']], flank_jns, on=['flanks'], how='inner')
flank_jns_group = flank_jns_group[flank_jns_group['mean_dpsi_per_lsv_junction'] == flank_jns_group['max_dPSI']]

# FILTER 4: If flank has 1+ junctions with same max dPSI value, keep one
flank_jns_group.drop_duplicates(subset='flanks', keep='first', inplace=True)

# # bookkeeping
del(flank_jns_group['max_dPSI'])
flank_jns_group = flank_jns_group[['gene_id', 'lsv_id', 'seqid', 'junction0', 'mean_dpsi_per_lsv_junction',
        'probability_changing', 'flanks', 'start', 'stop', 'strand']]


# Get all filtered flanks
flank_jns_group.drop_duplicates().to_csv('0_Files/all_flanks.csv', sep='\t', index=False)

# Get new annotation file with flanks (to annotate MANorm peaks)
flank_jns_group[['seqid', 'start', 'stop']].drop_duplicates().to_csv('0_Files/filtered_flanks.bed', index=False, sep='\t', header=False)

# Get dPSI values of filtered flanks
flank_jns_group[['flanks', 'gene_id', 'mean_dpsi_per_lsv_junction']].drop_duplicates().to_csv('0_Files/Filtered_dPSI.csv', index=False, sep='\t')
