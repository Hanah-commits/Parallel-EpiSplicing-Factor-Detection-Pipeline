import pandas as pd


ref = pd.read_csv('gencode.v39.annotation.gff3', delimiter='\t', skiprows=7, header=None)
ref = ref[ref[2] == 'exon']

# ref['diff_start'] = ref[3].diff()
# ref['diff_stop'] = ref[4].diff()
# diff = [1, -1]
# ref = ref[~((ref['diff_start'].isin(diff)) & (ref['diff_stop'].isin(diff)))]

start_df = ref.drop(ref.columns[4], 1)
stop_df = ref.drop(ref.columns[3], 1)

start_df[4] = pd.to_numeric(start_df[3]) + 200  # extend 200bp 5' direction
start_df[3] = pd.to_numeric(start_df[3]) - 200  # extend 200bp 3' direction

stop_df[3] = pd.to_numeric(stop_df[4]) - 200  # extend 200bp 5' direction
stop_df[4] = pd.to_numeric(stop_df[4]) + 200  # extend 200bp 3' direction

new_ref = pd.concat([start_df, stop_df]).sort_index(kind='merge')

new_ref[3] = new_ref[3].astype(int)
new_ref[4] = new_ref[4].astype(int)

new_ref[[0, 3, 4]].to_csv('flanks' + str(200) + '.bed', sep='\t', index=False, header=False)

