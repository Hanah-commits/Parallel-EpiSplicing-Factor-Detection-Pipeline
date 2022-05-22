import pandas as pd
import json

with open('paths.json') as f:
    d = json.load(f)

ref = d['Reference genome']
ref = pd.read_csv(ref, delimiter='\t', skiprows=7, header=None)
ref = ref[ref[2] == 'exon']

flanks = [50, 100, 200]

for flank in flanks:

    start_df = ref.drop(ref.columns[4], 1)
    stop_df = ref.drop(ref.columns[3], 1)

    start_df[4] = pd.to_numeric(start_df[3]) + flank  # extend 200bp 5' direction
    start_df[3] = pd.to_numeric(start_df[3]) - flank  # extend 200bp 3' direction

    stop_df[3] = pd.to_numeric(stop_df[4]) - flank  # extend 200bp 5' direction
    stop_df[4] = pd.to_numeric(stop_df[4]) + flank  # extend 200bp 3' direction

    new_ref = pd.concat([start_df, stop_df]).sort_index(kind='merge')

    new_ref[3] = new_ref[3].astype(int)
    new_ref[4] = new_ref[4].astype(int)

    new_ref[[0, 3, 4]].to_csv('0_Files/flanks' + str(flank) + '.bed', sep='\t', index=False, header=False)


