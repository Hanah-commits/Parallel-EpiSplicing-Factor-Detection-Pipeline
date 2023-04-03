import pandas as pd
import json

with open('paths.json') as f:
    d = json.load(f)

ref = d['Reference genome']

ref = pd.read_csv(ref, delimiter='\t',header=None)[[0,1,2]]

flanks = [50, 100, 200]

for flank in flanks:

    start_df = ref.drop(ref.columns[2], 1)
    stop_df = ref.drop(ref.columns[1], 1)


    start_df[2] = pd.to_numeric(start_df[1]) + flank  # extend 200bp 5' direction
    start_df[1] = pd.to_numeric(start_df[1]) - flank  # extend 200bp 3' direction


    stop_df[1] = pd.to_numeric(stop_df[2]) - flank  # extend 200bp 5' direction
    stop_df[2] = pd.to_numeric(stop_df[2]) + flank  # extend 200bp 3' direction

    new_ref = pd.concat([start_df, stop_df]).sort_index(kind='merge')

    new_ref[1] = new_ref[1].astype(int)
    new_ref[2] = new_ref[2].astype(int)

    print(new_ref)
    new_ref.to_csv('0_Files/flanks' + str(flank) + '.bed', sep='\t', index=False, header=False)