import pandas as pd
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("output_dir")
p.add_argument("weights")
p.add_argument("--process", "-p",
    help="The name of the process")
args = p.parse_args()

proc = args.process

prefix = f'{proc}_0_Files/'
unscaled = ['features_epi.csv', 'features_nonepi.csv', 'all_features.csv']
scaled = ['features_scaled_epi.csv', 'features_scaled_nonepi.csv', 'all_features_scaled.csv']

if args.weights == 'True':
    types = [unscaled, scaled]
else:
    types = [unscaled]

for file_list in types:
    epi_features = pd.read_csv(prefix+file_list[0], delimiter='\t')
    nonepi_features = pd.read_csv(prefix+file_list[1], delimiter='\t')

    epi_features['label'] = 'epigene'
    nonepi_features['label'] = 'non-epigene'

    all_features = pd.concat([epi_features, nonepi_features], axis=0)
    all_features.drop(['mean_dpsi_per_lsv_junction', 'gene', 'flanks'], axis=1, inplace=True)
    col = all_features.pop("label")
    all_features.insert(0, col.name, col)

    all_features.to_csv(prefix+file_list[2], sep='\t', index=False)


