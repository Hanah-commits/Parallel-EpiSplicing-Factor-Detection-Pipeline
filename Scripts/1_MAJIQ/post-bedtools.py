import pandas as pd
import os
import json
from argparse import ArgumentParser

# Get the process name, use it in the output directory

print('STEP 5 is in progress')

p = ArgumentParser()
p.add_argument("--process", "-p",
help="The name of the process")
args = p.parse_args()
proc = args.process

tmp_out_dir = proc + '_0_Files'


flank_lens = [50, 100, 200]
junctions = pd.read_csv(f'{tmp_out_dir}/majiq_junctions.csv', delimiter='\t')
with open('paths.json') as f:
    data = json.load(f)
d = data[proc]

fasta = d['Reference fasta']
ref_genome= fasta+".fai"

annot = []
for length in flank_lens:

        # adjust flank length to 200bp +- without exceeding chromosome bounds
        if length < 200:
                file = "majiq_flanks" + str(length) + ".bed"
                adjust_size = str(200 -length)

                # separate start,stop flank coords
                os.system(f"cut -f 1-3 {tmp_out_dir}/{file} > {tmp_out_dir}/coords.bed")

                # adjust flank boundaries
                os.system(f"bedtools slop -i {tmp_out_dir}/coords.bed" + " -g " + ref_genome + " -b " + adjust_size + f" > {tmp_out_dir}/coords_adjusted.bed")

                # replace flank coords with adjusted coords
                os.system("awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' " + f"{tmp_out_dir}/coords_adjusted.bed {tmp_out_dir}/" + file + f" > {tmp_out_dir}/adjusted_flanks.bed")
                os.system("awk 'FNR==NR{a[NR]=$3;next}{$3=a[FNR]}1' " + f"{tmp_out_dir}/coords_adjusted.bed {tmp_out_dir}/adjusted_flanks.bed > {tmp_out_dir}/adjusted.bed")
                os.system(f"sed 's/ /\t/g' {tmp_out_dir}/adjusted.bed > {tmp_out_dir}/{file}")

                # remove intermediate files
                os.system(f"rm {tmp_out_dir}/coords*.bed")
                os.system(f"rm  {tmp_out_dir}/adjusted*.bed")


        flanks = pd.read_csv(f'{tmp_out_dir}/majiq_flanks' + str(length) + '.bed', delimiter='\t', header=None)

        # drop flanks that have no junction
        ##chrY    13359417        13360117        Exon    .       -       chrY    13359767        13359768        flank   .       - 
        ##chr10   100041843       100042543       Exon    .       -       .       -1      -1      .       -1      . 

        flanks = flanks[flanks[8] != -1]

        # merge the flanks df with the jns df
        flanks[12] = flanks[[1, 2]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        flanks.drop([3, 4, 5, 6, 8, 9, 10, 11], axis=1, inplace=True)
        flanks.set_axis(['seqid', 'start', 'stop', 'junction0', 'flanks'], axis=1, inplace=True)

        junctions['index'] = junctions.index
        flank_jns = pd.merge(junctions, flanks, on=['junction0', 'seqid'])

        flank_jns.sort_values(['index'], inplace=True)
        flank_jns = flank_jns.reset_index(drop=True)
        del (flank_jns['index'])
        flank_jns['dpsi_'+str(length)] = pd.to_numeric(flank_jns['mean_dpsi_per_lsv_junction'])

        # # get all junctions that belong to each flank
        flank_jns_group = flank_jns.groupby(['flanks'])['dpsi_'+str(length)] \
                .apply(lambda val: ','.join(str(v) for v in val)).reset_index()
        annot.append(flank_jns_group)


for a in annot:
        a.set_index('flanks',inplace=True)

df = pd.concat(annot,axis=1,sort=False).reset_index()
df.dpsi_50.fillna(df.dpsi_100, inplace=True)
df.dpsi_50.fillna(df.dpsi_200, inplace=True)
flank_jns_group = df.drop(['dpsi_100', 'dpsi_200'], axis=1)
flank_jns_group.columns = ['flanks', 'mean_dpsi_per_lsv_junction']

# FILTER 3: if flank has 1+ junctions, keep junction with highest dPSI value
flank_jns_group['max_dPSI'] = flank_jns_group['mean_dpsi_per_lsv_junction'].str.split(',')\
        .apply(lambda x: max(map(float, x)))  # string -> list of strings -> list of floats -> max float

# # get the corresponding junction for each flank's max dPSI value
flank_jns_group = pd.merge(flank_jns_group[['flanks', 'max_dPSI']], flank_jns, on=['flanks'], how='inner')
flank_jns_group = flank_jns_group[flank_jns_group['mean_dpsi_per_lsv_junction'] == flank_jns_group['max_dPSI']]

# FILTER 4: If flank has 1+ junctions with same max dPSI value, keep one
flank_jns_group.drop_duplicates(subset='flanks', keep='first', inplace=True)

# # bookkeeping
del(flank_jns_group['max_dPSI'])
flank_jns_group = flank_jns_group[['gene_id', 'lsv_id', 'seqid', 'junction0', 'mean_dpsi_per_lsv_junction',
        'probability_changing', 'flanks', 'start', 'stop', 'strand']]

# Get all filtered flanks
flank_jns_group.drop_duplicates().to_csv(f'{tmp_out_dir}/all_flanks.csv', sep='\t', index=False)

# Get new annotation file with flanks (to annotate MANorm peaks)
flank_jns_group[['seqid', 'start', 'stop']].drop_duplicates().to_csv(f'{tmp_out_dir}/filtered_flanks.bed', index=False, sep='\t', header=False)

# Get dPSI values of filtered flanks
flank_jns_group[['flanks', 'gene_id', 'mean_dpsi_per_lsv_junction']].drop_duplicates().to_csv(f'{tmp_out_dir}/Filtered_dPSI.csv', index=False, sep='\t')
