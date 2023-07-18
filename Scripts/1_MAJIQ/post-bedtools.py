import pandas as pd
import os
import json


flank_lens = [50, 100, 200]
junctions = pd.read_csv('0_Files/majiq_junctions.csv', delimiter='\t')
with open('paths.json') as f:
    d = json.load(f)

fasta = d['Reference fasta']
ref_genome= fasta+".fai"

annot = []
for length in flank_lens:

        # adjust flank length to 200bp +- without exceeding chromosome bounds
        if length < 200:
                file = "majiq_flanks" + str(length) + ".bed"
                adjust_size = str(200 -length)

                # separate start,stop flank coords
                os.system("cut -f 1-3 0_Files/"+  file +" > 0_Files/coords.bed")

                # adjust flank boundaries
                os.system("bedtools slop -i 0_Files/coords.bed" + " -g " + ref_genome + " -b " + adjust_size + " > 0_Files/coords_adjusted.bed")

                # replace flank coords with adjusted coords
                os.system("awk 'FNR==NR{a[NR]=$2;next}{$2=a[FNR]}1' 0_Files/coords_adjusted.bed 0_Files/" + file + " > 0_Files/adjusted_flanks.bed")
                os.system("awk 'FNR==NR{a[NR]=$3;next}{$3=a[FNR]}1' 0_Files/coords_adjusted.bed 0_Files/adjusted_flanks.bed > 0_Files/adjusted.bed")
                os.system("sed 's/ /\t/g' 0_Files/adjusted.bed > 0_Files/"+file)

                # remove intermediate files
                os.system("rm 0_Files/coords*.bed")
                os.system("rm  0_Files/adjusted*.bed")


        flanks = pd.read_csv('0_Files/majiq_flanks' + str(length) + '.bed', delimiter='\t', header=None)

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
flank_jns_group.drop_duplicates().to_csv('0_Files/all_flanks.csv', sep='\t', index=False)

# Get new annotation file with flanks (to annotate MANorm peaks)
flank_jns_group[['seqid', 'start', 'stop']].drop_duplicates().to_csv('0_Files/filtered_flanks.bed', index=False, sep='\t', header=False)

# Get dPSI values of filtered flanks
flank_jns_group[['flanks', 'gene_id', 'mean_dpsi_per_lsv_junction']].drop_duplicates().to_csv('0_Files/Filtered_dPSI.csv', index=False, sep='\t')
