import os
import json

with open('paths.json') as f:
    d = json.load(f)

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
hms = d["Histone modifications"]

for hm in hms:
    
    file = hm + '_' + tissue2 + '_peak_vs_' + hm + '_' + tissue2 +  '_peak_all_MAvalues.xls'

    os.system('bedtools intersect -loj -a 0_Files/filtered_flanks.bed -b ' + file + ' | sort | uniq > ' + hm + '_flanks.bed')