import os
import sys
import json

with open('paths.json') as f:
    d = json.load(f)

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
hms = d["Histone modifications"]

prefix = sys.argv[1]+ 'MANorm/'

for hm in hms:
    
    input = prefix + hm + '_' + tissue1 + '_peak_vs_' + hm + '_' + tissue2 +  '_peak_all_MAvalues.xls'
    output = prefix + hm + '_flanks.bed'

    os.system('bedtools intersect -loj -a 0_Files/filtered_flanks.bed -b ' + input + ' | sort | uniq > ' + output)