import os
import json
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("output_dir")
p.add_argument("--process", "-p",
    help="The name of the process")

args = p.parse_args()
proc = args.process

with open('paths.json') as f:
    data = json.load(f)
d = data[proc]

tissue1 = d["tissue1"]
tissue2 = d["tissue2"]
hms = d["Histone modifications"]

prefix = args.output_dir + 'MANorm/'

for hm in hms:
    
    input = prefix + hm + '_' + tissue1 + '_peak_vs_' + hm + '_' + tissue2 +  '_peak_all_MAvalues.xls'
    output = prefix + hm + '_flanks.bed'

    os.system(f'bedtools intersect -loj -a {proc}_0_Files/filtered_flanks.bed -b ' + input + ' | sort | uniq > ' + output)