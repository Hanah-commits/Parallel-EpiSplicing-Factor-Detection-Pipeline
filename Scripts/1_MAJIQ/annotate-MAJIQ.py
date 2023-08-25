import os
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("--process", "-p",
    help="The name of the process")
args = p.parse_args()
proc = args.process

tmp_out_dir = proc + '_0_Files'

flanks = [50,100,200]

for flank in flanks:

    os.system('bedtools intersect -loj -s -a ' + f'{tmp_out_dir}/flanks' + str(flank) +f'.bed -b {tmp_out_dir}/majiq.bed | sort | uniq > {tmp_out_dir}/majiq_flanks' + str(flank) + '.bed')
