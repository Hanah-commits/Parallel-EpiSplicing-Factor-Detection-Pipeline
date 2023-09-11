import json
import os

from argparse import ArgumentParser

# Get the process name, use it in the output directory
p = ArgumentParser()
p.add_argument("--process", "-p",
    help="The name of the process")
args = p.parse_args()
proc = args.process

tmp_out_dir = proc + '_0_Files'


with open('paths.json') as f:
    data = json.load(f)
d = data[proc]

ref = d['Reference genome']
fasta = d['Reference fasta']

# exon coordinates
os.system(f'bash PreProcessing/get_exons.sh {proc}')

# genome file
ref_genome= fasta+".fai"

# flanks
flanks = ["50", "100", "200"]

for flank in flanks:

    # exon boundary external flanks
    os.system(f"bedtools flank -i {tmp_out_dir}/exon_coords.bed -g " + ref_genome + " -b "+ flank + f" > {tmp_out_dir}/flanks.bed" )

    # separate start,stop flank coords
    os.system(f"sed -n 'n;p' {tmp_out_dir}/flanks.bed > {tmp_out_dir}/stop.bed")
    os.system(f"sed -n 'p;n' {tmp_out_dir}/flanks.bed > {tmp_out_dir}/start.bed")

    # exon boundary internal flanks
    os.system(f"bedtools slop -i {tmp_out_dir}/start.bed -g " + ref_genome + " -l 0 -r " + flank + f" > {tmp_out_dir}/start_flanks.bed")
    os.system(f"bedtools slop -i {tmp_out_dir}/stop.bed -g " + ref_genome +" -l " + flank + f" -r 0 > {tmp_out_dir}/stop_flanks.bed")

    # combine start,stop flank coords
    os.system(f"paste -d'\n' {tmp_out_dir}/start_flanks.bed {tmp_out_dir}/stop_flanks.bed > {tmp_out_dir}/flanks" + flank + ".bed")

    # remove intermediate files
    os.system(f"rm {tmp_out_dir}/start*.bed")
    os.system(f"rm  {tmp_out_dir}/stop*.bed")
    os.system(f"rm {tmp_out_dir}/flanks.bed")
