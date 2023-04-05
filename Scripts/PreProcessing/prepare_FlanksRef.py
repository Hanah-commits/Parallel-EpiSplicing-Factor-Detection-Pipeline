
import json
import os

with open('paths.json') as f:
    d = json.load(f)

ref = d['Reference genome']
fasta = d['Reference fasta']

# exon coordinates
os.system('Rscript PreProcessing/get_exons.R')

# genome file
os.system("samtools faidx " + fasta)
ref_genome= fasta+".fai"

# flanks
flanks = ["50", "100", "200"]

for flank in flanks:

    # exon boundary external flanks
    os.system("bedtools flank -i exon_coords.bed -g " + ref_genome + " -b "+ flank + " > 0_Files/flanks.bed" )

    # separate start,stop flank coords
    os.system("sed -n 'n;p' flanks.bed > 0_Files/stop.bed")
    os.system("sed -n 'p;n' flanks.bed > 0_Files/start.bed")

    # exon boundary internal flanks
    os.system("bedtools slop -i start.bed -g " + ref_genome + " -l 0 -r " + flank + " > 0_Files/start_flanks.bed")
    os.system("bedtools slop -i stop.bed -g " + ref_genome +" -l " + flank + " -r 0 > 0_Files/stop_flanks.bed")

    # combine start,stop flank coords
    os.system("paste -d'\n' start_flanks.bed stop_flanks.bed > 0_Files/flanks" + flank + ".bed")

    # remove intermediate files
    os.system("rm 0_Files/start*.bed")
    os.system("rm  0_Files/stop*.bed")
    os.system("rm 0_Files/flanks.bed")