import os

flanks = [50,100,200]

for flank in flanks:

    os.system('bedtools intersect -loj -a ' + '0_Files/flanks' + str(flank) +'.bed -b 0_Files/majiq.bed | sort | uniq > 0_Files/majiq_flanks' + str(flank) + '.bed')