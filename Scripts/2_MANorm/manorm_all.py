import os


data_dir = '/home/hanah/data/bam_dir/histone_data'
op_dir = '/home/hanah/data/bam_dir/histone_data/manorm'


hms = ['H3K27ac', 'H3K27me3', 'H3K4me3', 'H3K9me3']

os.chdir(data_dir)

for hm in hms:

    os.system('manorm --p1 ' + hm + '_H1_peak.bed --p2 ' + hm + '_ectodermalcell_peak.bed --r1 ' + hm + '_H1_alignment.bed --r2 ' + hm + '_ectodermalcell_alignment.bed -o ' + op_dir)