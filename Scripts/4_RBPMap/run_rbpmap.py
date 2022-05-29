import shutil
import os
import glob
import json


# wdir = os.path.dirname(os.path.realpath(__file__))

# # STEP 0: Delete output directories from previous run
# path = '../Input Files/RBPmap'
# result_dirs = [x[0] for x in os.walk(path)]
# result_dirs = [p for p in result_dirs if 'results' in p]
# for dirpath in result_dirs:
#     shutil.rmtree(dirpath)

# STEP 1: Run RBPmap

with open('paths.json') as f:
            d = json.load(f)

rbp_path = d["RBPmap directory"]

currdir = os.getcwd()
file = os.getcwd() + '/0_Files/input.txt'
os.chdir(rbp_path)

sfs = ['HNRNPL', 'HNRNPH1', 'HNRNPK', 'SFPQ', 'HNRNPA1', 'HNRNPA2B1', 'PTBP1', 'HNRNPF', 'HNRNPH2', 'HNRNPM', 'FUS',
          'YBX1', 'PCBP1', 'HNRNPC', 'HuR', 'TARDBP', 'HNRNPU', 'PCBP2', 'SRSF9', 'SRSF1', 'SRSF7', 'U2AF2', 'SRSF10',
          'SRSF2', 'RALY', 'MBNL1', 'RBM4', 'PABPN1', 'RBM3', 'TIA1', 'KHDRBS1', 'RBM28', 'PABPC1', 'RBM5', 'SART3',
          'SNRNP70', 'FXR2', 'FXR1', 'ESRP2', 'HNRNPA1L2', 'ZNF638', 'SNRPA', 'RBM8A', 'FMR1', 'DAZAP1', 'RBM42',
          'ZCRB1', 'KHDRBS3', 'QKI', 'KHDRBS2', 'ZC3H10', 'RBM24', 'RBFOX1', 'BRUNOL5', 'BRUNOL4', 'BRUNOL6']

os.system("parallel --joblog parallel_log -a " + file + " perl RBPmap_EpiSplicing.pl -input {1} -genome 'human' -db 'hg38' -db_motifs " + ",".join(sfs))

os.chdir(currdir)






