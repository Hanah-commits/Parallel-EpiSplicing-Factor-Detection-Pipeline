import os
import json
from argparse import ArgumentParser

# Get the process name, use it in the output directory

p = ArgumentParser()
p.add_argument("--process", "-p",
    help="The name of the process")

args = p.parse_args()
proc = args.process

# STEP 1: Run RBPmap

with open('paths.json') as f:
            data = json.load(f)
d = data[proc]

rbp_path = d["RBPmap directory"]

currdir = os.getcwd()
file = os.getcwd() + f'/../{proc}_RBPmap/input.txt'
os.chdir(rbp_path)

sfs = ['BRUNOL4', 'BRUNOL5', 'BRUNOL6', 'DAZAP1', 'ESRP2', 'FMR1', 'FUS', 'FXR1', 'FXR2', 'HNRNPA1', 'HNRNPA1L2', 'HNRNPA2B1', 'HNRNPC', 'HNRNPF', 'HNRNPH1', 'HNRNPH2', 'HNRNPK', 'HNRNPL', 'HNRNPM', 'HNRNPU', 'HuR', 'KHDRBS1', 'KHDRBS2', 'KHDRBS3', 'MBNL1', 'PABPC1', 'PABPN1', 'PCBP1', 'PCBP2', 'PTB3', 'QKI', 'RALY', 'RBFOX1', 'RBM24', 'RBM28', 'RBM3', 'RBM4', 'RBM42', 'RBM5', 'RBM8A', 'SART3', 'SFPQ', 'SNRNP70', 'SNRPA', 'SRSF1', 'SRSF10', 'SRSF2', 'SRSF7', 'SRSF9', 'TARDBP', 'TIA1', 'U2AF2', 'YBX1', 'ZC3H10', 'ZCRB1', 'ZNF638']
          

os.system("parallel --joblog parallel_log -a " + file + f" perl RBPmap_EpiSplicing.pl {proc}" + " -input {1} -genome 'human' -db 'hg38' -db_motifs " + ",".join(sfs))

os.chdir(currdir)