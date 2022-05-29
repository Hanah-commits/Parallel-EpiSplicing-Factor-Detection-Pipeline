import os
import json
import sys
from pathlib import Path

if __name__ == "__main__":

    with open('paths.json') as f:
            d = json.load(f)

    manorm_files = d['ChIPSeq files']
    hms = d["Histone modifications"]
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]

    currdir = os.getcwd()
    opdir = sys.argv[1] + 'MANorm'
    Path(opdir).mkdir(parents=True, exist_ok=True)

    os.chdir(manorm_files)
   
    for hm in hms:

        os.system('manorm --p1 ' + hm + '_'+tissue1+'_peak.bed --p2 ' + hm + '_'+tissue2+'_peak.bed --r1 ' + hm + '_'+tissue1+'_alignment.bed --r2 ' + hm + '_'+tissue2+'_alignment.bed -o ' + opdir)

    os.chdir(currdir)