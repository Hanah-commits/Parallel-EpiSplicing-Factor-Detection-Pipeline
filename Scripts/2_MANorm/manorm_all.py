import os
import json

if __name__ == "__main__":

    with open('paths.json') as f:
            d = json.load(f)

    manorm_files = d['ChIPSeq files']
    hms = d["Histone modifications"]
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]

    currdir = os.getcwd()
    opdir = currdir + '/../Input_Files/MANorm'
   
    for hm in hms:

        os.system('manorm --p1 ' + hm + '_'+tissue1+'_peak.bed --p2 ' + hm + '_'+tissue2+'_peak.bed --r1 ' + hm + '_'+tissue1+'_alignment.bed --r2 ' + hm + '_'+tissue2+'_alignment.bed -o ' + opdir)