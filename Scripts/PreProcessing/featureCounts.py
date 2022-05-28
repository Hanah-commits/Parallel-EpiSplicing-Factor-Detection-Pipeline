import os
import json

if __name__ == "__main__":

    with open('paths.json') as f:
        d = json.load(f)

    bam_files = d['RNASeq files']
    ref = d['Reference genome']
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]
    threads = d['threads']

    currdir = os.getcwd()
    opdir = currdir + '/0_Files/' + tissue1+'_'+tissue2+'_counts'

    # -T -> specify cores
    # -a -> required option for specifying path to GTF
    # -o -> required option for specifying path to, and name of the text output (count matrix)


    os.system('featureCounts -p -T '+ threads +' -a ' + ref + ' -o ' + opdir + ' ' + bam_files+ '/*.bam')
