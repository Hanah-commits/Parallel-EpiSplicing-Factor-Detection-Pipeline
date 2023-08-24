
import os
import json
from argparse import ArgumentParser

# Get the process name, use it in the output directory
def get_argument_parser():
    p = ArgumentParser()
    p.add_argument("--process", "-p",
        help="The name of the process")
    return p

def main(args):
    proc = args.process

    with open('paths.json') as f:
        data = json.load(f)

    d = data[proc]

    bam_files = d['RNASeq files']
    ref = d['Reference genome']
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]
    threads = d['threads']
    strand = d['strandedness']

    currdir = os.getcwd()
    opdir = currdir + f'/{proc}_0_Files/' + tissue1+'_'+tissue2+'_counts'

    # -T -> specify cores
    # -a -> required option for specifying path to GTF
    # -o -> required option for specifying path to, and name of the text output (count matrix)


    os.system('featureCounts -p -s ' + strand + ' -T '+ threads +' -a ' + ref + ' -o ' + opdir + ' ' + bam_files+ '/*.bam')

if __name__ == "__main__":
    p = get_argument_parser()
    args = p.parse_args()
    main(args)
