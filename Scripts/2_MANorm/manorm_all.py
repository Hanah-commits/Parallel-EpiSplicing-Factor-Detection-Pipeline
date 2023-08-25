import os
import json
from pathlib import Path
from argparse import ArgumentParser

# Get the process name, use it in the output directory
def get_argument_parser():
    p = ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("--process", "-p",
        help="The name of the process")
    return p


def main(args):
    proc = args.process

    with open('paths.json') as f:
            data = json.load(f)
    d = data[proc]

    manorm_files = d['ChIPSeq files']
    hms = d["Histone modifications"]
    tissue1 = d["tissue1"]
    tissue2 = d["tissue2"]

    currdir = os.getcwd()
    opdir = args.output_dir + 'MANorm/'
    Path(opdir).mkdir(parents=True, exist_ok=True)

    os.chdir(manorm_files)
   
    for hm in hms:

        os.system('manorm --p1 ' + hm + '_'+tissue1+'_peak.bed --p2 ' + hm + '_'+tissue2+'_peak.bed --r1 ' + hm + '_'+tissue1+'_alignment.bed --r2 ' + hm + '_'+tissue2+'_alignment.bed -o ' + opdir)

    os.chdir(currdir)


if __name__ == "__main__":    
      p = get_argument_parser()
      args = p.parse_args()
      main(args)
