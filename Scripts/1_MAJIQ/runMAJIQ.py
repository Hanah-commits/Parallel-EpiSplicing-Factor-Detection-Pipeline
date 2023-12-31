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

      majiq_files = d['RNASeq files']
      ref = d['Reference genome']
      config = d["MAJIQ config"]
      threads = d["threads"]
      tissue1 = d["tissue1"]
      tissue2 = d["tissue2"]


      currdir = os.getcwd()
      
      output = args.output_dir + 'MAJIQ/'
      build_output = output+'build/'
      tissue1_output = output+'psi_tissue1/'
      tissue2_output = output+'psi_tissue2/'
      deltapsi_output = output+ 'deltapsi/'

      dirs = [output, build_output, tissue1_output, tissue2_output, deltapsi_output]
      for dir in dirs:
            Path(dir).mkdir(parents=True, exist_ok=True)

      # switch to the directory with the bam files
      os.chdir(majiq_files)

      # Detect Splice Variations
      os.system('majiq build ' + ref + ' -j ' + threads + ' -o '+build_output + ' -c ' + config + ' --disable-ir')

      # collect .majiq files
      all_files = []
      for file in os.listdir(build_output):
            if file.endswith(".majiq"):
                  all_files.append(file)      

      tissue1_files = [os.path.join(build_output, file) for file in all_files if tissue1 in file]
      tissue2_files = [os.path.join(build_output, file) for file in all_files if tissue2 in file]


      tissue1_count = len(tissue1_files)
      tissue2_count = len(tissue2_files)

      tissue1_files = " ".join(tissue1_files)
      tissue2_files = " ".join(tissue2_files)

      # Quantify Splice Variation
      os.system('majiq psi ' + tissue1_files + ' -j '+ str(tissue1_count) + ' -o '+ tissue1_output + ' -n ' + tissue1)
      os.system('majiq psi '+ tissue2_files + ' -j '+ str(tissue2_count) + ' -o '+ tissue2_output + ' -n ' + tissue2)

      # Quantify Differential Splice Variation
      os.system('majiq deltapsi -grp1 ' + tissue1_files + ' -grp2 ' + tissue2_files + ' -j ' + threads + ' -o ' + deltapsi_output + ' -n ' + tissue1 + ' ' + tissue2)

      os.chdir(currdir)

      # Obtain tsv file
      os.system('voila tsv ' + build_output + 'splicegraph.sql ' + deltapsi_output + tissue1+'-'+tissue2 +'.deltapsi.voila' +' -f ' + output + 'majiq_output')


if __name__ == "__main__":    
      p = get_argument_parser()
      args = p.parse_args()
      main(args)
