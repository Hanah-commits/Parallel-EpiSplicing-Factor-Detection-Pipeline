import shutil
import os


wdir = os.path.dirname(os.path.realpath(__file__))

# STEP 0: Delete output directories from previous run
path = '../Input Files/RBPmap'
result_dirs = [x[0] for x in os.walk(path)]
result_dirs = [p for p in result_dirs if 'results' in p]
for dirpath in result_dirs:
    shutil.rmtree(dirpath)

# STEP 1: Run RBPmap
#TODO





