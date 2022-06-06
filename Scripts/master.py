import os
import shutil
import sys
from HelperFunctions.check_args import check_args

# Check input arguments from paths.json
output_dir = check_args()

#check command line arguments
weights = False
if len(sys.argv) > 1 and sys.argv[1] == "-w":
    weights = True

# STEP 0: Preprocessing

# Differential Expression Analysis
if weights:
    try:
        exec(open("PreProcessing/featureCounts.py").read())
        os.system("Rscript PreProcessing/Limma.R")
    except Exception as ex:
        print(ex)
        sys.exit(1)

# Prepare flank reference : 50, 100, 200 bp
try:
    os.system("python PreProcessing/prepare_FlanksRef.py")
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 1: Execute MAJIQ - Differential Exon Usage
try:
    os.system("python 1_MAJIQ/runMAJIQ.py " + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 2: Execute MANorm -  Differential Histone Modifications
try:
    os.system("python 2_MANorm/manorm_all.py " + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 3: Process MAJIQ output
try:
    os.system("python 1_MAJIQ/post-MAJIQ.py " + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 4: BEDTools - Annotate exon flanks with MAJIQ junctions
try:
    exec(open("1_MAJIQ/annotate-MAJIQ.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 5: Process BEDTools output
try:
    exec(open("1_MAJIQ/post-bedtools.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 6: BEDTools - Annotate exon flanks with MANorm peaks
try:
    os.system('python 2_MANorm/annotate-MANorm.py ' + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 7: Process peak-annotated flanks
try:
    os.system("python 2_MANorm/post-manorm.py " + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 8: DEU - DHM Correlation
try:
    exec(open("3_Episplicing/correlation.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 9: Prepare RBPmap input
try:
    exec(open("4_RBPMap/pre-rbp.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 10: Execute RBPmap
try:
    exec(open("4_RBPMap/run_rbpmap.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 11: Process RBPMap output
try:
    exec(open("4_RBPMap/post-rbp.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

try:    
    exec(open("5_Classification & Enrichment/rbp_pvals.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 12: Add logFC weights to binding scores from RBPMap
if weights:
    try:
        exec(open("4_RBPMap/rbp-weights.py").read())
    except Exception as ex:
        print(ex)
        sys.exit(1)

# STEP 13: Prep Feature Matrix
try:
    exec(open("5_Classification & Enrichment/features.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

try: 
    exec(open("5_Classification & Enrichment/classifier_features.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)

try:
    exec(open("5_Classification & Enrichment/rbp_pvals.py").read())
except Exception as ex:
    print(ex)
    sys.exit(1)
    
# STEP 14: Binary Classification
try:
    os.system("python 5_Classification&Enrichment/classifier.py " + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 15: Enrichment
try:
    os.system("python 5_Classification&Enrichment/enrichment.py " + output_dir)
except Exception as ex:
    print(ex)
    sys.exit(1)

# STEP 15: Move files generated from current pipeline run to
shutil.move('0_Files/', output_dir)
shutil.move('../RBPmap/', output_dir)