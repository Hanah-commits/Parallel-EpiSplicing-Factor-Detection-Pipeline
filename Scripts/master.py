import os
import shutil
import sys
from HelperFunctions.check_args import check_args

# Check input arguments from paths.json
output_dir = check_args()

#check command line arguments
weights = False
if len(sys.argv) > 0 and sys.argv[1] == "-w":
    weights = True

# STEP 0: Preprocessing

# Differential Expression Analysis
if weights:
    exec(open("PreProcessing/featureCounts.py").read())
    os.system("Rscript PreProcessing/Limma.R")

# Prepare flank reference : 50, 100, 200 bp
os.system("python PreProcessing/prepare_FlanksRef.py")

# STEP 1: Execute MAJIQ - Differential Exon Usage
os.system("python 1_MAJIQ/runMAJIQ.py " + output_dir)

# STEP 2: Execute MANorm -  Differential Histone Modifications
os.system("python 2_MANorm/manorm_all.py " + output_dir)

# STEP 3: Process MAJIQ output
exec(open("1_MAJIQ/post-MAJIQ.py").read())

# STEP 4: BEDTools - Annotate exon flanks with MAJIQ junctions
exec(open("1_MAJIQ/annotate-MAJIQ.py").read())

# STEP 5: Process BEDTools output
exec(open("1_MAJIQ/post-bedtools.py").read())

# STEP 6: BEDTools - Annotate exon flanks with MANorm peaks
exec(open("2_MANorm/annotate-MANorm.py").read())

# STEP 7: Process peak-annotated flanks
exec(open("2_MANorm/post-manorm.py").read())

# STEP 8: DEU - DHM Correlation
exec(open("3_Episplicing/correlation.py").read())

# STEP 9: Prepare RBPmap input
exec(open("4_RBPMap/pre-rbp.py").read())

# STEP 10: Execute RBPmap
exec(open("4_RBPMap/run_rbpmap.py").read())

# STEP 11: Process RBPMap output
exec(open("4_RBPMap/post-rbp.py").read())
exec(open("5_Classification & Enrichment/rbp_pvals.py").read())

# STEP 12: Add logFC weights to binding scores from RBPMap
if weights:
    exec(open("4_RBPMap/rbp-weights.py").read())

# STEP 13: Prep Feature Matrix
exec(open("5_Classification & Enrichment/features.py").read())
exec(open("5_Classification & Enrichment/classifier_features.py").read())
exec(open("5_Classification & Enrichment/rbp_pvals.py").read())

# STEP 14: Binary Classification
os.system("python 5_Classification&Enrichment/classifier.py " + output_dir)

# STEP 15: Enrichment
os.system("python 5_Classification&Enrichment/enrichment.py " + output_dir)

# STEP 15: Move files generated from current pipeline run to
shutil.move('0_Files/', output_dir)
shutil.move('../RBPmap/', output_dir)