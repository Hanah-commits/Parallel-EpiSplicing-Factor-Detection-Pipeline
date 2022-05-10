import shutil
import os
from glob import glob
import natsort
import pathlib

#TODO: MAJIQ, MANorm, Feature Counts, prepare flank files, delete existing result folders RBPMap (all in single script)


# STEP 1: Execute MAJIQ - Differential Exon Usage

# STEP 2: Execute MANorm -  Differential Histone Modifications

# STEP 3: Process MAJIQ output
exec(open("1_MAJIQ/post-MAJIQ.py").read())

# STEP 4: BEDTools - Annotate exon flanks with MAJIQ junctions
#TODO

# STEP 5: Process BEDTools output
exec(open("1_MAJIQ/post-bedtools.py").read())

# STEP 6: BEDTools - Annotate exon flanks with MANorm peaks
#TODO

# STEP 7: Process peak-annotated flanks
exec(open("2_MANorm/post-manorm.py").read())

# STEP 8: DEU - DHM Correlation
exec(open("3_Episplicing/correlation.py").read())

# STEP 9: Prepare RBPmap input
exec(open("4_RBPMap/pre-rbp.py").read())

# STEP 10: Execute RBPmap
# TODO

# STEP 11: Process RBPMap output
# exec(open("4_RBPMap/post-rbp.py").read())
# exec(open("5_Classification & Enrchment/rbp_pvals.py").read())

# STEP 12: Add logFC weights to binding scores from RBPMap
exec(open("4_RBPMap/rbp-weights.py").read())

# STEP 13: Prep Feature Matrix
exec(open("5_Classification & Enrichment/features.py").read())
exec(open("5_Classification & Enrichment/classifier_features.py").read())
exec(open("5_Classification & Enrichment/rbp_pvals.py").read())

# STEP 14: Binary Classification
exec(open("5_Classification & Enrichment/classifier.py").read())

# STEP 15: Enrichment
exec(open("5_Classification & Enrichment/enrichment.py").read())

