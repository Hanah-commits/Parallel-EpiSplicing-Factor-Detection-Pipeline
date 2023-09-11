import os
import sys
from HelperFunctions.check_args import check_args, move_dirs

from multiprocessing import Process

# Check input arguments from paths.json
list_of_processes, output_dirs = check_args()

# The function that will be executed by all processes (semantics is not changed)
def master_function(proc, output_dir):

    #check command line arguments
    weights = False
    if len(sys.argv) > 1 and sys.argv[1] == "-w":
        weights = True

    # STEP 0: Preprocessing

    # Differential Expression Analysis
    if weights:
        try:
            os.system(f"PreProcessing/featureCounts.py -p {proc}")
            os.system(f"Rscript PreProcessing/Limma.R {proc}")
        except Exception as ex:
            print(ex)
            move_dirs(output_dir, proc)
            sys.exit(1)

    # Prepare flank reference : 50, 100, 200 bp
    try:
        os.system(f"python PreProcessing/prepare_FlanksRef.py -p {proc}")
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    # # STEP 1: Execute MAJIQ - Differential Exon Usage
    try:
        os.system(f"python 1_MAJIQ/runMAJIQ.py {output_dir} -p {proc}")
    except Exception as ex:
        print(ex)
        move_dirs(output_dir, proc)
        sys.exit(1)

    # # STEP 2: Execute MANorm -  Differential Histone Modifications
    # try:
    #     os.system(f"python 2_MANorm/manorm_all.py {output_dir} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 3: Process MAJIQ output
    # try:
    #     os.system(f"python 1_MAJIQ/post-MAJIQ.py {output_dir} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 4: BEDTools - Annotate exon flanks with MAJIQ junctions
    # try:
    #     os.system(f"1_MAJIQ/annotate-MAJIQ.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 5: Process BEDTools output
    # try:
    #     os.system(f"1_MAJIQ/post-bedtools.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 6: BEDTools - Annotate exon flanks with MANorm peaks
    # try:
    #     os.system(f'python 2_MANorm/annotate-MANorm.py {output_dir} -p {proc}')
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 7: Process peak-annotated flanks
    # try:
    #     os.system(f"python 2_MANorm/post-manorm.py {output_dir} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 8: DEU - DHM Correlation
    # try:
    #     os.system(f"3_Episplicing/correlation.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 9: Prepare RBPmap input
    # try:
    #     os.system(f"4_RBPMap/pre-rbp.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 10: Execute RBPmap
    # try:
    #     os.system(f"4_RBPMap/run_rbpmap.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 11: Process RBPMap output
    # try:
    #     os.system(f"4_RBPMap/post-rbp.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # try:    
    #     os.system(f"5_Classification/rbp_pvals.py -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 12: Add logFC weights to binding scores from RBPMap
    # if weights:
    #     try:
    #         os.system(f"PreProcessing/get_weights.py -p {proc}")
    #         os.system(f"4_RBPMap/rbp-weights.py -p {proc}")
    #     except Exception as ex:
    #         print(ex)
    #         move_dirs(output_dir, proc)
    #         sys.exit(1)

    # # STEP 13: Prep Feature Matrix
    # try:
    #     os.system(f"python 5_Classification/features.py {output_dir} {str(weights)} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # try: 
    #     os.system(f"python 5_Classification/classifier_features.py {output_dir} {str(weights)} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 14: Binary Classification
    # try:
    #     os.system(f"python 5_Classification/classifier.py {output_dir} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # # STEP 15: Enrichment
    # try:
    #     os.system(f"python 6_Enrichment/enrichment.py {output_dir} -p {proc}")
    # except Exception as ex:
    #     print(ex)
    #     move_dirs(output_dir, proc)
    #     sys.exit(1)

    # STEP 15: Move files generated from current pipeline run to
    move_dirs(output_dir, proc)

active_processes = []

# Start the processes
for pr, out_dir in zip(list_of_processes, output_dirs):
    proc = Process(target=master_function, args=(pr, out_dir,))
    print('Proc is created and launched!')
    active_processes.append(proc)
    proc.start()

# Wait for the processes to complete and free resources
for proc in active_processes:
        proc.join()
