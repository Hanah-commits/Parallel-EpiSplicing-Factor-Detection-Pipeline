import time
import os
import json
from pathlib import Path
import shutil


def check_args():

    with open('paths.json') as f:
        data = json.load(f)

    procs = data["list_of_processes"]
    output_dirs = []

    new_paths_file = {"list_of_processes":data["list_of_processes"]}

    for proc in procs:

        d = data[proc]

        args = {
        "RNASeq files" : d['RNASeq files'],
        "Reference genome" : d['Reference genome'],
        "Reference fasta": d['Reference fasta'],
        "tissue1" : d["tissue1"],
        "tissue2" : d["tissue2"],
        "Histone modifications" : d["Histone modifications"],
        "ChIPSeq files" : d["ChIPSeq files"],
        "MAJIQ config" : d["MAJIQ config"],
        "RBPmap directory" : d["RBPmap directory"],
        "threads" : d['threads'],
        "strandedness": int(d['strandedness'])
        }

        #check histone modifications & tissue names
        args["Histone modifications"] = [x for x in args["Histone modifications"] if x]
        if len(args["Histone modifications"]) == 0:
            raise ValueError('No Histone Modifications')

        if len(args["tissue1"].strip()) == 0 or len(args["tissue2"].strip()) == 0:
            raise ValueError('No Tissue Name(s)')
        
        #check datatypes
        try:
            int(args["threads"])
        except:
            raise ValueError('Invalid Input ' + args["threads"])


        # check path validity of directories
        dirs = ["RNASeq files", "ChIPSeq files", "RBPmap directory"]
        dir_paths = []
        for dir in dirs:

            if args[dir][-1] != '/':
                args[dir] += '/'

            dir_paths.append(args[dir])
        

        for path in dir_paths:
        
            if not os.path.exists(os.path.dirname(path)):
                raise ValueError('Path does not exist ' + path)

        # check path validity of files
        file_paths = [args["Reference genome"], args["Reference fasta"], args["MAJIQ config"]]
        for file in file_paths:
        
            if not os.path.isfile(file):
                raise ValueError('File does not exist ' + file)

        
        # check if temp directories already exist
        temp_dirs = [f'{proc}_0_Files/', f'../{proc}_RBPmap/']
        for dir in temp_dirs:
            if os.path.exists(dir):
                # temp dir not empty
                if len(os.listdir(dir)) != 0:
                    raise ValueError('Delete or move directory to another location ' + dir)
            else: 
                #create temp dir
                Path(dir).mkdir(parents=True, exist_ok=True)

        # check if log files already exist
        log_file_name = f'{proc}_output.log'
        try:
            f = open(log_file_name, 'x')
        except FileExistsError:
            raise ValueError('Delete or move the log file to another location ' + log_file_name)

        # create custome output directory process_tissue1_tissue2_timestamp
        output_dir = str(Path(os.getcwd()).parent.absolute()) + "/Output/"+ proc + "_" +  args["tissue1"]+ "_" + args["tissue2"]+ "_" + str(time.time()) +"/"
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        output_dirs.append(output_dir)

        new_paths_file[proc] = args

    
    with open('paths.json', 'w') as fp:
            json.dump(new_paths_file, fp)

    for out_dir in output_dirs:
        # copy input arguments (paths.json) to output_dir
        shutil.copyfile('paths.json', out_dir+'paths.json')

    return procs, output_dirs


def move_dirs(output_dir, proc):
    shutil.move(f'{proc}_0_Files/', output_dir)
    shutil.move(f'../{proc}_RBPmap/', output_dir)
    shutil.move(f'{proc}_output.log', output_dir)
