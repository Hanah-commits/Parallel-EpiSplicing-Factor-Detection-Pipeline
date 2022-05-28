from ast import arg
import os
import json


def check_args():

    with open('paths.json') as f:
        d = json.load(f)

    args = {
    "RNASeq files" : d['RNASeq files'],
    "Reference genome" : d['Reference genome'],
    "tissue1" : d["tissue1"],
    "tissue2" : d["tissue2"],
    "Histone modifications" : d["Histone modifications"],
    "ChIPSeq files" : d["ChIPSeq files"],
    "MAJIQ config" : d["MAJIQ config"],
    "RBPmap directory" : d["RBPmap directory"],
    "threads" : d['threads']
    }

    #check histone modifications & tissue names
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

    dirs = ["RNASeq files", "ChIPSeq_files", "RBPmap directory"]
    dir_paths = []
    for dir in dirs:

        if args[dir][-1] != '/':
            args[dir] += '/'

        dir_paths.append(args[dir])
        

    for path in dir_paths:
        
        if not os.path.exists(os.path.dirname(path)):
            raise ValueError('Path does not exist ' + path)

    # check path validity of files
    file_paths = [args["Reference genome"], args["MAJIQ config"]]
    for file in file_paths:
        
        if not os.path.isfile(file):
            raise ValueError('File does not exist ' + file)


    with open('paths.json', 'w') as fp:
        json.dump(args, fp)

