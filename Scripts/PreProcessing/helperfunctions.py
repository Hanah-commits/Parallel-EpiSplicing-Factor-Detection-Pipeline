from ast import arg
import os
import json


def check_args():

    with open('paths.json') as f:
        d = json.load(f)

    args = {
    "bam_files" : d['RNASeq files'],
    "ref" : d['Reference genome'],
    "tissue1" : d["tissue1"],
    "tissue2" : d["tissue2"],
    "hms" : d["Histone modifications"],
    "chipseq_files" : d["ChIPSeq files"],
    "config" : d["MAJIQ config"],
    "rbpmap" : d["RBPmap directory"],
    "threads" : d['threads']

    }

    #check histone modifications & tissue names
    if len(args["hms"]) == 0:
        raise ValueError('No Histone Modifications')

    if len(args["tissue1"].strip()) == 0 or len(args["tissue2"].strip()) == 0:
        raise ValueError('No Tissue Name(s)')
        
    #check datatypes
    try:
        int(args["threads"])
    except:
        raise ValueError('Invalid Input ' + args["threads"])


    # check path validity of directories

    dirs = ["bam_files", "chipseq_files", "rbpmap"]
    dir_paths = []
    for dir in dirs:

        if args[dir][-1] != '/':
            args[dir] += '/'

        dir_paths.append(args[dir])
        

    for path in dir_paths:
        
        if not os.path.exists(os.path.dirname(path)):
            raise ValueError('Path does not exist ' + path)

    # check path validity of files
    file_paths = [args["ref"], args["config"]]
    for file in file_paths:
        
        if not os.path.isfile(file):
            raise ValueError('File does not exist ' + file)


    return args

