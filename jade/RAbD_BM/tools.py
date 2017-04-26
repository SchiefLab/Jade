import datetime
import glob
import re
import os
import sys

def get_pdb_paths(in_dir, exp_name, match_name = None, use_ensemble = False):
    exp_name = os.path.basename(in_dir)
    print exp_name
    if use_ensemble:
        files = glob.glob(in_dir+"/"+"*.pdb*")
    else:
        files = glob.glob(in_dir+"/"+exp_name+"*.pdb*")

    #print files
    #Match any names like relax
    if match_name:
        new_files = []
        for file in files:
            if re.search(match_name, file):
                new_files.append(file)
            else:
                continue
        files = new_files


    #Get full paths of the files
    new_files = []
    for file in files:
        if re.search("initial", file):
            continue
        new_files.append(in_dir+"/"+os.path.basename(file))
    files = new_files

    return files

def get_lambda_kappa_pdb_ids(dataset,pdb_type, root_dataset_dir = "datasets/pdblists" ):
    """
    Get two lists: lambda and kappa pdbids


    :param dataset: str
    :param root_dataset_dir: str
    :rtype: [str],[str]
    """
    def get_ids(path):
        ids = []
        if not os.path.exists(path):
            return ids

        HANDLE = open(path, "r")
        for line in HANDLE:
            line = line.strip()
            if not line or line.startswith("#"):continue
            ids.append(line)
        HANDLE.close()
        return ids

    lambda_path = os.path.join(root_dataset_dir, ".".join([dataset,pdb_type,"lambda", "PDBLIST.txt"  ]))
    kappa_path = os.path.join(root_dataset_dir, ".".join([dataset,pdb_type,"kappa", "PDBLIST.txt"  ]))

    if not os.path.exists(lambda_path) or not os.path.exists(kappa_path):
        sys.exit("Lambda and Kappa PDBLIST paths must exist for RAbD Benchmarks!!")

    return get_ids(lambda_path),get_ids(kappa_path)