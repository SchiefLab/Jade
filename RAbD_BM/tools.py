import datetime
import glob
import re
import os

def get_today():
    return datetime.date.today().strftime("%Y/%m/%d")


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

