import os
import re
import gzip
import glob

#A Collection of pathing/file/dir functions


def make_dir_if_not_exists(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_make_get_dirs(root, dirs):
    """
    Recursively make dirs and return the final path
    :param root:
    :param dirs:
    :return:
    """
    make_dir_if_not_exists(root)

    outpath = root
    for dir in dirs:
        outpath = outpath+"/"+dir
        make_dir_if_not_exists(outpath)

    return outpath

def get_file_paths(pattern, dir, ext = ".pdb"):
    files = glob.glob(dir+"/"+pattern+ext)
    file_paths = [dir+"/"+os.path.basename(x) for x in files]
    return file_paths

def open_file(file):
    if file.split(".")[-1] =="gz":
        #print "opening gzipped file"
        INFILE = gzip.open(file, 'rb')
    else:
        INFILE = open(file, 'r')

    return INFILE
