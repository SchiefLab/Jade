import os
import re
import gzip
import glob

#A Collection of pathing/file/dir functions


def get_database_path():
    p = os.path.split(os.path.abspath(__file__))[0]+"/../database"

    return p

def get_rosetta_features_root():
    return os.getenv('ROSETTA3_DB')+"/../tests/features"

def get_feat_input_path():
    p  = os.path.split(os.path.abspath(__file__))[0]+"/../rosetta_general/features_input"
    return p

def get_xml_scripts_path():
    p  = os.path.split(os.path.abspath(__file__))[0]+"/../rosetta_general/xml_scripts"
    return p


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

def open_file(file, mode = 'r'):
    if file.split(".")[-1] =="gz":
        #print "opening gzipped file"
        INFILE = gzip.open(file, mode+'b')
    else:
        INFILE = open(file, mode)

    return INFILE
