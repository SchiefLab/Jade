import os
import re
import gzip
import glob

#A Collection of pathing/file/dir functions


def get_database_path():
    p = os.path.split(os.path.abspath(__file__))[0]+"/../database"

    return p

def get_bin_path():
    return os.path.split(os.path.abspath(__file__))[0]+"/../bin"

def get_rosetta_features_root():
    return os.getenv('ROSETTA3_DB')+"/../tests/features"

def get_feat_input_path():
    p  = os.path.split(os.path.abspath(__file__))[0]+"/../rosetta_gen/features_input"
    return p

def get_pdb_path(decoy, alternate_paths = None):
    return get_decoy_path(decoy, alternate_paths)

def get_decoy_name(decoy):
    """
    Get the decoy name from path or name, whether .pdb, .pdb.gz or no extension.
    :param decoy:
    :rtype:str
    """

    name = os.path.basename(decoy)

    if re.search(".pdb.gz", name):
        return '.'.join(name.split(".")[0:-2])
    elif re.search(".pdb", name):
        return '.'.join(name.split('.')[0:-1])
    else:
        return name

def get_decoy_path(decoy, alternate_paths = None):
    """
    Search no extensions or with .pdb or .pdb.gz.
    Search alternative search paths.
    Return found path or NONE.

    :param decoy:
    :param alternate_paths:
    :rtype:str
    """

    #This is a hack due to wierd issues with the score file vs pdb file and an extra '_'

    decoy = decoy.replace('pre_model_1_', 'pre_model_1__')

    if alternate_paths:
        for dir in alternate_paths:
            decoy = dir +"/"+decoy
        if os.path.exists(decoy):
            return decoy
        elif os.path.exists(decoy+".pdb"):
            return decoy+".pdb"
        elif os.path.exists(decoy+".pdb.gz"):
            return decoy+".pdb.gz"
        else:
            return None
    else:
        if os.path.exists(decoy):
            return decoy
        elif os.path.exists(decoy+".pdb"):
            return decoy+".pdb"
        elif os.path.exists(decoy+".pdb.gz"):
            return decoy+".pdb.gz"
        else:
            return None

def get_xml_scripts_path():
    p  = os.path.split(os.path.abspath(__file__))[0]+"/../rosetta_gen/xml_scripts"
    return p


def make_dir_if_not_exists(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def get_make_get_dirs(root, dirs):
    """
    Recursively make dirs and return the final path
    :param root:
    :param dirs:
    :rtype: str
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
