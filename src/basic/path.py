import os
import re
import gzip
import glob



extensions = [".pdb", ".cif", ".xml"]
compressions = [".gz", ".tar.gz", ""]


#A Collection of pathing/file/dir functions

def get_Jade_root():
    """
    Get the root path of Jade directory.
    :rtype: str
    """
    return os.path.split(os.path.abspath(__file__))[0]+"/../.."

def get_database_path():
    """
    Get the path to the Jade Database
    :rtype: str
    """
    return get_Jade_root()+"/database"

def get_bin_path():
    """
    Get the path to the Jade apps directory
    :rtype: str
    """
    return get_Jade_root()+"/apps"

def get_testing_path():
    """
    Get the path to the Jade testing directory
    :rtype: str
    """
    return get_Jade_root()+"/testing"

def get_testing_inputs_path():
    """
    Get the path to testing inputs (PDBs,fasta,etc.)
    :rtype:str
    """
    return get_testing_path()+"/inputs"


def get_rosetta_features_root():
    """
    Get the path to Rosetta features directory through set ROSETTA3_DB env variable.
    :rtype: str
    """
    return os.getenv('ROSETTA3_DB')+"/../tests/features"

def get_rosetta_features_run_script():
    """
    Get the path to Rosetta features script dir through the set ROSETTA3_DB env variable.
    :rtype: str
    """
    return os.getenv('ROSETTA3_DB')+"/../source/scripts/features/run_features.R"

def get_xml_scripts_path():
    """
    Get the path to the Rosetta xml script directory.  Useful for variable substitutions.
    :rtype: str
    """
    return get_bin_path()+"/xml_scripts"

def get_rosetta_flags_path():
    return get_database_path()+"/rosetta/flags"

def get_rosetta_json_run_path():
    return get_database_path()+"/rosetta/jsons"

def get_rosetta_features_json_path():
    return get_database_path()+"/rosetta/features"

###########################################################################


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
    Search no extensions or with .pdb, .pdb.gz, .cif, .cif.gz, .xml, .xml.gz
    Search alternative search paths.
    Return found path or NONE.

    :param decoy:
    :param alternate_paths:
    :rtype:str
    """

    #This is a hack due to wierd issues with the score file vs pdb file and an extra '_'
    #decoy = decoy.replace('pre_model_1_', 'pre_model_1__')

    def find_decoy(decoy):
        for ext in extensions:
            for comp in compressions:
                search_path = decoy+ext+comp
                if os.path.exists(search_path):
                    return search_path

        return None

    if alternate_paths:
        for dir in alternate_paths:
            new_decoy_path = dir +"/"+decoy
            return find_decoy(new_decoy_path)

    else:
        return find_decoy(decoy)

def get_decoy_extension(decoy):
    """
    Return the extension of the decoy.  .pdb, .pdb.gz, .cif, .cif.gz, etc.
    :param decoy: str
    :rtype: str
    """

    for ext in extensions:
        if re.search(ext, decoy):
            if re.search(".gz", decoy):
                return ext+".gz"
            else:
                return ext

    return ""


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

def get_directories_recursively(inpath):
    """
    Get a list of directories recursively in a path.  Skips hidden directories.
    :param inpath: str
    :rtype: list
    """

    all_dirs = []
    for root, dirs, files in os.walk(inpath):
        all_dirs.extend([root+"/"+d for d in dirs if d[0] != '.'])
    return all_dirs


