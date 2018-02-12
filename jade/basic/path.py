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
    return get_Jade_root()+"/jade/tests"

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
    return os.getenv('ROSETTA3')+"/../tests/features"

def get_rosetta_features_run_script():
    """
    Get the path to Rosetta features script dir through the set ROSETTA3_DB env variable.
    :rtype: str
    """
    return get_database_path()+"/rosetta/features/run_features.R"

def get_xml_scripts_path():
    """
    Get the path to the Rosetta xml script directory.  Useful for variable substitutions.
    :rtype: str
    """
    return get_database_path()+"/rosetta/xml_scripts"

def get_rosetta_flags_path():
    return get_database_path()+"/rosetta/flags"

def get_rosetta_json_run_path():
    return get_database_path()+"/rosetta/jsons"

def get_rosetta_features_json_path():
    return get_database_path()+"/rosetta/features"

def get_nnk_database_path():
    return get_database_path()+"/nnk"

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

    if re.search(".pdb.gz", name) or re.search(".cif.gz", name):
        return '.'.join(name.split(".")[0:-2])
    elif re.search(".pdb", name) or re.search(".cif", name):
        return '.'.join(name.split('.')[0:-1])
    else:
        return name

def get_decoy_path(decoy, alternate_paths = None):
    """
    Search .pdb, .pdb.gz, .cif, .cif.gz, .xml, .xml.gz
    In addition, Search alternative search paths.
    Return found path or NONE.

    :param decoy:
    :param alternate_paths:
    :rtype:str
    """

    #This is a hack due to wierd issues with the score file vs pdb file and an extra '_'
    #decoy = decoy.replace('pre_model_1_', 'pre_model_1__')

    def find_decoy(f):
        found_ext = get_decoy_extension(f)
        if found_ext:
            #print "Found extension "+found_ext
            #print "Replacing "+f
            f = f.replace(found_ext, "")

        for ext in extensions:
            for comp in compressions:
                search_path = f+ext+comp
                #print search_path
                if os.path.exists(search_path):
                    return search_path


    if alternate_paths:
        for dir in alternate_paths:
            new_decoy_path = dir +"/"+decoy
            return find_decoy(new_decoy_path)

    else:
        return find_decoy(decoy)

    return None

def get_decoy_extension(decoy):
    """
    Return the extension of the decoy.  .pdb, .pdb.gz, .cif, .cif.gz, etc.
    :param decoy: str
    :rtype: str
    """
    #print "finding extension for "+ decoy
    for ext in extensions:
        for comp in compressions:
            extension = ext+comp
            #print extension
            if re.search(extension, decoy):
                return extension

    return None


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
    """
    Get file paths matching the exact pattern and extension.
    :param pattern: 
    :param dir: 
    :param ext: 
    :return: 
    """
    files = glob.glob(dir+"/"+pattern+ext)
    file_paths = [dir+"/"+os.path.basename(x) for x in files]
    return file_paths

def get_matching_pdbs(directory, pattern, ext='.pdb'):
    """
    Get pdbs in a directory matching a pattern.
    :param directory: 
    :param pattern: 
    :param ext: 
    :return: 
    """
    files = glob.glob(directory+"/"+'*'+pattern+'*'+ext)
    return [os.path.basename(f) for f in files]

def get_all_pdbs(directory, ext='.pdb'):
    return get_matching_pdbs(directory, '', ext)

def get_all_pdb_paths(directory, ext='.pdb'):
    pdbs = get_all_pdbs(directory, ext)
    return [directory+'/'+pdb for pdb in pdbs]

def open_file(file_path, mode = 'r'):
    if file_path.split(".")[-1] =="gz":
        #print "opening gzipped file"
        INFILE = gzip.open(file_path, mode+'b')
    else:
        INFILE = open(file_path, mode)

    return INFILE

def parse_contents(file_path):
    """
    Return a list of (stripped) content, skipping empty lines and comments.
    :param file: string
    :rtype: list 
    """
    return [x.strip() for x in open_file(file_path, 'r').readlines() if not x.startswith('#') and x]

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


###############

def get_database_testing_path():
    """
    Get the path to the database testing file.
    :return: 
    """
    return get_database_path()+'/test_file.txt'


