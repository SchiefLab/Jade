import sys
import os
from collections import defaultdict

class BenchmarkInfo:
    """
    Simple Class for holding info for a particular benchmark.
    Parses the Run_Settings.txt file in the decoy directory.  This file is output by RunRosettaBenchmarks.

        The settings dictionary then holds key/value pairs.  Here is an example of this file for RAbD:

        CDR = ALL
        DATASET = bm2_ten
        DOCK = False
        INNER_CYCLE_ROUNDS = 1
        INPUT_PDB_TYPE = pareto
        L_CHAIN = kappa
        MINTYPE = relax
        OUTER_CYCLE_ROUNDS = 100
        PAPER_AB_DB = True
        PROTOCOL = even_cluster_mc
        RANDOM_START = True
        REMOVE_ANTIGEN = True
        SEPARATE_CDRS = False


    """
    def __init__(self, decoy_path, full_name, final_name, scorefunction = "talaris2014"):
        self._initialize_data(decoy_path, full_name, final_name,  scorefunction)

    def _initialize_data(self, decoy_path, full_name, final_name, scorefunction):

        self.decoy_dir = decoy_path
        self.full_name = full_name
        self.final_name = final_name
        self.scorefunction = scorefunction

        self.settings = get_run_settings(self.decoy_dir)


    def get_full_name(self):
        """
        Get the full name of the benchmark.
        :rtype: str
        """
        return self.full_name

    def get_final_name(self):
        """
        Get the final name of the benchmark (used mainly for features dbs or comparisons between benchmarks.)
        :rtype: str
        """
        return self.final_name

    def get_dataset(self):
        """
        Get the dataset used for benchmarking.
        :rtype: str
        """
        return self.settings["dataset"]

    def get_input_pdb_type(self):
        """
        Get the input pdb type used, ex: native vs pareto
        :rtype: str
        """
        return self.settings["input_pdb_type"]

    def get_scorefunction_name(self):
        """
        Get the scorefunction name set in this info instance.
        :rtype: str
        """
        return self.scorefunction

    def get_decoy_path(self):
        """
        Get the directory of all of the decoys for this benchmark.
        :rtype:
        """
        return self.decoy_dir

    ### RAbD ###
    def has_log_path(self):
        return False



def get_run_settings(dir, fname="RUN_SETTINGS.txt"):
    """
    Gets a dict of the settings used to run the benchmark in the directory.

    The settings file looks like this, and is output by RunRosettaBenchmarks into the decoy directory:

        CDR = ALL
        DATASET = bm2_ten
        DOCK = False
        INNER_CYCLE_ROUNDS = 1
        INPUT_PDB_TYPE = pareto
        L_CHAIN = kappa
        MINTYPE = relax
        OUTER_CYCLE_ROUNDS = 100
        PAPER_AB_DB = True
        PROTOCOL = even_cluster_mc
        RANDOM_START = True
        REMOVE_ANTIGEN = True
        SEPARATE_CDRS = False


    :param dir: str
    :rtype: defaultdict
    """

    fpath = dir+"/"+fname
    print "attempting to open "+fpath
    settings = defaultdict()
    if not os.path.exists(fpath):
        print fpath+" does not exist"
        return settings

    FILE = open(fpath)
    for line in FILE:
        line = line.strip()
        if not line or line.startswith("#"): continue
        lineSP = line.split('=')
        lineSP = [i.strip() for i in lineSP]
        settings[lineSP[0]] = lineSP[1]
        settings[lineSP[0].lower()] = lineSP[1]
    FILE.close()
    return settings
