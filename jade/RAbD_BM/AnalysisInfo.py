import json, os, sys

import jade.rosetta_jade.BenchmarkInfo as rosetta_bms
from jade.basic.path import get_decoy_name


class AnalysisInfo:
    """
    Simple class that parses a json file which defines (USING RELATIVE PATHS):

    1) exp - The name of the experiment - whatever you want it to be.
    2) decoy_dir - the directory of the decoys.
    3) features_db - the db where the features reporters have been run.

    The class will store this information, and parse the benchmark info in the decoy dir, storing a BenchmarkInfo object.
    Benchmark classes and scripts will take lists to these analysis files and use them to generate plots and data.
    """
    def __init__(self, json_path):
        print "Loading "+json_path
        self.json_path = json_path
        JSON_FILE = open(self.json_path, 'r')
        self.data = json.load(JSON_FILE)
        JSON_FILE.close()

        self.bm_info = rosetta_bms.BenchmarkInfo(self.get_decoy_dir(), os.path.basename(self.get_decoy_dir()), self.get_exp())

    def get_bm_info(self):
        """
        Get the benchmark info
        :rtype: rosetta_bms.BenchmarkInfo
        """
        return self.bm_info

    def get_exp(self):
        return self.data["short_name"]

    def get_decoy_dir(self):
        return self.data["decoy_dir"]

    def get_features_db(self):
        return self.data["features_db"]

class NativeInfo:
    """
    Simple class to hold native information.
    """
    def __init__(self, dataset, input_pdb_type, root_dataset_dir="datasets"):
        def load_pdbids(path):
            """
            Read PDBLISTS and return an array of PDB names
            :param path: str
            :rtype: [str]
            """
            pdbids = []
            print "reading from "+path
            if not os.path.exists(path):
                return pdbids

            INFILE = open(path, 'r')
            for line in INFILE:
                line = line.strip()
                if not line or line.startswith("#"):continue
                pdbids.append(get_decoy_name(line))
            INFILE.close()
            return pdbids

        self.dataset = dataset
        self.input_pdb_type = input_pdb_type
        self.root_dataset_dir = root_dataset_dir

        self.lambda_pdbids = self.pdbids = load_pdbids(os.path.join(root_dataset_dir, "pdblists",".".join([self.dataset, "lambda.PDBLIST.txt"])))
        self.kappa_pdbids = self.pdbids = load_pdbids(os.path.join(root_dataset_dir, "pdblists", ".".join([self.dataset, "kappa.PDBLIST.txt"])))

        self.pdbids = self.lambda_pdbids
        self.pdbids.extend(self.kappa_pdbids)


        self.db_path = os.path.join(self.root_dataset_dir, "databases", ".".join([self.input_pdb_type, "native_ab_features.db"]))
        if not os.path.exists(self.db_path):
            sys.exit("DB Path, "+self.db_path+" , does not exist. Please use the generate_rabd_features_dbs script to generate.")
        self.decoy_path = os.path.join(self.root_dataset_dir, self.input_pdb_type)

    def get_features_db(self):
        return self.db_path

    def get_decoy_dir(self):
        return self.decoy_path