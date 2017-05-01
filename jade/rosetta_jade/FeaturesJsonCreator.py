#!/usr/bin/env python

# Author: Jared Adolf-Bryfogle

#This script will create either cluster features or antibody features json for use in Features R script.
#I am just really sick of doing this by hand.
#Example Cmd-line:  python create_features_json.py --database databases/baseline_comparison.txt --scripts cluster

import json
import shutil
from collections import defaultdict

from jade.basic.path import *
from jade.basic.threading.Threader import *

########################################################################################################################
####### Class based implementation
########################################################################################################################

class JsonCreator:
    """
    Basic implementation of a simple JsonCreator to create Jsons.  Could be expanded to not load jsons with pre-set scripts.
    A nicer implementation would be a GUI for running the FeaturesReporter scripts.
    """
    def __init__(self, out_path, script_type):
        self.out_path = out_path
        self.script_types = ["antibody", "interface", "cluster", "antibody_minimal","antibody_minimal_hbond_analysis", "antibody_minimal_no_hbond_analysis"]


        if script_type == "antibody_minimal_min_hbond_analysis":
            script_type = "antibody_minimal_hbond_analysis"

        if not script_type in self.script_types:
            sys.exit(script_type +" unrecognized.  Available JSON script types are: "+repr(self.script_types))

        self.json_dict = initialize_json_dict(out_path)
        setup_baseline_scripts_and_formats(self.json_dict, script_type)

    def add_sample_source_info(self, db_path, id, ref = False):
        info = _create_json_info_dic(db_path, id, ref)
        add_sample_source(self.json_dict, info)

    def add_features_script(self, rel_script_path):
        """
        Add a features script to run.
        """
        self.json_dict["analysis_scripts"].append(rel_script_path)

    def add_output_method(self, output_method):
        """
        Add an output method
        """
        self.json_dict["output_formats"].append(output_method)




    def save_json(self, out_path = "local_json.txt"):
        OUTFILE = open(out_path, 'w')


        json.dump(self.json_dict, OUTFILE, indent=1)
        OUTFILE.close()
        self.json_path = out_path

        print "Json written to: "+out_path
        print "Json path set to JsonCreator"

    def run_json(self, backround = False):
        run_features_json(self.json_path, backround, self.out_path)



###################################################################################################################


def setup_baseline_scripts_and_formats(json_dict, type):
    feat_path = get_rosetta_features_json_path()

    FILE = open(feat_path+"/"+type+"_features.json")
    data = json.load(FILE)
    FILE.close()
    append_scripts_formats_to_json_dict(data, json_dict)

    FILE = open(feat_path+"/output_formats.json")
    data = json.load(FILE)
    FILE.close()
    json_dict["output_formats"] = data["output_formats"]

def append_scripts_formats_to_json_dict(data, json_dict):
    if data.has_key("output_formats"):
        json_dict["output_formats"].extend(data["output_formats"])
    if data.has_key("analysis_scripts"):
        json_dict["analysis_scripts"].extend(data["analysis_scripts"])

def initialize_json_dict(out_dir):

    json_dict = defaultdict()
    json_dict["output_formats"] = []
    json_dict["analysis_scripts"] = []
    json_dict["sample_sources"] = []
    json_dict["output_dir"] = out_dir

    return json_dict

'''
def add_sample_source_comparisons(current_json_dict, json_file_path = None):
    if json_file_path:
        if not os.path.exists(json_file_path):
            sys.exit("Passed Json file path does not exist!!")
        FILE = open(json_file_path, 'r')
        json_dict = json.load(FILE)
        FILE.close()
        current_json_dict["sample_source_comparisons"].extend(json_dict["sample_source_comparisons"])
'''

def add_sample_source(json_dict, sample_source_dict):
    json_dict["sample_sources"].append(sample_source_dict)


def write_json_for_single_recovery_experiment(db_path_exp, db_path_natives, exp_id, out_path):
    """
    Create a JSON file for recovery of a single experiment.
    """
    if not os.path.exists(db_path_exp):
        sys.exit("SQLITE3 database path does not exist!!")

    json_dict = initialize_json_dict( out_path )
    setup_baseline_scripts_and_formats(json_dict, "cluster")

    #Ref HAS to go first here!
    add_sample_source(json_dict, _create_json_info_dic(db_path_natives, "ref", True))
    add_sample_source(json_dict, _create_json_info_dic(db_path_exp, exp_id, False))


    if not os.path.exists(out_path):
        os.mkdir(out_path)
    if not os.path.exists(out_path+"/jsons"):
        os.mkdir(out_path+"/jsons")
    OUTFILE = open(out_path+"/jsons/"+"cluster_features."+exp_id+".json", 'w')
    json.dump(json_dict, OUTFILE, indent=1)
    OUTFILE.close()

    print "Complete"


########################################################################################################################
def _create_json_info_dic(db_path, id, ref=False):
    """
    Create a simple dict that holds info to create the json.
    """
    databases = defaultdict()
    databases["database_path"] = db_path
    databases["id"] = id
    databases["reference"] = ref

    return databases




def run_features_json(json_path, backround = False, outpath = ""):
    """
    Convenience function
    Outputs an R script for running a JSON file, and runs it.
    Works with the new Library structure of the Features Reporter Framework.
    """
    def move_delete_build():
        """
        Used if outpath is not set in the json.
        :return:
        """
        if os.path.exists("build"):
            c = glob.glob("build/*")
            for d in c:
                os.system("cp -r "+d+" "+outpath)
        if os.path.exists("build"):
            shutil.rmtree("build")



    r_cmd = get_rosetta_features_run_script()+" "+json_path

    print "Running: "+r_cmd

    if backround:
        thread = Threader()

        f = []
        f.append(lambda: os.system(r_cmd))
        f.append(lambda: move_delete_build())
        thread.run_functions(f)

    else:
        os.system(r_cmd)
        move_delete_build()

def run_features_json_old(json_path, backround = False, outpath = ""):
    """
    Convenience function
    Run compare_sample_sources with json path.
    """
    def move_delete_build():
        if os.path.exists("build"):
            c = glob.glob("build/*")
            for d in c:
                os.system("cp -r "+d+" "+outpath)
        if os.path.exists("build"):
            shutil.rmtree("build")

    r_cmd = get_rosetta_features_root()+"/compare_sample_sources.R --config "+json_path
    print "Running: "+r_cmd

    if backround:
        thread = Threader()

        f = []
        f.append(lambda: os.system(r_cmd))
        f.append(lambda: move_delete_build())
        thread.run_functions(f)

    else:
        os.system(r_cmd)
        move_delete_build()

