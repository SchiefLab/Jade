#!/usr/bin/env python

# Author: Jared Adolf-Bryfogle

#This script will create either cluster features or antibody features json for use in Features R script.
#I am just really sick of doing this by hand.
#Example Cmd-line:  python create_features_json.py --database databases/baseline_comparison.txt --scripts cluster

import subprocess
from argparse import ArgumentParser
import os
import sys
import json
import shutil
from collections import defaultdict
from tools.path import *
from tools.Threader import *


def setup_baseline_scripts_and_formats(json_dict, type):
    db_path = get_database_path()

    FILE = open(db_path+"/"+type+"_features.json")
    data = json.load(FILE)
    FILE.close()
    append_scripts_formats_to_json_dict(data, json_dict)

def append_scripts_formats_to_json_dict(data, json_dict):
    json_dict["sample_source_comparisons"][0]["output_formats"].extend(data["sample_source_comparisons"][0]["output_formats"])
    json_dict["sample_source_comparisons"][0]["analysis_scripts"].extend(data["sample_source_comparisons"][0]["analysis_scripts"])

def initialize_json_dict(out_dir):

    formats = defaultdict()
    formats["output_formats"] = []
    formats["analysis_scripts"] = []
    formats["sample_sources"] = []

    json_dict = defaultdict()
    json_dict["sample_source_comparisons"] = []
    json_dict["sample_source_comparisons"].append(formats)
    json_dict["output_dir"] = out_dir

    return json_dict


def add_sample_source_comparisons(current_json_dict, json_file_path = None):
    if json_file_path:
        if not os.path.exists(json_file_path):
            sys.exit("Passed Json file path does not exist!!")
        FILE = open(json_file_path, 'r')
        json_dict = json.load(FILE)
        FILE.close()
        current_json_dict["sample_source_comparisons"].extend(json_dict["sample_source_comparisons"])

def add_sample_source(json_dict, sample_source_dict):
    json_dict["sample_source_comparisons"][0]["sample_sources"].append(sample_source_dict)


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
        self.script_types = ["antibody", "interface", "cluster", "antibody_minimal","antibody_min_hbond_analysis", "antibody_no_hbond_analysis", "antibody_minimal_min_hbond_analysis"]
        if not script_type in self.script_types:
            sys.exit(script_type +" unrecognized.  Available JSON script types are: "+repr(self.script_types))

        self.json_dict = initialize_json_dict(out_path)
        setup_baseline_scripts_and_formats(self.json_dict, script_type)

    def add_sample_source_info(self, db_path, id, ref = False):
        info = _create_json_info_dic(db_path, id, ref)
        add_sample_source(self.json_dict, info)

    def add_current_sample_sources_to_json_dict(self, json_file_path):
        """
        Combine a JSON with data held in this class
        """
        if os.path.exists(json_file_path):
            add_sample_source_comparisons(self.json_dict, json_file_path)

    def save_json(self, out_path = "local_json.txt"):
        OUTFILE = open(out_path, 'w')


        json.dump(self.json_dict, OUTFILE, indent=1)
        OUTFILE.close()
        self.json_path = out_path

        print "Json written to: "+out_path
        print "Json path set to JsonCreator"

    def run_json(self, backround = False):
        run_features_json(self.json_path, backround, self.out_path)


def run_features_json(json_path, backround = False, outpath = ""):
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

if __name__ == "__main__":


    parser = ArgumentParser()

    parser.add_argument("--databases", "-l",
        help = "List of dbs: db_name,short_name,ref keyword if the reference databaseSeparated by white space. ",
        default = [],
        nargs = "*",


    )

    parser.add_argument("--script", "-s",
        help = "Script type.  Will setup the appropriate output formats and R scripts",
        default = "antibody_minimal",
        choices = ["cluster", "antibody", "interface", "antibody_minimal"]
    )

    parser.add_argument("--db_path", "-p",
        help = "Path to databases.  Default is pwd/databases",
        default=os.getcwd()+"/databases")

    parser.add_argument("--outdir", "-o",
        help = "Where to put the result of the analysis scripts.  Currently unsupported by the features framework.",
        default=os.getcwd()+"/plots")

    parser.add_argument("--outname", "-n",
        help = "Output file name of json file",
        default = "local_json_compare_ss.json" )

    parser.add_argument("--add_comparison_to_this_json", "-a",
        help = "Add all this data to this json as more sample sources.")

    parser.add_argument("--run", "-r",
        help = "Go ahead and run compare_sample_sources.R.  Must be in path!!",
        default = False,
        action = "store_true")

    options = parser.parse_args()

    if not options.databases:
        sys.exit("No databases given!!")

    if not os.path.exists(options.db_path):
        options.db_path = os.getcwd()

    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)

    json_creator = JsonCreator(options.outdir, options.script)

    for db in options.databases:
        dbSP = db.split(',')

        ref = False
        if len(dbSP) == 3 and dbSP[2] == "ref": ref = True
        db_path = options.db_path+"/"+dbSP[0]; id= dbSP[1]
        json_creator.add_sample_source_info(db_path, id, ref)

    if options.add_comparison_to_this_json:
        json_creator.add_current_sample_sources_to_json_dict(options.add_comparison_to_this_json)

    if not re.search("json", options.outname): options.outname = options.outname+".json"

    json_creator.save_json(options.outdir+"/"+options.outname)

    if options.run:
        json_creator.run_json()

    print "Complete"
