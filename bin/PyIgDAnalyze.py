#!/usr/bin/env python

import os
import sys
import glob
import re
import sqlite3
from collections import defaultdict
from optparse import OptionParser

#Rosetta Imports
from rosetta import *

#Module DB Imports
from structure.Structure import *
from structure.PythonPDB2 import *
from rosetta_general.DesignBreakdown import *
from sequence.SequenceStats import SequenceStats
from pymol.PyMolScriptWriter import *
import create_features_json as json_creator


rosetta.init("-ignore_unrecognized_res -ignore_zero_occupancy false -ex1 -ex2 -use_input_sc")



def main():
    ####################################################################################################################
    ##                                                  OPTIONS
    ####################################################################################################################


    parser = OptionParser()

    ############################
    ## Required Options
    ############################
    parser.add_option("--PDBLIST",
                      help = "Analyze a PDBLIST of the strategy.  Relative or full paths",
                      default=None)

    parser.add_option("--indir",
                      help = "Input directory used for either PDBLIST path (if PDBLIST and this is given or a path full of PDBs to analyze",
                      default = None)

    parser.add_option("--out_name",
                      help = "Name of the output Features db and other out names")


    parser.add_option("--out_db_batch",
                      help="Batch name for databases")


    ###########################
    ## Analysis Types
    ###########################

    parser.add_option("--do_all",
                      help = "Run all available analysis",
                      default = False,
                      action = "store_true")

    parser.add_option("--do_run_antibody_features_all",
                      help = "Run the AntibodyFeatures to create databases for the strategy on all PDBs",
                      default = False,
                      action = "store_true")

    parser.add_option("--do_run_cluster_features_all",
                      help = "Run the ClusterFeatures reporter on all PDBs in list or directory",
                      default = False,
                      action = "store_true")


    ############################
    ## Optional
    ############################

    parser.add_option("--outdir",
                      help = "\nOutput directory",
                      default = "antibody_design_analysis_results")

    parser.add_option("--score_weights",
                      help = "Weights to use during rescore or FeaturesReporters",
                      default = "talaris2013")

    parser.add_option("--rosetta_extension",
                      default = "linuxclangrelease")

    parser.add_option("--use_present_dbs",
                      default = False,
                      help = "Do not attempt to delete features databases present",
                      action = "store_true")

    (options, args) = parser.parse_args(sys.argv)

    if not options.PDBLIST and not options.indir:
        sys.exit("Cannot analyze strategy without a PDBLIST or an input directory full of PDBs")
    if options.indir and not os.path.exists(options.indir):sys.exit(options.indir+" does not exist - cannot continue")




    ##### Get or make the path to the PDBLIST ######

    if not options.PDBLIST and options.indir:
        pdb_list = make_pdblist(options.indir)
    elif options.PDBLIST and options.indir:
        PDBLIST = open(options.PDBLIST, 'r')
        NEW_PDBLIST = open(options.indir+"/FULLPATH_PDBLIST.txt")

        for line in PDBLIST:
            new_line = options.indir+"/"+line.strip()
            if not os.path.exists(new_line):
                sys.exit(new_line+" does not exist - cannot continue")

            NEW_PDBLIST.write(new_line+"\n")
        NEW_PDBLIST.close()
        PDBLIST.close()
        pdb_list = options.indir+"/FULLPATH_PDBLIST.txt"
    else:
        pdb_list = options.PDBLIST



    if not os.path.exists(pdb_list): sys.exit(pdb_list+" does not exist!")

    os.system('mkdir '+options.outdir)


    do_all = options.do_all


    ####################################################################################################################
    ##                     Analysis Components.  One option for each type of analysis.
    ####################################################################################################################


    if not options.out_name:
        sys.exit("Output name required")


    #### Create Features Databases ###
    if do_all or options.do_run_antibody_features_all:

        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        create_features_db(pdb_list, fdir,
                       'antibody_features',  options.rosetta_extension,
                       options.score_weights,options.out_name,
                       options.out_db_batch, options.outdir, options.use_present_dbs)

    if do_all or options.do_run_cluster_features_all:

        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        create_features_db(pdb_list, fdir,
                       'cluster_features',  options.rosetta_extension,
                       options.score_weights,options.out_name,
                       options.out_db_batch, options.outdir, options.use_present_dbs)


    print "Complete"

def rm_features_db(outdir, out_name, score_weights, xml_name):
    db_name = outdir+'/'+out_name+'.'+xml_name+'.'+score_weights+".db3"
    if os.path.exists(db_name):
        os.remove(db_name)

def create_features_db(pdb_list, features_dir, xml_name, extension, score_weights, out_db_name, out_db_batch, outdir, use_present_dbs, indir = ""):

    """
    old_db_name = outdir+'/'+out_db_name+'.'+score_weights+".db3"
    new_db_name = outdir+'/'+out_db_name+'.'+xml_name+'.'+score_weights+".db3"
    if os.path.exists(old_db_name):
        os.system('mv '+old_db_name+' '+new_db_name)
        print "Old db name already exists.  Moving."
        return
    """

    outdir = outdir+"/databases"
    if not os.path.exists(outdir): os.mkdir(outdir)

    if not out_db_batch:
        sys.exit("Must have output db name and db batch to run_features")

    if not use_present_dbs:
        rm_features_db(outdir, out_db_name, score_weights, xml_name)

    features_command = 'rosetta_scripts.'+extension +' -parser:protocol '+features_dir+\
                       '/'+xml_name+'.xml @ '+features_dir+'/features.flag -l '+pdb_list+ \
                       ' -parser:script_vars name='+out_db_name+'.'+xml_name+ \
                       ' score='+score_weights+' batch='+out_db_batch

    if indir:
        features_command = features_command+' in:path:pdb '+indir

    os.system(features_command)
    os.system('mv '+out_db_name+'.'+xml_name+'.'+score_weights+'.db3'+' '+outdir)

def get_data_from_features_db(db_path):
    command = """
                SELECT
                    structures.input_tag,
                    interfaces.dG as dG,
                    structure_scores.score_value as total_score,
                    interfaces.dG_cross as dG_cross,
                    interfaces.delta_unsatHbonds as delta_unsatHbonds,
                    interfaces.hbond_E_fraction as hbond_E_fraction,
                    interfaces.dSASA as dSASA,
                    interfaces.interface as interface
                FROM
                    interfaces,
                    score_types,
                    structure_scores,
                    structures
                WHERE
                    score_types.score_type_name='total_score' AND
                    structure_scores.score_type_id = score_types.score_type_id AND
                    structure_scores.struct_id = interfaces.struct_id AND
                    structure_scores.struct_id = structures.struct_id
                ORDER BY total_score
            """

    db = sqlite3.connect(db_path)
    cur = db.cursor()


    for row in cur.execute(command):
        pass


def make_pdblist(in_dir):
    pdb_list = glob.glob(in_dir+"/*.pdb*")
    OUTFILE = open(in_dir+"/PDBLIST.txt", 'w')
    for pdb in pdb_list:
        if re.search("initial_benchmark_perturbation", pdb) or re.search("excn", pdb):
            continue
        OUTFILE.write(pdb+"\n")
    OUTFILE.close()

    return in_dir+"/PDBLIST.txt"



if __name__ == "__main__":
    main()