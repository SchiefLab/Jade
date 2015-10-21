#!/usr/bin/env python

import os
import sys
import glob
import re
import sqlite3
from collections import defaultdict
from argparse import ArgumentParser

#Module DB Imports

from rosetta_general.features import *
from rosetta_general.RunRosetta import RunRosetta

from pymol.PyMolScriptWriter import *
from tools.general import *

def main():

    ####################################################################################################################
    ##                                                  OPTIONS
    ####################################################################################################################


    parser = ArgumentParser("Creates Features Databases for antibody design using MPI.  "
                            "This uses RunRosettaMPI, so that it can be run locally or on a cluster.")

    analysis_options = parser.add_argument_group("RAbD Analyze", "Options specific for Analysis of RAbD Output Structures")

    ############################
    ## Required Options
    ############################
    #parser.add_argument("--PDBLIST",
    #                  help = "Analyze a PDBLIST of the strategy.  Relative or full paths",
    #                  default=None)

    analysis_options.add_argument("--indir",
                      help = "Input directory used for either PDBLIST path (if PDBLIST and this is given) or a path full of PDBs to analyze",
                      default = None)

    analysis_options.add_argument("--analysis",
                        help = "Analysis to run on PDBs",
                        default = "all",
                        choices = ["all", "cluster_features", "antibody_features"])

    analysis_options.add_argument("--use_present_dbs",
                      default = False,
                      help = "Do not attempt to delete features databases present",
                      action = "store_true")

    analysis_options.add_argument("--db_prefix",
                      help = "Prefix to use for output databases.  Recommended to use the design and strategy name",
                      default = "")

    run_mpi_rosetta = RunRosetta(program = "rosetta_scripts", parser = parser)

    options = parser.parse_args()


    if not options.l and not options.indir:
        sys.exit("Cannot analyze strategy without a PDBLIST or an input directory full of PDBs")
    if options.indir and not os.path.exists(options.indir):sys.exit(options.indir+" does not exist - cannot continue")

    if options.outdir == "decoys":
        run_mpi_rosetta.set_outdir("databases")


    ##### Get or make the path to the PDBLIST ######

    pdb_lists = []
    if not options.l and options.indir:
        pdb_lists = make_pdblists(options.indir)
    elif options.l and options.indir:
        PDBLIST = open(options.l, 'r')
        NEW_PDBLIST = open(options.indir+"/FULLPATH_PDBLIST.txt")

        for line in PDBLIST:
            new_line = options.indir+"/"+line.strip()
            if not os.path.exists(new_line):
                sys.exit(new_line+" does not exist - cannot continue")

            NEW_PDBLIST.write(new_line+"\n")
        NEW_PDBLIST.close()
        PDBLIST.close()
        pdb_list = options.indir+"/FULLPATH_PDBLIST.txt"
        pdb_lists.append(["", pdb_list])
    else:
        pdb_lists.append(["",options.l])


    for pdb_list in pdb_lists:
        if not os.path.exists(pdb_list[1]): sys.exit(pdb_list+" does not exist!")


    ####################################################################################################################
    ##                     Analysis Components.  One option for each type of analysis.
    ####################################################################################################################

    possible_names = [x + "_"+run_mpi_rosetta.options.db_name for x in ["ab", "cl", "rela_ab", "rela_cl", "norm_ab", "norm"]]

    db_suffix = run_mpi_rosetta.options.db_name

    print repr(possible_names)
    if not options.use_present_dbs:
        rm_features_dbs(run_mpi_rosetta.options.outdir, possible_names)

    #### Create Features Databases ###

    for pdb_list in pdb_lists:
        print "PDBList "+ repr(pdb_list)
        if options.analysis == "all" or options.analysis == "antibody_features":
            run_mpi_rosetta.set_json_run("antibody_features.json")
            run_mpi_rosetta.options.l = pdb_list[1]
            run_mpi_rosetta.options.db_name = options.db_prefix+"."+pdb_list[0]+"ab_"+db_suffix
            run_mpi_rosetta.options.job_name = run_mpi_rosetta.options.job_name+"_"+pdb_list[0]
            print "DB Name " + run_mpi_rosetta.options.db_name
            run_mpi_rosetta.run()

            '''
            create_features_db(pdb_list,
                           'antibody_features',  options.compiler,
                           options.score_weights, options.out_name,
                           options.out_db_batch, options.outdir, options.use_present_dbs, "", options.mpi, options.np)
            '''

        if options.analysis == "all" or options.analysis == "cluster_features":
            run_mpi_rosetta.set_json_run("cluster_features.json")
            run_mpi_rosetta.options.l = pdb_list[1]
            run_mpi_rosetta.options.db_name = options.db_prefix+"."+pdb_list[0]+"cl_"+db_suffix
            run_mpi_rosetta.options.job_name = run_mpi_rosetta.options.job_name+"_"+pdb_list[0]
            run_mpi_rosetta.run()

            '''
            create_features_db(pdb_list,
                           'cluster_features',  options.compiler,
                           options.score_weights, options.out_name,
                           options.out_db_batch, options.outdir, options.use_present_dbs, "", options.mpi, options.np)
            '''



def make_pdblists(in_dir):

    pdblists = [] #List of PDBLists.

    pdb_list = []
    rel_list = []

    pdbs = glob.glob(in_dir+"/*.pdb*")

    relaxed_models = False

    for pdb in pdbs:
        if re.search("pre_model", pdb): relaxed_models = True

    for pdb in pdbs:
        if re.search("pre_model", pdb):
            pdb_list.append(pdb)
        elif re.search("ds_rel", pdb):
            rel_list.append(pdb)
        elif relaxed_models:
            rel_list.append(pdb)
        else:
            pdb_list.append(pdb)


    OUTFILE = open(in_dir+"/PDBLIST.txt", 'w')

    for pdb in pdb_list:
        if re.search("initial_benchmark_perturbation", pdb) or re.search("excn", pdb):
            continue
        OUTFILE.write(pdb+"\n")
    OUTFILE.close()

    pdblists.append(["norm_", in_dir+"/PDBLIST.txt"])


    if len(rel_list) > 0:
        REL_OUTFILE = open(in_dir+"/REL_PDBLIST.txt", 'w')
        for pdb in rel_list:


            if re.search("initial_benchmark_perturbation", pdb) or re.search("excn", pdb):
                continue
            REL_OUTFILE.write(pdb+"\n")
        REL_OUTFILE.close()
        pdblists.append(["rela_", in_dir+"/REL_PDBLIST.txt"])

    return pdblists
if __name__ == "__main__":
    main()