#!/usr/bin/env python

import os
import sys
import glob
import re
import sqlite3
from collections import defaultdict
from argparse import ArgumentParser

#Module DB Imports

from rosetta_general.DesignBreakdown import *
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

    ############################
    ## Required Options
    ############################
    parser.add_argument("--PDBLIST",
                      help = "Analyze a PDBLIST of the strategy.  Relative or full paths",
                      default=None)

    parser.add_argument("--indir",
                      help = "Input directory used for either PDBLIST path (if PDBLIST and this is given or a path full of PDBs to analyze",
                      default = None)

    parser.add_argument("--out_name",
                      help = "Name of the output Features db and other out names")


    parser.add_argument("--out_db_batch",
                      help="Batch name for databases")


    ###########################
    ## Analysis Types
    ###########################

    parser.add_argument("--analysis",
                        help = "Analysis to run on PDBs",
                        default = "all",
                        choices = ["all", "cluster_features", "antibody_features"])


    ############################
    ## Optional
    ############################

    #parser.add_argument("--outdir",
    #                  help = "\nOutput directory",
    #                  default = "antibody_design_analysis_results")

    parser.add_argument("--score_weights",
                      help = "Weights to use during FeaturesReporters",
                      default = "talaris2013")

    #parser.add_argument("--compiler",
    #                    help = "Compiler.  If on mac will automatically choose clang.",
    #                    default = "gcc")

    parser.add_argument("--use_present_dbs",
                      default = False,
                      help = "Do not attempt to delete features databases present",
                      action = "store_true")



    #parser.add_argument("--np",
    #                    default = 5,
    #                    help = "Number of processors for MPI")


    run_mpi_rosetta = RunRosetta(program = "rosetta_scripts", parser = parser);


    options = parser.parse_args()

    options.use_mpi = True

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




    ####################################################################################################################
    ##                     Analysis Components.  One option for each type of analysis.
    ####################################################################################################################


    if not options.out_name:
        sys.exit("Output name required")


    #### Create Features Databases ###
    if options.analysis == "all" or options.analysis == "antibody_features":

        create_features_db(pdb_list,
                       'antibody_features',  options.compiler,
                       options.score_weights, options.out_name,
                       options.out_db_batch, options.outdir, options.use_present_dbs, "", options.mpi, options.np)

    if options.analysis == "all" or options.analysis == "cluster_features":

        create_features_db(pdb_list,
                       'cluster_features',  options.compiler,
                       options.score_weights, options.out_name,
                       options.out_db_batch, options.outdir, options.use_present_dbs, "", options.mpi, options.np)


    print "Complete"





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