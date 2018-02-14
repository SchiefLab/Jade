#!/usr/bin/env python

import os
import sys
import sqlite3
from collections import defaultdict
from argparse import ArgumentParser
import glob
import re

from jade.pymol_jade.PyMolScriptWriter import PyMolScriptWriter

def get_parser():
    parser = ArgumentParser(description="This App aims to make pymol alignments using the PyIgClassify database and structures, matching specific criterion.")

    required = parser.add_argument_group("Required Arguments")

    required.add_argument("--db", "-d",
                        help = " Database to use from PyIgClassify.",
                        required = True)

    required.add_argument("--ab_dir", "-b",
                        help = "Directory with renumbered antibody PDBs (Full or CDRs-only)",
                        required = True)

    required.add_argument("--where", "-w",
                        help = "Your where clause for the db in quotes.  Not including WHERE. Use ' ' for string matches",
                        required = True)

    optional = parser.add_argument_group("Other Arguments")

    optional.add_argument("--outdir", "-o",
                        help = "Output directory.",
                        default = os.getcwd())

    optional.add_argument("--prefix", "-p",
                        help = "Output prefix")

    optional.add_argument("--cdr", "-c",
                        help = "Optionally load the CDR PDBs of the given type in the ab_dir. If this option is set, the ab_dir should be of CDRs only from PyIgClassify.")

    optional.add_argument("--native", "-n",
                          help = "Align everything to this PDB, the native or something you are interested in. ")

    return parser

if __name__ == "__main__":



    parser = get_parser()
    options = parser.parse_args()




    ###########
    ab_db = sqlite3.connect(options.db)


    ### Get Matching PDBs ###

    matching_files = defaultdict(int)
    if options.cdr:
        for row in ab_db.execute("select PDB, original_chain FROM cdr_data WHERE "+options.where +" and CDR='"+options.cdr+"'"):
            pdb = row[0].lower()+row[1]+"_"+options.cdr.upper()
            matching_files[pdb] += 1
    else:
        query = "select PDB, original_chain FROM cdr_data WHERE "+options.where
        print query
        for row in ab_db.execute(query):
            pdb = row[0].lower()+row[1]
            matching_files[pdb ] += 1


    print "Globbing: "+options.ab_dir+"/*.pdb*"
    all_pdb_files = glob.glob(options.ab_dir+"/*.pdb*")
    pdb_paths = []

    PDBLIST = open(options.outdir+"/"+options.prefix+"PDBLIST.txt", 'w')

    print "Matching PDB Files:"
    for found_pdb in matching_files:
        for globed_pdb in all_pdb_files:
            if re.search(found_pdb, globed_pdb):
                pdb_paths.append(globed_pdb)
                print globed_pdb
                PDBLIST.write(globed_pdb+"\n")
                break

    PDBLIST.close()
    if len(pdb_paths) == 0:
        sys.exit("No matching pdb paths found!")

    ### Make the PyMol session and save the script ###
    scripter = PyMolScriptWriter(options.outdir)

    if options.native:
        scripter.add_load_pdb(options.native)

    scripter.add_load_pdbs(pdb_paths)

    if options.native:
        scripter.add_align_all_to(options.native, limit_to_bb=False)
    else:
        scripter.add_align_all_to(pdb_paths[0], limit_to_bb=False)

    scripter.add_show("cartoon")
    scripter.add_center()
    scripter.add_antibody_script()
    scripter.add_save_session(options.outdir+"/"+options.prefix+"ab_session.pse")
    scripter.run_script()
