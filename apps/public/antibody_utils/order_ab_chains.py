#!/usr/bin/env python

#Jared Adolf-Bryfogle
#Reorders PDBFiles in a dirctory according to A_LH in order for Rosetta Antibody Design benchmarking. Removes HetAtm!!!

import sys
import os

from jade.basic.structure.BasicPose import *
from argparse import ArgumentParser

def reorder_and_save_chains(in_path, out_path, remove_het = False):
    blank_pdb = BasicPose()
    full_pdb = BasicPose(in_path)


    blank_pdb.copy_all_but_chains_into_pdb_map(full_pdb, ["L", "H"])
    blank_pdb.copy_chain_into_pdb_map(full_pdb, "L")
    blank_pdb.copy_chain_into_pdb_map(full_pdb, "H")
    if remove_het:
        blank_pdb.remove_hetatm_atoms()
        blank_pdb.remove_waters()

    blank_pdb.save_PDB(out_path)

def get_parser():
    parser = ArgumentParser(description="Reorders PDBFiles in a dirctory according to A_LH in order for Rosetta Antibody Design benchmarking. Removes HetAtm")

    parser.add_argument("--in_dir", "-i",
                        default = os.getcwd(),
                        help = "Input Directory of PDB files listed in any passed PDBLIST. Default=PWD")

    parser.add_argument("--in_pdblist", "-l",
                        help = "Input PDBList file. Assumes PDBList has no paths and requires an input directory as if we run Rosetta.",
                        default = "")

    parser.add_argument("--in_single", "-s",
                        help = "Path to Input PDB File, instead of list.",
                        default = "")

    parser.add_argument("--out_dir", "-d",
                        help = "Output Directory. Resultant PDB files will go here.",
                        default = "reordered")

    return parser
if __name__ == "__main__":


    parser = get_parser()
    options = parser.parse_args()

    in_dir = options.in_dir
    in_pdblist = options.in_pdblist
    out_dir = options.out_dir
    in_single = options.in_single

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    if in_pdblist:
        PDBLIST = open(in_pdblist, 'r')
        for line in PDBLIST:
            line = line.strip()
            line = in_dir+"/"+line

            print "Reordering "+line
            outpath = out_dir+"/"+os.path.basename(line)
            reorder_and_save_chains(line, outpath)
        PDBLIST.close()

    elif in_single:
        outpath = os.path.join(out_dir, in_single)
        reorder_and_save_chains(in_single, outpath)
    else:
        sys.exit("Must pass either -s or -l ")

    print "Done!!"
