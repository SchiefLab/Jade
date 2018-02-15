#!/usr/bin/env python

from __future__ import print_function

import os,sys,re
from argparse import ArgumentParser

from jade.basic.path import *
from collections import defaultdict

def get_parser():
    parser = ArgumentParser(description="This extremely simple script checks nstruct of the input files and outputs which nstruct number is missing.")


    parser.add_argument('-n', "--nstruct", default=1000)

    parser.add_argument("--pdb_files", help="Path to PDB files we will be checking.", nargs="*")

    parser.add_argument("--pdblist", "-l", help = "Optional INPUT PDBLIST (without 00s, etc. for which to check")
    parser.add_argument("--dir", help = "The Directory to check. As opposed to a list of pdb files.")

    return parser

if __name__ == "__main__":


    parser = get_parser()
    options = parser.parse_args()

    print(options)

    n = int(options.nstruct)
    missing_nstruct = range(1, int(options.nstruct) + 1)

    pdb_files = []
    if options.pdb_files:
        pdb_files = options.pdb_files
    elif options.dir:
        pdb_files.extend(get_all_pdb_paths(options.dir, ".pdb"))
        pdb_files.extend(get_all_pdb_paths(options.dir, ".pdb.gz"))
    elif options.pdblist:
        dir = "."
        if os.path.dirname(options.pdblist):
            dir = os.path.dirname(options.pdblist)

        pdb_files.extend(get_all_pdb_paths(dir, ".pdb"))
        pdb_files.extend(get_all_pdb_paths(dir, ".pdb.gz"))

    if len(pdb_files) == 0:
        sys.exit("Could not find any PDB Files")

    if options.pdblist:
        decoy_counts = defaultdict()
        origin_pdbs = [get_decoy_name(x) for x in parse_contents(options.pdblist)]
        for x in origin_pdbs:
            decoy_counts[x] = 0

        for f in pdb_files:
            name = get_decoy_name(f)
            for origin_pdb in origin_pdbs:
                if re.search(origin_pdb, name):
                    decoy_counts[origin_pdb]+=1

        for origin_pdb in sorted(decoy_counts.keys()):
            if decoy_counts[origin_pdb] != n:
                print("Missing",n-decoy_counts[origin_pdb], "pdbs for", origin_pdb)

    else :
        for f in pdb_files:
            name = get_decoy_name(f)
            decoy_number = int(name[-4:])
            if decoy_number in missing_nstruct:
                missing_nstruct.remove(decoy_number)

        print("\nMissing nstructs: "+repr(missing_nstruct))

