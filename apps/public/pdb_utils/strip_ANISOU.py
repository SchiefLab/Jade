#!/usr/bin/env python

from __future__ import print_function
from argparse import ArgumentParser
import re

def get_parser():
    parser = ArgumentParser(description=" Strips ANISOU lines out of PDBs.")

    parser.add_argument("pdb_files", help = "Path to PDB file we will be stripping.", nargs="*")

    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()

    for pdb_file in options.pdb_files:
        FILE = open(pdb_file, 'r')
        lines = FILE.readlines()
        FILE.close()



        OUTFILE = open(pdb_file, 'w')
        for line in lines:
            if not re.search("ANISOU", line):
                OUTFILE.write(line)
        OUTFILE.close()
