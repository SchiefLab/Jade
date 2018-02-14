#!/usr/bin/env python

from __future__ import print_function
import os,sys,re
from argparse import ArgumentParser

def get_parser():
    parser = ArgumentParser(description= "This script places ters between ATOM/HETATM columns.  This is currently needed to reload symmetrized glycan poses"
                            "created by the god aweful make_symm_file.pl Rosetta script. USE: place_TERs.py my_pdb - Does it in place. ")

    parser.add_argument("pdb_files", help="Path to PDB files we will be stripping.", nargs="*")

    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()


    for pdb_file in options.pdb_files:
        print("Working on pdb_file")
        FILE = open(pdb_file, 'r')
        lines = []
        for line in FILE:
            line = line.strip()
            lines.append(line)
        FILE.close()



        OUTFILE = open(pdb_file, 'w')
        outputlines = []
        outputlines.append(lines[0])
        for i in range(1, len(lines) ):
            prev_line = lines[i-1]
            current_line = lines[i]
            outputlines.append(current_line)
            #print(i, len(lines))
            if i == len(lines) - 1:
                if not re.search("TER", current_line):
                    outputlines.append("TER   ")

                break

            #Not at the end.  Not at the beginning.
            next_line = lines[i+1]
            if re.search("ATOM", current_line) and re.search("HETATM", next_line) and not re.search("TER", current_line):
                print("Adding TER between ATOM and HETATM")
                outputlines.append("TER   ")

        for line in outputlines:
            OUTFILE.write(line+"\n")
        OUTFILE.close()
        print("COMPLETE")








