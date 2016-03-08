#!/usr/bin/env python

from argparse import ArgumentParser
import os
import sys
import re

parser = ArgumentParser(" This simple script strips ters out of a PDB file and overwrites the input.  PyMol places ters "
                        "when th numbering is not 1-1.  And then Rosetta will F your Shit up.")

parser.add_argument("pdb_file", help = "Path to PDB file we will be stripping.")


options = parser.parse_args()

FILE = open(options.pdb_file, 'r')
lines = FILE.readlines()
FILE.close()



OUTFILE = open(options.pdb_file, 'w')
for line in lines:
    if not re.search("TER", line):
        OUTFILE.write(line)
OUTFILE.close()
