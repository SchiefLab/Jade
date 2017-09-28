#!/usr/bin/env python

from __future__ import print_function

import os,sys,re
from argparse import ArgumentParser

from jade.basic.path import *

if __name__ == "__main__":

    parser = ArgumentParser("This extremely simple script checks nstruct of the input files and outputs which nstruct number is missing.")


    parser.add_argument('-n', "--nstruct", default=1000)

    parser.add_argument("pdb_files", help="Path to PDB files we will be checking.", nargs="*")

    options = parser.parse_args()

    n = int(options.nstruct)
    missing_nstruct = range(1, int(options.nstruct) + 1)

    for f in options.pdb_files:
        name = get_decoy_name(f)
        decoy_number = int(name[-4:])
        if decoy_number in missing_nstruct:
            missing_nstruct.remove(decoy_number)

    print("\nMissing nstructs: "+repr(missing_nstruct))

