#!/usr/bin/env python


import re
import sys
import os

from argparse import ArgumentParser

from jade.basic import path

def get_parser():
    parser = ArgumentParser(description="Renames original files to new names for design ordering.  Copy all models going to be ordered into a single directory first. Run from directory with pdb files already copied in!")
    parser.add_argument("-i", "--new_names",
                        help = "File with new to old names.  Example line: new_name  *  filename.  Can have lines that don't have all three.  Will only rename if it has a star in the second column.",
                        required = True)

    return parser



if __name__ == "__main__":

    parser = get_parser()

    options = parser.parse_args()


    INFILE = open(options.new_names, 'r')
    for line in INFILE:
        line = line.strip()
        lineSP = line.split()
        if not line or len(lineSP) != 3:
            continue

        new_name = lineSP[0].strip()
        create = lineSP[1].strip()
        old_name = lineSP[2].strip()

        if create != '*': continue

        print "Renaming "+old_name+" to "+new_name+path.get_decoy_extension(old_name)

        os.system("mv "+old_name+" "+new_name+path.get_decoy_extension(old_name))

    INFILE.close()



