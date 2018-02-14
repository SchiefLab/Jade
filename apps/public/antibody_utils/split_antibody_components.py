#!/usr/bin/env python

# Author: Jared Adolf-Bryfogle
# Description: Script for splitting AHO renumbered antibodies into FAB, Fv, Fc, and linker regions

from jade.antibody.split_structure import *
from argparse import ArgumentParser

def get_parser():
    parser = ArgumentParser(description="Script for splitting AHO renumbered antibodies into Fv, Fc, and linker regions")

    parser.add_argument("--any_structure",
                        default=False,
                        action = "store_true",
                        help = "Be default, we only output structures with both L/H.  Pass this option to split structures that are L or H only.")

    parser.add_argument("--ab_dir", "-a",
                        help = "Antibody Directory with AHO-renumbered structures to split. Can be .pdb, or .pdb.gz",
                        required = True)

    parser.add_argument("--output_dir", "-o",
                        help = "Output Directory for antibody structures.",
                        required = True)
    return parser

if __name__ == '__main__':

    parser = get_parser()
    options = parser.parse_args()

    run_main(options.ab_dir, options.output_dir, not options.any_structure)