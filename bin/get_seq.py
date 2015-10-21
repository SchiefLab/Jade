#!/usr/bin/env python

#from sequence import fasta
from tools import biopython_util

from collections import defaultdict
import os
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Uses Biopython to print sequence information.  Example:\n"
                                                 "get_seq.py --pdb 2j88_A.pdb --format fasta --outpath test.txt")

    parser.add_argument("--pdb", "-s",
                        help = "Input PDB")

    parser.add_argument("--chain", "-c",
                        help = "A specific chain to output",
                        default = "")

    parser.add_argument("--format", "-f",
                        help = "The output format requried.",
                        default = "fasta",
                        choices = ["basic", "fasta"])

    parser.add_argument("--outpath", "-o",
                        help = "Output path.  If none is specified it will write to screen.")

    parser.add_argument("--prefix", "-t",
                        help = "Tag to add before chain")

    parser.add_argument("--region",
                        help = "specify a particular region, start:end:chain")

    options = parser.parse_args()

    #if not options.pdb:
    #    options.pdb = sys.argv[1]


    biostructure = biopython_util.get_biopython_structure(options.pdb)

    sequences = defaultdict()

    if options.chain:
        seq = biopython_util.get_seq_from_biostructure(biostructure, options.chain)
        sequences[options.chain] = seq
    else:
        for biochain in biostructure[0]:
            if biopython_util.get_chain_length(biochain) == 0:
                continue
            seq = biopython_util.get_seq_from_biochain(biochain)
            sequences[biochain.id] = seq

    outlines = []
    for chain in sequences:
        if options.format == "basic":
            outlines.append( options.prefix+chain+" : "+sequences[chain]+"\n")
        elif options.format == "fasta":
            outlines.append("> "+options.prefix+chain)
            outlines.append(sequences[chain]+"\n")


    if options.outpath:
        OUT = open(options.outpath, "w")
        for line in outlines:
            OUT.write(line+"\n")
        OUT.close()
    else:
        for line in outlines:
            print(line)

