#!/usr/bin/env python

#Deletes chains from PDB.  Not CIF, nothing fancy here for now.  Just single char chain IDs

#22
from __future__ import print_function
from argparse import ArgumentParser
import os
from jade.basic.path import open_file, get_decoy_name


def conditional_keep_lines(keep_chains, delete_chains, current_chain, outlines):
    if (len(keep_chains) > 0) and current_chain in keep_chains:
        outlines.append(line)
    elif (len(delete_chains) > 0) and current_chain not in delete_chains:
        outlines.append(line)

if __name__ == "__main__":

    parser = ArgumentParser("Delete chains from a PDB (not cif) file. Only works with ATOM records.  ")

    parser.add_argument("-s", help = "PDB file", required = True)

    parser.add_argument("-c", "--chains", help = "Comma-separated list of chains to remove")

    parser.add_argument("-k", "--keep_chains", help = "Chains to keep.  Use this instead of the chains option")

    parser.add_argument("--suffix", help = "Suffix to add to output file", default = "_del")

    options = parser.parse_args()



    INFILE = open_file(options.s, 'r')
    outlines = []

    delete_chains = []
    keep_chains = []
    if options.chains:
        delete_chains = options.chains.split(',')
    if options.keep_chains:
        keep_chains = options.keep_chains.split(',')

    for line in INFILE:
        current_chain = ""
        linked_chain = ""
        if line.startswith("ATOM"):
            current_chain = line[21]
            #print(line)
            #print(current_chain)

            conditional_keep_lines(keep_chains, delete_chains, current_chain, outlines)
            continue

        elif line.startswith("LINK"):
            current_chain = line[21]
            linked_chain = line[51]

            conditional_keep_lines(keep_chains, delete_chains, current_chain, outlines)
            conditional_keep_lines(keep_chains, delete_chains, linked_chain, outlines)
            continue

        elif line.startswith("HETATM"):
            current_chain = line[21]
            conditional_keep_lines(keep_chains, delete_chains, current_chain, outlines)
            continue

        elif line.startswith("TER"):
            current_chain = line[21]
            conditional_keep_lines(keep_chains, delete_chains, current_chain, outlines)
            continue
        elif line.startswith("SSBOND"):
            current_chain = line[15]
            linked_chain = line[29]

            conditional_keep_lines(keep_chains, delete_chains, current_chain, outlines)
            conditional_keep_lines(keep_chains, delete_chains, linked_chain, outlines)
            continue

        else:
            outlines.append(line)

    INFILE.close()
    OUTFILE = open(get_decoy_name(options.s)+options.suffix+".pdb", 'w')
    OUTFILE.writelines(outlines)
    print("Done")




