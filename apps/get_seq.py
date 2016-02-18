#!/usr/bin/env python

#from sequence import fasta
from basic import biopython_util
from basic import path

from collections import defaultdict
import re
import os
import sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Uses Biopython to print sequence information.  Example:\n"
                                                 "get_seq.py --pdb 2j88_A.pdb --format fasta --outpath test.txt")

    parser.add_argument("--pdb", "-s",
                        help = "Input PDB path")

    parser.add_argument("--pdblist", "-l",
                        help = "Input PDB List")

    parser.add_argument("--pdblist_input_dir", "-i",
                        help = "Input directory if needed for PDB list")

    parser.add_argument("--chain", "-c",
                        help = "A specific chain to output",
                        default = "")

    parser.add_argument("--format", "-f",
                        help = "The output format requried.",
                        default = "fasta",
                        choices = ["basic", "fasta", "ab_order"])

    parser.add_argument("--outpath", "-o",
                        help = "Output path.  If none is specified it will write to screen.")

    parser.add_argument("--prefix", "-t",
                        default = "",
                        help = "Tag to add before chain")

    parser.add_argument("--region",
                        help = "specify a particular region, start:end:chain")

    options = parser.parse_args()

    #if not options.pdb:
    #    options.pdb = sys.argv[1]




    sequences = defaultdict()
    ordered_ids = []

    pdbs = []
    if options.pdb:
        pdbs.append(options.pdb)

    if options.pdblist:
        INFILE = open(options.pdblist, 'r')
        for line in INFILE:
            line = line.strip()
            if not line or re.search("PDBLIST", line) or line[0]=='#':
                continue

            if options.pdblist_input_dir:
                pdb_path = options.pdblist_input_dir+"/"+line
            else:
                pdb_path = line
            pdbs.append(pdb_path)
        INFILE.close()

    for pdb in pdbs:
        #print "#Reading "+pdb
        biostructure = biopython_util.get_biopython_structure(pdb)
        if options.chain:
            seq = biopython_util.get_seq_from_biostructure(biostructure, options.chain)
            sequences[path.get_decoy_name(pdb)+"_"+options.chain] = seq
            ordered_ids.append(path.get_decoy_name(pdb)+"_"+options.chain)
        else:
            for biochain in biostructure[0]:
                if biopython_util.get_chain_length(biochain) == 0:
                    continue
                seq = biopython_util.get_seq_from_biochain(biochain)
                #b = os.path.basename(pdb).replace('.pdb.gz', '')
                sequences[path.get_decoy_name(pdb)+"_"+biochain.id] = seq
                ordered_ids.append(path.get_decoy_name(pdb)+"_"+biochain.id)

    outlines = []
    i = 1


    for name_chain in ordered_ids:

        if options.format == "basic":
            outlines.append( options.prefix+name_chain+" : "+sequences[name_chain]+"\n")
        elif options.format == "fasta":
            outlines.append("> "+options.prefix+name_chain)
            outlines.append(sequences[name_chain]+"\n")
        elif options.format == "ab_order":
            name = "IgG_"+name_chain+"_pFUSE"
            outlines.append(str(i)+". "+name+ " "+sequences[name_chain]+"\n")

        i+=1

    if options.outpath:
        OUT = open(options.outpath, "w")
        for line in outlines:
            OUT.write(line+"\n")
        OUT.close()
    else:
        for line in outlines:
            print(line)


