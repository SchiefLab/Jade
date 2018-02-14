#!/usr/bin/env python

import os
import sys
import re

from argparse import ArgumentParser

#This script converts an IMGT output file (5_AA-seqs.csv) to a FASTA.  All Framework and CDRs are concatonated.  * is skipped.
# The FASTA file can then be used by PyIgClassify.


#INPUTS:
# 1) IMGT file path
# 2) FASTA Out file path


#IMGT HEADER:
#Sequence number	Sequence ID	Functionality	V-GENE and allele	J-GENE and allele	D-GENE and allele	V-D-J-REGION	V-J-REGION	V-REGION	FR1-IMGT	CDR1-IMGT	FR2-IMGT	CDR2-IMGT	FR3-IMGT	CDR3-IMGT	JUNCTION	J-REGION	FR4-IMGT

def get_parser():
    parser = ArgumentParser(description="This script converts an IMGT output file (5_AA-seqs.csv) to a FASTA.  All Framework and CDRs are concatonated.  * is skipped.\n"
                            "The FASTA file can then be used by PyIgClassify.")

    parser.add_argument("--inpath", "-i",
                        help = "Input IMGT file path",
                        required = True)

    parser.add_argument("--outpath", "-o",
                        help = "Output Fasta outfile path.",
                        required = True)

    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()

    inpath = options.inpath
    outpath = options.outpath


    INFILE = open(inpath, 'r')
    INFILE.readline() ; #Read the first uncommented line.

    FASTA = open(outpath, 'w')

    i = 1

    for line in INFILE:
        line = line.strip()
        if line.startswith('#') or not line:continue
        lineSP = line.split()

        tag = lineSP[2]
        if tag != "productive": continue

        seq_id = lineSP[1]

        full_sequence = "".join(lineSP[-9:])
        full_sequence = full_sequence.replace('*', '')

        print repr(i)

        FASTA.write(">"+repr(i)+" "+seq_id+"\n")
        FASTA.write(full_sequence+"\n\n")
        i+=1

    INFILE.close()
    FASTA.close()
