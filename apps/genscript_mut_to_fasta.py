#/usr/bin/env python

#Jared Adolf-Bryfogle
#Written for Yumi.

#Format we read from genscript is:

#Mutagenesis Instruction:
#Variant name: DU422_SOSIP_MD39_N276D_mC

#Variant sequence: ACCGGTAACCTGGACCTGTGGGTGACCGTGTACTATGGCGTGCCAGTGTGGAAGGAGGCCAAGACCACAC

import sys,os
from argparse import ArgumentParser
from collections import defaultdict

import basic.sequence.fasta as fasta

def get_DNA(mutant_position, all_lines):
    """
    Get the DNA lines as a list knowning the mutant position.

    :param all_lines = [str]
    :param mutant_position: int
    :rtype: [str]
    """

    dna_lines = []
    for i in range(mutant_position + 2, len(all_lines)+1):
        line = all_lines[ i ]
        if not line:
            break

        line = line.replace("Variant sequence: ", "")
        dna_lines.append(line)
    return dna_lines


if __name__ == "__main__":
    parser = ArgumentParser("This script outputs fasta files from a genscript mutagenesis format."
                            "Ex: python genscript_mut_to_fasta.py MutagenesisFormatUadf")

    parser.add_argument("infile", help = "The mutagenesis format file.", required = True)

    options = parser.parse_args()

    lines = []
    mutant_positions = []
    INFILE = open(options.infile, 'r')
    i = 0
    for line in INFILE:
        line = line.strip()
        #if not line or line.startswith("#"): continue
        lines.append(line)


        lineSP = line.split(':')
        lineSP = [f.strip() for f in lineSP]
        if lineSP[0] == "Variant name":
            mutant_positions.append(i)

        i+=1

    INFILE.close()

    OUTFILE = open(os.path.basename(options.infile)+".fasta", "w")
    for position in mutant_positions:
        mutant_name_line = lines[position]
        mutantSP = mutant_name_line.split(':')
        mutantSP = [f.strip() for f in mutantSP]
        variant_name = mutantSP[1]
        dna_lines = get_DNA(position, lines)
        fasta.write_fasta("\n".join(dna_lines), variant_name, OUTFILE)
        OUTFILE.write("\n\n")

    OUTFILE.close()
    print "complete"



