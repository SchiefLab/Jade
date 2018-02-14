#!/usr/bin/env python

#Jared Adolf-Bryfogle
#Written for Yumi.

#Format we read from genscript is:

#Mutagenesis Instruction:
#Variant name: DU422_SOSIP_MD39_N276D_mC

#Variant sequence: ACCGGTAACCTGGACCTGTGGGTGACCGTGTACTATGGCGTGCCAGTGTGGAAGGAGGCCAAGACCACAC

import sys, os
from argparse import ArgumentParser


def get_DNA_seq(mutant_position, all_lines):
    """
    Get the DNA lines as a list knowning the mutant position.

    :param all_lines = [str]
    :param mutant_position: int
    :rtype: [str]
    """

    dna_lines = []
    for i in range(mutant_position + 2, len(all_lines)):
        line = all_lines[ i ]
        if not line:
            break

        line = line.replace("Variant sequence: ", "")
        dna_lines.append(line)
    return dna_lines


def write_fasta(sequence, label, HANDLE):
    """
    Writes a fasta with a sequence, chain, and open FILE handle.
    FULL Sequence on one line seems to be fine with HMMER3.
    """
    HANDLE.write(">"+label+"\n")
    HANDLE.write(sequence + "\n")

def get_parser():
    parser = ArgumentParser(description="This script outputs fasta files from a genscript format.  Pass the --format option to control which genscript format as input"
                            "  ~~~ Ex: python genscript_mut_to_fasta.py --format mutagenesis MutagenesisFormatU68  ~~~")

    parser.add_argument("infile", help = "The mutagenesis format file.")

    parser.add_argument("--format", help = "The genscript file format",
                        choices = ["mutagenesis", "GeneSynth"],
                        required = True)

    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()

    INFILE = open(options.infile, 'r')
    OUTFILE = open(os.path.basename(options.infile)+".fasta", "w")


    lines = []
    mutant_positions = []

    i = 0
    for line in INFILE:
        line = line.strip()
        #if not line or line.startswith("#"): continue
        lines.append(line)


        lineSP = line.split(':')
        lineSP = [f.strip() for f in lineSP]
        if lineSP[0] == "Variant name" or lineSP[0] == "Gene Name":
            mutant_positions.append(i)

        i+=1



    if options.format == "mutagenesis":
        for position in mutant_positions:
            mutant_name_line = lines[position]
            mutantSP = [f.strip() for f in mutant_name_line.split(':')]
            variant_name = mutantSP[1]
            dna_lines = get_DNA_seq(position, lines)
            write_fasta("\n".join(dna_lines), variant_name, OUTFILE)
            OUTFILE.write("\n\n")

    elif options.format == "GeneSynth":
        for position in mutant_positions:
            mutant_name_line = lines[position]
            mutantSP = [f.strip() for f in mutant_name_line.split(':')]
            mutantSP2 = [f.strip() for f in mutantSP[1].split(',')]
            variant_name = mutantSP2[0]
            dna_lines = get_DNA_seq(position, lines)
            write_fasta("\n".join(dna_lines), variant_name, OUTFILE)
            OUTFILE.write("\n\n")
    else:
        sys.exit("Unknown format: "+options.format)



    INFILE.close()
    OUTFILE.close()
    print "complete"



