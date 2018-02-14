#!/usr/bin/env python
# coding: utf-8

#from sequence import fasta
import argparse
import re
from collections import defaultdict

from jade.basic import path
from jade.basic.structure import util
from jade.basic import general
from jade.basic.structure.Structure import AntibodyStructure
from jade.basic.structure.BioPose import BioPose

begin_schief_order = """
Scripps PO C426550C

all genes are:
synthesis
codon optimization for mammalian (293) cells
various vectors are used
"""

end_schief_ab_order = """

-------- order nothing below this line --------------

###################################################################################################################
Notes for Bill:
GCACTTGTCACGAATTCG--[Antibody-Heavy-FV]--GCTAGCACCAAGGGCCCATC translates to ALVTNS--[Antibody-Heavy-FV]--ASTKGP
GCACTTGTCACGAATTCA--[Antibody-Kappa-Fv]--CGTACGGTGGCTGCACCA translates to ALVTNS--[Antibody-Kappa-Fv]--RTVAAP
(*** so kappa gene insert should end in K ***)
GCACTTGTCACGAATTCG--[Antibody-Lambda-FV]—-CTAGGTCAGCCCAAGGCTGCCCC translates to ALVTNS--[Antibody-Lambda-FV]—-LGQPKAA
(*** so lambda gene insert should normally end in V and should NOT include the last L ***)
###################################################################################################################
"""

igg_types = ["IgG_order", "IgG_order_lambda", "IgG_order_kappa", "IgG_order_heavy"]
format_types = ["basic", "fasta", "general_order"]
format_types.extend(igg_types)

def get_parser():



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

    parser.add_argument("--cdr",
                        help = "Pass a specific CDR to output alignments of.",
                        default = "")

    parser.add_argument("--format",
                        help = "The output format requried.",
                        default = "fasta",
                        choices = format_types)

    parser.add_argument("--outpath", "-o",
                        help = "Output path.  If none is specified it will write to screen.")

    parser.add_argument("--prefix", "-t",
                        default = "",
                        help = "Tag to add before chain")

    parser.add_argument("--region",
                        help = "specify a particular region, start:end:chain")

    #parser.add_argument("--strip_n_term",
    #                    help = "Strip this sequence off the N-term of resulting sequences. (Useful for antibodies")

    parser.add_argument("--strip_c_term",
                        help = "Strip this sequence off the C-term of resulting sequences. (Useful for antibodies")

    parser.add_argument("--pad_c_term",
                        help = "Pad this sequence with some C-term (Useful for antibodies")

    parser.add_argument("--output_original_seq",
                        default = False,
                        action = "store_true",
                        help = "Output the original sequence and the striped seqeunce if stripped.  Default FALSE.")


    return parser

if __name__ == "__main__":



    parser = get_parser()
    options = parser.parse_args()

    if options.format in igg_types:
        igg_type_format = True
    else:
        igg_type_format = False


    options.format = options.format.lower()

    #if not options.pdb:
    #    options.pdb = sys.argv[1]

    if (options.format == "igg_order_lambda" or options.format == "igg_order_kappa") and not options.chain:
        options.chain = "L"

    if options.format == "igg_order_heavy" and not options.chain:
        options.chain = "H"


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
        #print "Reading "+pdb
        biopose = BioPose(pdb)
        biostructure = biopose.structure()
        ab_info = AntibodyStructure()
        if options.cdr:
            seq = ab_info.get_cdr_seq(biopose, options.cdr.upper())
            sequences[path.get_decoy_name(pdb)+"_"+options.cdr] = seq
            ordered_ids.append(path.get_decoy_name(pdb)+"_"+options.cdr)
        elif options.chain:
            seq = biopose.get_sequence(options.chain)
            sequences[path.get_decoy_name(pdb)+"_"+options.chain] = seq
            ordered_ids.append(path.get_decoy_name(pdb)+"_"+options.chain)
        else:
            for biochain in biostructure[0]:
                if util.get_chain_length(biochain) == 0:
                    continue
                seq = util.get_seq_from_biochain(biochain)
                #b = os.path.basename(pdb).replace('.pdb.gz', '')
                sequences[path.get_decoy_name(pdb)+"_"+biochain.id] = seq
                ordered_ids.append(path.get_decoy_name(pdb)+"_"+biochain.id)

    outlines = []
    i = 1

    if options.format == "general_order" or igg_type_format:
        outlines.append(begin_schief_order)

        print "DOUBLE CHECK FOR MISSING DENSITY IN THE STRUCTURE!  MAKE SURE SEQUENCES MATCH!"

    if options.format == "igg_order_heavy":
        outlines.append("###################################################################################################################################################")
        outlines.append("Use pFUSEss-CHIg-hG1 (human heavy) vector and clone between the two flanking regions: "
                        "GCACTTGTCACGAATTCG--[Antibody-Heavy-FV]--GCTAGCACCAAGGGCCCATC")
        outlines.append("###################################################################################################################################################\n")

        print "Sequence should flank: ANY - WGqgtlVtvss\n"
        print ">Heavy example:\nEVQLVESGGGLVKPGGSLRLSCSASGFDFDNAAMSWVRQPPGKGLEWVGTITGPGEGWSVDYAAPVEGRFTISRLNSINFLYLEMNNLRMEDSGLYFCARGEWEFRNGETYSALTYWGRGTLVTVSS"

    elif options.format == "igg_order_lambda":
        outlines.append("#################################################################################################################################################")
        outlines.append("Use pFUSEss-CLIg-hL2 (human lambda) vector and clone between the two flanking regions: "
                        "GCACTTGTCACGAATTCG [Antibody-Lambda-Fv] CTAGGTCAGCCCAAGGCT")
        outlines.append("#################################################################################################################################################\n")

        print "Sequence should flank: ANY - gsgtqvTV\n"
        print ">Lambda example:\nSYELTQETGVSVALGRTVTITCRGDSLRYHYASWYQKKPGQAPILLFYGKNNRPSGVPDRFSGSASGNRASLTISGAQAEDDAEYYCMSAAKPGSWTRTFGGGTKLTVL"

    elif options.format == "igg_order_kappa":
        outlines.append("##########################################################################################################################################################")
        outlines.append("Use pFUSEss-CLIg-hk (human kappa) vector and clone between the two flanking regions: "
                        "GCACTTGTCACGAATTCA--[Antibody-Kappa-Fv]--CGTACGGTGGCTGCACCA")
        outlines.append("##########################################################################################################################################################\n")

        print "Sequence should flank: ANY - FGGGtkveik\n"
        print ">Kappa example:\nDIQMTQSPASLSASVGETVTITCRASENIYSYLTWYQQKQGKSPQLLVYNAKTLAEGVPSRFSGSGSGTQFSLKISSLQPEDFGNYYCQHHYGTRTFGGGTRLEIK"

    print "\n"
    for name_chain in ordered_ids:
        original_seq = sequences[name_chain]
        seq = sequences[name_chain]

        if options.strip_c_term:
            seq = general.strip_right(seq, options.strip_c_term)

        if options.pad_c_term:
            seq = seq+options.pad_c_term

        if options.format == "basic":
            if options.output_original_seq:
                outlines.append(options.prefix+name_chain+" : "+original_seq)
            outlines.append( options.prefix+name_chain+" : "+seq+"\n")

        elif options.format == "fasta":
            outlines.append("> "+options.prefix+name_chain)
            if options.output_original_seq:
                outlines.append(original_seq+"\n")

            outlines.append(seq+"\n")

        elif igg_type_format:
            name = "IgG_"+name_chain+"_pFUSE"
            if options.output_original_seq:
                outlines.append(str(i)+". "+name+ " "+original_seq+"\n")
            outlines.append(str(i)+". "+name+ " "+seq+"\n")

        elif options.format == "general_order":
            outlines.append(str(i)+". "+"vec"+ " "+seq+"\n")

        i+=1

    if igg_type_format:
        outlines.append(end_schief_ab_order)

    if options.outpath:
        OUT = open(options.outpath, "w")
        for line in outlines:
            OUT.write(line+"\n")
        OUT.close()
    else:
        for line in outlines:
            print(line)


    print "Done."
