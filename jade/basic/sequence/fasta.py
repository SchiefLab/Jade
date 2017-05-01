#Author Jared Adolf-Bryfogle

import random
from collections import defaultdict

from weblogolib import *

from .ClustalRunner import ClustalRunner
from jade.basic.structure.util import *


########  Fasta Input  ######

def get_label_from_fasta(fasta_path):
    """
    Gets the first chainID found - Should be a single chain fasta file.
    """
    FILE = open(fasta_path, 'r')
    for line in FILE:
        if line[0]==">":
            line = line.strip()
            label = line[1:]
            break
    FILE.close()
    return label

def get_sequence_from_fasta(fasta_path, label):
    FILE = open(fasta_path, 'r')
    lines = FILE.readlines()

    sequence = ""
    start = 0
    end = 0
    for i in range(0, len(lines)):
        lines[i] = lines[i].strip()
        if not lines[i]:continue
        if lines[i]==">"+label:
            start = i+1
            break

    for i in range(start, len(lines)):
        if not lines[i]:break
        elif lines[i][0]=="#":break
        elif lines[i][0]==">":break
        else:
            sequence = sequence+lines[i].strip()
    return sequence

def read_header_data_from_fasta(fasta_path):
    """
    Reads > from fasta (PDBAA) and returns a defaultdict of
    pdb_chain: [method, residues, resolution, R factor]
    """
    result = defaultdict()
    FILE = open(fasta_path, 'r')
    for line in FILE:
        if line[0]==">":
            line = line.strip()
            line = line[1:]
            lineSP = line.split()

            lis = [lineSP[2], lineSP[1], lineSP[3], lineSP[4]]
            result[lineSP[0]] = lis

    return result


########  Fasta Output #########

def write_fasta(sequence, label, HANDLE):
    """
    Writes a fasta with a sequence, chain, and open FILE handle.
    FULL Sequence on one line seems to be fine with HMMER3.
    """
    HANDLE.write(">"+label+"\n")
    HANDLE.write(sequence + "\n")

def fasta_from_sequences(sequences, outdir, outname):
    """
    Output a general fasta, with tag being 1_outname etc.
    Use write_fasta for more control.
    Returns path to Fasta File written
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outpath = outdir+"/"+outname+".fasta"
    OUTFILE = open(outpath, 'w')
    i= 1
    for sequence in sequences:
        name = repr(i)+"_"+outname
        write_fasta(sequence, name, OUTFILE)
        i+=1
    OUTFILE.close()
    return outpath

def output_fasta_from_pdbs_biopython(path_header_dict, out_path, native_path = None, native_label = "native", is_camelid = False):
    """
    Used only for L and H chains!  Concatonates the L and H in order if present, otherwise assumes camelid at H.
    """
    if not os.path.exists(os.path.dirname(out_path)):
        os.mkdir(os.path.dirname(out_path))

    OUTFILE = open(out_path, 'w')
    if native_path:
        structure =  get_biopython_structure(native_path)
        if is_camelid or len(structure[0]) == 1:
            seq = get_seq_from_biostructure(structure, 'H')
        else:
            seq = get_seq_from_biostructure(structure, "L")+get_seq_from_biostructure(structure, "H")
        write_fasta(seq, native_label, OUTFILE)

    for path in sorted(path_header_dict, key=path_header_dict.get):
        if not os.path.exists(path):
            print "File not found: "+path


        seq = ""
        print "Adding: "+os.path.basename(path)
        structure = get_biopython_structure(path)

        if is_camelid or len(structure[0]) == 1:
            seq = get_seq_from_biostructure(structure, 'H')
        else:
            seq = get_seq_from_biostructure(structure, "L")+get_seq_from_biostructure(structure, "H")

        write_fasta(seq, path_header_dict[path], OUTFILE)
    OUTFILE.close()

def fasta_from_pose(pose, fasta_label, outname, outdir):
    """
    Creates a fasta from the pose.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    seq = pose.sequence()
    FILE = open(outdir+"/"+outname+".fasta", 'w')
    write_fasta(seq, fasta_label, FILE)
    FILE.close()

    print "FASTA from pose written"

def chain_fasta_files_from_pose(pose, prefix, outdir):
    """
    Creates fasta for each chain in the pose.  Returns a list of paths for each fasta.
    """
    from rosetta.core.pose import get_chain_from_chain_id

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    num_chains = pose.conformation().num_chains()
    fasta_files = []
    for i in range(1, num_chains+1):
        seq = pose.chain_sequence(i)
        chain = get_chain_from_chain_id(i, pose)
        outname = prefix+"_"+chain+".fasta"
        FILE = open(outdir+"/"+outname, 'w')
        fasta_files.append(outdir+"/"+outname)
        write_fasta(seq, chain, FILE)
        FILE.close()

    print "FASTAs for chains written"
    return fasta_files

def chain_fasta_files_from_biostructure(structure, prefix, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    model = structure[0]
    num_chains = len(model)
    fasta_files = []
    for biochain in model:
        if get_chain_length(biochain) == 0:
            continue
        seq = get_biochain_sequence(biochain)
        chain = biochain.id
        outname = prefix+"_"+chain+".fasta"
        FILE = open(outdir+"/"+outname, 'w')
        fasta_files.append(outdir+"/"+outname)
        write_fasta(seq, chain, FILE)
        FILE.close()

    print "FASTAs for chains written"
    return fasta_files

def chain_fasta_from_biostructure(structure, outname, outdir):
    """
    Creates a single fasta from biopython structure, split by individual chains.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    model = structure[0]
    num_chains = len(model)
    FILE = open(outdir+"/"+outname+".fasta", 'w')
    for chain in model:
        seq = get_biochain_sequence(chain)
        write_fasta(seq, chain.id, FILE)
        FILE.write("\n")

    FILE.close()

    print "FASTA for chains written"

def chain_fasta_from_pose(pose, outname, outdir):
    """
    Creates a single fasta from pose, split by individual chains.
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    num_chains = pose.conformation().num_chains()
    FILE = open(outdir+"/"+outname+".fasta", 'w')
    for i in range(1, num_chains+1):
        seq = pose.chain_sequence(i)
        chain = get_chain_from_chain_id(i)


        write_fasta(seq, chain, FILE)
        FILE.write("\n")

    FILE.close()

    print "FASTA for chains written"

def split_fasta_from_fasta(fasta_path, prefix, outdir):
    """
    If we have a multiple fasta sequence, we split it into idividual files.  Makes analysis easier.
    Returns list of paths for each fasta
    """
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    INFILE = open(fasta_path, 'r')
    new_fasta = False
    chains = dict()
    found = False
    id = random.getrandbits(128)
    chain_num = 0
    n = ""
    for line in INFILE:
        if line[0]=='>':
            chain_num +=1
            #chain = line[1:]
            #chain = chain.strip().replace(' ', '_')
            #chains[chain] = ""
            #found = True
            #print chain
            found = True
            n = "chain_"+repr(id)+"_"+repr(chain_num)

            head = "_".join(line.split()) +"\n"
            chains[n] = head
            continue
        elif line[0] == '#':
            continue
        if found:
            chains[n] = chains[n]+line

    fasta_files = []
    for chain in chains:
        outname = prefix+"_"+chain+".fasta"
        #print outname
        OUTFILE = open(outdir+"/"+outname, 'w')
        OUTFILE.write(chains[chain])
        OUTFILE.close()
        fasta_files.append(outdir+"/"+outname)

    print "FASTA file split."
    INFILE.close()

    return fasta_files

def get_biochain_sequence(bio_chain):
    seq = ""
    d = RestypeDefinitions()

    for res in bio_chain:
        if res.id[0]==' ':
            aa = d.get_one_letter_from_three(res.resname)
            if not aa:
                print "Skipping non-canonical resname: "+res.resname
                print "This could pose a problem!"
                continue
            seq = seq+aa
    return seq

########  Etc ###########
def output_weblogo_for_sequences(sequences, outdir, outname, tag = "Dunbrack Lab - Antibody Database Team"):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fasta_outdir = outdir+"/fastas"

    fasta_path = fasta_from_sequences(sequences, fasta_outdir, outname)
    clustal_runner = ClustalRunner(fasta_path)
    clustal_runner.set_extra_options("--force")
    clustal_runner.output_alignment(outdir+"/clustal", outname+".clu")

    clustal_path = outdir+"/clustal/"+outname+".clu"

    weblogo_dir = outdir+"/weblogos"

    output_weblogo(clustal_path, weblogo_dir, outname, tag)

def output_weblogo(alignment_path, outdir, outname, tag = "Dunbrack Lab - Antibody Database Team"):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    MSA_IN = open(alignment_path)
    seqs = read_seq_data(MSA_IN)
    data = LogoData.from_seqs(seqs)
    options = LogoOptions()

    #Weblogo Options
    options.logo_title = outname
    #options.fineprint = datetime.datetime.now().strftime('%b-%d-%G')
    options.creator_text = tag
    options.show_fineprint = False
    options.resolution = 900
    options.number_interval = 1
    options.scale_width = True
    options.unit_name='probability'
    options.color_scheme=std_color_schemes["charge"]
    #options.color_scheme=std_color_schemes["hydrophobicity"]

    format = LogoFormat(data, options)
    LOGO_OUT = open(outdir+"/"+outname+"_weblogo.png", 'w')
    png_formatter( data, format, LOGO_OUT)

    MSA_IN.close()
    LOGO_OUT.close()


