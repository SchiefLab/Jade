import os, sys, re
from argparse import ArgumentParser
from collections import defaultdict, namedtuple

from jade.nnk.NNKIndex import TemplateAbNNKIndex, TestAbNNKIndex
from jade.nnk.NNKEnrichments import NNKEnrichments
from jade.nnk.NNKScoreFunction import AntibodyNNKScoreFunction
from jade.nnk.NNKAbMaturation import GetNNKData
from jade.basic.sequence import fasta

NNKSortIndex = namedtuple('NNKSortIndex', 'antibody antigen sort')
_sorts_ = ['S1', 'S2', 'S3']

if __name__ == "__main__":

    parser = ArgumentParser("This app scores putitative broadly neutralizing antibodies based on pre-processed NNK sorting data")

    parser.add_argument('--nnk_dir', '-d', default = "sort_data",         help = "Full or relative path to pre-processed NNK data")
    parser.add_argument('--out_dir', '-o', default = 'scores',            help = "Directory for all output files")
    parser.add_argument('--abs', '-a', default = ['CHA31'], nargs = '*',  help = "Antibodies to use for the scoring function. ")
    parser.add_argument('--gl_antigen', '-g', default = 'GT81',           help = 'The antigen used for Germline-scoring')
    parser.add_argument('--mat_antigens', '-m',                           help = "A list of mature antigens to use.  We will report both score as the means and them split.")
    parser.add_argument('--input', '-s',                                  help = "A fasta or PDB file")
    parser.add_argument('--input_seq_list', '-l',                         help = "A file with the ID in column 1 and the sequence in column 2. Space or tab delimited.  Header must have comments. Sequence can have '-' in them. ")

    options = parser.parse_args()



    index = TemplateAbNNKIndex()
    if not os.path.exists(options.nnk_dir):
        sys.exit('NNK Directory not found!  Please pass a proper directory!')

    if not os.path.exists(options.out_dir):
        os.mkdir(options.out_dir)

    #Create any FASTA files if we need them.

    base_name = ""
    if options.input_seq_list:
        base_name = ".".join(os.path.basename(options.input_seq_list).split('.')[0:-1])
        options.input = options.out_dir+'/'+base_name+".fasta"
        OUTFILE = open(options.input, 'w')
        INFILE = open(options.input_seq_list, 'r')
        for line in INFILE:
            line = line.strip()
            if not line or line.startswith('#'):continue
            lineSP = line.split()
            id = lineSP[0]
            seq = lineSP[1].strip('-')
            seq = seq.replace('-', '')

            fasta.write_fasta(seq, id, OUTFILE)
        OUTFILE.close()
        INFILE.close()

    else:
        base_name = ".".join((os.path.basename(options.input)).split('.')[0:-1])

    #Setup Enrichments and Scoring Functions
    combined_mature_antigens = []
    antigen_dict = defaultdict()

    #If we are using multiple antibodies, make sure we only use the subset antigens for combined antibody data.
    # We will also need to combine all other data as well...
    counts = defaultdict(int)
    for antibody in options.abs:
        mat_meta = GetNNKData(options.data_dir, 'mat'+antibody)
        if options.mat_antigens:
            antigens = list(set(mat_meta.antigens).intersection(options.mat_antigens))
        else:
            #Use all antigens
            antigens = mat_meta.antigens

        antigen_dict[antibody] = antigens
        for antigen in antigens:
            if not counts.has_key(antigen): counts[antigen] = 0
            else: counts[antigen]+=1

    for key in counts:
        if counts[key] == len(options.abs):
            combined_mature_antigens.append(key)

    #Create Scoring Functions
    gl_scoring_functions = defaultdict() #Key Tuple of Antibody, Antigen, Sort; Value is the scoring function.
    mat_scoring_functions = defaultdict()

    #Create Germline Scoring Functions
    for sort in _sorts_:
        es = []
        for antibody in options.abs:
            enrichments = NNKEnrichments(options.nnk_dir, 'VRC01', 'gl'+antibody, options.gl_antigen, sort)
            score = AntibodyNNKScoreFunction([enrichments], index)
            sort_index = NNKSortIndex('gl'+antibody, options.gl_antigen, sort)

            gl_scoring_functions[sort_index] = score
            es.append(score)
        if len(options.abs) > 1:
            combined_score = AntibodyNNKScoreFunction([es], index)
            sort_index = NNKSortIndex('glCombined', options.gl_antigen, sort)
            gl_scoring_functions[sort_index] = combined_score



    antibody_dict = defaultdict()
    for key, value in antigen_dict.iteritems():
        antigen_dict[value] = key

    #Create Mature Scoring Functions
    for sort in _sorts_:
        for antibody in options.abs:
            es = []
            for antigen in antigen_dict[antibody]:

                sort_index = NNKSortIndex('mat'+antibody, antigen, sort)
                enrichments = NNKEnrichments(options.nnk_dir, 'VRC01', 'mat'+antibody, antigen, sort)

                score = AntibodyNNKScoreFunction([enrichments], index)
                mat_scoring_functions[sort_index] = score

                if antigen in combined_mature_antigens:
                    es.append(enrichments)

            if len(antigen_dict[antibody]) > 1:
                combined_score = AntibodyNNKScoreFunction(es, index)
                sort_index = NNKSortIndex(antibody, 'AntigensCombined', sort)
                mat_scoring_functions[sort_index] = combined_score

        if len(options.abs) > 1:
            all = []
            for antigen in combined_mature_antigens:
                es []
                for antibody in options.ab:
                    enrichments = NNKEnrichments(options.nnk_dir, 'VRC01', 'mat'+antibody, antigen, sort)
                    es.append(enrichments)

                score = AntibodyNNKScoreFunction(es, index)
                sort_index = NNKSortIndex('matCombined', antigen, sort)
                mat_scoring_functions[sort_index] = score
            score = AntibodyNNKScoreFunction(all, index)
            sort_index = NNKSortIndex('matCombined', 'AntigensUnion', sort)
            mat_scoring_functions[sort_index] = score



    #For now, this is default HEAVY as we have none with LIGHT chains
    SCOREOUT = open(options.out_dir+'/'+base_name+'_scores.tsv')
    SCOREOUT.write("#\tid\tantibody\tantigen\tsort\n")

    test_data = TestAbNNKIndex(options.input)
    for classified_ab in test_data.classified_list():
        id = classified_ab.ig_chain.id
        print id
        for score_index, score in gl_scoring_functions.iteritems():

            s = score.relative_score_classified_ab(classified_ab)
            out = "\t".join([id, score_index.antibody, score_index.antigen, score_index.sort, "%.3f" % s])
            print out
            SCOREOUT.write(out+"\n")

    SCOREOUT.close()
    print "Done"

