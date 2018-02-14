#!/usr/bin/env python

import os, sys, re
from argparse import ArgumentParser
from collections import defaultdict, namedtuple

from jade.nnk.NNKIndex import TemplateAbNNKIndex, TestAbNNKIndex, NumericalIndexOptions
from jade.nnk.NNKEnrichments import NNKEnrichments
from jade.nnk.NNKScoreFunction import AntibodyNNKScoreFunction, MetaAntibodyNNKScoreFunction
from jade.nnk.NNKAbMaturation import GetNNKData
from jade.basic.sequence import fasta
from jade.basic import general

NNKSortIndex = namedtuple('NNKSortIndex', 'antibody antigen sort')
_sorts_ = ['S1', 'S2', 'S3']

def get_parser():
    parser = ArgumentParser("This app scores putitative broadly neutralizing antibodies based on pre-processed NNK sorting data")

    parser.add_argument('--nnk_dir', '-d', default = "sort_data",         help = "Full or relative path to pre-processed NNK data")
    parser.add_argument('--out_dir', '-o', default = 'scores',            help = "Directory for all output files")
    parser.add_argument('--gl_abs', '-a', default = ['CHA31'], nargs = '*',  help = "Antibodies to use for the germline scoring function. ")
    parser.add_argument('--mat_abs', default = ['CHA31'], nargs = '*', help = "Antibodies to use for the mature scoring functions.")
    parser.add_argument('--gl_antigen', '-g', default = 'GT81',           help = 'The antigen used for Germline-scoring')
    parser.add_argument('--mat_antigens', '-m',                           help = "A list of mature antigens to use.  We will report both score as the means and them split.")
    parser.add_argument('--skip_mat_antigens', '-k', default = ['GT8-1'], nargs = '*',
                                                                          help = "A list of mature antigens to skip. ")

    parser.add_argument('--input', '-s',                                  help = "A fasta or PDB file")
    parser.add_argument('--input_seq_list', '-l',                         help = "A file with the ID in column 1 and the sequence in column 2. Space or tab delimited.  Header must have comments. Sequence can have '-' in them. ")
    parser.add_argument('--zeros', '-z', default = -2,                    help = "The number to use for any zero top game enrichment.  Too high, and you may penalize these too much.  Too low, and they won't hold their weight. "
                                                                                 "-2 represents an erichment of .11; -2.5 at .08, and -3.0 at about .05.  -2 is about the last standard deviation in the erichment curve that resembles a gaussian.")
    parser.add_argument('--output_individual', '-c', default = False, action = 'store_true',
                                                                          help = "Should we output individual scores or only combined?")

    parser.add_argument('--use_factors', '-f', default = False, action = 'store_true',
                                                                          help = "Should we score the data as 'factor'?")
    parser.add_argument('--use_global_matches', default = False, action = 'store_true',
                                                                          help = "use global matches based on the input NNK sequences")
    parser.add_argument('--use_group_matches', default = False, action = 'store_true',
                                                                          help = "use group matches based on the input NNK sequences")

    #parser.add_argument('--plot', '-p', default=False, action = 'store_true',
    #                                                                      help = "Should we output plots of the data?")

    parser.add_argument('--award_max_only', default = False, action = 'store_true',
                                                                          help = "Award scores for only maximum matches.")

    parser.add_argument('--award_max_and_conserved_only', default=False, action = 'store_true',
                                                                          help = "Award a score if it is maximum or conserved.")

    parser.add_argument('--exp', '-e', required = True, help = "Name of the experiment.  Used for concatonating different experiments.")

    parser.add_argument('--additive_combine','-t', default = False, action = 'store_true',
                                                                          help = "Combine antigens/antibodies using addition instead of means.  ")

    parser.add_argument('--overwrite', default = False, action = "store_true",
                                                                          help = "Overwrite any output")
    parser.add_argument('--prefix', default = "", help = "Any prefix for output data.")

    parser.add_argument('--include_cdrs', default = ['L1','L2','L3','H1','H2'], nargs="*", help = "CDRs to include.  By default we skip H3")

    return parser

if __name__ == "__main__":

    parser = get_parser()
    options = parser.parse_args()


    score_function_options = defaultdict()
    score_function_options['additive_combine'] = options.additive_combine
    score_function_options['award_max_only'] = options.award_max_only
    score_function_options['award_max_and_conserved_only'] = options.award_max_and_conserved_only
    score_function_options['use_factors'] = options.use_factors
    score_function_options['use_global_matches'] = options.use_global_matches
    score_function_options['use_group_matches'] = options.use_group_matches

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
    for antibody in options.mat_abs:
        print "Getting meta data for mature "+antibody
        mat_meta = GetNNKData(options.nnk_dir, 'mat'+antibody)
        if options.mat_antigens:
            antigens = list(set(mat_meta.antigens).intersection(options.mat_antigens))
        elif options.skip_mat_antigens:
            antigens = mat_meta.antigens
            for a in options.skip_mat_antigens:
                if a in antigens:
                    antigens.remove(a)

        else:
            #Use all antigens
            antigens = mat_meta.antigens

        antigen_dict[antibody] = antigens
        for antigen in antigens:
            print antigen
            counts[antigen]+=1

    for key in counts:
        print key
        print repr(counts[key])
        if counts[key] == len(options.mat_abs):
            combined_mature_antigens.append(key)

    #### Create Scoring Functions
    gl_scoring_functions = defaultdict() #Key Tuple of Antibody, Antigen, Sort; Value is the scoring function.
    mat_scoring_functions = defaultdict()

    #Create Indexes
    indexes = defaultdict()
    skip_cdrs = []
    skip_residues_tuple = []
    skip_residues_number =[]
    regions = ['FRAME', 'CDR']

    _cdrs_ = ['L1','L2','L3','H1','H2','H3']
    for cdr in options.include_cdrs:
        if cdr not in _cdrs_:
            skip_cdrs.append(cdr)

    index_options = NumericalIndexOptions(skip_cdrs, skip_residues_tuple, skip_residues_number, regions)

    for antibody in ['gl'+ab for ab in options.gl_abs] + ['mat'+ab for ab in options.mat_abs]:
        indexes[antibody] = TemplateAbNNKIndex(antibody, options=index_options)

    #Create Germline Scoring Functions
    sort_index_order = []


    # Order groups first:
    for sort in _sorts_:
        scores = []
        for antibody in options.gl_abs:

            print "\nLoading Germline enrichment data for score function: "+antibody+' '+options.gl_antigen+" "+sort
            enrichments = NNKEnrichments(options.nnk_dir, options.zeros, 'VRC01', 'gl'+antibody, options.gl_antigen, sort)
            score = AntibodyNNKScoreFunction([enrichments], indexes['gl'+antibody], **score_function_options)
            sort_index = NNKSortIndex('gl'+antibody, options.gl_antigen, sort)

            gl_scoring_functions[sort_index] = score
            sort_index_order.append(sort_index)

            scores.append(score)
        if len(options.gl_abs) > 1:
            print "Combining Germline Antibodies"
            combined_score = MetaAntibodyNNKScoreFunction(scores, options.additive_combine)

            sort_index = NNKSortIndex('glCombined', options.gl_antigen, sort)
            sort_index_order.append(sort_index)
            gl_scoring_functions[sort_index] = combined_score



    antibody_dict = defaultdict()
    for key, value in antigen_dict.iteritems():
        antibody_dict[tuple(value)] = key

    #Create Mature Scoring Functions
    mature_scores = defaultdict()

    for sort in _sorts_:
        scores = []
        for antibody in options.mat_abs:
            sort_index_order.append(NNKSortIndex('mat'+antibody, 'AntigensCombined', sort))

            for antigen in antigen_dict[antibody]:
                print "\nLoading Mature enrichment data for score function: " + antibody + ' ' + antigen + " " + sort
                sort_index = NNKSortIndex('mat'+antibody, antigen, sort)
                enrichments = NNKEnrichments(options.nnk_dir, options.zeros, 'VRC01', 'mat'+antibody, antigen, sort)

                score = AntibodyNNKScoreFunction([enrichments], indexes['mat'+antibody], **score_function_options)
                mature_scores[(antibody, antigen)] = score
                mat_scoring_functions[sort_index] = score
                sort_index_order.append(sort_index)

                if antigen in combined_mature_antigens:
                    scores.append(score)

            if len(antigen_dict[antibody]) > 1:
                print "Combining mature antigens for sort: "+sort
                combined_score = MetaAntibodyNNKScoreFunction(scores, options.additive_combine)
                sort_index = NNKSortIndex('mat'+antibody, 'AntigensCombined', sort)
                mat_scoring_functions[sort_index] = combined_score

        if len(options.mat_abs) > 1:
            print "\nCombining enrichments"
            all = []
            for antigen in combined_mature_antigens:
                print repr(antigen)
                local_scores = []
                for antibody in options.ab:
                    #enrichments = NNKEnrichments(options.nnk_dir, options.zeros, 'VRC01', 'mat'+antibody, antigen, sort)
                    scores = mature_scores[(antibody, antigen)]
                    local_scores.append(score)
                    all.append(score)

                score = MetaAntibodyNNKScoreFunction(scores, options.additive_combine)
                sort_index = NNKSortIndex('matCombined', antigen, sort)
                mat_scoring_functions[sort_index] = score

            score = MetaAntibodyNNKScoreFunction(all, options.additive_combine)
            sort_index = NNKSortIndex('matCombined', 'AntigensCombined', sort)
            mat_scoring_functions[sort_index] = score

    #For now, this is default HEAVY as we have none with LIGHT chains
    score_outname = options.out_dir+'/'+options.prefix+'all_scores.tsv'

    how = 'w'
    if os.path.exists(score_outname) and not options.overwrite:
        how = 'a'
    SCOREOUT = open(score_outname, how)

    header = "exp\tid\tzeros\ts\tantigen\tsort\tscore\tmax_score\trelative_score\n"
    if how == 'w':
        SCOREOUT.write(header)

    print header

    #GROUPS = open(options.out_dir+'/'+'score_function_info.tsv', 'w')
    #GROUPS.write('#exp\tzeros\tantibody\tantigen\tsort\tmin\tmin_sequence\tmax\tmax_sequence\n')

    all_scoring_functions = general.merge_dicts(gl_scoring_functions, mat_scoring_functions)
    test_data = TestAbNNKIndex(options.input)
    print "Test Chains"+repr(len(test_data.classified_list))
    for classified_ab in test_data.classified_list:
        id = classified_ab.ig_chain.id
        print id
        for score_index in sort_index_order:
            if not all_scoring_functions.has_key(score_index): continue
            score = all_scoring_functions[score_index]

            s = score.relative_score_classified_ab(classified_ab)
            s_all = score.score_classified_ab(classified_ab)

            index_score = score.nnk_index_score

            if re.search('gl', score_index.antibody) or options.output_individual or score_index.antigen == 'AntigensCombined':
                out = "\t".join([ options.exp, id, str(options.zeros), score_index.antibody, score_index.antigen, score_index.sort, "%.3f" % s_all, "%.3f" % index_score, "%.3f"%s])
                print out
                SCOREOUT.write(out+"\n")


                #out2 = [options.exp, str(options.zeros), score_index.antibody, score_index.antigen, score_index.sort, "%.3f"%score.min_score, score.min_sequence, "%.3f"%score.max_score, score.max_sequence]
                #GROUPS.write("\t".join(out2)+"\n")

    SCOREOUT.close()
    #GROUPS.close()
    print "Done"

