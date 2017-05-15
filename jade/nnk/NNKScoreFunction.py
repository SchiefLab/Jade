import sys, numpy, pandas
from collections import defaultdict

from jade.basic.RestypeDefinitions import RestypeDefinitions, ResTypeSergey
from jade.nnk.NNKEnrichments import NNKEnrichments, combine_enrichments
from jade.nnk.NNKIndex import TemplateAbNNKIndex, TestAbNNKIndex
from pyigclassify.modules.chains.util import ClassifiedAb
from pyigclassify.modules.Structure import PDBResInfo, Residue, AntibodyResidue
from jade.basic.numeric import linear_rescale

class AntibodyNNKScoreFunction(object):
    """
    A class that takes a list of NNKEnrichment instances to score a potential antibody.
    """
    def __init__(self, list_of_enrichments, template_index,
                 award_max_only = False, use_factors = False, award_max_and_conserved_only = False,
                 use_global_matches = False, use_group_matches = False, additive_combine = False, ):
        """
        For now, we calculate the score up to H3.  If we knew the length of H3 for each NNK,
          we could do the other side of the tail.

        :type list_of_enrichments: [NNKEnrichments]
        :type template_index: TemplateAbNNKIndex
        """

        #Options:
        self.award_max_only = award_max_only
        self.use_factors = use_factors
        self.award_max_conservative_only = award_max_and_conserved_only
        self.use_global_matches = use_global_matches
        self.use_group_matches = use_group_matches

        self.template_index = template_index
        self.enrichments = list_of_enrichments
        self.additive_combine = additive_combine
        self.enrich = combine_enrichments(list_of_enrichments, additive_combine)

        if use_factors:
            print "Calculating factors"
            self.factors = self.enrich.calculate_factors()

        self.res = RestypeDefinitions()
        self.groups = ResTypeSergey()

        if not isinstance( self.enrich, NNKEnrichments ): sys.exit()

        #Calculate Min and Maxes based on the germline indexes.
        self.global_matches, self.group_matches = self.__calculate_global_and_group_matches()
        self.__calculate_min_max()




    def score_whole_antibody(self, antibody):
        """
        The Antibody can be made of multiple light/heavy chains.
        We will score all of them that match the template - so if you have multiple Heavy chains and your
        template is heavy, then you will get scores that do NOT make sense!
        -> Prune the classified_list of the antibody to fix this or prune your PDB/FASTA

        :type antibody: TestAbNNKIndex
        :rtype: float
        """
        #For now, we are using 0 to 1, but what we really need to do is match the numbers.


        score_sum = 0
        for classified_ab in antibody.classified_list():
            score_sum += self.score_classified_ab(classified_ab)

        return score_sum

    def score_classified_ab(self, classified_ab):
        """

        :type classified_ab: ClassifiedAb
        :rtype: float
        """
        score_sum = 0
        if not isinstance(classified_ab, ClassifiedAb): sys.exit()
        res_info = classified_ab.ab_chain.get_pdb_res_info()
        if not isinstance(res_info, PDBResInfo): sys.exit()

        for i in range(1, res_info.total_residue() + 1):
            tup = (res_info.get_residue(i).tuple())
            res = res_info.get_residue(i)
            if not self.template_index.index.has_key(tup):
                continue
            else:
                # Get the energy here:
                position = self.template_index.index[tup]
                score = self.__score(position, self.res.get_three_letter_from_one(res.get_aa()))
                score_sum += score

        return score_sum

    def __score(self, position, three_letter):

        template_aa = self.res.get_three_letter_from_one(self.template_index.res_info.get_residue(position).get_aa())

        v = self.enrich.value(position, three_letter )
        #print repr(v)
        v_template = self.enrich.value(position, template_aa)


        if self.use_factors:
            v = self.factors.get_value( three_letter, str(position) )
            #print "F: "+repr(v)
            v_template = self.factors.get_value( template_aa, str(position))

        max_value, max_aa = self.enrich.max(position)




        if self.award_max_only or self.award_max_conservative_only:
            if three_letter == max_aa:
                return v

            elif self.award_max_conservative_only and self.res.is_conserved(three_letter, max_aa):
                return v

        elif self.use_global_matches or self.use_group_matches:
            if self.use_global_matches and self.global_matches[position]:
                #print "Global Match: "+str(position)+" "+repr(v)
                if three_letter == max_aa:
                    return v
                else:
                    #Subtract points for it being a global match, but not actually matching.
                    if v < 0:
                        return v
                    else:
                        return -v

            elif self.use_group_matches and self.group_matches[position]:
                #print repr(v)
                if self.additive_combine:
                    return v +v_template
                else:
                    v = (v + v_template)/2.0

                #print "Group Match: " + str(position) + " " + repr(v)
                #Return the mean of the grouped values.
                return v

        else:
            #print "None max"
            return v

        return 0

    def relative_score_whole_antibody(self, antibody):
        """

        :type antibody: TestAbNNKIndex
        :rtype: float
        """
        score = self.score_whole_antibody(antibody)

        return linear_rescale(0, self.max_score, score)

    def relative_score_classified_ab(self, classified_ab):
        """

        :type classified_ab: ClassifiedAb
        :rtype: float
        """

        score = self.score_classified_ab(classified_ab)
        return linear_rescale(0, self.nnk_index_score, score)


    def max_sequence(self):
        """
        Return the amino acid sequence of the maximum score.
        :rtype: str
        """
        return self.max_sequence

    def min_sequence(self):
        """
        Return the amino acid sequence of the minimum score.
        :rtype: str
        """
        return self.min_sequence


    def __max_score(self):
        """
        Return the maximum score you can obtain from this function.
        :rtype: float
        """
        return self.max_score

    def __min_score(self):
        """
        Return the minimum score you can obtain from this function.
        :rtype: float
        """
        return self.min_score

    def __index_score(self):
        """
        Return the score of the NNKIndex used here.  This now corresponds to the maximumn score of the scorefunction.
        
        :return: 
        """
        return self.nnk_index_score

    def __calculate_min_max(self):

        self.min_score = 0; self.max_score = 0; self.nnk_index_score = 0
        self.min_sequence = ""; self.max_sequence = ""

        for position in sorted(self.template_index.reverse_index.keys()):

            nnk_three_letter = self.res.get_three_letter_from_one( self.template_index.res_info.get_residue(position).get_aa() )

            min_sc, min_seq = self.enrich.min(position)
            max_sc, max_seq = self.enrich.max(position)

            self.min_score += min_sc; self.max_score += max_sc

            self.min_sequence += self.res.get_one_letter_from_three(min_seq)
            self.max_sequence += self.res.get_one_letter_from_three(max_seq)

            #print position, nnk_three_letter
            #print repr(self.__score(position, nnk_three_letter))

            self.nnk_index_score += self.__score(position, nnk_three_letter)

        print repr(self.nnk_index_score)

    def __calculate_global_and_group_matches(self):
        group_matches = defaultdict()
        global_matches = defaultdict()

        for position in sorted(self.template_index.reverse_index.keys()):
            group_matches[position] = False
            global_matches[position] = False

            aa = self.res.get_three_letter_from_one( self.template_index.res_info.get_residue(position).get_aa() )
            max_sc, max_seq = self.enrich.max(position)
            if aa == max_seq:
                global_matches[position] = True
                continue

            if self.groups.has_common_group(aa, max_seq):
                group_matches[position] = True

        return global_matches, group_matches




class MetaAntibodyNNKScoreFunction(object):
    def __init__(self, list_of_scorefunctions, additive_combine = False):
        """
        
        :rtype list_of_scorefunctions: [AntibodyNNKScoreFunction]
        """
        self.scorefxns = list_of_scorefunctions
        self.additive_combine = additive_combine

        #print repr(numpy.array([scorefxn.nnk_index_score for scorefxn in self.scorefxns]))

        if self.additive_combine:
            self.nnk_index_score = numpy.sum(numpy.array([scorefxn.nnk_index_score for scorefxn in self.scorefxns]))
        else:
            self.nnk_index_score = numpy.mean(numpy.array([scorefxn.nnk_index_score for scorefxn in self.scorefxns]))

    def relative_score_classified_ab(self, classified_ab):
        """

        :type classified_ab: ClassifiedAb
        :rtype: float
        """

        s = self.score_classified_ab(classified_ab)

        if self.additive_combine:
            totals = numpy.sum(numpy.array([scorefxn.nnk_index_score for scorefxn in self.scorefxns]))
            return linear_rescale(0, totals, s)
        else:
            #print repr(numpy.array([scorefxn.nnk_index_score for scorefxn in self.scorefxns]))
            #print repr(numpy.array([scorefxn.relative_score_classified_ab(classified_ab) for scorefxn in self.scorefxns]))
            return linear_rescale(0, self.nnk_index_score, self.score_classified_ab(classified_ab))

    def score_classified_ab(self, classified_ab):
        scores = [scorefxn.score_classified_ab(classified_ab) for scorefxn in self.scorefxns]

        if self.additive_combine:
            return numpy.sum(numpy.array(scores))
        else:
            return numpy.mean(numpy.array(scores))

    def __score(self, position, three_letter):
       scores = [scorefxn.__score(position, three_letter) for scorefxn in self.scorefxns]
       if self.additive_combine:
           return numpy.sum(numpy.array(scores))
       else:
           return numpy.mean(numpy.array(scores))