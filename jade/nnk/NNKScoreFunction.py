import sys, numpy, pandas

from jade.basic.RestypeDefinitions import RestypeDefinitions
from jade.nnk.NNKEnrichments import NNKEnrichments, combine_enrichments
from jade.nnk.NNKIndex import TemplateAbNNKIndex, TestAbNNKIndex
from pyigclassify.modules.chains.util import ClassifiedAb
from pyigclassify.modules.Structure import PDBResInfo, Residue, AntibodyResidue
from jade.basic.numeric import linear_rescale

class AntibodyNNKScoreFunction(object):
    """
    A class that takes a list of NNKEnrichment instances to score a potential antibody.
    """
    def __init__(self, list_of_enrichments, template_index):
        """
        For now, we calculate the score up to H3.  If we knew the length of H3 for each NNK,
          we could do the other side of the tail.

        :type list_of_enrichments: [NNKEnrichments]
        :type template_index: TemplateAbNNKIndex
        """
        self.template_index = template_index
        self.enrichments = list_of_enrichments
        self.enrich = combine_enrichments(list_of_enrichments)

        self.res = RestypeDefinitions()

        if not isinstance( self.enrich, NNKEnrichments ): sys.exit()

        #Calculate Min and Maxes based on the germline indexes.
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
                score = self.enrich.value(position, self.res.get_three_letter_from_one( res.get_aa() ))
                score_sum += score

        return score_sum

    def relative_score_whole_antibody(self, antibody):
        """

        :type antibody: TestAbNNKIndex
        :rtype: float
        """
        score = self.score_whole_antibody(antibody)

        return linear_rescale(self.min_score, self.max_score, score)

    def relative_score_classified_ab(self, classified_ab):
        """

        :type classified_ab: ClassifiedAb
        :rtype: float
        """

        score = self.score_classified_ab(classified_ab)
        return linear_rescale(self.min_score, self.max_score, score)


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

    def __calculate_min_max(self):

        s = len(self.template_index.get_seq_up_to_h3())

        self.min_score = 0; self.max_score = 0
        self.min_sequence = ""; self.max_sequence = 0

        for i in range(1, s+1):
            min_sc, min_seq = self.enrich.min(i)
            max_sc, max_seq = self.enrich.max(i)

            self.min_score += min_sc; self.max_score += max_sc

            self.min_sequence += self.res.get_one_letter_from_three(min_seq)
            self.max_sequence += self.res.get_one_letter_from_three(max_seq)

