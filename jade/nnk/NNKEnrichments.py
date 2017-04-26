import sys, os, numpy
from collections import defaultdict
import pandas

from jade.nnk import NNKAbMaturation

class NNKEnrichments(object):
    """
    Simple class that holds all the enrichment data for a particular class, antibody, and antigen.
    """
    def __init__(self, data_dir, class_type = 'VRC01', antibody = 'glVRC01', antigen='GT81', sort = 'S1'):

        data_loader = NNKAbMaturation.GetNNKData(data_dir, antibody)

        self.df = data_loader.get_2D_data_freq_nnk_data(antigen = antigen, sort=sort)
        self.df = self.df.applymap(numpy.log)
        if not isinstance(self.df, pandas.DataFrame):sys.exit()

    def max(self, position):
        """

        Get the maximum enrichment at a particular position.

        :param position:
        :return:
        """

        return numpy.max(self.df[str(position)])

    def value(self, position, three_letter_code):
        """

        Get the enrichment value of a particular position and code.

        :param position:
        :param three_letter_code:
        :return:
        """
        return self.df.get_value(three_letter_code, position)