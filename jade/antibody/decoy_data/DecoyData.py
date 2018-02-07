import math
import os

import pandas

from jade.basic.sql.StatementCreator import *


########################################################################################################################
###   DecoyData
########################################################################################################################
class DecoyDataTriple(object):
    """
    Struct for holding data instead of a tupple
    """
    def __init__(self, strategy, struct_id, decoy, score, out_name, raw_name):
        self.strategy = strategy
        self.struct_id = struct_id
        self.decoy = decoy
        self.score = score
        self.score_type = out_name
        self.raw_name = raw_name

class DecoyData(object):
    def __init__(self, name, has_real_values = True, reverse_top = False):


        self.name = name
        self.interface = 'LH_A'

        self._has_real_values = has_real_values
        self.reverse_top = reverse_top

        self.all_data = defaultdict()
        #all_data is dict of [ strategy ][ input_tag (decoy) ][ DecoyDataTriple ]

        self.filters = None
        self.filter_name = None

    def get_pandas_dataframe(self, top_n = None, drop_dir_prfix = False):
        """
        Gets all data as a pandas dataframe.  Uses the set name as the score.
        You can then order, or select specific ones using the data frame.
        :return: pandas.DataFrame
        """

        temp_dict = defaultdict(list)
        for strategy in sorted(self.all_data):

            decoys = self.get_ordered_decoy_list(strategy, top_n)

            for decoy in decoys:
                triple = self.all_data[strategy][decoy]
                if isinstance(triple, DecoyDataTriple): pass

                temp_dict["strategy"].append(strategy)
                if drop_dir_prfix:
                    temp_dict["decoy"].append(os.path.basename(decoy))
                else:
                    temp_dict["decoy"].append(decoy)

                temp_dict[self.name].append(float(triple.score))

        columns = ["strategy", "decoy", self.name]

        df = pandas.DataFrame(temp_dict, columns=columns)
        df.index = df["decoy"]
        del df["decoy"]
        df = df.convert_objects(convert_numeric=True)

        return df


    def set_interface(self, interface):
        """
        Set the Antibody-Antigen interface - used mainly for H_A vs LH_A
        """
        self.interface = interface

    def add_filters(self, filters, filter_name):
        self.filters = filters
        self.filter_name = filter_name

    def get_outname(self):
        if not self.filters: return self.name
        else: return self.name+"_"+self.filter_name

    def has_real_values(self):
        return self._has_real_values

    def add_data(self, strategy, con):
        """
        Baseclass method - needs to be overridden in subclass
        :param strategy: Strategy to which we are adding data.
        :param con: Sqlite3 Connection object

        """
        pass

    def _get_add_data(self, strategy, stmt_creator, con):

        if isinstance(stmt_creator, StatementCreator): pass

        if self.filters:
            for filter in self.filters:
                stmt_creator.add_data_filter(filter)

        print "Getting data for: "+self.name

        stmt = stmt_creator.create_statement()

        data = defaultdict()
        cur = con.cursor()
        for row in cur.execute(stmt):
            #print repr(row)
            triple = DecoyDataTriple(strategy, row[0], str(row[1]), row[2], self.get_outname(), self.name)
            data[str(row[1])] = triple
        self._add_data(strategy, data)

    def _add_data(self, strategy, decoy_data_map):
        """
        Add data in the form of a dict of decoy:DataTriple
        """
        if not self.all_data.has_key(strategy):
            self.all_data[strategy] = defaultdict()
        self.all_data[strategy] = decoy_data_map

    def get_data_for_decoy(self, strategy, decoy):
        """
        Get the held data for the decoy
        :param strategy: Strategy Name
        :param decoy: Decoy name (with dir and suffix)
        :rtype: DecoyDataTriple
        """
        return self.all_data[ strategy ][ decoy ]

    def get_top_x_percent_cutoff_value(self, strategy, top_percent):

        if top_percent < 1.0: top_percent = top_percent*100

        decoys = self.get_ordered_decoy_list(strategy)
        last_entry = math.ceil(len(decoys)/float(top_percent))
        print "Using top "+repr(last_entry)+" entries as subset"
        return self.get_data_for_decoy(strategy, decoys[ int(last_entry) -1]).score


    #############################################################################
    def get_strategy_data(self, strategy, by_score_tuple = False):
        """
        For a particular strategy:
        Return a dictionary of decoy:DataTriple
        or if by_score_tuple:
            [score, decoy] = DataTriple
        """
        #print self.name
        if by_score_tuple:
            out_dict = defaultdict()
            for decoy in self.all_data[strategy]:
                triple = self.all_data[strategy][decoy]
                out_dict[(triple.score, triple.decoy)] = triple
            return out_dict
        else:
            return self.all_data[strategy]

    def get_top_strategy_data(self, strategy, top_n, by_score_tuple = False):
        """
        For a particular strategy:
        Return a dictionary of decoy:DataTriple
        or if by_score_tuple:
            [score, decoy] = DataTriple

        For only the top scoring decoys
        """
        top_n = int(top_n)
        data = self.get_strategy_data(strategy, True)

        reverse = False
        if self.reverse_top: reverse = True

        top_points = sorted(data.keys(), reverse = reverse)[0:top_n]

        top_data = defaultdict()
        for point in top_points:
            if by_score_tuple:
                top_data[point] = data[point]
            else:
                top_data[point[1]] = data[point]

        return top_data

    def get_top_all_data(self, top_n, by_score_tuple = False):
        """
        Over all the strategies:
        Return a dictionary of decoy:DataTriple
        or if by_score_tuple:
            [score, decoy] = DataTriple
        """
        top_n = int(top_n)
        data = self.get_concatonated_map(True)

        reverse = False
        if self.reverse_top: reverse = True

        top_points = sorted(data.keys(), reverse = reverse)[0:top_n]

        top_data = defaultdict()
        for point in top_points:
            if by_score_tuple:
                top_data[point] = data[point]
            else:
                top_data[point[1]] = data[point]

        return top_data

    def get_ordered_decoy_list(self, strategy, top_n = None):
        """
        Get an ordered array of decoy names by energy for a particular strategy
        :rtype: list of str
        """
        if self.reverse_top: reverse = True
        else: reverse = False

        decoys = []
        if top_n:
            data = self.get_top_strategy_data(strategy, top_n, True)
        else:
            data = self.get_strategy_data(strategy, True)
        #print repr(data)
        points = sorted(data.keys(), reverse = reverse)
        #print points
        for tup in points:
            decoys.append(data[tup].decoy)

        return decoys

    def get_ordered_decoy_list_all(self, top_n = None):
        """
        Get an ordered array of decoy names by energy over all the strategies
        :rtype: list of str
        """
        if self.reverse_top: reverse = True
        else: reverse = False

        if top_n:
            all_data = self.get_top_all_data(top_n, True)
        else:
            all_data = self.get_concatonated_map(True)

        decoys = []
        points = sorted(all_data.keys(), reverse = reverse)[0:top_n]
        for tup in points:
            decoys.append(all_data[tup].decoy)

        return decoys

    def get_concatonated_map(self, by_score_tuple = False):
        """
        Returns a defaultDic:
        Default:
            decoy: DecoyDataTriple

        by_score_tuple (for sorting on score and having possible redundancy)
            [score, decoy]: DecoyDataTriple
        """

        result_data = defaultdict()
        for strategy in self.all_data:
            for decoy in self.all_data[strategy]:
                triple = self.all_data[strategy][decoy]
                if isinstance(triple, DecoyDataTriple): pass

                if by_score_tuple:
                    result_data[(triple.score, decoy)] = triple
                else:
                    result_data[decoy] = triple
        return result_data
