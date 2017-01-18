import pandas
import os
from collections import defaultdict


class CDRDataInfo(object):
    """
    Simple class for holding and accessing Cluster and Length data for a particular Decoy.
    """
    def __init__(self, name, strategy, decoy):
        self.name = name
        self.strategy = strategy
        self.decoy = decoy

        self.data = defaultdict()
        self.cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]

    def __repr__(self):
        line = ""
        for cdr in self.cdrs:
            if self.data.has_key(cdr):
                line = line + cdr+":"+repr(self.get_value_for_cdr(cdr))
        line = line+"\n"
        #print line
        return line


    def get_data(self):
        return self.data

    def get_value_for_cdr(self, cdr):
        return self.data[cdr]

    def set_value_for_cdr(self, cdr, value):
        self.data[cdr] = value

    def set_data(self, data):
        """
        Dictionary for each CDR: L1, L2, L3, H1, H2, H3
        """
        self.data = data

    def set_value(self, cdr, value):
        self.data[cdr] = value

    def has_data(self, cdr):
        if self.data.has_key(cdr):
            return True
        else:
            return False

    def is_camelid(self):
        """
        Return True if missing light chain data
        """
        if (not self.has_data('L1')) and (not self.has_data('L2')) and (not self.has_data('L3')):
            return True
        else:
            return False

class CDRData(object):
    """
    Class holding cluster and length data from cluster or antibody features database.
    """
    def __init__(self, name, native_path, is_camelid = False):
        self.is_camelid = is_camelid
        self.name = name

        self.all_data = defaultdict()
        if is_camelid:
            self.cdrs = ["H1", "H2", "H3"]
        else:
            self.cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]


        self.native_data = None
        self._setup_native_data(native_path)

#Public

    def get_pandas_dataframe(self, cdrs=None, drop_dir_prefix = False):
        """
        Gets all data as a pandas dataframe.  Uses the set name as the score.
        You can then order, or select specific ones using the data frame.
        :rtype: pandas.DataFrame
        """

        if not cdrs:
            cdrs = self.cdrs


        temp_dict = defaultdict(list)

        for strategy in self.all_data:
            for decoy in self.all_data[strategy]:
                triple = self.all_data[strategy][decoy]
                if isinstance(triple, CDRDataInfo): pass

                temp_dict["strategy"].append(strategy)

                if drop_dir_prefix:
                    temp_dict["decoy"].append(os.path.basename(decoy))
                else:
                    temp_dict["decoy"].append(decoy)
                for cdr in cdrs:
                    temp_dict["_".join([cdr, self.name])].append(triple.get_value_for_cdr(cdr))

        columns = ["strategy", "decoy"]
        columns = columns.extend(["_".join([cdr, self.name]) for cdr in self.cdrs])
        df = pandas.DataFrame(temp_dict, columns=columns)
        df.index = df["decoy"]
        del df["decoy"]
        df = df.convert_objects(convert_numeric=True)

        return df

    def add_data(self, strategy, con):
        """
        Function to add data to the class.  Needs to be defined in subclass.
        :param strategy: Strategy for which we are adding data
        :param con: Sqlite3 Connection object

        """
        pass

    def get_strategy_data(self, strategy):
        """
        Get data for each decoy
        :param strategy: Strategy string
        :rtype: list of CDRDataInfo
        """
        return self.all_data[strategy]

    def get_strategy_data_for_decoy(self, strategy, decoy):
        """
        Get the data of the decoy
        :param strategy: Strategy string
        :param decoy: decoy including path and suffix
        :rtype: CDRDataInfo
        """
        return self.all_data[strategy][decoy]

    def get_native_data(self):
        return self.native_data

    def set_native_data_input_tag(self, con, input_tag):
        self._set_native_data_input_tag(con, input_tag, self.column_name)

    def get_concatonated_map(self, cdr = None, decoy_list = None):
        """
        Returns a defaultDic:
        Default:
            decoy: CDRDataInfo
          If CDR != None:
             [value, cdr] = CDRDataTriple
            #->CDR to get back cdr_value, decoy for sorting on cdr_value
        :rtype: defaultdict
        """

        result_data = defaultdict()
        for strategy in self.all_data:
            for decoy in self.all_data[strategy]:
                if decoy_list and decoy not in decoy_list:
                    continue
                triple = self.all_data[strategy][decoy]
                if isinstance(triple, CDRDataInfo): pass

                if cdr:
                    result_data[(triple.get_value_for_cdr(cdr), decoy)] = triple
                else:
                    result_data[decoy] = triple

        return result_data


#Private

    def _get_stmt(self, column_name):
        if not self.is_camelid:
            stmt = "SELECT "+ \
                        "structures.input_tag as decoy,"+ \
                        "H1."+column_name+" as H1," +\
                        "H2."+column_name+" as H2," +\
                        "H3."+column_name+" as H3," +\
                        "L1."+column_name+" as L1," +\
                        "L2."+column_name+" as L2," +\
                        "L3."+column_name+" as L3"

            stmt = stmt + """
                    FROM
                        cdr_clusters as L1,
                        cdr_clusters as L2,
                        cdr_clusters as L3,
                        cdr_clusters as H1,
                        cdr_clusters as H2,
                        cdr_clusters as H3,
                        structures
                    WHERE
                        structures.struct_id = L1.struct_id AND
                        structures.struct_id = L2.struct_id AND
                        structures.struct_id = L3.struct_id AND
                        structures.struct_id = H1.struct_id AND
                        structures.struct_id = H2.struct_id AND
                        structures.struct_id = H3.struct_id AND
                        L1.CDR = 'L1' AND
                        L2.CDR = 'L2' AND
                        L3.CDR = 'L3' AND
                        H1.CDR = 'H1' AND
                        H2.CDR = 'H2' AND
                        H3.CDR = 'H3'
                    """
        else:
            stmt = "SELECT "+ \
                        "structures.input_tag as decoy,"+ \
                        "H1."+column_name+" as H1," +\
                        "H2."+column_name+" as H2," +\
                        "H3."+column_name+" as H3"

            stmt = stmt + """
                    FROM
                        cdr_clusters as H1,
                        cdr_clusters as H2,
                        cdr_clusters as H3,
                        structures
                    WHERE
                        structures.struct_id = H1.struct_id AND
                        structures.struct_id = H2.struct_id AND
                        structures.struct_id = H3.struct_id AND
                        H1.CDR = 'H1' AND
                        H2.CDR = 'H2' AND
                        H3.CDR = 'H3'
                    """
        #print stmt
        return stmt

    def _get_add_data(self, strategy, con, column_name):
        self.column_name = column_name

        data = defaultdict()
        cur = con.cursor()
        print "Adding "+column_name+" data for "+strategy
        for row in cur.execute(self._get_stmt(column_name)):
            #print repr(row)
            d = CDRDataInfo(self.name, strategy, row[0])
            d.set_value('H1', row[1])
            d.set_value('H2', row[2])
            d.set_value('H3', row[3])

            if not self.is_camelid:
                d.set_value('L1', row[4])
                d.set_value('L2', row[5])
                d.set_value('L3', row[6])
            data[row[0]] = d
            #print self.name+" type: "+str(type(row[1]))
            #print repr(d)

        self._add_data(strategy, data)

    def _add_data(self, strategy, decoy_data_map):
        """
        Add data in the form of a dict of decoy:DataTriple
        """
        if not self.all_data.has_key(strategy):
            self.all_data[strategy] = defaultdict()
        self.all_data[strategy] = decoy_data_map
        #print repr(decoy_data_map)

    def _setup_native_data(self, pdb_path):
        if not pdb_path: return None

    def _set_native_data_from_biopose(self, pdb_path):
        pass

    def _set_native_data(self, data):
        self.native_data = data

    def _set_native_data_input_tag(self, con, input_tag, column_name):
        pass