#!/usr/bin/python


#This Script parses the log of MPI-run antibody_design benchmarks to get the frequency of clusters and lengths against a native database.
#It also gets DB frequences from our main databases and will soon parse the resultant length and cluster recoveries to determine the log-odds of each experiment for conclusive statistics.
#I really did not want to write this.

#exp_list is the same thing we will use for Json skipping ref.  full names on left, reference name on right.

#This was also written way before pandas was even on the radar, and all the plotting functionality within Python was laughable.

import gzip
import os
import sys

try:
    import sqlite3
    import numpy
except ImportError:
    print "Cannot import SQLITE3 or NUMPY.  Pleas install to use pool data"

from collections import defaultdict
from optparse import OptionParser
import re
import glob

import math

import RAbD_BM.tools as tools
from rosetta_jade.BenchmarkInfo import BenchmarkInfo
from RAbD_BM.AnalysisInfo import AnalysisInfo
from RAbD_BM.AnalysisInfo import NativeInfo
from basic.numeric import *

class AnalyzeRecovery:
    """
    Pools Recovery and RR data, outputs to DB
    """
    def __init__(self, pyig_design_db_path, analysis_info, native_info):
        """

        :param pyig_design_dbpath: str
        :param analysis_info: AnalysisInfo
        :param native_info: NativeInfo
        :return:
        """



        if not isinstance(analysis_info, AnalysisInfo): sys.exit()
        if not isinstance(native_info, NativeInfo): sys.exit()

        self.pyig_design_db = sqlite3.connect(pyig_design_db_path)
        self.analysis_info = analysis_info
        self.native_info = native_info
        self.structures_lam_kap = defaultdict()

        self.structures_lam_kap["lambda"], self.structures_lam_kap["kappa"] = tools.get_lambda_kappa_pdb_ids(self.analysis_info.bm_info.get_dataset(),
                                                                                               self.analysis_info.bm_info.get_input_pdb_type())

        self.recovery_calculator = RecoveryCalculator(self.native_info.get_features_db(), self.analysis_info.get_features_db())

        if self.analysis_info.bm_info.settings["CDR"] == "ALL":
            self.cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]
        else:
            self.cdrs = [self.analysis_info.bm_info.settings["CDR"]]

        self.heavy_cdrs = ["H1", "H2", "H3"]
        self.recovery_types = ["length", "cluster"]

        self.native_lengths = defaultdict(lambda: defaultdict())
        self.native_clusters = defaultdict(lambda: defaultdict())


        self.exp_totals = defaultdict(lambda: defaultdict(dict))
        self.ab_db_cdr_totals = defaultdict(lambda: defaultdict(dict))

        self.exp_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.ab_db_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.natives = defaultdict(lambda: defaultdict(dict))

        self.initialize()

    def initialize(self):
        native_db = self.native_info.get_features_db()

        if not os.path.exists(native_db):
            sys.exit("native_db path not good: "+native_db)


        con = sqlite3.connect(native_db)
        structures = []
        for row in con.execute("SELECT input_tag from structures"):
            structures.append(str(row[0]))

        """
        Initializes all data dictionaries to zero

        """
        exp = self.analysis_info.get_exp()
        for structure in structures:
            for cdr in self.cdrs:
                self.exp_totals[exp][structure][cdr] = 0
                self.ab_db_cdr_totals[exp][structure][cdr] = 0


        for type in self.recovery_types:
            exp = self.analysis_info.get_exp()
            for structure in structures:
                for cdr in self.cdrs:
                    self.exp_data[type][exp][structure][cdr] = 0
                    self.ab_db_data[type][exp][structure][cdr] = 0

        """
        Calculates Native totals

        """
        for structure in structures:
            for cdr in self.cdrs:

                for row in con.execute("SELECT cdr_clusters.length " +
                                     "FROM cdr_clusters, structures " +
                                     "WHERE structures.input_tag = ? AND " +
                                     "structures.struct_id == cdr_clusters.struct_id AND " +
                                     "cdr_clusters.CDR=?", (structure, cdr)):
                    self.native_lengths[structure][cdr] = str(row[0])
                    self.natives["length"][structure][cdr] = str(row[0])

                for row in con.execute("SELECT cdr_clusters.fullcluster " +
                                     "FROM cdr_clusters, structures " +
                                     "WHERE structures.input_tag = ? AND " +
                                     "structures.struct_id == cdr_clusters.struct_id AND " +
                                     "cdr_clusters.CDR=?", (structure, cdr)):
                    self.native_clusters[structure][cdr] = str(row[0])
                    self.natives["cluster"][structure][cdr] = str(row[0])

        con.close()


        """
        Loads native counts of the clusters and lengths of each CDR from the PyIgClassify Database.
        This data goes into the self.ab_db_data dictionary.

        This is what we will use to calculate the Risk Ratios, as this is what we are sampling on.
        """
        con = self.pyig_design_db
        def __add_data(data_dict, structure, cdr, n):
            exp = self.analysis_info.get_exp()
            data_dict[exp][structure][cdr] = n

        for lam_kap_type in self.structures_lam_kap.keys():
            for structure in self.structures_lam_kap[lam_kap_type]:
                for cdr in self.cdrs:
                    n_clusters = 0
                    n_lengths = 0

                    length = self.native_lengths[structure][cdr]
                    cluster = self.native_clusters[structure][cdr]
                    gene = lam_kap_type
                    if cdr in self.heavy_cdrs:
                        gene = "heavy"

                    #TOTAL of NATIVE CLUSTER in PyIgClassify DB
                    for row in con.execute("SELECT * from cdr_data WHERE CDR=? AND fullcluster=? AND gene = ? AND datatag != ?", (cdr, cluster, gene, 'loopKeyNotInPaper')):
                        n_clusters+=1
                    #__add_data(self.ab_db_totals_cluster, structure, cdr, n_clusters)
                    __add_data(self.ab_db_data["cluster"], structure, cdr, n_clusters)


                    #TOTAL of NATIVE LENGTH in PyIgClassify DB
                    for row in con.execute("SELECT * from cdr_data WHERE CDR=? AND length=? AND gene=? AND datatag != ?", (cdr, int(length), gene, 'loopKeyNotInPaper')):
                        n_lengths+=1
                    __add_data(self.ab_db_data["length"], structure, cdr, n_lengths)

            #HERE, we calculate the TOTAL NUMBER OF ENTRIES in the db for each CDR for a particular gene.
            for cdr in self.cdrs:
                gene = lam_kap_type
                if cdr in self.heavy_cdrs:
                    gene = "heavy"
                exp = self.analysis_info.get_exp()
                for row in con.execute("SELECT * from cdr_data WHERE CDR=? AND gene=? AND datatag != ?",(cdr, gene, 'loopKeyNotInPaper')):

                    for structure in self.structures_lam_kap[lam_kap_type]:
                        self.ab_db_cdr_totals[exp][structure][cdr]+=1 #TOTAL NUMBER OF CDR ENTRIES IN DATABASE

        con.close()


    def apply(self, db_path, drop_tables = False):
        """
        Calculate and Output all the data

        :param db_path: str
        :param drop_tables: bool
        :return:
        """
        db_dir = os.path.dirname(db_path)
        db_name = os.path.basename(db_path)
        filenames = get_filenames(self.analysis_info.get_decoy_dir(), os.path.basename(self.analysis_info.get_decoy_dir()))

        if len(filenames) == 0:
            sys.exit("No MPI Log filenames were found!")
        for name in filenames:
            #print name
            print "Parsing "+name
            self.read_pdb_graft_log_data(name, self.analysis_info.get_exp())



        def _setup_output_db(drop_tables = False):

            for type in ["length", "cluster"]:

                table = "all_"+type
                if drop_tables:
                    self.db.execute("DROP TABLE IF EXISTS "+table)


                if type == "length":
                    self.db.execute("CREATE TABLE IF NOT EXISTS "+table+"(" +
                                        "id INT, native TEXT, CDR TEXT, exp TEXT, exp_group TEXT, exp_type TEXT, length INT," +
                                        "chosen_perc REAL, db_perc REAL," +
                                        "chosen_freq INT, graft_total INT, db_freq INT, db_total INT," +
                                        "top_rec REAL, top_rr REAL, topx_rec REAL, topx_rr REAL, x INT)")
                if type == "cluster":
                    self.db.execute("CREATE TABLE IF NOT EXISTS "+table+"(" +
                                        "id INT, native TEXT, CDR TEXT, exp TEXT, exp_group TEXT, exp_type TEXT, cluster TEXT," +
                                        "chosen_perc REAL, db_perc REAL," +
                                        "chosen_freq INT, graft_total INT, db_freq INT, db_total INT," +
                                        "top_rec REAL, top_rr REAL, topx_rec REAL, topx_rr REAL, x INT)")



                table = "cdr_"+type
                if drop_tables:
                    self.db.execute("DROP TABLE IF EXISTS "+table)


                self.db.execute("CREATE TABLE IF NOT EXISTS "+table+"(" +
                                        "id INT, CDR TEXT, exp TEXT, exp_group TEXT, exp_type TEXT," +
                                        "chosen_perc REAL, db_perc REAL, avg_chosen_perc REAL, avg_db_perc REAL, "+
                                        "chosen_freq INT, graft_total INT, db_freq INT, db_total INT, " +
                                        "top_rec REAL, top_rr REAL, top_rr_avg REAL, top_rr_log REAL, "+
                                        "topx_rec REAL, topx_rr REAL, topx_rr_avg REAL, topx_rr_log REAL, x INT)")

                table = "exp_"+type
                if drop_tables:
                    self.db.execute("DROP TABLE IF EXISTS "+table)


                self.db.execute("CREATE TABLE IF NOT EXISTS "+table+"(" +
                                        "id INT, exp TEXT, exp_group TEXT, exp_type TEXT, h3_present INT, " +
                                        "chosen_perc REAL, db_perc REAL, avg_chosen_perc REAL, avg_db_perc REAL, "+
                                        "chosen_freq INT, graft_total INT, db_freq INT, db_total INT, " +
                                        "top_rec REAL, top_rr REAL, top_rr_avg REAL, top_rr_log REAL, "+
                                        "topx_rec REAL, topx_rr REAL, topx_rr_avg REAL, topx_rr_log REAL, x INT)")
        ####################################

        def _output_data(db_dir, name):

            if not os.path.exists(db_dir):
                os.mkdir(db_dir)

            data_dict = self.exp_data[name]

            ab_db_dict = self.ab_db_data[name]
            native_dict = self.natives[name]

            rec_freq = self.recovery_calculator.get_all_freq(name)
            rec_totals = self.recovery_calculator.get_all_totals(name)
            self.exp_avgs = defaultdict(lambda: defaultdict(dict))
            self.ab_db_avgs = defaultdict(lambda: defaultdict(dict))
            self.risk_ratios = defaultdict(lambda: defaultdict(dict))

            #Per Input tag + CDR
            print "Outputting all " + name
            #OUT = open(outname, "w")
            MISSING = open(db_dir+"/missing_data.txt", 'w')

            #Add Log ODDS to the end if file given!!


            #OUT.write("NATIVE\tCDR\tEXP\tNATIVE_"+name.upper()+"\tChosen_Percent\t\tDb_Percent\tChosen_Freq\tGraft_Total\tDb_Freq\tDb_Total\ttop_rec top_rr topx_rec, topx_rr, x\n")

            i = 0
            for structure in sorted(self.native_clusters.keys()):
                for cdr in self.cdrs:

                    exp = self.analysis_info.get_exp()
                    exp_freq = data_dict[exp][structure][cdr]
                    exp_total= self.exp_totals[exp][structure][cdr]
                    if exp_total == 0:

                        print("CDR Not present in data!!!??")
                        MISSING.write(name+" "+exp+" "+structure+" "+cdr+"\n")
                        #sys.exit("")
                        continue

                    i+=1
                    exp_perc = get_perc(exp_freq, exp_total)
                    ab_db_freq = ab_db_dict[exp][structure][cdr]
                    ab_db_total = self.ab_db_cdr_totals[exp][structure][cdr]
                    ab_db_perc = get_perc(ab_db_freq, ab_db_total)

                    self.exp_avgs[exp][structure][cdr] = exp_perc
                    self.ab_db_avgs[exp][structure][cdr] = ab_db_perc


                    native = native_dict[structure][cdr]
                    if name == 'length':
                        native = int(native)

                    exp_type, exp_group = self.get_exp_split(exp)

                    data = [i, structure, cdr, exp, exp_group, exp_type, native, exp_perc, ab_db_perc, exp_freq, exp_total, ab_db_freq, ab_db_total]
                    rec_data = self.get_recoveries(rec_freq, rec_totals, exp_perc, exp, structure, cdr)
                    data.extend(rec_data)

                    exec_string = "INSERT INTO all_"+name+" VALUES "+get_question_mark_str(len(data) )
                    #print repr(exec_string)
                    #print repr(data)
                    self.db.execute(exec_string, data)
                    self.db.commit()


            #OUT.close()
            MISSING.close()

            ####### Per CDR ##########
            self.output_cdr_data(name)

            ###### Per Exp No H3 ######
            cdrs = ["L1", "L2", "L3", "H1", "H2"]
            self.output_exp_data(name, cdrs)

            ###### Per Exp ############
            self.output_exp_data(name, self.cdrs)
            print "Complete"
        ###############################################################################################################
        db_path = db_dir+"/"+db_name
        if not os.path.exists(db_dir):
            os.mkdir(db_dir)

        #self._check_data()
        self.db = sqlite3.connect(db_path)
        _setup_output_db(drop_tables)
        _output_data(db_dir, "length")
        _output_data(db_dir, "cluster")
        self.db.close()

    def calc_per_cdr_data(self, data_dict, total_dict, avg_dict=None):
        """
        Calculate the per-cdr data and return the result.

        result[exp][cdr][param]

        :param data_dict:
        :param total_dict:
        :param avg_dict:
        :rtype: defaultdict(dict)
        """
        result = defaultdict(lambda: defaultdict(dict))
        for exp in data_dict.keys():
            for cdr in self.cdrs:
                n_perc = 0
                for structure in data_dict[exp].keys():

                    if not result[exp].has_key(cdr):
                        result[exp][cdr]["totals"] = 0
                        result[exp][cdr]["freq"] = 0
                        result[exp][cdr]["perc_total"] = 0

                    result[exp][cdr]["totals"]+=total_dict[exp][structure][cdr]
                    result[exp][cdr]["freq"]+=data_dict[exp][structure][cdr]
                    if avg_dict and avg_dict[exp].has_key(structure) and avg_dict[exp][structure].has_key(cdr):
                        result[exp][cdr]["perc_total"]+=avg_dict[exp][structure][cdr]

                    n_perc+=1


                if avg_dict:
                    result[exp][cdr]["avg_perc"] = result[exp][cdr]["perc_total"]/n_perc

        #This means we are calculating recoveries
        if not avg_dict:
            for exp in data_dict.keys():
                for cdr in self.cdrs:
                    rr = []
                    logrr = []
                    for structure in data_dict[exp].keys():
                        try:
                            if self.risk_ratios[exp][structure][cdr] != None:
                                rr.append(self.risk_ratios[exp][structure][cdr])
                                if self.risk_ratios[exp][structure][cdr] != 0:
                                    logrr.append(math.log(self.risk_ratios[exp][structure][cdr]))
                        except KeyError:
                            continue
                    if not rr: continue
                    result[exp][cdr]['rr_avg'] = numpy.mean(rr)
                    result[exp][cdr]['rr_log'] = numpy.mean(logrr)**math.e

        return result

    def calc_per_exp_data(self, cdrs, data_dict, total_dict, avg_dict=None):
        """
        Calculate the per-exp data when we used to do multiple experiments.

        result[exp][cdr][param]

        :param cdrs:
        :param data_dict:
        :param total_dict:
        :param avg_dict:
        :rtype: defaultdict(dict)
        """
        result = defaultdict(dict)
        for exp in data_dict.keys():
            n_perc = 0
            for structure in data_dict[exp].keys():
                for cdr in cdrs:
                    if not result.has_key(exp):
                        result[exp]["totals"] = 0
                        result[exp]["freq"] = 0
                        result[exp]["perc_total"] = 0
                    result[exp]["totals"]+=total_dict[exp][structure][cdr]
                    result[exp]["freq"]+=data_dict[exp][structure][cdr]
                    if avg_dict and avg_dict[exp].has_key(structure) and avg_dict[exp][structure].has_key(cdr):
                        result[exp]["perc_total"]+=avg_dict[exp][structure][cdr]


                    n_perc+=1
            if avg_dict:
                result[exp]["avg_perc"] = result[exp]["perc_total"]/n_perc

        #This means we are calculating recoveries
        if not avg_dict:
            for exp in data_dict.keys():
                rr = []
                logrr = []
                for cdr in self.cdrs:
                    for structure in data_dict[exp].keys():
                        try:
                            if self.risk_ratios[exp][structure][cdr] != None:
                                rr.append(self.risk_ratios[exp][structure][cdr])
                                if self.risk_ratios[exp][structure][cdr] != 0:
                                    logrr.append(math.log(self.risk_ratios[exp][structure][cdr]))
                        except KeyError:
                            continue

                result[exp]['rr_avg'] = numpy.mean(rr)
                result[exp]['rr_log'] = numpy.mean(logrr)**math.e

        return result

    def output_cdr_data(self, name):
        """
        Commit per-cdr data to the current open database at self.db

        :param name:
        :return:
        """
        data_dict = self.exp_data[name]
        ab_db_dict = self.ab_db_data[name]
        native_dict = self.natives[name]
        rec_freq = self.recovery_calculator.get_all_freq(name)
        rec_totals = self.recovery_calculator.get_all_totals(name)

        print "Outputting cdr "+name

        exp_result = self.calc_per_cdr_data(data_dict, self.exp_totals, self.exp_avgs)
        ab_db_result = self.calc_per_cdr_data(ab_db_dict, self.ab_db_cdr_totals, self.ab_db_avgs)
        rec_result = self.calc_per_cdr_data(rec_freq, rec_totals)
        #OUT = open(outname, "w")
        #OUT.write("CDR\tEXP\tAvg_Chosen_Percent\tAvg_Db_Percent\tChosen_Percent\tDb_Percent\tChosen_Freq\tGraft_Total \tDb_Freq Db_Total\t"+
        #          "top_rec top_rr top_rr_avg top_rr_log topx_rec topx_rr topx_rr_avg topx_rr_log x\n")

        i = 0
        cur = self.db.cursor()
        for cdr in self.cdrs:
            exp = self.analysis_info.get_exp()
            i+=1
            exp_perc =  get_perc(exp_result[exp][cdr]["freq"], exp_result[exp][cdr]["totals"])
            ab_db_perc = get_perc(ab_db_result[exp][cdr]["freq"], ab_db_result[exp][cdr]["totals"])
            avg_perc = exp_result[exp][cdr]["avg_perc"]
            avg_db_perc = ab_db_result[exp][cdr]["avg_perc"]
            exp_freq = exp_result[exp][cdr]["freq"]
            exp_total = exp_result[exp][cdr]["totals"]
            db_freq = ab_db_result[exp][cdr]["freq"]
            db_total = ab_db_result[exp][cdr]["totals"]

            exp_type, exp_group = self.get_exp_split(exp)

            data = [i, cdr, exp, exp_group, exp_type, avg_perc, avg_db_perc, exp_perc, ab_db_perc, exp_freq, exp_total,db_freq,db_total ]
            #print data
            rec_data = get_recoveries_avg(rec_result, exp, exp_perc, cdr)

            data.extend(rec_data)

            exec_string = "INSERT INTO cdr_"+name+" VALUES"+get_question_mark_str(len(data))

            #print repr(exec_string)
            #print repr(data)

            cur.execute(exec_string, data)

        cur.close()
        self.db.commit()
        #OUT.close()

    def output_exp_data(self, name, cdrs):
        """
        Commit the per-exp recovery and risk ratio data to the currently open database.

        :param name: str
        :param cdrs: [str]
        :return:
        """
        data_dict = self.exp_data[name]
        ab_db_dict = self.ab_db_data[name]
        native_dict = self.natives[name]
        rec_freq = self.recovery_calculator.get_all_freq(name)
        rec_totals = self.recovery_calculator.get_all_totals(name)

        print "Outputting exp "+name
        exp_result = self.calc_per_exp_data(cdrs, data_dict, self.exp_totals, self.exp_avgs)
        ab_db_result = self.calc_per_exp_data(cdrs, ab_db_dict, self.ab_db_cdr_totals, self.ab_db_avgs)
        rec_result = self.calc_per_exp_data(cdrs, rec_freq, rec_totals)

        #OUT = open(outname, "w")
        #OUT.write("EXP h3_present, Avg_Chosen_Percent Avg_Db_Percent Chosen_Percent Db_Percent Chosen_Freq Graft_Total  Db_Freq Db_Total\t"+
        #    "top_rec top_rr top_rr_avg top_rr_log topx_rec topx_rr topx_rr_avg topx_rr_log x\n")

        i = 0
        cur = self.db.cursor()


        #Bad assumption here!
        h3_present = True
        if len(cdrs) != 6:
            h3_present = False

        exp = self.analysis_info.get_exp()
        i+=1
        exp_perc =  get_perc(exp_result[exp]["freq"], exp_result[exp]["totals"])
        ab_db_perc = get_perc(ab_db_result[exp]["freq"], ab_db_result[exp]["totals"])
        avg_perc = exp_result[exp]["avg_perc"]
        avg_db_perc = ab_db_result[exp]["avg_perc"]
        exp_freq = exp_result[exp]["freq"]
        exp_total = exp_result[exp]["totals"]
        db_freq = ab_db_result[exp]["freq"]
        db_total = ab_db_result[exp]["totals"]

        exp_type, exp_group = self.get_exp_split(exp)

        data = [i, exp, exp_group, exp_type, h3_present, avg_perc, avg_db_perc, exp_perc, ab_db_perc, exp_freq, exp_total,db_freq,db_total ]
        rec_data = get_recoveries_avg(rec_result, exp, exp_perc)

        data.extend(rec_data)

        exec_string = "INSERT INTO exp_"+name+" VALUES"+get_question_mark_str(len(data))
        cur.execute(exec_string, data)

        cur.close()
        self.db.commit()
        #OUT.close()

    def get_recoveries(self, rec_freq, rec_totals, exp_perc, exp, structure, cdr):
        """
        Get the recovery array from the risk ratios dictionary.

        :param rec_freq:
        :param rec_totals:
        :param exp_perc:
        :param exp:
        :param structure:
        :param cdr:
        :return:
        """
        rec_array = []
        x = ""
        rec_exp = exp

        print rec_exp+" "+structure+" "+cdr
        print rec_freq[rec_exp][structure][cdr]
        print rec_totals[rec_exp][structure][cdr]
        rec_perc = get_perc(rec_freq[rec_exp][structure][cdr], rec_totals[rec_exp][structure][cdr])

        if exp_perc == 0.0:
            rr = None
        else:
            rr = rec_perc/exp_perc
        self.risk_ratios[rec_exp][structure][cdr] = rr
        rec_array.append(rec_perc)
        rec_array.append(rr)

        #Due to removing topx - Remove this for New DBs
        rec_array.append(0)
        rec_array.append(0)

        rec_array.append(x)
        return rec_array

    def read_pdb_graft_log_data(self, filename, ref_name):
        """
        Reads a RAbD PDB file and adds the grafted data from it to help calculate Risk Ratios.

        :param filename:
        :param ref_name:
        :return:
        """
        if os.path.basename(filename).split('.')[-1] == "gz":
            INFILE = gzip.open(filename, 'rb')
        else:
            INFILE = open(filename, 'r')
        for line in INFILE:
            if re.search(" DATA GRAFT_CLOSURE ", line):
                lineSP = line.strip().split()
                data_index = 0
                for i in range(0, len(lineSP)):
                    if lineSP[i]=="DATA":
                        data_index = i
                        break

                input_tag = lineSP[data_index+3]
                cluster = lineSP[data_index+5]
                cdr = cluster.split("-")[0]
                length = cluster.split("-")[1]


                # Now, Add the data to self.exp_data
                print "Adding data for ref_name "+ref_name
                print "Adding data for input_tag "+input_tag
                print "Adding data for cdr " +cdr

                self.exp_totals[ref_name][input_tag][cdr]+=1

                if length == self.native_lengths[input_tag][cdr]:
                    #self.exp_lengths[ref_name][input_tag][cdr] +=1
                    self.exp_data["length"][ref_name][input_tag][cdr]+=1

                if cluster == self.native_clusters[input_tag][cdr]:
                    #self.exp_clusters[ref_name][input_tag][cdr] +=1
                    self.exp_data["cluster"][ref_name][input_tag][cdr]+=1


        INFILE.close()

    def get_exp_split(self, exp):
        """
        This NO longer works and needs to be deleted!
        :param exp:
        :return:
        """
        expSP = exp.split(".")
        exp_type = None
        exp_group = None
        if len(expSP) == 2:
            exp_group = expSP[0]
            exp_type = ".".join(expSP[1:])
        else:
            exp_group = ".".join(expSP[0:2])
            exp_type = expSP[-1]

        return exp_type, exp_group


class RecoveryCalculator:
    def __init__(self, native_db_path,bm_db_path ):
        self.native_db = sqlite3.connect(native_db_path)
        self.bm_db = sqlite3.connect(bm_db_path)

        self.all_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))

    def get_all_freq(self, type):

        return self.all_data[type]["freq"]

    def get_all_totals(self, type):

        return self.all_data[type]["total"]

    def apply(self):
        """
        Calculate length and cluster recoveries.  Store them the same way we used to for the recovery parser.
        :return:
        """
        pass

"""
class RecoveryParser:
    def __init__(self):
        self.all_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))

    def parse_data(self, features_output_path):
        self.in_path = features_output_path
        if not os.path.exists(features_output_path):
            sys.exit("Could not find features path: "+features_output_path)

        self._parse_all_data("length")
        self._parse_all_data("cluster")

    def _parse_all_data(self, type):
        filename = "cdr_recovery_by_pdb_"+type+"*"
        glob_str = self.in_path+"/"+filename
        print glob_str
        filenames = glob.glob(glob_str)

        assert len(filenames) == 1
        filename = filenames[0]
        if not os.path.exists(filename):
            sys.exit("Features file does not exist "+filename)

        INFILE = open(filename, 'r')
        INFILE.readline() #Skip first line

        for line in INFILE:
            line = line.strip()
            if not line: continue
            if line[0] == "#": continue
            lineSP = line.split()
            cdr = lineSP[0]
            structure = lineSP[1]+".pdb"
            exp = lineSP[2]
            freq = int(lineSP[4])
            total = int(lineSP[5])
            self.all_data[type]["freq"][exp][structure][cdr] = freq
            self.all_data[type]["total"][exp][structure][cdr] = total



    def get_all_freq(self, type):

        return self.all_data[type]["freq"]

    def get_all_totals(self, type):

        return self.all_data[type]["total"]

"""

def run_R_scripts(db_path, out_path, db_info_array):

    print db_path
    print out_path
    if not os.path.exists(out_path):
            os.mkdir(out_path)

    exp_names = [ db.final_name for db in db_info_array ]

    exp_str = " ".join(exp_names)

    print exp_str
    script = "pooled_data/plot_pooled_data.R"

    base_cmd = script+" "+db_path+" "+out_path
    cmd = base_cmd+" length "+exp_str
    print cmd
    os.system(cmd)

    cmd = base_cmd+" cluster "+exp_str
    print cmd
    os.system(cmd)

def get_question_mark_str(length):

    s = "("
    for i in range(0, length-1):
        s = s+"?,"
    s = s + "?)"
    return s

def get_filenames(input_dir, tag):
    """
    Use GLOB to Match on tag for file names in the input dir.
    This should skip all the extra PDBs like excn, initial, relax, etc.
    :param input_dir: str
    :param tag: str
    """
    if not os.path.exists(input_dir):
        sys.exit("MPI Log Dir does not exist!")
    search_name = input_dir+"/"+tag+"*"
    #print search_name
    return glob.glob(search_name)

def convert_types_to_line(array_of_types, line):
    for entry in array_of_types:

        if type(entry) == float:
            line = line + " " + get_n_s(entry)
        elif type(entry) == int:
            line = line + " "+ repr(entry)
        elif type(entry) == None:
            line = line + " "+ repr(entry)
        elif type(entry) == str:
            line = line + " "+entry
        else:
            line = line + " "+repr(entry)

    return line

def get_recoveries_avg(rec_data, exp, exp_perc, cdr=None):
    """
    I really have no idea what the hell this does.

    :param rec_data:
    :param exp:
    :param exp_perc:
    :param cdr:
    :return:
    """
    x = ""
    rec_array = []
    rec_exp = exp
    rec_perc = ""
    rr_avg =""
    if cdr:
        print rec_data[rec_exp][cdr]
        if not rec_data[rec_exp][cdr].has_key('rr_avg'):
            rec_data[rec_exp][cdr]['rr_avg'] = 0
        if not rec_data[rec_exp][cdr].has_key('rr_log'):
            rec_data[rec_exp][cdr]['rr_log'] = 0

        print "EXP "+repr(rec_data[rec_exp][cdr]['freq']) + " "+ repr(rec_data[rec_exp][cdr]['totals'])
        rec_perc = get_perc(rec_data[rec_exp][cdr]['freq'], rec_data[rec_exp][cdr]['totals'])
        rr_avg = rec_data[rec_exp][cdr]['rr_avg']
        rr_log = rec_data[rec_exp][cdr]['rr_log']
        rr = get_perc(rec_perc, exp_perc)/100
    else:
        if not rec_data[rec_exp].has_key('rr_avg'):
            rec_data[rec_exp]['rr_avg'] = 0
        if not rec_data[rec_exp].has_key('rr_log'):
            rec_data[rec_exp]['rr_log'] = 0

        rec_perc = get_perc(rec_data[rec_exp]['freq'], rec_data[rec_exp]['totals'])
        rr_avg = rec_data[rec_exp]['rr_avg']
        rr_log = rec_data[rec_exp]['rr_log']
        rr = get_perc(rec_perc, exp_perc)/100

    rec_array.append(rec_perc)
    rec_array.append(rr)
    rec_array.append(rr_avg)
    rec_array.append(rr_log)

    #Due to removing topx - Remove this for New DBs
    rec_array.append(0)
    rec_array.append(0)
    rec_array.append(0)
    rec_array.append(0)

    rec_array.append(x)
    return rec_array







