#!/usr/bin/python


#This Script parses the log of MPI-run antibody_design benchmarks to get the frequency of clusters and lengths against a native database.
#It also gets DB frequences from our main databases and will soon parse the resultant length and cluster recoveries to determine the log-odds of each experiment for conclusive statistics.
#I really did not want to write this.

#exp_list is the same thing we will use for Json skipping ref.  full names on left, reference name on right.

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

class PoolData:
    """
    Pools Recovery and RR data, outputs to DB
    """
    def __init__(self, lambda_db, kappa_db, ab_db, features_base_dir, benchmark_info):
        """

        :param lambda_db: str
        :param kappa_db: str
        :param ab_db: str
        :param features_base_dir: str
        :param benchmark_info: BenchmarkInfo
        """

        self.recovery_parser = RecoveryParser()
        self.features_base_dir = features_base_dir
        self.recovery_base_dir = self.features_base_dir+"/recovery"

        if not isinstance(benchmark_info, BenchmarkInfo):sys.exit()

        self.benchmark_info = benchmark_info

        if self.benchmark_info.settings["CDR"] == "ALL":
            self.cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]
        else:
            self.cdrs = [self.benchmark_info.settings["CDR"]]

        self.heavy = ["H1", "H2", "H3"]
        self.types = ["length", "cluster"]

        self.native_lengths = defaultdict(lambda: defaultdict())
        self.native_clusters = defaultdict(lambda: defaultdict())


        self.exp_totals = defaultdict(lambda: defaultdict(dict))
        self.ab_db_cdr_totals = defaultdict(lambda: defaultdict(dict))

        self.exp_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.ab_db_data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        self.natives = defaultdict(lambda: defaultdict(dict))

        self.structures_lam_kap = defaultdict()

        print "LAMBDA:" +lambda_db+" KAPPA: "+kappa_db
        self._parse_recovery_tables()
        self._parse_native_db(lambda_db, "lambda")
        self._parse_native_db(kappa_db, "kappa")
        self._parse_ab_design_db(ab_db)

    def apply(self, db_out_dir, db_out_name, append_database, old_style = False):
        """
        Pool data.
        Output Database to db_out_dir using name given.
        Optionally, append the database in that directory with the given name.

        Old_style is deprecated.
        If the DatabaseInfo object has a log dir, use that.  Else, we parse it from the PDB files.
        """


        exp_path = self.benchmark_info.full_name
        if self.benchmark_info.has_log_path():
            filenames = self.get_filenames(self.benchmark_info.log_dir, exp_path)
        else:
            filenames = self.get_filenames(self.benchmark_info.decoy_dir, exp_path)

        if len(filenames) == 0:
            sys.exit("No MPI Log filenames were found!")
        for name in filenames:
            #print name
            print "Parsing "+name
            if old_style:
                self.parse_single_file_old_style(name, self.benchmark_info.final_name)
            else:
                self.parse_single_file_new_style(name, self.benchmark_info.final_name)

        self.output_data(append_database, db_out_dir, db_out_name)

    def _parse_recovery_tables(self):
        #print repr(self.db_info_array)
        print repr(self.benchmark_info)
        final_name = self.benchmark_info.final_name
        features_path = self.recovery_base_dir+"/"+final_name+"/output_tab_delimited_table"
        self.recovery_parser.parse_data(features_path)

    def _check_data(self):
        print "EXP LENGTHS"
        print repr(self.exp_lengths)

        print "EXP CLUSTERS"
        print repr(self.exp_clusters)

        print "ABDB LENGTHS"
        print repr(self.ab_db_totals_length)

        print "ABDB CLUSTERS"
        print repr(self.ab_db_totals_cluster)

        print "EXP Totals"
        print repr(self.exp_totals)

        print "ABDB Totals"
        print repr(self.ab_db_cdr_totals)

    def _initialize_data(self, structures):

        exp = self.benchmark_info.final_name
        for structure in structures:
            for cdr in self.cdrs:
                self.exp_totals[exp][structure][cdr] = 0
                #self.exp_clusters[exp][structure][cdr] = 0
                #self.exp_lengths[exp][structure][cdr] = 0

                self.ab_db_cdr_totals[exp][structure][cdr] = 0
                #self.ab_db_totals_cluster[exp][structure][cdr] = 0
                #self.ab_db_totals_length[exp][structure][cdr] = 0

        for type in self.types:
            exp = self.benchmark_info.final_name
            for structure in structures:
                for cdr in self.cdrs:
                    self.exp_data[type][exp][structure][cdr] = 0
                    self.ab_db_data[type][exp][structure][cdr] = 0

    def _parse_native_db(self, native_db, lam_kap_type):
        if not os.path.exists(native_db):
            sys.exit("native_db path not good: "+native_db)


        con = sqlite3.connect(native_db)
        structures = []
        for row in con.execute("SELECT input_tag from structures"):
            structures.append(str(row[0]))
        self._initialize_data(structures)

        for structure in structures:
            self.structures_lam_kap[structure] = lam_kap_type
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

    def _parse_ab_design_db(self, ab_db):
        if not os.path.exists(ab_db):
            sys.exit("ab_design_db path not good: "+ab_db)
        self.__parse_ab_design_db(ab_db, "lambda")
        self.__parse_ab_design_db(ab_db, "kappa")

    def __parse_ab_design_db(self, ab_db, lam_kap_type):
        def __add_data(data_dict, structure, cdr, n):
            exp = self.benchmark_info.final_name
            data_dict[exp][structure][cdr] = n

        con = sqlite3.connect(ab_db)
        for structure in self.structures_lam_kap.keys():
            if not self.structures_lam_kap[structure] == lam_kap_type: continue
            for cdr in self.cdrs:
                n_clusters = 0
                n_lengths = 0

                length = self.native_lengths[structure][cdr]
                cluster = self.native_clusters[structure][cdr]
                gene = lam_kap_type
                if cdr in self.heavy:
                    gene = "heavy"

                #Clusters
                for row in con.execute("SELECT * from cdr_data WHERE CDR=? AND fullcluster=? AND gene = ? AND datatag != ?", (cdr, cluster, gene, 'loopKeyNotInPaper')):
                    n_clusters+=1
                #__add_data(self.ab_db_totals_cluster, structure, cdr, n_clusters)
                __add_data(self.ab_db_data["cluster"], structure, cdr, n_clusters)

                #self.ab_db_totals_cluster[structure][cdr] = n_clusters

                #Lengths

                for row in con.execute("SELECT * from cdr_data WHERE CDR=? AND length=? AND gene=? AND datatag != ?", (cdr, int(length), gene, 'loopKeyNotInPaper')):
                    n_lengths+=1
                #__add_data(self.ab_db_totals_length, structure, cdr, n_lengths)
                __add_data(self.ab_db_data["length"], structure, cdr, n_lengths)


        for cdr in self.cdrs:
            gene = lam_kap_type
            if cdr in self.heavy:
                gene = "heavy"
            for row in con.execute("SELECT * from cdr_data WHERE CDR=? AND gene=? AND datatag != ?",(cdr, gene, 'loopKeyNotInPaper')):
                exp = self.benchmark_info.final_name
                for structure in self.structures_lam_kap.keys():
                    if not self.structures_lam_kap[structure] == lam_kap_type: continue
                    self.ab_db_cdr_totals[exp][structure][cdr]+=1
            #print cdr+" "+" "+gene+" "+repr(self.ab_db_cdr_totals[exp][structure][cdr])
        con.close()

    def get_filenames(self, input_dir, tag):

        if not os.path.exists(input_dir):
            sys.exit("MPI Log Dir does not exist!")
        search_name = input_dir+"/"+tag+"*"
        #print search_name
        return glob.glob(search_name)

    def add_data(self, ref_name, input_tag, cdr, length, cluster):

        #print repr(self.exp_totals[ref_name])

        #print repr(self.exp_totals.keys())

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

    def parse_single_file_old_style(self, filename, ref_name):

        if os.path.basename(filename).split('.')[-1] == "gz":
            INFILE = gzip.open(filename, 'rb')
        else:
            INFILE = open(filename, 'r')
        input_tag = ""

        for line in INFILE:
            if re.search("filling pose from PDB", line):
                input_tag = line.strip().split()[-1]
            if re.search(" DATA ", line):
                #print line
                if re.search(" DATA ", line) and re.search("false", line):
                    cluster = line.strip().split()[-2]
                else:
                    cluster = line.strip().split()[-1]

                cdr = cluster.split("-")[0]
                #print input_tag+" "+ref_name+" "+cdr+" "+cluster

                length = cluster.split("-")[1]

                self.add_data(ref_name, input_tag, cdr, length, cluster)
            else:
                continue

        INFILE.close()

    def parse_single_file_new_style(self, filename, ref_name):
        #print filename

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
                #print line
                input_tag = lineSP[data_index+3]
                cluster = lineSP[data_index+5]
                cdr = cluster.split("-")[0]
                length = cluster.split("-")[1]
                self.add_data(ref_name, input_tag, cdr, length, cluster)

        INFILE.close()
        #print self.exp_data

    def output_data(self, append_database, db_dir, db_name, drop_tables = False):

        db_path = db_dir+"/"+db_name

        if not os.path.exists(db_dir):
            os.mkdir(db_dir)

        if not append_database and os.path.exists(db_path):
            os.remove(db_path)

        #self._check_data()

        self.db = sqlite3.connect(db_path)

        self._setup_output_db(drop_tables)

        self._output_data(db_dir, "length")
        self._output_data(db_dir, "cluster")
        self.db.close()

    def _setup_output_db(self, drop_tables = False):

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

    def _output_data(self, db_dir, name):

        if not os.path.exists(db_dir):
            os.mkdir(db_dir)

        data_dict = self.exp_data[name]

        ab_db_dict = self.ab_db_data[name]
        native_dict = self.natives[name]

        rec_freq = self.recovery_parser.get_all_freq(name)
        rec_totals = self.recovery_parser.get_all_totals(name)
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

                exp = self.benchmark_info.final_name
                exp_freq = data_dict[exp][structure][cdr]
                exp_total= self.exp_totals[exp][structure][cdr]
                if exp_total == 0:

                    print("CDR Not present in data!!!??")
                    MISSING.write(name+" "+exp+" "+structure+" "+cdr+"\n")
                    #sys.exit("")
                    continue

                i+=1
                exp_perc = self._get_perc(exp_freq, exp_total)
                ab_db_freq = ab_db_dict[exp][structure][cdr]
                ab_db_total = self.ab_db_cdr_totals[exp][structure][cdr]
                ab_db_perc = self._get_perc(ab_db_freq, ab_db_total)

                self.exp_avgs[exp][structure][cdr] = exp_perc
                self.ab_db_avgs[exp][structure][cdr] = ab_db_perc


                native = native_dict[structure][cdr]
                if name == 'length':
                    native = int(native)

                exp_type, exp_group = self.get_exp_split(exp)

                data = [i, structure, cdr, exp, exp_group, exp_type, native, exp_perc, ab_db_perc, exp_freq, exp_total, ab_db_freq, ab_db_total]
                rec_data = self._get_recoveries(rec_freq, rec_totals, exp_perc, exp, structure, cdr)
                data.extend(rec_data)



                #print line

                #OUT.write(self._convert_types_to_line(data, "")+"\n")

                exec_string = "INSERT INTO all_"+name+" VALUES "+self._get_question_mark_str(len(data) )
                #print repr(exec_string)
                #print repr(data)
                self.db.execute(exec_string, data)
                self.db.commit()


        #OUT.close()
        MISSING.close()

        ####### Per CDR ##########
        self._output_cdr_data(name)

        ###### Per Exp No H3 ######
        cdrs = ["L1", "L2", "L3", "H1", "H2"]
        self._output_exp_data(name, cdrs)

        ###### Per Exp ############
        self._output_exp_data(name, self.cdrs)
        print "Complete"

    def get_exp_split(self, exp):
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

    def _output_cdr_data(self, name):
        data_dict = self.exp_data[name]
        ab_db_dict = self.ab_db_data[name]
        native_dict = self.natives[name]
        rec_freq = self.recovery_parser.get_all_freq(name)
        rec_totals = self.recovery_parser.get_all_totals(name)

        print "Outputting cdr "+name

        exp_result = self._calc_per_cdr_data(data_dict, self.exp_totals, self.exp_avgs)
        ab_db_result = self._calc_per_cdr_data(ab_db_dict, self.ab_db_cdr_totals, self.ab_db_avgs)
        rec_result = self._calc_per_cdr_data(rec_freq, rec_totals)
        #OUT = open(outname, "w")
        #OUT.write("CDR\tEXP\tAvg_Chosen_Percent\tAvg_Db_Percent\tChosen_Percent\tDb_Percent\tChosen_Freq\tGraft_Total \tDb_Freq Db_Total\t"+
        #          "top_rec top_rr top_rr_avg top_rr_log topx_rec topx_rr topx_rr_avg topx_rr_log x\n")

        i = 0
        cur = self.db.cursor()
        for cdr in self.cdrs:
            for db in self.db_info_array:
                exp = db.final_name
                i+=1
                exp_perc =  self._get_perc(exp_result[exp][cdr]["freq"], exp_result[exp][cdr]["totals"])
                ab_db_perc = self._get_perc(ab_db_result[exp][cdr]["freq"], ab_db_result[exp][cdr]["totals"])
                avg_perc = exp_result[exp][cdr]["avg_perc"]
                avg_db_perc = ab_db_result[exp][cdr]["avg_perc"]
                exp_freq = exp_result[exp][cdr]["freq"]
                exp_total = exp_result[exp][cdr]["totals"]
                db_freq = ab_db_result[exp][cdr]["freq"]
                db_total = ab_db_result[exp][cdr]["totals"]

                exp_type, exp_group = self.get_exp_split(exp)

                data = [i, cdr, exp, exp_group, exp_type, avg_perc, avg_db_perc, exp_perc, ab_db_perc, exp_freq, exp_total,db_freq,db_total ]
                #print data
                rec_data = self._get_recoveries_avg(rec_result, exp, exp_perc, cdr)

                data.extend(rec_data)

                exec_string = "INSERT INTO cdr_"+name+" VALUES"+self._get_question_mark_str(len(data))

                #print repr(exec_string)
                #print repr(data)

                cur.execute(exec_string, data)

        cur.close()
        self.db.commit()
        #OUT.close()

    def _output_exp_data(self, name, cdrs):
        data_dict = self.exp_data[name]
        ab_db_dict = self.ab_db_data[name]
        native_dict = self.natives[name]
        rec_freq = self.recovery_parser.get_all_freq(name)
        rec_totals = self.recovery_parser.get_all_totals(name)

        print "Outputting exp "+name
        exp_result = self._calc_per_exp_data(cdrs, data_dict, self.exp_totals, self.exp_avgs)
        ab_db_result = self._calc_per_exp_data(cdrs, ab_db_dict, self.ab_db_cdr_totals, self.ab_db_avgs)
        rec_result = self._calc_per_exp_data(cdrs, rec_freq, rec_totals)

        #OUT = open(outname, "w")
        #OUT.write("EXP h3_present, Avg_Chosen_Percent Avg_Db_Percent Chosen_Percent Db_Percent Chosen_Freq Graft_Total  Db_Freq Db_Total\t"+
        #    "top_rec top_rr top_rr_avg top_rr_log topx_rec topx_rr topx_rr_avg topx_rr_log x\n")

        i = 0
        cur = self.db.cursor()


        #Bad assumption here!
        h3_present = True
        if len(cdrs) != 6:
            h3_present = False

        for db in self.db_info_array:
            exp = db.final_name
            i+=1
            exp_perc =  self._get_perc(exp_result[exp]["freq"], exp_result[exp]["totals"])
            ab_db_perc = self._get_perc(ab_db_result[exp]["freq"], ab_db_result[exp]["totals"])
            avg_perc = exp_result[exp]["avg_perc"]
            avg_db_perc = ab_db_result[exp]["avg_perc"]
            exp_freq = exp_result[exp]["freq"]
            exp_total = exp_result[exp]["totals"]
            db_freq = ab_db_result[exp]["freq"]
            db_total = ab_db_result[exp]["totals"]

            exp_type, exp_group = self.get_exp_split(exp)

            data = [i, exp, exp_group, exp_type, h3_present, avg_perc, avg_db_perc, exp_perc, ab_db_perc, exp_freq, exp_total,db_freq,db_total ]
            rec_data = self._get_recoveries_avg(rec_result, exp, exp_perc)

            data.extend(rec_data)
            #OUT.write(self._convert_types_to_line(rec_data, "")+"\n")

            exec_string = "INSERT INTO exp_"+name+" VALUES"+self._get_question_mark_str(len(data))
            cur.execute(exec_string, data)

        cur.close()
        self.db.commit()
        #OUT.close()

    ####################################################################################################################

    def _convert_types_to_line(self, array_of_types, line):
        for entry in array_of_types:

            if type(entry) == float:
                line = line + " " + self._get_n_s(entry)
            elif type(entry) == int:
                line = line + " "+ repr(entry)
            elif type(entry) == None:
                line = line + " "+ repr(entry)
            elif type(entry) == str:
                line = line + " "+entry
            else:
                line = line + " "+repr(entry)

        return line

    def _get_recoveries(self, rec_freq, rec_totals, exp_perc, exp, structure, cdr):

        rec_array = []
        x = ""
        rec_exp = exp

        print rec_exp+" "+structure+" "+cdr
        print rec_freq[rec_exp][structure][cdr]
        print rec_totals[rec_exp][structure][cdr]
        rec_perc = self._get_perc(rec_freq[rec_exp][structure][cdr], rec_totals[rec_exp][structure][cdr])

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

    def _get_recoveries_avg(self, rec_data, exp, exp_perc, cdr=None):

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
            rec_perc = self._get_perc(rec_data[rec_exp][cdr]['freq'], rec_data[rec_exp][cdr]['totals'])
            rr_avg = rec_data[rec_exp][cdr]['rr_avg']
            rr_log = rec_data[rec_exp][cdr]['rr_log']
            rr = self._get_perc(rec_perc, exp_perc)/100
        else:
            if not rec_data[rec_exp].has_key('rr_avg'):
                rec_data[rec_exp]['rr_avg'] = 0
            if not rec_data[rec_exp].has_key('rr_log'):
                rec_data[rec_exp]['rr_log'] = 0

            rec_perc = self._get_perc(rec_data[rec_exp]['freq'], rec_data[rec_exp]['totals'])
            rr_avg = rec_data[rec_exp]['rr_avg']
            rr_log = rec_data[rec_exp]['rr_log']
            rr = self._get_perc(rec_perc, exp_perc)/100

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

    def _calc_per_cdr_data(self, data_dict, total_dict, avg_dict=None):
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

    def _calc_per_exp_data(self, cdrs, data_dict, total_dict, avg_dict=None):
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


    ####################################################################################################################
    def _get_perc(self, freq, total):
        """
        Get percent
        """
        freq = int(freq)
        total = float(total)

        if freq==0 and total==0:
            return 1000
        if total==0:
            sys.exit("cannot calculate percent as total is 0!")
        return freq/total *100


    def _get_s_perc(self, freq, total):
        """
        Get string of percent
        """
        return self._get_n_s(self._get_perc(freq, total))

    def _get_n_s(self, num):
        """
        Get a string for a float at .2
        """
        if num == None:
            return 'None'
        return "%.2f"%num

    ####################################################################################################################

    def _get_question_mark_str(self, length):

        s = "("
        for i in range(0, length-1):
            s = s+"?,"
        s = s + "?)"
        return s

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














