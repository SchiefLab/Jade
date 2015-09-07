#!/usr/bin/env python

import sqlite3
import os
import sys
import re
import numpy
import matplotlib.pyplot as plot
from pylab import polyfit, poly1d

from collections import defaultdict

class NativeCDRData:
    def __init__(self, datatype, native_path, data_table = "cdr_metrics"):
        self.native_path = native_path
        self.db = sqlite3.connect(native_path)
        self.datatype = datatype
        self.data_table = data_table

        self.data = defaultdict()
        self.setup_data(datatype)

    def setup_data(self, datatype):
        cur = self.db.cursor()
        query = "SELECT input_tag, CDR, "+datatype+" FROM "+self.data_table+", structures WHERE structures.struct_id = "+self.data_table+".struct_id"
        print query
        for row in cur.execute(query):
            pdbid = row[0]; cdr = row[1]; data = row[2]

            if not self.data.has_key(pdbid):
                self.data[pdbid] = defaultdict()
            self.data[pdbid][cdr] = float(data)

    def get_all_data(self):
        return self.data

    def get_data(self, pdbid, cdr):
        return self.data[pdbid][cdr]


class RecoveryCDRData:
    def __init__(self, db_paths, type = "length"):
        self.type = type
        self.table = "all_"+type
        self.db_paths = db_paths
        self.recovery_data = defaultdict()
        self.risk_ratio_data = defaultdict()
        self.setup_data()

    def setup_data(self):
        self.connections = []
        for db_path in self.db_paths:
            nameSP = os.path.basename(db_path).split(".")
            #exp_name = nameSP[0]+"_"+nameSP[1]+"_"+nameSP[3]
            exp_name = nameSP[1]+"_"+nameSP[3]
            print exp_name
            self.recovery_data[exp_name] = defaultdict()
            self.risk_ratio_data[exp_name] = defaultdict()

            con = sqlite3.connect(db_path)
            cur = con.cursor()
            query = "SELECT native, CDR, top_rec, top_rr FROM "+self.table
            print query
            for row in cur.execute(query):
                pdb = row[0]; cdr = row[1]; rec = row[2]; rr = row[3]
                if rr == None:
                    rr = 0
                #print repr(row)
                if not self.recovery_data[exp_name].has_key(pdb):
                    self.recovery_data[exp_name][pdb] = defaultdict()
                    self.risk_ratio_data[exp_name][pdb] = defaultdict()
                self.recovery_data[exp_name][pdb][cdr] = float(rec)
                self.risk_ratio_data[exp_name][pdb][cdr] = float(rr)


class PlotData:
    def __init__(self,native_data, rec_data):
        if isinstance(native_data, NativeCDRData):pass
        if isinstance(rec_data, RecoveryCDRData): pass

        self.native_data = native_data
        self.rec_data = rec_data

    def get_xy_of_exp(self, exp, rec = True, skip_H3 = True):
        x = []
        y = []
        for pdb in self.native_data.data.keys():
            for cdr in self.native_data.data[pdb].keys():
                if cdr == "H3" and skip_H3: continue;
                pdb_cdr = pdb.split(".")[0]+"_"+cdr
                x0 = self.native_data.data[pdb][cdr]

                if rec:
                    y0 = self.rec_data.recovery_data[exp][pdb][cdr]
                else:
                    y0 = self.rec_data.risk_ratio_data[exp][pdb][cdr]

                x.append(x0)
                y.append(y0)
        return x, y

    def plot_data(self, outname, rec = True):
        for exp in self.rec_data.recovery_data.keys():
            x, y = self.get_xy_of_exp(exp, rec)

            plot.plot(x, y, 'o', label = exp)


            fitx = polyfit(x,y,1)
            fit_fnx=poly1d(fitx)

            #plot.plot(x,fit_fnx(x),'b-')  # Regr line (x,fit_fn(x))
        ylab = self.rec_data.type+" "
        if rec:
            ylab = ylab+"Recovery"
        else:
            ylab = ylab+"Risk Ratio"

        plot.ylabel(ylab)
        plot.xlabel(self.native_data.datatype)

        plot.legend()
        plot.show()

if __name__ == "__main__":
    outname = "testing"
    type = "ag_ab_contacts_total"
    #type = "ag_ab_dSASA"
    native_db = sys.argv[1]
    dbs = sys.argv[2:]

    native_data = NativeCDRData(type, native_db)
    rec_data = RecoveryCDRData(dbs, "cluster")
    plotter = PlotData(native_data, rec_data)
    plotter.plot_data(outname, False)