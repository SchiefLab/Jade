#!/usr/bin/env python

import os
import sys
from collections import defaultdict
from RAbD_BM import tools as bm_tools
from tools import path as path_tools
import json
import re
import numpy
from plotting.MakeFigure import MakeFigure

import matplotlib.pyplot as plot
import matplotlib as mpl
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, WeekdayLocator

import argparse

class BMInfo:
    def __init__(self, json_path):
        self.json_path = json_path
        JSON_FILE = open(self.json_path, 'r')
        self.data = json.load(JSON_FILE)
        JSON_FILE.close()

    def get_exp(self):
        return self.data["exp"]

    def get_dir(self):
        return self.data["dir"]

    def get_name_match(self):
        return self.data["name_match"]

    def get_dataset(self):
        return self.data["dataset"]

class AnalyzeMCAcceptance:
    def __init__(self):
        self.parse_args()
        self.bm_info = []
        if not self.options.jsons:
            sys.exit("JSON files need to be set for analysis to occur")

        for js in self.options.jsons:
            self.bm_info.append(BMInfo(js))

        self.all_data = defaultdict()

    def analyze(self):
        for bm_info in self.bm_info:
            bm_data = MCData(self.options.dataset_dir, bm_info)
            self.all_data[bm_info.get_exp()] = bm_data

        self.plot_cum_acceptance()
        self.plot_cycle_by_e()

    def parse_args(self):
        parser = argparse.ArgumentParser(description= "Calculates and plots monte carlo acceptance values for antibody design benchmarking.")


        ############################
        ## Required Options
        ############################
        parser.add_argument("--jsons","-j",
                            help = "Analysis JSONs to use",
                            nargs = "*")

        parser.add_argument("--data_outdir","-o",
                            help = "Path to outfile",
                            default = "data")

        parser.add_argument("--plot_outdir", "-p",
                            help = "DIR for plots",
                            default = "plots/mc_benchmarks")

        parser.add_argument("--dataset_dir",
                            help = "List of PDBIds to use for individual PDB output.",
                            default = "datasets/all")

        self.options = parser.parse_args()


    def plot_cycle_by_e(self):
        def add_data(maker, all_models, include_native, final_e, range_only = None):
                i = 0
                for cy_data in all_models:
                    i+=1
                    assert isinstance(cy_data, CycleData)


                    y = []

                    if range_only:
                        x = [z for z in range(range_only[0], range_only[1] + 1)]
                    elif include_native:
                        x = [z for z in range(-1, 152)]
                        y.append(cy_data.native_e)
                        y.append(cy_data.start_e)
                    else:
                        x = [z for z in range(0, 152)]
                        y.append(cy_data.start_e)



                    if range_only:
                        r = range(range_only[0], range_only[1] )
                        #print repr(r)
                    else:
                        r = cy_data.values.keys()

                    if final_e:
                        y.extend([cy_data.get_final_e(ii) for ii in r])
                    else:
                        y.extend([cy_data.get_mc_e(ii) for ii in r])

                    #print str(len(x))
                    #print str(len(y))

                    y.append(cy_data.get_protocol_final_e())

                    maker.add_data(x, y, "model_"+repr(i))

        for exp in self.all_data:
            root_exp_dir = path_tools.get_make_get_dirs(self.options.plot_outdir, [exp, "bm_pdbs"])
            mc_data = self.all_data[exp]
            assert isinstance(mc_data, MCData)


            for pdb_id in mc_data.data.keys():

                all_models = mc_data.get(pdb_id)
                make_figure = MakeFigure(2, 1)

                #root_pdbid_dir = path_tools.get_make_get_dirs(root_exp_dir, [pdb_id])


                #### Plot not including antigen ####

                ###### Pre MC ######

                add_data(make_figure, all_models, include_native=True, final_e=False)
                make_figure.fill_subplot("Pre MC E "+pdb_id.split("_")[0], make_figure.labels, y_axis_label= "REU" , add_legend=False)



                ###### Post MC ######
                add_data(make_figure, all_models, include_native=True, final_e=True)
                make_figure.fill_subplot("Post MC E "+pdb_id.split("_")[0], make_figure.labels, x_axis_label="Cycle", y_axis_label="REU", add_legend=False)




                make_figure.set_y_scale('symlog')
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0]+".cycle-vs-mc_e_log.with_native.pdf")


                #### Plot including antigen ####
                make_figure = MakeFigure(2, 1)

                ###### Pre MC ######
                add_data(make_figure, all_models, include_native=False, final_e=False)

                make_figure.fill_subplot("Pre MC E "+pdb_id.split("_")[0], make_figure.labels, y_axis_label= "REU" , add_legend=False)



                ###### Post MC ######
                add_data(make_figure, all_models, include_native=False, final_e=True)

                #Make first subplot
                make_figure.fill_subplot("Post MC E "+pdb_id.split("_")[0], make_figure.labels, x_axis_label="Cycle", y_axis_label="REU", add_legend=False)
                make_figure.set_y_scale('symlog')
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0] +".cycle-vs-mc_e_log.no_native.pdf")

                make_figure = MakeFigure(1, 1)
                add_data(make_figure, all_models, include_native=False, final_e=True, range_only=[50, 150])
                make_figure.fill_subplot("Post MC E "+ pdb_id.split("_")[0] , make_figure.labels, x_axis_label="Cycle", y_axis_label="REU", add_legend=False)
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0] +".cycle-vs-mc_e_From50.pdf")


    def plot_cum_acceptance(self):
        def add_data(maker, all_models):
                i = 0
                for cy_data in all_models:
                    i+=1
                    assert isinstance(cy_data, CycleData)

                    y = []

                    x = [z for z in range(1, len(cy_data.values.keys()) +1)]
                    #print repr(x)
                    r = cy_data.values.keys()
                    y.extend([cy_data.get_cum_acceptance(ii) for ii in r])

                    maker.add_data(x, y, "model_"+repr(i))


        make_exp_figure = MakeFigure(1, 1)


        x_values = None
        for bm_info in self.bm_info:
            exp = bm_info.get_exp()

            make_all_pdbids_figure = MakeFigure(1, 1)
            root_exp_dir = path_tools.get_make_get_dirs(self.options.plot_outdir, [exp, "bm_pdbs"])
            mc_data = self.all_data[exp]
            all_y_values = []

            for pdb_id in mc_data.data.keys():

                all_models = mc_data.get(pdb_id)
                make_figure = MakeFigure(1, 1)

                add_data(make_figure, all_models)
                make_figure.fill_subplot("Cumulative MC Acceptance "+pdb_id.split("_")[0], make_figure.labels, x_axis_label="Cycle", y_axis_label="Cummulative Accepts", add_legend=False)
                make_figure.add_grid()
                make_figure.save_plot(root_exp_dir+"/"+pdb_id.split("_")[0] + ".cumulative_acceptance.pdf")
                if not x_values:
                    x_values = make_figure.get_x_data(make_figure.labels[0])

                make_all_pdbids_figure.add_data(x_values, numpy.mean(make_figure.get_y_as_list(make_figure.labels), axis = 0), pdb_id)
                all_y_values.extend(make_figure.get_y_as_list(make_figure.labels))

            make_exp_figure.add_data(x_values, numpy.mean(all_y_values, axis=0), exp)

            make_all_pdbids_figure.fill_subplot("Cumulative MC Acceptance", make_all_pdbids_figure.labels, x_axis_label="Cycle", y_axis_label="Avg Cumulative Accepts")
            make_all_pdbids_figure.save_plot(self.options.plot_outdir+"/"+"all_pdbs_avg_cumulative_acceptance_"+exp+".pdf")

        make_exp_figure.fill_subplot("Cumulative MC Acceptance", make_exp_figure.labels, x_axis_label="Cycle", y_axis_label="Avg Cumulative Accepts", add_legend=True)
        make_exp_figure.save_plot(self.options.plot_outdir+"/"+"avg_cumulative_acceptance_"+"_".join([exp for exp in self.all_data])+".pdf")

class MCData:
    def __init__(self, dataset_dir, bm_info):
        if isinstance(bm_info, BMInfo): pass
        self.dataset_dir = dataset_dir
        self.data = defaultdict(lambda: [])
        self.bm_info = bm_info

        self.data = defaultdict()

        self.native_pdblist = self.dataset_dir+"/"+self.bm_info.get_dataset()+".PDBLIST.txt"


        if not os.path.exists(self.native_pdblist):
            sys.exit("Native PDBList path, "+self.native_pdblist+" does not exist.  Make sure to concatonate both lambda and kappa abs")

        FILE = open(self.native_pdblist)
        self.pdb_ids = [line.strip() for line in FILE]
        FILE.close()

        self._read_pdb_files()
        self._calculate_cum_accepts()
        self._calculate_delta_e()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.get(item)

    def __setitem__(self, key, value):
        self.data[key] = value

    def get(self, pdb_id):
        return self.data[pdb_id]

    def get_pdb_id_of_path(self, file_path):
        for pdb_id in self.pdb_ids:
            pdb_id = pdb_id.split(".")[0]
            if re.search(pdb_id, file_path):
                return pdb_id

        sys.exit("unknown PDB ID of resultant benchmark PDB")

    def _calculate_cum_accepts(self):
        for pdb_id in self.data:
            print pdb_id
            for cy_data in self.data[pdb_id]:
                cum = 0
                for n in cy_data.values:
                    print repr(n)
                    cum = cum + cy_data.get_acceptance(n)
                    cy_data.set_cum_acceptance(n, cum)

    def _calculate_delta_e(self):
        for pdb_id in self.data:
            for cy_data in self.data[pdb_id]:
                for n in cy_data.values:
                    if n == 1:
                        cy_data.set_delta_e(n, 0)
                    else:
                        cy_data.set_delta_e(n, cy_data.get_final_e(n) - cy_data.get_final_e(n-1))

    def _read_pdb_files(self):
        file_paths = bm_tools.get_pdb_paths("decoys/"+self.bm_info.get_dir(), self.bm_info.get_exp(), self.bm_info.get_name_match())
        for file_path in file_paths:
            self._parse_file_data(file_path)

    def _parse_file_data(self, file_path):
        print "Reading "+file_path
        INFILE = path_tools.open_file(file_path)
        pdb_id = self.get_pdb_id_of_path(file_path)
        cy_data = CycleData(file_path, pdb_id)

        for line in INFILE:
            line = line.strip()
            lineSP = line.split()

            if not lineSP:                continue
            if not lineSP[0] == "ACCEPT": continue
            if not lineSP[1] == "LOG":    continue

            if lineSP[3] == "FINAL":
                if lineSP[4] == "END":
                    cy_data.set_protocol_final_e(lineSP[5])
                else:
                    cy_data.set_final_energy(lineSP[4], lineSP[5])
            else:
                if lineSP[3] == "0":
                    cy_data.set_protocol_start_e(lineSP[4])
                elif lineSP[3] == "-1":
                    cy_data.set_protocol_native_e(lineSP[4])
                else:
                    cy_data.set_mc_data(lineSP[3], lineSP[4], lineSP[5])

        INFILE.close()
        if not self.data.has_key(pdb_id):
            self.data[pdb_id] = []

        #print cy_data.values.keys()
        #for i in cy_data.values:
        #    print repr(i)+" "+repr(cy_data[i])
        #    print repr(i)+" "+repr(cy_data.get_final_e(i))
        self.data[pdb_id].append(cy_data)

class CycleData:
    def __init__(self, pdb_path, pdb_id):
        self.pdb_path = pdb_path
        self.pdb_id = pdb_id
        self.values = defaultdict(dict)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, item):
        return self.get_cycle(int(item))

    def __repr__(self):
        return repr(self.values)

    def get_cycle(self, n):
        return self.values[int(n)]


    def get_acceptance(self, n):
        return self[n]["mc_accept"]

    def get_mc_e(self, n):
        return self[n]["mc_e"]

    def get_final_e(self,n):
        return self[n]["final_e"]

    def get_protocol_final_e(self):
        return self.final_e

    def get_native_e(self):
        return self.get_protocol_final_e()

    def set_protocol_final_e(self, e):
        self.final_e = float(e)

    def set_protocol_start_e(self, e):
        self.start_e = float(e)

    def set_protocol_native_e(self, e):
        self.native_e = float(e)


    def set_final_energy(self, n, e):
        n = int(n)
        e = float(e)
        self.values[n]["final_e"] = e

    def set_mc_data(self, n, e, acceptance):
        n = int(n)
        e = float(e)
        #acceptance = bool(acceptance)
        self.values[n]["mc_e"] = e
        self.values[n]["mc_accept"] = int(acceptance)

    def set_cum_acceptance(self, n, cum_accepts):
        self.values[n]["cum_accepts"] = cum_accepts
        print "Setting "+repr(n)+" "+repr(cum_accepts)

    def get_cum_acceptance(self, n):
        return self.values[n]["cum_accepts"]

    def set_delta_e(self, n, delta_e):
        self.values[n]["delta_e"] = delta_e

    def get_delta_e(self, n):
        return self.values[n]["delta_e"]

if __name__ == "__main__":

    analyzer = AnalyzeMCAcceptance()
    analyzer.analyze()


