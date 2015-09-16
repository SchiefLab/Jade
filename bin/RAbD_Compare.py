#!/usr/bin/env python

#Yes, the imports are all over.  Basically everything I've coded over the past few years.
import sys
import os
import re
import glob
import sqlite3
import numpy
import math
import multiprocessing
import shutil
import copy

from optparse import OptionParser
from collections import defaultdict

#from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB import PDBIO

#TkInter
from Tkinter import *
from tkFont import *
import tkFileDialog
import tkSimpleDialog
import tkMessageBox

p = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath

#PyIgD
import bin.RAbD_Analyze as analyze_strat
from sequence.ClustalRunner import *
from tools.StatementCreator import *
from antibody.cdr_data.CDRData import *
from antibody.cdr_data.CDRDataTypes import *
from antibody.decoy_data.DecoyData import *
from antibody.decoy_data.DecoyDataTypes import *
from tools.filters.DataFilter import *
from tools.filters.DataFilters import *
from tools.filters.FilterSettings import *

from RAbD.window_modules.FilterSettingsWindow import *
from sequence import fasta


#Rosetta Tools
import create_features_json as json_creator





p = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(p); #Allows all modules to use all other modules, without needing to update pythonpath

from rosetta import *

#Rosetta is only used for set_native_data_from_rosetta function for clusters.
rosetta.init(" -ignore_unrecognized_res -ignore_zero_occupancy false -ex1 -ex2 -use_input_sc"
             " -antibody:numbering_scheme AHO_Scheme "
             " -antibody:cdr_definition North")

#We assume that the names of each decoy is different.
#  This may be a problem if people do not use out:prefix or out:suffix to name their decoys...

### Setup Enums ###
#length = 1; cluster = 2; sequence = 3;



def main():

    main_dir = ""
    out_dir_name = ""
    strategies = []

    if len(sys.argv) >= 2: main_dir      =  sys.argv[1]
    if len(sys.argv) >= 3: out_dir_name  =  sys.argv[2]
    if len(sys.argv) >= 4: strategies    =  sys.argv[3:]


    GUI = CompareAntibodyDesignStrategies_GUI(Tk(), main_dir, out_dir_name, strategies)
    GUI.run()

class Listbox(Listbox):
    def autowidth(self,maxwidth, list = None):
        f = Font(font=self.cget("font"))
        pixels = 0
        if not list:
            for item in self.get(0, "end"):
                pixels = max(pixels, f.measure(item))
        else:
            for item in list:
                pixels = max(pixels, f.measure(item))

        # bump listbox size until all entries fit
        pixels = pixels + 10
        width = int(self.cget("width"))
        for w in range(0, maxwidth+1, 5):
            if self.winfo_reqwidth() >= pixels:
                break
            self.config(width=width+w)

class CompareAntibodyDesignStrategies_GUI:
    def __init__(self, main , main_dir = "", out_dir_name = "", strategies = []):

        self._tk_ = main
        self._tk_.title("PyIgDesign Compare")
        self.compare_designs = CompareAntibodyDesignStrategies(main_dir, out_dir_name, strategies)
        if self.compare_designs.main_dir and not self.compare_designs.strategies:
            self.compare_designs.set_strategies_from_databases()

        self.current_dir = os.path.split(os.path.abspath(__file__))[0]

        self.clustal_procs = IntVar(value = multiprocessing.cpu_count())
        self.clustal_output_format = StringVar(value = "clu")
        self.native_path = None
        self.clustal_output_formats = ['fasta', 'clustal', 'msf', 'phylip', 'selex', 'stockholm', 'vienna',
                                        'a2m', 'fa', 'clu', 'phy', 'st', 'vie']
        self.base_clustal_options = " -v --auto --force"

        ###Sub Windows ####
        self.filter_settings_window = FilterSettingsWindow(self.compare_designs.filter_settings)

    def run(self):
        self.set_tk()
        self.set_menu()
        self.sho_tk()
        self._tk_.mainloop()

    def set_tk(self):
        self.main_dir_entry = Entry(self._tk_, textvariable = self.compare_designs.main_dir, justify = CENTER)
        self.out_dir_entry = Entry(self._tk_, textvariable = self.compare_designs.out_dir_name, justify = CENTER)

        self.main_dir_label = Label(self._tk_, text = "Main Analysis Directory", justify = CENTER)
        self.out_dir_label = Label(self._tk_, text = "Root Name", justify = CENTER)

        self.all_strategies_listbox = Listbox(self._tk_)
        self.current_strategies_listbox = Listbox(self._tk_)

        self.ab_features_button = Button(self._tk_, text = "Run Antibody Features", command = lambda: self.run_features_reporter("antibody"), justify = CENTER)
        self.clus_features_button = Button(self._tk_, text = "Run Cluster Features", command = lambda: self.run_features_reporter("cluster"), justify = CENTER)

        self.ab_features_options_label = Label(self._tk_, text = "Antibody Features Options", justify = CENTER)

        self.normal_hbond_radio = Radiobutton(self._tk_, text = "All Hbond R Scripts", variable = self.compare_designs.features_hbond_set, value = 0)
        self.min_hbond_radio = Radiobutton(self._tk_, text = "Minimal Hbond R Scripts", variable = self.compare_designs.features_hbond_set, value = 1)
        self.no_hbond_radio = Radiobutton(self._tk_, text = "No Hbond R Scripts", variable = self.compare_designs.features_hbond_set, value = 2)
        #self.separator = Separator(self._tk_, orient = HORIZONTAL)

    def sho_tk(self, r = 0, c = 0):

        self.main_dir_label.grid(row = r+0, column = c+0, columnspan = 2, sticky = W+E, pady = 7)
        self.main_dir_entry.grid( row = r+1, column = c+0, columnspan = 2, sticky = W+E, padx = 5)

        self.all_strategies_listbox.grid(row = r+2, column = c+0, padx = 6, pady = 10)
        self.current_strategies_listbox.grid(row = r+2, column = c+1, padx = 6, pady = 10)



        #self.separator.grid(row = r+3, column = c, columnspan = 2, sticky = W+E, pady = 15)

        self.out_dir_label.grid(row = r+5, column = c+0, columnspan = 1, pady = 5)
        self.out_dir_entry.grid(row = r+5, column = c+1, columnspan = 1, padx = 5, pady = 5)

        self.all_strategies_listbox.bind("<Double-Button-1>", lambda event: self.add_to_current(self.all_strategies_listbox, self.current_strategies_listbox))
        self.all_strategies_listbox.bind("<Button-2>", lambda event: self.show_strat_items())

        self.ab_features_button.grid(row = r+6, column = c+0, columnspan = 1, pady = 3, sticky = W+E)
        self.clus_features_button.grid(row = r+6, column = c+1, columnspan = 1, pady = 3, sticky = W+E)

        #self.ab_features_options_label.grid(row = r+7, column = c+0, columnspan = 2, pady = 3, sticky = W+E)

        self.normal_hbond_radio.grid(row = r+8, column = c+0, columnspan = 1, pady = 1, sticky = W)
        self.min_hbond_radio.grid(row = r+9, column = c+0, columnspan = 1,  pady = 1, sticky = W)
        self.no_hbond_radio.grid(row = r+10, column = c+0, columnspan = 1, pady = 1, sticky = W)

        self.current_strategies_listbox.bind("<Double-Button-1>", lambda event: self.delete_current(self.current_strategies_listbox))


        self.populate_all_strategies()

    def set_menu(self):


        self.main_menu = Menu(self._tk_)

        ## File Menu ##
        self.file_menu = Menu(self.main_menu, tearoff=0)

        self.file_menu.add_checkbutton(label = "Camelid Antibody", variable = self.compare_designs.is_camelid)
        #self.dtypes_menu = Menu(self.main_menu, tearoff = 0)
        #self.dtypes_menu.add_checkbutton(label = "Group by dG", variable = self.compare_designs.group_dG)

        self.file_menu.add_command(label = "Filter Models", command = lambda: self.filter_settings_window.setup_sho_gui(Toplevel(self._tk_)))
        self.file_menu.add_command(label = "Read Strategies from Main", command = lambda: self.read_from_main_set_strategies())
        self.file_menu.add_command(label = "Add Strategy", command = lambda: self.add_main_strategy())
        self.file_menu.add_separator()
        self.file_menu.add_command(label = "Set Reference Native", command = lambda :self.set_native_path())
        self.file_menu.add_command(label = "Set Reference Database", command = lambda: self.set_reference_db())
        #self.file_menu.add_command(label = "Set Scorefunction", command = lambda: self.set_scorefunction())
        self.file_menu.add_command(label = "Set top N", command = lambda: self.set_top_n())
        self.file_menu.add_command(label = "Set top N For Combined", command = lambda: self.set_top_n_combined())
        self.file_menu.add_separator()
        self.file_menu.add_command(label = "Set Strategy DIR as Working Directory")
        self.file_menu.add_separator()
        self.file_menu.add_checkbutton(label = "Reload Query Data", variable = self.compare_designs.reload_scores)
        self.file_menu.add_separator()

        for name in sorted(self.compare_designs.scores_on.keys()):
            self.file_menu.add_checkbutton(label = name, variable = self.compare_designs.scores_on[name])

        self.main_menu.add_cascade(label = "File", menu = self.file_menu)


        ## Score Menu ##
        self.score_menu = Menu(self.main_menu, tearoff=0)
        self.score_menu.add_command(label = "Output All Score Stats", command = lambda: self.run_output_all_stats())
        self.score_menu.add_command(label = "Output Per-Strategy Score Data", command = lambda: self.run_output_strategy_stats())
        self.pymol_menu = Menu(self.main_menu, tearoff = 0)
        self.pymol_menu.add_command(label = "Top Per Strategy", command = lambda: self.compare_designs.copy_top_strategy(self.native_path))
        self.pymol_menu.add_command(label = "Top Combined", command = lambda: self.compare_designs.copy_top_combined(self.native_path))
        self.pymol_menu.add_command(label = "All Models", command = lambda: self.compare_designs.copy_all_models(self.native_path))

        self.score_menu.add_cascade(label = "Create PyMol Sessions", menu = self.pymol_menu)
        self.main_menu.add_cascade(label = "PyMol", menu = self.score_menu)

        ## Clustal Menu ##
        self.clustal_menu = Menu(self.main_menu, tearoff = 0)
        self.clustal_menu.add_command(label = "Set Max Processors", command = lambda: self.set_max_clustal_procs())
        self.clustal_menu.add_command(label = "Set Output format", command = lambda: self.set_clustal_output_format())
        self.clustal_menu.add_command(label = "Set Soft Wrap", command = lambda: self.set_clustal_soft_wrap())
        self.clustal_menu.add_separator()
        self.clustal_menu.add_command(label = "Run Clustal Omega on Strategies", command = lambda: self.run_clustal_on_strategies())
        self.clustal_menu.add_command(label = "Run Clustal Omega on Top Combined Decoys", command = lambda: self.run_clustal_on_top_combined())
        self.clustal_menu.add_command(label = "Run Clustal Omega on ALL Combined Decoys", command = lambda: self.run_clustal_on_all_combined())
        self.main_menu.add_cascade(label = "Clustal", menu = self.clustal_menu)

        ## Alignment ##
        self.alignment_menu = Menu(self.main_menu, tearoff = 0)

        self.alignment_menu.add_command(label = "Output Length Alignments", command = lambda: self.run_cdr_type_alignment("length"))
        self.alignment_menu.add_command(label = "Output Cluster Alignments", command = lambda: self.run_cdr_type_alignment("cluster"))
        self.alignment_menu.add_command(label = "Output CDR Sequence Alignments", command = lambda: self.run_cdr_type_alignment("aligned_sequence"))
        self.alignment_menu.add_separator()
        self.main_menu.add_cascade(label = "Alignment", menu = self.alignment_menu)


        ## Recovery ##
        self.recovery_menu = Menu(self.main_menu, tearoff = 0)
        self.recovery_menu.add_command(label = "Output Length Recovery", command = lambda: self.run_cdr_type_recovery("length"))
        self.recovery_menu.add_command(label = "Output Cluster Recovery", command = lambda: self.run_cdr_type_recovery("cluster"))
        self.main_menu.add_cascade(label = "Recovery", menu = self.recovery_menu)

        ## Enrichment ##
        self.enrichment_menu = Menu(self.main_menu, tearoff = 0)
        self.enrichment_menu.add_command(label = "Output Length Enrichments")
        self.enrichment_menu.add_command(label = "Output Cluster Enrichments")

        self.main_menu.add_cascade(label = "Enrichment", menu = self.enrichment_menu)

        self.subset_menu = Menu(self.main_menu, tearoff = 0)

        #Have to do this manually:
        #for score_name in :
        #    x = copy.deepcopy(score_name)
        #    self.subset_menu.add_command(label = "Create DB of Top Subset: "+score_name, command = lambda: self.create_subset_databases(x))

        score_names = self.compare_designs._get_score_names()
        self.subset_menu.add_command(label = "Create DB of Top Subset: "+score_names[0], command = lambda: self.create_subset_databases(score_names[0]))
        self.subset_menu.add_command(label = "Create DB of Top Subset: "+score_names[1], command = lambda: self.create_subset_databases(score_names[1]))
        self.subset_menu.add_command(label = "Create DB of Top Subset: "+score_names[2], command = lambda: self.create_subset_databases(score_names[2]))
        self.subset_menu.add_command(label = "Create DB of Top Subset: "+"Top N dG of Top Total", command = lambda: self.create_subset_databases(score_names[3]))

        self.main_menu.add_cascade(label = "Features Subsets", menu = self.subset_menu)

        self._tk_.config(menu = self.main_menu)



        ########### Callbacks ############

    def show_strat_items(self):
        item = self.all_strategies_listbox.get(self.all_strategies_listbox.curselection())
        items = glob.glob(self.compare_designs.main_dir.get()+"/*"+item+"*")
        #for i in items:
            #print i

        if os.path.exists(self.compare_designs.main_dir.get()+"/databases"):
            print "\n Databases:"
            dbs = glob.glob(self.compare_designs.main_dir.get()+"/databases/*"+item+"*")
            for db in dbs:
                print db

    def add_to_current(self, from_listbox, to_listbox):
        item = from_listbox.get(from_listbox.curselection())
        to_listbox.insert(END, item)
        strategies = self.get_full_strategy_list()
        self.compare_designs.set_strategies(strategies)

    def delete_current(self, listbox):
        listbox.delete(listbox.curselection())
        strategies = self.get_full_strategy_list()
        self.compare_designs.set_strategies(strategies)

    def populate_all_strategies(self):
        self.all_strategies_listbox.delete(0, END)
        for strategy in self.compare_designs.strategies:
            self.all_strategies_listbox.insert(END, strategy)

        self.all_strategies_listbox.autowidth(100)
        self.current_strategies_listbox.autowidth(100, self.compare_designs.strategies)


        ######### Auxilliary Functions ###########

    def add_main_strategy(self):
        strategy_name = tkSimpleDialog.askstring(title = "Strategy", prompt = "Strategy Name")
        if not strategy_name:
            return

        self.all_strategies_listbox.insert(END, strategy_name)

    def read_from_main_set_strategies(self):
        self.compare_designs.strategies = []

        if not self.compare_designs.main_dir.get():

            strat_dir = tkFileDialog.askdirectory(initialdir = self.current_dir, title = "Strategy Analysis Directory")
            if not strat_dir:
                return
            self.current_dir = strat_dir
            self.main_dir.set(strat_dir)
        self.compare_designs.set_strategies_from_databases()
        self.populate_all_strategies()

        strategies = self.get_full_strategy_list()
        self.compare_designs.set_strategies(strategies)

    def get_full_strategy_list(self):
        strategies = self.current_strategies_listbox.get(0, END)
        return strategies

    def set_reference_db(self):
        d = tkFileDialog.askopenfilename(title = "Reference DB", initialdir = self.current_dir)
        if not d: return
        self.current_dir = os.path.dirname(d)
        self.compare_designs.reference_db.set(d)

    def set_top_n(self):
        top = tkSimpleDialog.askinteger(title = "Top N", prompt = "Number of top scoring", initialvalue=self.compare_designs.top_n.get())
        if not top: return
        self.compare_designs.top_n.set(top)

    def set_top_n_combined(self):
        top = tkSimpleDialog.askinteger(title = "Top N Combined", prompt = "Number of top scoring Combined", initialvalue=self.compare_designs.top_n_combined.get())
        if not top: return
        self.compare_designs.top_n_combined.set(top)

    def set_max_clustal_procs(self):
        max = tkSimpleDialog.askinteger(title = "Max P", prompt = "Max NP.  Clustal by default uses all.", initialvalue=self.clustal_procs.get())
        if not max: return
        self.clustal_procs.set(max)

    def set_clustal_output_format(self):
        f = tkSimpleDialog.askstring(title = "Clustal output format", initialvalue = self.clustal_output_format.get())
        if not f:
            return
        if not f in self.clustal_output_formats:
            print "Format "+f+" not recognized.  Available formats are: \n"+repr(self.clustal_output_formats)
            return

        self.clustal_output_format.set(f)

    def set_clustal_soft_wrap(self):
        wrap = tkSimpleDialog.askinteger(title = "Wrap", prompt = "Set Soft Wrap", initialvalue = self.compare_designs.clustal_soft_wrap.get())
        if not wrap:
            return
        self.compare_designs.clustal_soft_wrap.set(wrap)

    #def set_scorefunction(self):
        #score = tkSimpleDialog.askstring(title="Score", prompt = "Set Scorefunction", initialvalue = self.compare_designs.scorefunction.get())
        #if not score:
            #return
        #self.compare_designs.scorefunction.set(score)

    def set_native_path(self):
        native_path = tkFileDialog.askopenfilename(title = "Native path", initialdir = self.current_dir)
        if not native_path:
            return
        else:
            self.current_dir = os.path.dirname(native_path)
            self.native_path = native_path


        ######## Main Analysis ############

    def run_features_reporter(self, type):
        strategies = self.get_full_strategy_list()
        if len(strategies) == 0:
            print "No strategies selected..."
            return

        self.compare_designs.set_strategies(strategies)
        self.compare_designs.run_features(type)

    def run_output_all_stats(self):
        #do_hb_analysis = tkMessageBox.askyesno(title = "Hbonds", message="Query Hbonds for ALL stats?  May take a few minutes...")
        #self.compare_designs.query_hbonds.set(do_hb_analysis)
        self.compare_designs.output_stats()

    def run_output_strategy_stats(self):
        #do_hb_analysis = tkMessageBox.askyesno(title = "Hbonds", message="Query Hbonds for ALL stats?  May take a few minutes...")
        #self.compare_designs.query_hbonds.set(do_hb_analysis)
        self.compare_designs.output_score_extra_stats()

    def run_copy_all(self):
        self.compare_designs.copy_all_models(self.native_path)

    def run_clustal_on_strategies(self):

        extra_options = tkSimpleDialog.askstring(title = "Extra Options", prompt="Clustal Extra Options", initialvalue = self.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega_on_strategies(self.clustal_procs.get(), self.clustal_output_format.get(), extra_options, self.native_path)

    def run_clustal_on_all_combined(self):

        extra_options = tkSimpleDialog.askstring(title = "Extra Options", prompt = "Clustal Extra Options", initialvalue = self.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega_on_all_combined(self.compare_designs.out_dir_name.get(), self.clustal_procs.get(), self.clustal_output_format.get(), extra_options, self.native_path)

    def run_clustal_on_top_combined(self):

        extra_options = tkSimpleDialog.askstring(title = "Extra Options", prompt = "Clustal Extra Options", initialvalue = self.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega_on_top_combined(self.clustal_procs.get(), self.clustal_output_format.get(), extra_options, self.native_path)

    def run_cdr_type_alignment(self, alignment_type):

        self.compare_designs.output_len_or_clus_alignment(alignment_type, 'antibody', self.compare_designs.is_camelid.get(), self.native_path)

    def run_cdr_type_recovery(self, alignment_type):

        self.compare_designs.output_len_or_clus_recovery(alignment_type, 'antibody', False, self.native_path)

    def create_subset_databases(self, score_name):
        rosetta_extension = tkSimpleDialog.askstring(title= "Rosetta Extension", prompt = "Please set the Rosetta Extension", initialvalue=self.compare_designs.rosetta_extension.get())
        if not rosetta_extension:
            return
        self.compare_designs.rosetta_extension.set(rosetta_extension)

        prefix = tkSimpleDialog.askstring(title = "Prefix", prompt = "Please set the prefix that will be used for the new databases", initialvalue = score_name)
        if not prefix:
            print "The prefix needs to be set"
            return

        self.compare_designs.create_score_subset_database(score_name, prefix)
        #self.read_from_main_set_strategies()

class CompareAntibodyDesignStrategies:
    """
    Class mainly for comparing different Antibody Design strategies using our Features Databases.
    """
    def __init__(self, main_analysis_dir, out_dir_name, strategies = []):

        self.main_dir = StringVar(value = main_analysis_dir)
        self.out_dir_name = StringVar(value = out_dir_name)
        self.reference_db = StringVar()
        self.clustal_soft_wrap = IntVar(value = 100)
        self.reload_scores = IntVar(value = 1)

        self.top_n = IntVar(value=10)
        self.top_n_combined = IntVar(value = 20)

        self.strategies = strategies
        self.strategy_scorefxns = defaultdict()

        self.filter_settings = FilterSettings()

        self.features_hbond_set = IntVar()
        self.features_hbond_sets = ["", "_min_hbond_analysis", "_no_hbond_analysis"]

        self.is_camelid = IntVar(); self.is_camelid.set(0)
        self.top_total_percent = IntVar(); self.top_total_percent.set(10)

        total_scores = TotalDecoyData(); dg_scores = dGDecoyData(); dsasa_scores = dSASADecoyData(); top10_by_10 = dGTotalScoreSubset()

        #hbond_counts = InterfaceHbondCountDecoyData();
        #hbond_energies = InterfaceHbondEnergyDecoyData();

        #self.scores = [total_scores, dg_scores, dsasa_scores, top10_by_10, hbond_counts, hbond_energies]

        self.scores = [total_scores, dg_scores, dsasa_scores, top10_by_10]

        self.rosetta_extension = StringVar(); self.rosetta_extension.set("linuxclangrelease")

        self.score_names = self._get_score_names()
        self.scores_on = defaultdict()
        for score_name in self.score_names:
            self.scores_on[score_name] = IntVar(value = 0)

        self.scores_on["dG"].set(1)
        self.scores_on["dG_top_Ptotal"].set(0)

        self.query_hbonds = IntVar(value = 0)

    def _get_score_names(self):

        names = []
        for score in self.scores:
            names.append(score.name)
        return names

    def _setup_outdir(self, subdirs = [], use_out_dir_name = True):
        """
        Sets up the main output dir in the main_analys_dir, and any subdirectories such as 'decoys' or decoys/combined_3
        Returns the final output directory
        """

        filters = self._setup_filters()

        if self.out_dir_name.get() and use_out_dir_name and filters:
            outdir = self.main_dir.get()+"/"+self.out_dir_name.get()+"_"+self.filter_settings.name.get()
        elif self.out_dir_name.get() and use_out_dir_name:
            outdir = self.main_dir.get()+"/"+self.out_dir_name.get()
        elif not self.main_dir.get():
            sys.exit("Main StrategyAnalysis DIR must exist.  Please use analyze_antibody_design_strategy to populate")
        elif filters:
            outdir = self.main_dir.get()+"/"+self.filter_settings.name.get()
        else:
            outdir = self.main_dir.get()
        if not os.path.exists(outdir): os.mkdir(outdir)

        for subdir in subdirs:
            if not subdir: continue
            outdir = outdir+"/"+subdir
            if not os.path.exists(outdir): os.mkdir(outdir)

        return outdir

    def _setup_outdir_individual(self, subdirs):

        return self._setup_outdir(["individual_data", self.out_dir_name.get()] + subdirs, False)

    def _setup_outdir_combined(self, subdirs):
        return self._setup_outdir(["combined_data", self.out_dir_name.get()]+subdirs, False)

    def _setup_scores(self, features_type = "antibody", use_all = False):
        """
        Setup the Score Classes.  If not use_all, will use only use those set.
        """

        scores = self.scores
        score_subset = []
        query_hbonds = self.query_hbonds.get()

        for score in self.scores:
            if score.name == "hbond_count" or score.name == "hbond_energy":

                if self.scores_on[score.name].get():
                    query_hbonds = True

        #Quickly analyze scores with the same settings...
        if not self.reload_scores.get():
            for score_class in scores:
                if self.scores_on[score_class.name].get() or use_all:
                    score_subset.append(copy.deepcopy(score_class))
            return score_subset




        hb_loader = InterfaceHBondDecoyDataLoader();


        filters = self._setup_filters()

        if self.is_camelid.get():
            for score in scores:
                score.set_interface('H_A')

        if filters:
            scores.append(CombinedStrDecoyData(filters, self.filter_settings.name.get()))

        for strategy in self.strategies:
            db_path = self.get_db_path(strategy, features_type)
            if not os.path.exists(db_path):
                sys.exit("DB path does not exist.  Please Run Features reporter with StructureScores and ScoreTypes for this strategy\n"+db_path)
            print db_path
            con = sqlite3.connect(db_path)

            if query_hbonds:
                if filters:
                    hb_loader.add_filters(filters, self.filter_settings.name.get())
                hb_loader.add_data(strategy, con)

            for score_class in scores:



                if filters:
                    score_class.add_filters(filters, self.filter_settings.name.get())

                #We set these up later
                if score_class.name == "hbond_count" or score_class.name == "hbond_energy":
                    continue

                if score_class.name == "dG_top_Ptotal":
                    score_class.add_data(strategy, con, self.top_total_percent.get())
                else:
                    score_class.add_data(strategy, con)

        #Setup each Hbond Data class
        if query_hbonds:
            for score_class in scores:
                if score_class.name == "hbond_count" or score_class.name == "hbond_energy":
                    score_class.setup_from_loader(hb_loader)

        for score_class in scores:
            if self.scores_on[score_class.name].get() or use_all:
                score_subset.append(copy.deepcopy(score_class))

        self.scores = scores

        return score_subset

    def _get_score(self, score_name):
        for score in self.scores:
            if score.name == score_name:
                return score
        return None


    def _setup_filters(self):
        filters = []

        if self.filter_settings.h3_filter.get():
            filters.append(H3ExtendedFilter())

        if self.filter_settings.extra_required_where:
            custom_filter = DataFilter("custom_filter", type = "unknown")
            custom_filter.required_tables = self.filter_settings.extra_required_tables
            custom_filter.required_wheres = self.filter_settings.extra_required_where
            filters.append(custom_filter)

        energy_filters = defaultdict()
        energy_filters["dG"] = dGCutoffFilter(0)
        energy_filters["dSASA"] = dSASACutoffFilter(0)
        energy_filters["total"] = TotalScoreCutoffFilter(0)

        for energy_type in self.filter_settings.energy_types:
            if not self.filter_settings.get_energy_enabled(energy_type):continue

            filter = energy_filters[energy_type]
            filter.set_value(self.filter_settings.get_energy_cutoff(energy_type))
            filters.append(filter)

        return filters

    def _setup_cdr_types(self, data_type, is_camelid, features_type, native_path = None):
        if data_type == "length":
            data_class = CDRLengthData(native_path, is_camelid)

        elif data_type == "cluster":
            data_class = CDRClusterData(native_path, is_camelid)

        elif data_type == "sequence":
            data_class = CDRSequenceData(native_path, is_camelid)

        elif data_type == "aligned_sequence":
            data_class = CDRAlignedSequenceData(self._setup_outdir_individual(['clustal']), self._setup_outdir_combined(['clustal']), native_path, is_camelid)

        else:
            sys.exit("data_type not understood!")


        for strategy in self.strategies:
            db_path = self.get_db_path(strategy, features_type)
            if not os.path.exists(db_path):
                sys.exit("DB path does not exist.  Please Run Features reporter with StructureScores and ScoreTypes for this strategy\n"+db_path)
            print db_path
            con = sqlite3.connect(db_path)
            data_class.add_data(strategy, con)

        return data_class

    def set_strategies_from_main_dir_top_dir(self):
        dirs = glob.glob(self.main_dir.get()+"/TOP_*")
        non_redun = defaultdict()
        for d in sorted(dirs):
            #print d
            non_redun ["_".join( d.split( "_" )[ 3: ] ) ] = None


        for x in sorted(non_redun.keys()):
            #print x
            self.strategies.append(x)

    def set_strategies_from_databases(self):
        """
        Set the strategies from the main_dir/databases directory
        :return:
        """
        if not os.path.exists(self.main_dir.get()+"/databases"):
            print "Could not find any databases to use for analysis.  Please make sure databases are in "+self.main_dir.get()+"/databases"
            print "Databases should end with .db or .db3 extension.  Naming format of databases: strategy.features_type.scorefunction.db"
            return

        dbs = glob.glob(self.main_dir.get()+"/databases/*.db*")
        if len(dbs) == 0:
            print "Could not find any databases to use for analysis.  Please make sure databases are in "+self.main_dir.get()+"/databases"
            print "Databases should end with .db or .db3 extension.  Naming format of databases: strategy.features_type.scorefunction.db"
            return

        nr_dbs = defaultdict()
        for db in sorted(dbs):
            strategy = os.path.basename(".".join(db.split('.')[:-3]))
            if not nr_dbs.has_key(strategy):

                self.strategies.append(strategy)
                self.strategy_scorefxns[strategy] = db.split('.')[-2]
                nr_dbs[strategy] = " "

    def set_strategies(self, strategies):
        self.strategies = strategies

    def get_strategies(self):
        return self.strategies

    def get_db_path(self, strategy, features_type = 'antibody'):
        db_dir = self.main_dir.get()+"/databases"
        db_path = db_dir+"/"+strategy+"."+features_type+"_features."+self.strategy_scorefxns[strategy]+".db3"; #This may need to change later
        return db_path

    def get_full_features_type(self, type):
        if type == "cluster_features":
            return type
        else:
            return type+self.features_hbond_sets[ self.features_hbond_set.get() ]

    ################ Main Functions #######################
    def run_features(self, type, plot_name = ""):

        if not plot_name:
            plot_name = self.out_dir_name.get()

        if not plot_name:
            print "No root name set!"
            return

        outdir = self._setup_outdir(["features_plots", plot_name], False)
        db_dir = self.main_dir.get()+"/databases"

        if len(self.strategies) == 0:
            return

        if not os.path.exists(db_dir): sys.exit("Please run Features reporter and copy databases to main Strategy Analysis DIR /databases.")

        fulltype = self.get_full_features_type(type)
        creator = json_creator.JsonCreator(outdir, fulltype)

        if self.reference_db.get() and os.path.exists(self.reference_db.get()):
            creator.add_sample_source_info(self.reference_db.get(), "ref", True)

        for strategy in self.strategies:
            id = strategy
            id = id.replace("talaris2013", "talaris")

            #Fix up the name so it is not too long.  Will add options to do this manually later:
            """
            id = id.replace("cluster", "clus")
            id = id.replace("exclude", "excl")
            id = id.replace("include", "incl")

            """

            db_path = self.get_db_path(strategy, type)
            if not os.path.exists(db_path):
                sys.exit(db_path +" does not exist!")

            creator.add_sample_source_info(db_path, id)

        creator.save_json(outdir+"/"+type+"_"+plot_name+".json")
        creator.run_json()

        #pwd = os.getcwd()
        #os.chdir("build")
        #os.system("cp -r *../"+outdir)
        #os.system("mv build build_old")
        #os.chdir(pwd)

        #This is because outdir still doesn't work in Features.  I swear I have fixed it like 3 times already.
        dirs = glob.glob("build/*")
        for d in dirs:
            os.system('cp -r '+d+' '+outdir)
        shutil.rmtree("build")
        print "Plots ignore any set filters.  To plot with filters, create new databases through query..."
        print "Complete..."

    def output_stats(self):

        if len(self.strategies) == 0:
            print "No strategies set..."
            return

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        outdir_stats = self._setup_outdir_combined(["stats"])
        outdir_lists = self._setup_outdir_combined(["ordered_pdb_lists"])

        scores = self._setup_scores("antibody", True)
        top_n = self.top_n.get()
        top_n_combined = self.top_n_combined.get()

        #Output Strategy data:
        for score_class in scores:
            if isinstance(score_class, DecoyData): pass

            if score_class.name == "dSASA": reverse = True
            else: reverse = False

            for strategy in self.strategies:
                OUTFILE = open(outdir_lists+"/"+strategy+"_ORDERED_"+score_class.get_outname()+".txt", 'w')
                data = score_class.get_strategy_data(strategy, True)

                for tup in sorted(data.keys(), reverse = reverse):
                    triple = data[tup]
                    OUTFILE.write(os.path.basename(triple.decoy)+"\t"+get_str(triple.score)+"\t"+triple.strategy+"\n")
                OUTFILE.close()


        #This is broken for some unkown reason. - Ends up skipping pretty much everything.

        #Output Combined Data:
        COMBINED = open(outdir_stats+"/combined_selection_data.txt", 'w')
        COMBINED.write("#strategy decoy "+" ".join([score_class.get_outname() for score_class in scores])+"\n")
        total_scores = scores[0]
        all_data = total_scores.get_concatonated_map()
        for decoy in all_data:
            strategy = all_data[decoy].strategy
            line = strategy+" "+os.path.basename(decoy)
            #print strategy
            for score_class in scores:
                if score_class.name == "combined_str_score" or score_class.name == "dG_top_Ptotal":continue

                if not score_class.all_data[strategy].has_key(decoy):
                    print "Skipping "+strategy+" "+decoy
                    continue

                line = line+" "+get_str(score_class.all_data[strategy][decoy].score)
            COMBINED.write(line+"\n")
        COMBINED.close()


        #Output Stats:
        STATS = open(outdir_stats+"/combined_selection_data_stats.txt", 'w')
        STATS.write("#strategy "+" ".join([score_class.get_outname()+"_avg"+" "+score_class.get_outname()+"_sd " for score_class in scores])+"\n")

        for strategy in self.strategies:
            line = strategy
            for score_class in scores:
                if not score_class.has_real_values() or score_class.name == "combined_str_score":continue
                raw_scores_tuple = score_class.get_strategy_data(strategy, True)
                raw_scores = [s[0] for s in raw_scores_tuple]
                m = numpy.mean(raw_scores)
                sd = numpy.std(raw_scores)
                line = line +" %.3f"%m+" "+"%.3f"%sd
            STATS.write(line+"\n")
        STATS.close()

        #Output Top Stats:
        STATS = open(outdir_stats+"/combined_selection_data_top_"+repr(top_n)+"_stats.txt", 'w')
        STATS.write("#strategy "+" ".join([score_class.get_outname()+"_avg"+" "+score_class.get_outname()+"_sd " for score_class in scores])+"\n")
        for strategy in self.strategies:
            line = strategy
            for score_class in scores:
                if not score_class.has_real_values() or score_class.name == "combined_str_score":continue
                raw_scores_tuple = score_class.get_top_strategy_data(strategy, top_n, True)
                raw_scores = [s[0] for s in raw_scores_tuple]
                m = numpy.mean(raw_scores)
                sd = numpy.std(raw_scores)
                line = line +" %.3f"%m+" "+"%.3f"%sd
            STATS.write(line+"\n")
        STATS.close()
        print "Complete"

    def output_score_extra_stats(self):
        top_n = self.top_n.get()
        main_scores = self._setup_scores()
        all_scores = self._setup_scores("antibody", True)

        out_dir = self._setup_outdir_individual(["raw_score_data"])

        #Top, then all
        for strategy in self.strategies:
            for score in main_scores:
                if score.name == "combined_str_score":
                    continue

                if isinstance(score, DecoyData): pass

                #Get TopN of that particular score
                decoy_list = score.get_ordered_decoy_list(strategy, top_n)

                outfile = out_dir+"/scores_of_top_"+score.name+"_"+strategy+".txt"
                print "writing "+outfile
                OUT = open(outfile, 'w')

                if not score.name == "dG_top_Ptotal":
                    header = "#decoy\t"+score.name
                else:
                    header = "#decoy"

                for a_score in all_scores:
                    if a_score.name == "combined_str_score" or a_score.name == "dG_top_Ptotal": continue
                    if a_score.name == score.name: continue

                    header = header+"\t"+a_score.name

                OUT.write(header+"\n")

                for decoy in decoy_list:

                    if not score.name == "dG_top_Ptotal":
                        line = os.path.basename(decoy)+"\t"+get_str(score.get_score_for_decoy(strategy, decoy))
                    else:
                        line = os.path.basename(decoy)

                    for a_score in all_scores:
                        if a_score.name == "combined_str_score": continue
                        if a_score.name == score.name: continue

                        line = line+ "\t"+get_str(a_score.get_score_for_decoy(strategy, decoy))
                    OUT.write(line+"\n")
                OUT.close()

        #All decoys

        score = main_scores[0]

        for strategy in self.strategies:
            all_decoys = score.get_ordered_decoy_list(strategy)

            out_file = out_dir+"/all_scores_"+strategy+".txt"
            print "writing "+out_file

            OUT = open(out_file,'w')

            if not score.name == "dG_top_Ptotal":
                header = "#decoy\t"+score.name
            else:
                header = "#decoy"

            for a_score in all_scores:
                if a_score.name == "combined_str_score" or a_score.name == "dG_top_Ptotal": continue
                if a_score.name == score.name: continue

                header = header+"\t"+a_score.name
            OUT.write(header+"\n")

            for decoy in all_decoys:
                if not score.name == "dG_top_Ptotal":

                    line = os.path.basename(decoy)+"\t"+get_str(score.get_score_for_decoy(strategy, decoy))
                else:
                    line = os.path.basename(decoy)
                for a_score in all_scores:
                    if a_score.name == "combined_str_score": continue
                    if a_score.name == score.name: continue

                    line = line+ "\t"+get_str(a_score.get_score_for_decoy(strategy, decoy))
                OUT.write(line+"\n")
            OUT.close()

        print "Complete"

    def copy_top_strategy(self, native_path = None):

        top_n = self.top_n.get()
        scores = self._setup_scores()

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        print "Copying Top Models.."
        #Each Strategy Top Scoring (Skip total score here for now):
        for strategy in self.strategies:
            for score in scores:
                if score.name == "combined_str_score":
                    continue

                out_dir = self._setup_outdir_individual(["pdbs_sessions", "top_"+repr(top_n)+"_"+score.get_outname()+"_"+strategy])
                SCORELIST = open(out_dir+"/MODELS.txt", 'w')
                print "Copying "+strategy+" "+score.get_outname()+" to: "+out_dir
                if isinstance(score, DecoyData): pass

                decoys = score.get_top_strategy_data(strategy, top_n)
                decoy_list = score.get_ordered_decoy_list(strategy, top_n)
                load_as = []
                i = 1
                for decoy in decoy_list:
                    load_as.append("model_"+repr(i)+"_"+score.get_outname()+"_"+get_str(decoys[decoy].score))
                    os.system('cp '+decoy+" "+out_dir+"/"+"top_"+repr(i)+"_"+os.path.basename(decoy))
                    SCORELIST.write(repr(i)+"\t"+get_str(decoys[decoy].score)+"\t"+os.path.basename(decoy)+"\n")
                    i+=1
                analyze_strat.make_pymol_session_on_top(out_dir, decoy_list, load_as, out_dir, score.get_outname(), top_n, native_path)
                SCORELIST.close()

        print "Complete"

    def copy_top_combined(self, native_path = None):
        """
        Outputs total_score,
        """
        if not self.out_dir_name.get():
            print "No root name set!"
            return

        top_n = self.top_n_combined.get()
        scores = self._setup_scores()

        #Overall Strategy:
        for score in scores:
            if score.name == "combined_str_score":
                continue

            if isinstance(score, DecoyData): pass
            outdir_top_pdbs = self._setup_outdir_combined(["top_structures", score.get_outname()])
            outdir_top_sessions = self._setup_outdir_combined(["top_sessions"])

            SCORELIST = open(outdir_top_pdbs+"/MODELS.txt", 'w')
            print "Copying "+score.get_outname()+" to: "+outdir_top_pdbs
            decoys = score.get_top_all_data(top_n)
            decoy_list = score.get_ordered_decoy_list_all(top_n)
            load_as = []
            i = 1
            for decoy in decoy_list:
                load_as.append("model_"+repr(i)+"_"+score.get_outname()+"_"+get_str(decoys[decoy].score))
                os.system('cp '+decoy+" "+outdir_top_pdbs+"/top_"+repr(i)+"_"+os.path.basename(decoy))
                SCORELIST.write(repr(i)+"\t"+get_str(decoys[decoy].score)+"\t"+os.path.basename(decoy)+"\n")
                i+=1

            analyze_strat.make_pymol_session_on_top(outdir_top_pdbs, decoy_list, load_as, outdir_top_sessions, score.get_outname(), top_n, native_path)
            SCORELIST.close()

        print "Complete"

    def copy_all_models(self, native_path = None):

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        if len(self.strategies) == 0:
            print "No strategies set..."
            return


        #outdir = self._setup_outdir()
        scores = self._setup_scores()

        #Each Strategy Top Scoring (Skip total score here for now):
        for strategy in self.strategies:
            for score in scores[1:]:
                out_dir = self._setup_outdir_individual(["all"+"_"+score.get_outname()+"_"+strategy], False)
                print "Copying "+strategy+" "+score.get_outname()+" to: "+out_dir
                if isinstance(score, DecoyData): pass

                decoys = score.get_strategy_data(strategy)
                decoy_list = score.get_ordered_decoy_list(strategy)
                load_as = []
                i = 1
                for decoy in decoy_list:
                    load_as.append("model_"+repr(i)+"_"+score.get_outname()+"_"+get_str(decoys[decoy].score))
                    os.system('cp '+decoy+" "+out_dir+"/top_"+repr(i)+"_"+os.path.basename(decoy))
                    i+=1
                analyze_strat.make_pymol_session_on_top(out_dir, decoy_list, load_as, out_dir, score.get_outname(), None, native_path)

        #Overall Strategy:
        for score in scores:
            if isinstance(score, DecoyData): pass
            out_dir = self._setup_outdir_combined(["all_structures", score.get_outname()])
            print "Copying "+score.get_outname()+" to: "+out_dir
            decoys = score.get_concatonated_map()
            decoy_list = score.get_ordered_decoy_list_all()
            load_as = []
            i = 1
            for decoy in decoy_list:
                load_as.append("model_"+repr(i)+"_"+score.get_outname()+"_"+get_str(decoys[decoy].score))
                os.system('cp '+decoy+" "+out_dir+"/top_"+repr(i)+"_"+os.path.basename(decoy))
                i+=1
            session_dir = out_dir = self._setup_outdir_combined(["all_sessions"])
            analyze_strat.make_pymol_session_on_top(out_dir, decoy_list, load_as, session_dir, score.get_outname(), None, native_path)

        print "Complete"

    def run_clustal_omega_on_strategies(self, processors, output_format = "fasta", extra_options = "", native_path = None):

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        scores = self._setup_scores()
        top_n = self.top_n.get()

        for strategy in self.strategies:
            print "\nRunning Clustal on: "+strategy
            score_zero = scores[0]
            if isinstance(score_zero, DecoyData): pass

            #Output per strategy ALL
            decoys = score_zero.get_strategy_data(strategy)
            decoy_header_dict = defaultdict()

            for decoy in decoys:
                decoy_header_dict[decoy]=os.path.basename(decoy)

            fasta_dir = self._setup_outdir_individual(["sequences"])
            fasta_path = fasta_dir+"/"+strategy+"_all.fasta"
            fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native", self.is_camelid.get())

            aln_dir = self._setup_outdir_individual(["clustal"])
            aln_name = strategy+"_all.clus"
            clustal_runner = ClustalRunner(fasta_path)
            clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
            clustal_runner.set_threads(processors)
            clustal_runner.set_extra_options(extra_options)
            clustal_runner.set_output_format(output_format)
            clustal_runner.output_alignment(aln_dir, aln_name)

            #Output on Top Scoring:

            for score in scores:
                print score
                if isinstance(score, DecoyData):pass

                basename = strategy+"_"+score.get_outname()+"_top_"+str(top_n)
                fasta_path = fasta_dir+"/"+basename+".fasta"
                aln_name = basename+".clus"

                decoy_header_dict = defaultdict()
                decoys = score.get_top_strategy_data(strategy, top_n)
                for decoy in decoys:
                    header = "v"+get_str(decoys[decoy].score)+"::"+os.path.basename(decoy)
                    decoy_header_dict[decoy]= header
                fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native", self.is_camelid.get())
                clustal_runner.set_fasta_path(fasta_path)
                clustal_runner.set_threads(processors)
                clustal_runner.set_extra_options(extra_options)
                clustal_runner.set_output_format(output_format)
                clustal_runner.output_alignment(aln_dir, aln_name)

        print "Complete"

    def run_clustal_omega_on_top_combined(self, processors, output_format, extra_options = "", native_path = None):

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        scores = self._setup_scores()
        top_n = self.top_n.get()

        for score in scores:
            print "Running Clustal Omega on "+score.get_outname()
            if isinstance(score, DecoyData): pass
            decoys = score.get_top_all_data(top_n, False)

            decoy_header_dict = defaultdict()
            for decoy in decoys:
                header = "v"+get_str(decoys[decoy].score)+"::"+os.path.basename(decoy)
                decoy_header_dict[decoy] = header

            root_name = self.out_dir_name.get()
            basename = root_name+"_"+score.get_outname()+"_top_"+repr(top_n)
            fasta_dir = self._setup_outdir_combined(["sequences"])
            fasta_path = fasta_dir+"/"+basename+".fasta"

            clustal_dir = self._setup_outdir_combined(["clustal"])
            clustal_name = basename+".aln"

            fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native", self.is_camelid.get())
            clustal_runner = ClustalRunner(fasta_path)
            clustal_runner.set_threads(processors)
            clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
            clustal_runner.set_extra_options(extra_options)
            clustal_runner.set_output_format(output_format)
            clustal_runner.output_alignment(clustal_dir, clustal_name)

        print "Complete"

    def run_clustal_omega_on_all_combined(self, processors, output_format, extra_options = "", native_path = None):

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        print "Running Clustal on All Combined"
        scores = self._setup_scores()
        score = scores[0]
        if isinstance(score, DecoyData):pass

        root_name = self.out_dir_name.get()
        fasta_dir = self._setup_outdir_combined(["sequences"])
        fasta_path = fasta_dir+"/"+root_name+"_all.fasta"

        clustal_dir = self._setup_outdir_combined(["clustal"])
        clustal_name = root_name+"_all.aln"

        all_data = score.get_concatonated_map(False)

        all_data_array = []
        for s in scores:
            if s.name == "combined_str_score":
                continue
            all_data_array.append(s.get_concatonated_map(False))

        decoy_header_dict = defaultdict()
        for decoy in all_data:
            a = all_data_array[0]
            header = "v"+get_str(a[decoy].score)

            for a in all_data_array[1:]:
                header = header+":"+get_str(a[decoy].score)
            header = header+":"+os.path.basename(decoy)
            decoy_header_dict[decoy] = header

        fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native", self.is_camelid.get())
        clustal_runner = ClustalRunner(fasta_path)
        clustal_runner.set_threads(processors)
        clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
        clustal_runner.set_extra_options(extra_options)
        clustal_runner.set_output_format(output_format)
        clustal_runner.output_alignment(clustal_dir, clustal_name)

        print "Complete"

    def output_len_or_clus_alignment(self, alignment_type, features_type = 'antibody', is_camelid = False, native_path = None):

        def _output_alignment(self, outdir, top_decoys, score, type_data, strategy = None, extra_name = "top"):

            if isinstance(type_data, CDRData):pass
            if isinstance(score, DecoyData):pass

            if not strategy:
                outname = "cdr_type_alignments_"+alignment_type+"_"+score.get_outname()+"_"+extra_name+".txt"
            else:
                outname = "cdr_type_alignments_"+alignment_type+"_"+score.get_outname()+"_"+strategy+"_"+extra_name+"_.txt"

            print "Outputting cdr type alignment: "+outdir+"/"+outname

            OUTFILE = open(outdir+"/"+outname, 'w')
            if native_path:
                header = "#decoy\tscore\tnative_matches"
            else:
                header = "#decoy\tscore"
            for cdr in type_data.cdrs:
                header = header+"\t"+cdr
            OUTFILE.write(header+"\n")

            all_data = score.get_concatonated_map()
            all_type_data = type_data.get_concatonated_map()

            if native_path:
                line = "Native::"+os.path.basename(native_path)+"\t...\tNA"
                if all_type_data.has_key((cdr, "native")):
                    native_info = type_data[(cdr, "native")]
                else:
                    native_info = type_data.get_native_data()

                if isinstance(native_info, CDRDataInfo): pass

                for cdr in type_data.cdrs:
                    line = line+"\t"+str(native_info.get_value_for_cdr(cdr)).rjust(10)
                OUTFILE.write(line+"\n")

                for decoy in top_decoys:
                    decoy_info = all_type_data[decoy]
                    score_info = all_data[decoy]
                    counts = count_native_matches(decoy_info, native_info, type_data.cdrs)
                    line = os.path.basename(decoy)+"\t"+get_str(score_info.score)+"\t"+repr(counts)
                    for cdr in type_data.cdrs:
                        line = line+"\t"+str(get_star_if_native(decoy_info, native_info, cdr))

                    OUTFILE.write(line+"\n")
            else:
                for decoy in top_decoys:
                    print decoy
                    decoy_info = all_type_data[decoy]
                    score_info = all_data[decoy]
                    print repr(score_info)
                    print repr(score_info.score)

                    line = os.path.basename(decoy)+"\t"+get_str(score_info.score)
                    for cdr in type_data.cdrs:
                        line = line+"\t"+str(decoy_info.get_value_for_cdr(cdr))
                    OUTFILE.write(line+"\n")

            OUTFILE.close()

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        top_n = self.top_n.get()
        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type, native_path)
        scores = self._setup_scores(features_type)
        if isinstance(data_class, CDRData):pass

        ### Top/All  each Strategy
        outdir = self._setup_outdir_individual(['cdr_alignments'])
        for strategy in self.strategies:

            for score in scores:
                if isinstance(score, DecoyData): pass
                top_decoys = score.get_ordered_decoy_list(strategy, top_n)
                all_decoys = score.get_ordered_decoy_list(strategy)
                #print "Top: "+repr(top_decoys)

                print "Working on : "+strategy+" "+score.name
                _output_alignment(self, outdir, top_decoys, score, data_class, strategy, "top")
                _output_alignment(self, outdir, all_decoys, score, data_class, strategy, "all")


        ### Top/All each combined
        outdir = self._setup_outdir_combined(['cdr_alignments'])
        for score in scores:
            top_decoys = score.get_ordered_decoy_list_all(top_n)
            all_decoys = score.get_ordered_decoy_list_all()

            print "Top: "+repr(len(top_decoys))+score.name
            _output_alignment(self, outdir, top_decoys, score, data_class, None, "top")
            _output_alignment(self, outdir, all_decoys, score, data_class, None, "all")

        print "Complete"

    def output_len_or_clus_recovery(self, alignment_type, features_type = 'antibody', is_camelid = False, native_path = None):
        if not native_path:
            print "Must pass select native path to calculate recoveries"
            return

        def _get_header(self, type_data):
            header = "#name\tavg"
            for cdr in type_data.cdrs:
                header = header+"\t"+cdr
            return header

        def _add_native_line(self, native_info, type_data, OUTFILE):
            line = "Native::"+os.path.basename(native_path)+"\tNA"

            if isinstance(native_info, CDRDataInfo): pass

            for cdr in type_data.cdrs:
                line = line+"\t"+str(native_info.get_value_for_cdr(cdr))
            OUTFILE.write(line+"\n")

        def _add_recovery_line(self, label, decoys, type_data, OUTFILE, strategy = None):

            if isinstance(type_data, CDRData):pass

            total = 0
            for cdr in type_data.cdrs:
                enrichment_data = calculate_recovery(type_data.get_native_data(), type_data, cdr, decoys)
                total = total + enrichment_data.get_perc_decimal()
            avg = total/6.0
            line = label+"\t%.3f"%avg
            for cdr in type_data.cdrs:
                enrichment_data = calculate_recovery(type_data.get_native_data(), type_data, cdr, decoys)

                line = line+"\t%.3f"%enrichment_data.get_perc_decimal()

            OUTFILE.write(line+"\n")


        if not self.out_dir_name.get():
            print "No root name set!"
            return

        top_n = self.top_n.get()

        data_class = self._setup_cdr_types(alignment_type, is_camelid, features_type, native_path)
        scores = self._setup_scores(features_type)
        if isinstance(data_class, CDRData):pass

        #Each strategy + Top Values
        outdir = self._setup_outdir_individual(['enrichment'])
        OUTFILE = open(outdir+"/"+"cdr_type_recoveries_"+alignment_type+"_.txt", 'w')
        OUTFILE.write(_get_header(self,data_class)+"\n")
        _add_native_line(self, data_class.get_native_data(), data_class, OUTFILE)
        for strategy in self.strategies:
            outname = strategy
            decoys = scores[0].get_strategy_data(strategy).keys()
            _add_recovery_line(self, outname, decoys, data_class, OUTFILE, strategy)
            for score in scores:
                if isinstance(score, DecoyData): pass
                if score.name == "combined_str_score":continue

                top_decoys = score.get_ordered_decoy_list(strategy, top_n)
                outname = strategy+"_"+score.get_outname()+"_top_"+repr(top_n)
                _add_recovery_line(self, outname, top_decoys, data_class, OUTFILE, strategy)
        OUTFILE.close()


        #Combined + Top Values
        outdir = self._setup_outdir_combined(['enrichment'])
        OUTFILE = open(outdir+"/"+"cdr_type_recoveries_all_"+alignment_type+"_.txt", 'w')
        OUTFILE.write(_get_header(self,data_class)+"\n")
        _add_native_line(self, data_class.get_native_data(), data_class, OUTFILE)
        outname = "combined_all"
        _add_recovery_line(self, outname, score.get_concatonated_map().keys(), data_class, OUTFILE)
        for score in scores:
            if score.name == "combined_str_score":continue
            top_decoys = score.get_ordered_decoy_list_all(top_n)
            outname = "combined_"+score.get_outname()+"_top_"+repr(top_n)
            #print "Top: "+repr(top_decoys)
            _add_recovery_line(self, outname, top_decoys, data_class, OUTFILE)
        OUTFILE.close()
        print "Complete"

    def output_len_or_clus_enrichment(self, alignment_type, top_n, features_type = 'antibody', is_camelid = False):
        if not self.out_dir_name.get():
            print "No root name set!"
            return

    def create_score_subset_database(self, score_name, prefix, features_type = 'antibody'):
        self._setup_scores()
        score = self._get_score(score_name)
        if not isinstance(score, DecoyData):
            print "Score type not found!"
            return

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        top_n = self.top_n.get()
        fdir= os.path.split(os.path.abspath(__file__))[0]+"/features_inputs"

        for strategy in self.strategies:
            temp_name = "temp_PDBLIST.txt"
            OUTFILE = open(temp_name, 'w')
            decoys = score.get_ordered_decoy_list(strategy, top_n)
            for decoy in decoys:
                OUTFILE.write(decoy+"\n")
            OUTFILE.close()

            out_db_name = prefix+"_"+strategy
            out_db_batch = "Subset"

            #This should be done manually by getting struct_id and copying all the data in a new database via python sqlite3
            #I do not have time to figure that out and get it working right now, so this will have to do.

            analyze_strat.create_features_db(temp_name, fdir, features_type+"_features", self.rosetta_extension.get(), self.strategy_scorefxns[ strategy ], out_db_name, out_db_batch, self.main_dir.get(), False)


            os.remove(temp_name)

########################################################################################################################
###   Window Modules
########################################################################################################################


########################################################################################################################
###   Filters
########################################################################################################################



########################################################################################################################
###   CDRData
########################################################################################################################



########################################################################################################################
### Helper Classes
########################################################################################################################


class EnrichmentInfo:
    """
    Simple class for holding enrichment/recovery information
    """
    def __init__(self, count, total):
        self.count = count
        self.total = total
        self._calc_per()

    def _calc_per(self):
        self.perc = self.count/float(self.total)

    def get_count(self):
        return self.count

    def get_total(self):
        return self.total

    def get_perc_decimal(self):
        return self.perc

    def get_perc_whole(self):
        return self.perc*100

    def get_formated_perc(self, perc):
        return "%.3f"%perc





########################################################################################################################
### Helper Functions
########################################################################################################################
def get_str(value):
    if type(value) == str:
        return value
    else:
        return "%.3f"%value

def calculate_recovery(native_data, all_decoy_data, cdr, decoy_list = None):
    """
    Calculate the recovery of some value to native
    Returns
    """
    if isinstance(native_data, CDRDataInfo): pass
    if isinstance(all_decoy_data, CDRData): pass

    raw_decoy_map = all_decoy_data.get_concatonated_map()
    if not decoy_list:

        decoy_list = raw_decoy_map.keys()

    count = 0
    for decoy in decoy_list:
        cdr_info = raw_decoy_map[decoy]
        if isinstance(cdr_info, CDRDataInfo): pass
        if cdr_info.get_value_for_cdr(cdr) == native_data.get_value_for_cdr(cdr):
            count+=1

    enrich_info = EnrichmentInfo(count, len(decoy_list))
    return enrich_info

def calculate_observed_value(value, all_decoy_data, cdr, decoy_list = None):
    """
    Calculate the enrichment of some value
    """
    if isinstance(all_decoy_data, CDRData): pass

    raw_decoy_map = all_decoy_data.get_concatonated_map()
    if not decoy_list:

        decoy_list = raw_decoy_map.keys()

    count = 0
    for decoy in decoy_list:
        cdr_info = raw_decoy_map[decoy]
        if isinstance(cdr_info, CDRDataInfo): pass
        if cdr_info.get_value_for_cdr(cdr) == value:
            count+=1

    enrich_info = EnrichmentInfo(count, len(decoy_list))
    return enrich_info

def count_native_matches(decoy_data, native_data, cdrs):
    if isinstance(decoy_data, CDRDataInfo): pass
    if isinstance(native_data, CDRDataInfo): pass

    count = 0
    for cdr in cdrs:
        if native_data.get_value_for_cdr(cdr) == decoy_data.get_value_for_cdr(cdr):
            count+=1
    return count

def get_star_if_native(decoy_data, native_data, cdr):
    if native_data.get_value_for_cdr(cdr) == decoy_data.get_value_for_cdr(cdr):
        return "*"
    else:
        return decoy_data.get_value_for_cdr(cdr)



if __name__ == "__main__":

    # 1) Main StrategyAnalysis DIR
    # 2) Output DIR
    # 3) Any Names of strategies to run.  Else will attempt load from Main StrategyAnalysisDIR

    main()