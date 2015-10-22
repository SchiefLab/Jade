#!/usr/bin/env python

# Yes, the imports are all over.  Basically everything I've coded over the past few years.
import sys
import os
import re
import glob
import sqlite3
import atexit
import shutil

from argparse import ArgumentParser
from collections import defaultdict

# TkInter
from Tkinter import *
from tkFont import *
import tkFileDialog
import tkSimpleDialog
import tkMessageBox

p = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(p);  # Allows all modules to use all other modules, without needing to update pythonpath

# PyIgD
import bin.RunRabD_Feat as analyze_strat
from pymol.PyMolScriptWriter import *
from antibody.cdr_data.CDRDataTypes import *
from antibody.decoy_data.DecoyDataTypes import *
from tools.filters.DataFilters import *
from tools.filters.FilterSettings import *
from tools.Threader import *
from RAbD.window_modules.FilterSettingsWindow import *
from sequence import fasta

# Rosetta Tools
import create_features_json as json_creator

p = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(p);  # Allows all modules to use all other modules, without needing to update pythonpath

from rosetta import *

# Rosetta is only used for set_native_data_from_rosetta function for clusters. (CDRClusterer - which uses a pose to get dihedrals.  needs refactoring.)
rosetta.init(" -ignore_unrecognized_res -ignore_zero_occupancy false -ex1 -ex2 -use_input_sc"
             " -antibody:numbering_scheme AHO_Scheme "
             " -antibody:cdr_definition North")


# We assume that the names of each decoy is different.  How can we remove this dependency?
#  This may be a problem if people do not use out:prefix or out:suffix to name their decoys...

### Setup Enums ###
# length = 1; cluster = 2; sequence = 3;



def main():
    parser = ArgumentParser()
    parser.add_argument("--db_dir",
                        help="Directory with databases to compare",
                        default="databases")

    parser.add_argument("--analysis_name",
                        help="Main directory to complete analysis",
                        default="prelim_analysis")

    parser.add_argument("--native",
                        help="Any native structure to compare to")

    parser.add_argument("--root_dir",
                        help="Root directory to run analysis from",
                        default=os.getcwd())

    parser.add_argument("--cdrs",
                        help="A list of CDRs for the analysis (Not used for Features Reporters)",
                        choices=["L1", "H1", "L1", "H2", "L3", "H3"],
                        nargs='*',
                        default=["L1", "L2", "L3", "H1", "H2", "H3"])

    options = parser.parse_args()

    if options.root_dir != os.getcwd():
        print "Changing to root."
        os.chdir(options.root_dir)

    GUI = CompareAntibodyDesignStrategies_GUI(Tk(), options.db_dir, options.analysis_name)

    # Set any values
    if options.native:
        GUI.native_path = options.native

    GUI.compare_designs.set_cdrs_from_list(options.cdrs)

    GUI.run()


class Listbox(Listbox):
    def autowidth(self, maxwidth, list=None):
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
        for w in range(0, maxwidth + 1, 5):
            if self.winfo_reqwidth() >= pixels:
                break
            self.config(width=width + w)


class CompareAntibodyDesignStrategies_GUI:
    def __init__(self, main, db_dir="", analysis_dir="", strategies=[]):

        self._tk_ = main
        self._tk_.title("PyIgDesign Compare")
        self.compare_designs = CompareAntibodyDesignStrategies(db_dir, analysis_dir, strategies)
        if self.compare_designs.db_dir and not self.compare_designs.strategies:
            self.compare_designs.set_strategies_from_databases()

        self.current_dir = os.getcwd()

        self.clustal_procs = IntVar(value=multiprocessing.cpu_count())
        self.clustal_output_format = StringVar(value="clu")
        self.native_path = None
        self.clustal_output_formats = ['fasta', 'clustal', 'msf', 'phylip', 'selex', 'stockholm', 'vienna',
                                       'a2m', 'fa', 'clu', 'phy', 'st', 'vie']
        self.base_clustal_options = " -v --auto --force"

        ###Sub Windows ####
        self.filter_settings_window = FilterSettingsWindow(self.compare_designs.filter_settings)

        ###Tracers###
        self.compare_designs.is_camelid.trace_variable('w', self.camelid_tracer)
        #atexit.register(self.exit)
        self._tk_.protocol("WM_DELETE_WINDOW", self.exit)

    def exit(self):
        if threads.n_alive() > 1:
            if tkMessageBox.askyesno(title = "Exit?", message = repr(threads.n_alive())+" Processes still running.  Exit?"):
                threads.kill_all()
                self._tk_.destroy()
            else:
                self._tk_.mainloop()
        else:
            self._tk_.destroy()


    def run(self):
        self.set_tk()
        self.set_menu()
        self.sho_tk()
        self._tk_.mainloop()

    def set_tk(self):
        # self.db_dir_entry = Entry(self._tk_, textvariable = self.compare_designs.out_dir_name, justify = CENTER)
        self.out_dir_entry = Entry(self._tk_, textvariable=self.compare_designs.out_dir_name, justify=CENTER)

        # self.root_dir_label = Label(self._tk_, text = "Root Directory", justify = CENTER)
        self.out_dir_label = Label(self._tk_, text="Analysis Name", justify=CENTER)

        self.all_strategies_listbox = Listbox(self._tk_)
        self.current_strategies_listbox = Listbox(self._tk_)

        self.ab_features_button = Button(self._tk_, text="Run Antibody Features",
                                         command=lambda: self.run_features_reporter("antibody"), justify=CENTER)
        self.clus_features_button = Button(self._tk_, text="Run Cluster Features",
                                           command=lambda: self.run_features_reporter("cluster"), justify=CENTER)

        self.ab_features_options_label = Label(self._tk_, text="Antibody Features Options", justify=CENTER)

        self.normal_hbond_radio = Radiobutton(self._tk_, text="All Hbond R Scripts",
                                              variable=self.compare_designs.features_hbond_set, value=0)
        self.min_hbond_radio = Radiobutton(self._tk_, text="Minimal Hbond R Scripts",
                                           variable=self.compare_designs.features_hbond_set, value=1)
        self.no_hbond_radio = Radiobutton(self._tk_, text="No Hbond R Scripts",
                                          variable=self.compare_designs.features_hbond_set, value=2)

        # Setup CDRs
        self.L_chain_buttons = []
        self.H_chain_buttons = []
        for cdr_name in ["L1", "L2", "L3", "H1", "H2", "H3"]:
            button = Checkbutton(self._tk_, text=cdr_name, variable=self.compare_designs.cdrs[cdr_name])
            if re.search("L", cdr_name):
                self.L_chain_buttons.append(button)
            else:
                self.H_chain_buttons.append(button)

                # self.separator = Separator(self._tk_, orient = HORIZONTAL)

        self.individual_analysis = Checkbutton(self._tk_, text = "Individual Analysis", variable=self.compare_designs.individual_analysis)
        self.combined_analysis = Checkbutton(self._tk_, text = "Combined Analysis", variable = self.compare_designs.combined_analysis)

    def sho_tk(self, r=0, c=0):

        # self.root_dir_label.grid(row = r+0, column = c+0, columnspan = 2, sticky = W+E, pady = 7)
        # self.db_dir_entry.grid( row = r+1, column = c+0, columnspan = 2, sticky = W+E, padx = 5)



        self.all_strategies_listbox.grid(row=r + 1, column=c + 0, columnspan=3, padx=6, pady=10)
        self.current_strategies_listbox.grid(row=r + 1, column=c + 3, columnspan=3, padx=6, pady=10)

        position = 0

        for cdr_button in self.L_chain_buttons:
            cdr_button.grid(row=r + 2, column=c + position)
            position += 1

        for cdr_button in self.H_chain_buttons:
            cdr_button.grid(row=r + 2, column=c + position)
            position += 1

        # self.separator.grid(row = r+3, column = c, columnspan = 2, sticky = W+E, pady = 15)

        self.individual_analysis.grid(row = r+3, column = c+3, columnspan = 3, pady=5, sticky=W+E)
        self.combined_analysis.grid(row = r+4, column = c+3, columnspan = 3, pady = 5, sticky= W+E)

        self.out_dir_label.grid(row=r + 5, column=c + 0, columnspan=3, pady=5)
        self.out_dir_entry.grid(row=r + 5, column=c + 3, columnspan=3, padx=5, pady=5)

        self.all_strategies_listbox.bind("<Double-Button-1>",
                                         lambda event: self.add_to_current(self.all_strategies_listbox,
                                                                           self.current_strategies_listbox))
        self.all_strategies_listbox.bind("<Button-2>", lambda event: self.show_strat_items())

        self.ab_features_button.grid(row=r + 6, column=c + 0, columnspan=3, pady=3, sticky=W + E)
        self.clus_features_button.grid(row=r + 6, column=c + 3, columnspan=3, pady=3, sticky=W + E)

        # self.ab_features_options_label.grid(row = r+7, column = c+0, columnspan = 2, pady = 3, sticky = W+E)

        self.normal_hbond_radio.grid(row=r + 8, column=c + 0, columnspan=3, pady=1, sticky=W)
        self.min_hbond_radio.grid(row=r + 9, column=c + 0, columnspan=3, pady=1, sticky=W)
        self.no_hbond_radio.grid(row=r + 10, column=c + 0, columnspan=3, pady=1, sticky=W)

        self.current_strategies_listbox.bind("<Double-Button-1>",
                                             lambda event: self.delete_current(self.current_strategies_listbox))

        self.populate_all_strategies()

    def set_menu(self):

        self.main_menu = Menu(self._tk_)

        ## File Menu ##
        self.file_menu = Menu(self.main_menu, tearoff=0)

        self.file_menu.add_checkbutton(label="Camelid Antibody", variable=self.compare_designs.is_camelid)
        # self.dtypes_menu = Menu(self.main_menu, tearoff = 0)
        # self.dtypes_menu.add_checkbutton(label = "Group by dG", variable = self.compare_designs.group_dG)

        self.file_menu.add_command(label="Filter Models",
                                   command=lambda: self.filter_settings_window.setup_sho_gui(Toplevel(self._tk_)))
        self.file_menu.add_command(label="Read Strategies from DB DIR",
                                   command=lambda: self.read_from_db_dir_set_strategies())
        self.file_menu.add_command(label="Add Strategy", command=lambda: self.add_main_strategy())
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Set Reference Native", command=lambda: self.set_native_path())
        self.file_menu.add_command(label="Set Reference Database", command=lambda: self.set_reference_db())
        # self.file_menu.add_command(label = "Set Scorefunction", command = lambda: self.set_scorefunction())
        self.file_menu.add_command(label="Set top N", command=lambda: self.set_top_n())
        self.file_menu.add_command(label="Set top N For Combined", command=lambda: self.set_top_n_combined())
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Change Root Dir", command=lambda: self.change_root())
        self.file_menu.add_command(label = "Print current threads", command = lambda: self.print_threads())
        self.file_menu.add_separator()
        self.file_menu.add_checkbutton(label="Reload Query Data", variable=self.compare_designs.reload_scores)
        self.file_menu.add_checkbutton(label="Backround Features", variable=self.compare_designs.backround_features)

        self.file_menu.add_separator()


        for name in sorted(self.compare_designs.scores_on.keys()):
            self.file_menu.add_checkbutton(label=name, variable=self.compare_designs.scores_on[name])

        self.main_menu.add_cascade(label="File", menu=self.file_menu)


        ## Score Menu ##
        self.score_menu = Menu(self.main_menu, tearoff=0)
        self.score_menu.add_command(label="Output Score Stats", command=lambda: self.compare_designs.output_stats())
        self.pymol_menu = Menu(self.main_menu, tearoff=0)
        self.pymol_menu.add_command(label="Top Models",
                                    command=lambda: self.compare_designs.copy_top(self.native_path))
        self.pymol_menu.add_command(label="All Models",
                                    command=lambda: self.compare_designs.copy_all_models(self.native_path))

        self.score_menu.add_cascade(label="Create PyMol Sessions", menu=self.pymol_menu)
        self.main_menu.add_cascade(label="PyMol", menu=self.score_menu)

        ## Clustal Menu ##
        self.clustal_menu = Menu(self.main_menu, tearoff=0)
        self.clustal_menu.add_command(label="Set Max Processors", command=lambda: self.set_max_clustal_procs())
        self.clustal_menu.add_command(label="Set Output format", command=lambda: self.set_clustal_output_format())
        self.clustal_menu.add_command(label="Set Soft Wrap", command=lambda: self.set_clustal_soft_wrap())
        self.clustal_menu.add_separator()
        self.clustal_menu.add_command(label="Run Clustal Omega on Top Decoys",
                                      command=lambda: self.run_clustal_omega())
        self.clustal_menu.add_command(label="Run Clustal Omega on ALL Combined Decoys",
                                      command=lambda: self.run_clustal_on_all_combined())
        self.main_menu.add_cascade(label="Clustal", menu=self.clustal_menu)

        ## Alignment ##
        self.alignment_menu = Menu(self.main_menu, tearoff=0)

        self.alignment_menu.add_command(label="Output Length Alignments",
                                        command=lambda: self.compare_designs.output_len_or_clus_alignment('length', 'antibody',native_path = self.native_path))
        self.alignment_menu.add_command(label="Output Cluster Alignments",
                                        command=lambda: self.compare_designs.output_len_or_clus_alignment('cluster', 'antibody',native_path = self.native_path))
        self.alignment_menu.add_command(label="Output CDR Sequence Alignments",
                                        command=lambda: self.compare_designs.output_len_or_clus_alignment('aligned_sequence', 'antibody',native_path = self.native_path))
        self.alignment_menu.add_separator()
        self.main_menu.add_cascade(label="Alignment", menu=self.alignment_menu)


        ## Recovery ##
        self.recovery_menu = Menu(self.main_menu, tearoff=0)
        self.recovery_menu.add_command(label="Output Length Recovery",
                                       command=lambda: self.compare_designs.output_len_or_clus_recovery('length', 'antibody', native_path = self.native_path))
        self.recovery_menu.add_command(label="Output Cluster Recovery",
                                       command=lambda: self.compare_designs.output_len_or_clus_recovery('cluster', 'antibody', native_path = self.native_path))
        self.main_menu.add_cascade(label="Recovery", menu=self.recovery_menu)

        ## Enrichment ##
        self.enrichment_menu = Menu(self.main_menu, tearoff=0)
        self.enrichment_menu.add_command(label="Output Length Enrichments",
                                         command = lambda: self.compare_designs.output_len_or_clus_enrichment("length", native_path = self.native_path))
        self.enrichment_menu.add_command(label="Output Cluster Enrichments",
                                         command = lambda: self.compare_designs.output_len_or_clus_enrichment("cluster", native_path = self.native_path))

        self.main_menu.add_cascade(label="Enrichment", menu=self.enrichment_menu)

        self.subset_menu = Menu(self.main_menu, tearoff=0)

        # Have to do this manually:
        # for score_name in :
        #    x = copy.deepcopy(score_name)
        #    self.subset_menu.add_command(label = "Create DB of Top Subset: "+score_name, command = lambda: self.create_subset_databases(x))

        score_names = self.compare_designs._get_score_names()
        self.subset_menu.add_command(label="Create DB of Top Subset: " + score_names[0],
                                     command=lambda: self.create_subset_databases(score_names[0]))
        self.subset_menu.add_command(label="Create DB of Top Subset: " + score_names[1],
                                     command=lambda: self.create_subset_databases(score_names[1]))
        self.subset_menu.add_command(label="Create DB of Top Subset: " + score_names[2],
                                     command=lambda: self.create_subset_databases(score_names[2]))
        self.subset_menu.add_command(label="Create DB of Top Subset: " + "Top N dG of Top Total",
                                     command=lambda: self.create_subset_databases(score_names[3]))

        self.main_menu.add_cascade(label="Features Subsets", menu=self.subset_menu)

        self._tk_.config(menu=self.main_menu)



        ########### Callbacks ############

    def print_threads(self):
        print "Total running threads: "+repr(threads.n_alive())
        for pid in range(0, len(threads)):
            if threads.is_alive(pid):
                print "Thread "+repr(pid)+" is alive."

    def show_strat_items(self):
        item = self.all_strategies_listbox.get(self.all_strategies_listbox.curselection())
        items = glob.glob(self.compare_designs.db_dir.get() + "/*" + item + "*")
        # for i in items:
        # print i

        if os.path.exists(self.compare_designs.db_dir.get() + "/databases"):
            print "\n Databases:"
            dbs = glob.glob(self.compare_designs.db_dir.get() + "/databases/*" + item + "*")
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

        ######### Tracers ##############

    def camelid_tracer(self, name, index, mode):
        varValue = self.compare_designs.is_camelid.get()
        if varValue == 1:
            self.compare_designs.cdrs["L1"].set(0)
            self.compare_designs.cdrs["L2"].set(0)
            self.compare_designs.cdrs["L3"].set(0)

        return


        ######### Auxilliary Functions ###########

    def add_main_strategy(self):
        strategy_name = tkSimpleDialog.askstring(title="Strategy", prompt="Strategy Name")
        if not strategy_name:
            return

        strategy_path = tkFileDialog.askdirectory(initialdir=self.current_dir, title="Strategy Path")
        self.compare_designs.strategies.append(strategy_name)
        self.compare_designs.db_paths[strategy_name] = strategy_path

        self.all_strategies_listbox.insert(END, strategy_name)

    def change_root(self):
        root = tkFileDialog.askdirectory(initialdir=self.current_dir, title="Root Directory")
        if not root:
            return
        self.current_dir = root
        os.chdir(root)
        print "Root directory changed to: " + root

    def read_from_db_dir_set_strategies(self):
        self.compare_designs.strategies = []

        if not self.compare_designs.db_dir.get() or not os.path.exists(self.compare_designs.db_dir.get()):

            db_dir = tkFileDialog.askdirectory(initialdir=self.current_dir, title="Database DIR")
            if not db_dir:
                return
            self.current_dir = db_dir
            self.compare_designs.db_dir.set(db_dir)
        self.compare_designs.set_strategies_from_databases()
        self.populate_all_strategies()

        strategies = self.get_full_strategy_list()
        self.compare_designs.set_strategies(strategies)

    def get_full_strategy_list(self):
        strategies = self.current_strategies_listbox.get(0, END)
        return strategies

    def set_reference_db(self):
        d = tkFileDialog.askopenfilename(title="Reference DB", initialdir=self.compare_designs.db_dir.get())
        if not d: return
        self.current_dir = os.path.dirname(d)
        self.compare_designs.reference_db.set(d)

    def set_top_n(self):
        top = tkSimpleDialog.askinteger(title="Top N", prompt="Number of top scoring",
                                        initialvalue=self.compare_designs.top_n.get())
        if not top: return
        self.compare_designs.top_n.set(top)

    def set_top_n_combined(self):
        top = tkSimpleDialog.askinteger(title="Top N Combined", prompt="Number of top scoring Combined",
                                        initialvalue=self.compare_designs.top_n_combined.get())
        if not top: return
        self.compare_designs.top_n_combined.set(top)

    def set_max_clustal_procs(self):
        max = tkSimpleDialog.askinteger(title="Max P", prompt="Max NP.  Clustal by default uses all.",
                                        initialvalue=self.clustal_procs.get())
        if not max: return
        self.clustal_procs.set(max)

    def set_clustal_output_format(self):
        f = tkSimpleDialog.askstring(title="Clustal output format", initialvalue=self.clustal_output_format.get())
        if not f:
            return
        if not f in self.clustal_output_formats:
            print "Format " + f + " not recognized.  Available formats are: \n" + repr(self.clustal_output_formats)
            return

        self.clustal_output_format.set(f)

    def set_clustal_soft_wrap(self):
        wrap = tkSimpleDialog.askinteger(title="Wrap", prompt="Set Soft Wrap",
                                         initialvalue=self.compare_designs.clustal_soft_wrap.get())
        if not wrap:
            return
        self.compare_designs.clustal_soft_wrap.set(wrap)

        # def set_scorefunction(self):
        # score = tkSimpleDialog.askstring(title="Score", prompt = "Set Scorefunction", initialvalue = self.compare_designs.scorefunction.get())
        # if not score:
        # return
        # self.compare_designs.scorefunction.set(score)

    def set_native_path(self):
        native_path = tkFileDialog.askopenfilename(title="Native path", initialdir=self.current_dir)
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

    def run_copy_all(self):
        self.compare_designs.copy_all_models(self.native_path)


    def run_clustal_on_all_combined(self):

        extra_options = tkSimpleDialog.askstring(title="Extra Options", prompt="Clustal Extra Options",
                                                 initialvalue=self.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega_on_all_combined(self.clustal_procs.get(),
                                                               self.clustal_output_format.get(),
                                                               extra_options=extra_options, native_path=self.native_path)

    def run_clustal_omega(self):

        extra_options = tkSimpleDialog.askstring(title="Extra Options", prompt="Clustal Extra Options",
                                                 initialvalue=self.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega(self.clustal_procs.get(),
                                                               self.clustal_output_format.get(),
                                                               extra_options=extra_options, native_path=self.native_path)

    def create_subset_databases(self, score_name):
        rosetta_extension = tkSimpleDialog.askstring(title="Rosetta Extension",
                                                     prompt="Please set the Rosetta Extension",
                                                     initialvalue=self.compare_designs.rosetta_extension.get())
        if not rosetta_extension:
            return
        self.compare_designs.rosetta_extension.set(rosetta_extension)

        prefix = tkSimpleDialog.askstring(title="Prefix",
                                          prompt="Please set the prefix that will be used for the new databases",
                                          initialvalue=score_name)
        if not prefix:
            print "The prefix needs to be set"
            return

        self.compare_designs.create_score_subset_database(score_name, prefix)
        # self.read_from_main_set_strategies()


class CompareAntibodyDesignStrategies:
    """
    Class mainly for comparing different Antibody Design strategies using our Features Databases.
    """

    def __init__(self, db_dir, out_dir_name, strategies=[]):

        #Init construction options
        self.db_dir = StringVar(value=db_dir)
        self.out_dir_name = StringVar(value=out_dir_name)

        self.strategies = strategies

        #Init Classes and data
        self.db_paths = defaultdict()
        self.filter_settings = FilterSettings()

        #Init Components
        self._init_default_options()
        self._init_scores()
        self._init_default_scores()
        self._init_cdrs()



    def _init_default_options(self):

        self.scorefxn = "talaris2013"

        self.main_dir = StringVar()

        self.reference_db = StringVar()
        self.clustal_soft_wrap = IntVar(value=100)
        self.reload_scores = IntVar(value=1)
        self.top_n = IntVar(value=10)
        self.top_n_combined = IntVar(value=15)

        self.features_hbond_set = IntVar()
        self.features_hbond_sets = ["", "_min_hbond_analysis", "_no_hbond_analysis"]
        self.features_hbond_set.set(1)
        self.query_hbonds = IntVar(value=0)

        self.is_camelid = IntVar();
        self.is_camelid.set(0)
        self.top_total_percent = IntVar();
        self.top_total_percent.set(10)
        self.backround_features = IntVar();
        self.backround_features.set(1)

        self.rosetta_extension = StringVar();
        self.rosetta_extension.set("linuxclangrelease")

        self.individual_analysis = IntVar(value = 1);
        self.combined_analysis = IntVar(value = 0);

    def _init_scores(self):
        total_scores = TotalDecoyData()
        dg_scores = dGDecoyData()
        dsasa_scores = dSASADecoyData()
        top10_by_10 = dGTotalScoreSubset()
        self.scores = [total_scores, dg_scores, dsasa_scores, top10_by_10]

    def _init_default_scores(self):
        #Setup the scores that are on
        self.score_names = self._get_score_names()
        self.scores_on = defaultdict()
        for score_name in self.score_names:
            self.scores_on[score_name] = IntVar(value=0)

        self.scores_on["dG"].set(1)
        self.scores_on["dG_top_Ptotal"].set(1)
        self.scores_on["total"].set(1)

    def _init_cdrs(self):
        self.cdrs = defaultdict()
        self.cdrs["L1"] = IntVar(value=1)
        self.cdrs["L2"] = IntVar(value=1)
        self.cdrs["L3"] = IntVar(value=1)
        self.cdrs["H1"] = IntVar(value=1)
        self.cdrs["H2"] = IntVar(value=1)
        self.cdrs["H3"] = IntVar(value=1)

    def _get_score_names(self):

        names = []
        for score in self.scores:
            names.append(score.name)
        return names

    def _setup_outdir(self, subdirs=[], use_out_dir_name=True):
        """
        Sets up the main output dir in the main_analys_dir, and any subdirectories such as 'decoys' or decoys/combined_3
        Returns the final output directory
        """

        filters = self._setup_filters()

        # if self.out_dir_name.get() and use_out_dir_name and filters:
        #    outdir = self.main_dir.get()+"/"+self.out_dir_name.get()+"_"+self.filter_settings.name.get()
        # elif self.out_dir_name.get() and use_out_dir_name:
        #    outdir = self.main_dir.get()+"/"+self.out_dir_name.get()
        # elif filters:
        #    outdir = self.out_dir_name.get()+"/"+self.filter_settings.name.get()
        # else:

        outdir = os.getcwd()

        if not os.path.exists(outdir): os.mkdir(outdir)

        for subdir in subdirs:
            if not subdir: continue
            outdir = outdir + "/" + subdir
            if not os.path.exists(outdir): os.mkdir(outdir)

        return outdir

    def _setup_outdir_individual(self, subdirs, use_outdir_name=False):

        return self._setup_outdir(["analysis_individual", self.out_dir_name.get()] + subdirs, use_outdir_name)

    def _setup_outdir_combined(self, subdirs):
        return self._setup_outdir(["analysis_combined", self.out_dir_name.get()] + subdirs, False)

    def _setup_cdrs(self):
        """
        Get a list of CDRs to process for [many] non-features related tasks.
        """

        cdrs = []
        camelid_cdrs = ["L1", "L2", "L3"]

        for cdr_name in self.cdrs:
            if self.is_camelid.get() and cdr_name in camelid_cdrs:
                continue

            if self.cdrs[cdr_name].get():
                cdrs.append(cdr_name)

        return cdrs

    def _setup_scores(self, features_type="antibody", use_all=False):
        """
        Setup the Score Classes.  If not use_all, will use only use those set.
        """

        self._init_scores()
        score_subset = []
        query_hbonds = self.query_hbonds.get()

        for score in self.scores:
            if score.name == "hbond_count" or score.name == "hbond_energy":

                if self.scores_on[score.name].get():
                    query_hbonds = True

        # Quickly analyze scores with the same settings...
        if not self.reload_scores.get():
            for score_class in self.scores:
                if self.scores_on[score_class.name].get() or use_all:
                    score_subset.append(copy.deepcopy(score_class))
            return score_subset

        hb_loader = InterfaceHBondDecoyDataLoader();

        filters = self._setup_filters()

        if self.is_camelid.get():
            for score in self.scores:
                score.set_interface('H_A')

        if filters:
            self.scores.append(CombinedStrDecoyData(filters, self.filter_settings.name.get()))

        for strategy in self.strategies:
            db_path = self.get_db_path(strategy, features_type)
            if not os.path.exists(db_path):
                sys.exit(
                    "DB path does not exist.  Please Run Features reporter with StructureScores and ScoreTypes for this strategy\n" + db_path)
            print db_path
            con = sqlite3.connect(db_path)

            if query_hbonds:
                if filters:
                    hb_loader.add_filters(filters, self.filter_settings.name.get())
                hb_loader.add_data(strategy, con)

            for score_class in self.scores:

                if filters:
                    score_class.add_filters(filters, self.filter_settings.name.get())

                # We set these up later
                if score_class.name == "hbond_count" or score_class.name == "hbond_energy":
                    continue

                if score_class.name == "dG_top_Ptotal":
                    score_class.add_data(strategy, con, self.top_total_percent.get())
                else:
                    score_class.add_data(strategy, con)

        # Setup each Hbond Data class
        if query_hbonds:
            for score_class in self.scores:
                if score_class.name == "hbond_count" or score_class.name == "hbond_energy":
                    score_class.setup_from_loader(hb_loader)

        for score_class in self.scores:
            if self.scores_on[score_class.name].get() or use_all:
                score_subset.append(copy.deepcopy(score_class))

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
            custom_filter = DataFilter("custom_filter", type="unknown")
            custom_filter.required_tables = self.filter_settings.extra_required_tables
            custom_filter.required_wheres = self.filter_settings.extra_required_where
            filters.append(custom_filter)

        energy_filters = defaultdict()
        energy_filters["dG"] = dGCutoffFilter(0)
        energy_filters["dSASA"] = dSASACutoffFilter(0)
        energy_filters["total"] = TotalScoreCutoffFilter(0)

        for energy_type in self.filter_settings.energy_types:
            if not self.filter_settings.get_energy_enabled(energy_type): continue

            filter = energy_filters[energy_type]
            filter.set_value(self.filter_settings.get_energy_cutoff(energy_type))
            filters.append(filter)

        return filters

    def _setup_cdr_types(self, data_type, is_camelid, features_type, native_path=None):
        if data_type == "length":
            data_class = CDRLengthData(native_path, is_camelid)

        elif data_type == "cluster":
            data_class = CDRClusterData(native_path, is_camelid)

        elif data_type == "sequence":
            data_class = CDRSequenceData(native_path, is_camelid)

        elif data_type == "aligned_sequence":
            data_class = CDRAlignedSequenceData(self._setup_outdir_individual(['clustal']),
                                                self._setup_outdir_combined(['clustal']), native_path, is_camelid)

        else:
            sys.exit("data_type not understood!")

        for strategy in self.strategies:
            db_path = self.get_db_path(strategy, features_type)
            if not os.path.exists(db_path):
                sys.exit(
                    "DB path does not exist.  Please Run Features reporter with StructureScores and ScoreTypes for this strategy\n" + db_path)
            print db_path
            con = sqlite3.connect(db_path)
            data_class.add_data(strategy, con)

        return data_class

    def set_cdrs_from_list(self, cdr_list):
        for name in self.cdrs:
            if name in cdr_list:
                self.cdrs[name].set(1)
            else:
                self.cdrs[name].set(0)

    def set_strategies_from_db_dir_top_dir(self):
        dirs = glob.glob(self.db_dir.get() + "/TOP_*")
        non_redun = defaultdict()
        for d in sorted(dirs):
            # print d
            non_redun["_".join(d.split("_")[3:])] = None

        for x in sorted(non_redun.keys()):
            # print x
            self.strategies.append(x)

    def set_strategies_from_databases(self):
        """
        Set the strategies from the db_dir/databases directory
        :return:
        """
        if not os.path.exists(self.db_dir.get()):
            print "Could not find any databases to use for analysis.  Please make sure databases are in " + self.db_dir.get()
            print "Databases should end with .db or .db3 extension.  Naming format of databases: strategy.features_type.scorefunction.db"
            print "WTF?"
            return

        dbs = glob.glob(self.db_dir.get() + "/*.db*")

        if len(dbs) == 0:
            print "Could not find any databases to use for analysis.  Please make sure databases are in " + self.db_dir.get()
            print "Databases should end with .db or .db3 extension.  Naming format of databases: strategy.features_type.scorefunction.db"
            return

        print repr(dbs)
        nr_dbs = defaultdict()
        for db in sorted(dbs):
            print db
            x = os.path.basename(db)

            ##Example naming convention: 'ch103_5_CDR_prelim.norm_ab_features.db'
            strategy = '.'.join(x.split('.')[:-2])
            print strategy
            if len(x.split('.')[-2].split('_')) > 2:
                strategy = strategy + '.' + '_'.join(x.split('.')[-2].split('_')[:-2])
                print strategy

            if not nr_dbs.has_key(strategy):

                self.strategies.append(strategy)
                nr_dbs[strategy] = " "
                self.db_paths[strategy] = []
                self.db_paths[strategy].append(db)
            else:
                self.db_paths[strategy].append(db)

    def set_strategies(self, strategies):
        self.strategies = strategies

    def get_strategies(self):
        return self.strategies

    def get_db_path(self, strategy, features_type='antibody'):

        names = {
            'antibody': ['ab', 'antibody'],
            'cluster': ['cl', 'cluster', 'ab', 'antibody']
        }

        for p in self.db_paths[strategy]:
            for match in names[features_type]:
                if re.search(match, p):
                    return p

        print "Matching database name not found for features type: " + features_type
        print "Database must have any of these names in them: " + repr(names[features_type])

    def get_full_features_type(self, type):
        if type == "cluster":
            return type
        else:
            return type + "_minimal" + self.features_hbond_sets[self.features_hbond_set.get()]

    ################ Main Functions #######################
    def run_features(self, type, plot_name=""):

        if not plot_name:
            plot_name = self.out_dir_name.get()

        if not plot_name:
            print "No root name set!"
            return

        os.system("rm build")
        outdir = self._setup_outdir(["plots", plot_name], False)
        db_dir = self.db_dir.get()

        if len(self.strategies) == 0:
            return

        if not os.path.exists(db_dir): sys.exit("Please run Features reporter and copy databases to db directory.")

        fulltype = self.get_full_features_type(type)
        creator = json_creator.JsonCreator(outdir, fulltype)

        if self.reference_db.get() and os.path.exists(self.reference_db.get()):
            creator.add_sample_source_info(self.reference_db.get(), "ref", True)

        for strategy in self.strategies:
            id = strategy
            id = id.replace("talaris2013", "talaris")

            # Fix up the name so it is not too long.  Will add options to do this manually later:
            """
            id = id.replace("cluster", "clus")
            id = id.replace("exclude", "excl")
            id = id.replace("include", "incl")

            """

            db_path = self.get_db_path(strategy, type)
            if not os.path.exists(db_path):
                sys.exit(db_path + " does not exist!")

            creator.add_sample_source_info(db_path, id)

        creator.save_json(outdir + "/" + type + "_" + plot_name + ".json")
        creator.run_json(self.backround_features.get())

        # pwd = os.getcwd()
        # os.chdir("build")
        # os.system("cp -r *../"+outdir)
        # os.system("mv build build_old")
        # os.chdir(pwd)

        print "Plots ignore any set filters.  To plot with filters, create new databases through query..."
        print "Complete..."

    def output_stats(self):

        def output_all_stats():

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

            # Output Strategy data:
            for score_class in scores:
                if isinstance(score_class, DecoyData): pass

                if score_class.name == "dSASA":
                    reverse = True
                else:
                    reverse = False

                for strategy in self.strategies:
                    OUTFILE = open(outdir_lists + "/" + strategy + "_ORDERED_" + score_class.get_outname() + ".txt", 'w')
                    data = score_class.get_strategy_data(strategy, True)

                    for tup in sorted(data.keys(), reverse=reverse):
                        triple = data[tup]
                        OUTFILE.write(
                            get_str(triple.score) + "\t" + triple.strategy + "\t" + os.path.basename(triple.decoy) + "\n")
                    OUTFILE.close()


            # This is broken for some unkown reason. - Ends up skipping pretty much everything.

            # Output Combined Data:
            COMBINED = open(outdir_stats + "/combined_selection_data.txt", 'w')
            COMBINED.write("#strategy decoy " + " ".join([score_class.get_outname() for score_class in scores]) + "\n")
            total_scores = scores[0]
            all_data = total_scores.get_concatonated_map()
            for decoy in all_data:
                strategy = all_data[decoy].strategy
                line = strategy + " " + os.path.basename(decoy)
                # print strategy
                for score_class in scores:
                    if score_class.name == "combined_str_score" or score_class.name == "dG_top_Ptotal": continue

                    if not score_class.all_data[strategy].has_key(decoy):
                        print "Skipping " + strategy + " " + decoy
                        continue

                    line = line + " " + get_str(score_class.all_data[strategy][decoy].score)
                COMBINED.write(line + "\n")
            COMBINED.close()


            # Output Stats:
            STATS = open(outdir_stats + "/combined_selection_data_stats.txt", 'w')
            STATS.write("#strategy " + " ".join(
                [score_class.get_outname() + "_avg" + " " + score_class.get_outname() + "_sd " for score_class in
                 scores]) + "\n")

            for strategy in self.strategies:
                line = strategy
                for score_class in scores:
                    if not score_class.has_real_values() or score_class.name == "combined_str_score": continue
                    raw_scores_tuple = score_class.get_strategy_data(strategy, True)
                    raw_scores = [s[0] for s in raw_scores_tuple]
                    m = numpy.mean(raw_scores)
                    sd = numpy.std(raw_scores)
                    line = line + " %.3f" % m + " " + "%.3f" % sd
                STATS.write(line + "\n")
            STATS.close()

            # Output Top Stats:
            STATS = open(outdir_stats + "/combined_selection_data_top_" + repr(top_n) + "_stats.txt", 'w')
            STATS.write("#strategy " + " ".join(
                [score_class.get_outname() + "_avg" + " " + score_class.get_outname() + "_sd " for score_class in
                 scores]) + "\n")
            for strategy in self.strategies:
                line = strategy
                for score_class in scores:
                    if not score_class.has_real_values() or score_class.name == "combined_str_score": continue
                    raw_scores_tuple = score_class.get_top_strategy_data(strategy, top_n, True)
                    raw_scores = [s[0] for s in raw_scores_tuple]
                    m = numpy.mean(raw_scores)
                    sd = numpy.std(raw_scores)
                    line = line + " %.3f" % m + " " + "%.3f" % sd
                STATS.write(line + "\n")
            STATS.close()
            print "Complete"

        def output_score_extra_stats():
            top_n = self.top_n.get()
            main_scores = self._setup_scores()
            all_scores = self._setup_scores("antibody", True)

            if self.individual_analysis.get():
                out_dir = self._setup_outdir_individual(["raw_score_data"])

                # Top, then all
                for strategy in self.strategies:
                    for score in main_scores:
                        if score.name == "combined_str_score":
                            continue

                        if isinstance(score, DecoyData): pass

                        # Get TopN of that particular score
                        decoy_list = score.get_ordered_decoy_list(strategy, top_n)

                        outfile = out_dir + "/scores_of_top_" + score.name + "_" + strategy + ".txt"
                        print "writing " + outfile
                        OUT = open(outfile, 'w')

                        if not score.name == "dG_top_Ptotal":
                            header = "#decoy\t" + score.name
                        else:
                            header = "#decoy"

                        for a_score in all_scores:
                            if a_score.name == "combined_str_score" or a_score.name == "dG_top_Ptotal": continue
                            if a_score.name == score.name: continue

                            header = header + "\t" + a_score.name

                        OUT.write(header + "\n")

                        for decoy in decoy_list:

                            if not score.name == "dG_top_Ptotal":
                                line = os.path.basename(decoy) + "\t" + get_str(score.get_score_for_decoy(strategy, decoy))
                            else:
                                line = os.path.basename(decoy)

                            for a_score in all_scores:
                                if a_score.name == "combined_str_score": continue
                                if a_score.name == score.name: continue

                                line = line + "\t" + get_str(a_score.get_score_for_decoy(strategy, decoy))
                            OUT.write(line + "\n")
                        OUT.close()

                # All decoys

                score = main_scores[0]

                for strategy in self.strategies:
                    all_decoys = score.get_ordered_decoy_list(strategy)

                    out_file = out_dir + "/all_scores_" + strategy + ".txt"
                    print "writing " + out_file

                    OUT = open(out_file, 'w')

                    if not score.name == "dG_top_Ptotal":
                        header = "#decoy\t" + score.name
                    else:
                        header = "#decoy"

                    for a_score in all_scores:
                        if a_score.name == "combined_str_score" or a_score.name == "dG_top_Ptotal": continue
                        if a_score.name == score.name: continue

                        header = header + "\t" + a_score.name
                    OUT.write(header + "\n")

                    for decoy in all_decoys:
                        if not score.name == "dG_top_Ptotal":

                            line = os.path.basename(decoy) + "\t" + get_str(score.get_score_for_decoy(strategy, decoy))
                        else:
                            line = os.path.basename(decoy)
                        for a_score in all_scores:
                            if a_score.name == "combined_str_score": continue
                            if a_score.name == score.name: continue

                            line = line + "\t" + get_str(a_score.get_score_for_decoy(strategy, decoy))
                        OUT.write(line + "\n")
                    OUT.close()

        if self.individual_analysis.get():
            print "Outputting individual stats"
            output_score_extra_stats()
        if self.combined_analysis.get():
            print "Outputting combined stats"
            output_all_stats()
        print "Complete"

    def copy_top(self, native_path=None):
        def copy_top_strategy(native_path=None):

            top_n = self.top_n.get()
            scores = self._setup_scores()

            if not self.out_dir_name.get():
                print "No root name set!"
                return

            print "Copying Top Models.."
            # Each Strategy Top Scoring (Skip total score here for now):
            for strategy in self.strategies:
                for score in scores:
                    if score.name == "combined_str_score":
                        continue

                    out_dir = self._setup_outdir_individual(
                        ["pdbs_sessions", "top_" + repr(top_n) + "_" + score.get_outname() + "_" + strategy])
                    SCORELIST = open(out_dir + "/MODELS.txt", 'w')
                    print "Copying " + strategy + " " + score.get_outname() + " to: " + out_dir
                    if isinstance(score, DecoyData): pass

                    decoys = score.get_top_strategy_data(strategy, top_n)
                    decoy_list = score.get_ordered_decoy_list(strategy, top_n)
                    load_as = []
                    i = 1
                    for decoy in decoy_list:
                        load_as.append("model_" + repr(i) + "_" + score.get_outname() + "_" + get_str(decoys[decoy].score))
                        os.system('cp ' + decoy + " " + out_dir + "/" + "top_" + repr(i) + "_" + os.path.basename(decoy))
                        SCORELIST.write(
                            repr(i) + "\t" + get_str(decoys[decoy].score) + "\t" + os.path.basename(decoy) + "\n")
                        i += 1
                    make_pymol_session_on_top(decoy_list, load_as, out_dir, out_dir, score.get_outname(), top_n,
                                              native_path)
                    SCORELIST.close()

        def copy_top_combined(native_path=None):
            """
            Outputs total_score,
            """
            if not self.out_dir_name.get():
                print "No root name set!"
                return

            top_n = self.top_n_combined.get()
            scores = self._setup_scores()

            # Overall Strategy:
            for score in scores:
                if score.name == "combined_str_score":
                    continue

                if isinstance(score, DecoyData): pass
                outdir_top_pdbs = self._setup_outdir_combined(["top_structures", score.get_outname()])
                outdir_top_sessions = self._setup_outdir_combined(["top_sessions"])

                SCORELIST = open(outdir_top_pdbs + "/MODELS.txt", 'w')
                print "Copying " + score.get_outname() + " to: " + outdir_top_pdbs
                decoys = score.get_top_all_data(top_n)
                decoy_list = score.get_ordered_decoy_list_all(top_n)
                load_as = []
                i = 1
                for decoy in decoy_list:
                    load_as.append("model_" + repr(i) + "_" + score.get_outname() + "_" + get_str(decoys[decoy].score))
                    os.system('cp ' + decoy + " " + outdir_top_pdbs + "/top_" + repr(i) + "_" + os.path.basename(decoy))
                    SCORELIST.write(repr(i) + "\t" + get_str(decoys[decoy].score) + "\t" + os.path.basename(decoy) + "\n")
                    i += 1

                make_pymol_session_on_top(decoy_list, load_as, outdir_top_pdbs, outdir_top_sessions, score.get_outname(),
                                          top_n,
                                          native_path)
                SCORELIST.close()

        if self.individual_analysis.get():
            print "Outputting individual sessions"
            copy_top_strategy(native_path)
        if self.combined_analysis.get():
            print "Outputting combined sessions"
            copy_top_combined(native_path)

    def copy_all_models(self, native_path=None):

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        if len(self.strategies) == 0:
            print "No strategies set..."
            return


        # outdir = self._setup_outdir()
        scores = self._setup_scores()

        # Each Strategy Top Scoring (Skip total score here for now):
        for strategy in self.strategies:
            for score in scores[1:]:
                out_dir = self._setup_outdir_individual(["all" + "_" + score.get_outname() + "_" + strategy], False)
                print "Copying " + strategy + " " + score.get_outname() + " to: " + out_dir
                if isinstance(score, DecoyData): pass

                decoys = score.get_strategy_data(strategy)
                decoy_list = score.get_ordered_decoy_list(strategy)
                load_as = []
                i = 1
                for decoy in decoy_list:
                    load_as.append("model_" + repr(i) + "_" + score.get_outname() + "_" + get_str(decoys[decoy].score))
                    os.system('cp ' + decoy + " " + out_dir + "/top_" + repr(i) + "_" + os.path.basename(decoy))
                    i += 1
                make_pymol_session_on_top(decoy_list, load_as, out_dir, out_dir, score.get_outname(), None, native_path)

        # Overall Strategy:
        for score in scores:
            if isinstance(score, DecoyData): pass
            out_dir = self._setup_outdir_combined(["all_structures", score.get_outname()])
            print "Copying " + score.get_outname() + " to: " + out_dir
            decoys = score.get_concatonated_map()
            decoy_list = score.get_ordered_decoy_list_all()
            load_as = []
            i = 1
            for decoy in decoy_list:
                load_as.append("model_" + repr(i) + "_" + score.get_outname() + "_" + get_str(decoys[decoy].score))
                os.system('cp ' + decoy + " " + out_dir + "/top_" + repr(i) + "_" + os.path.basename(decoy))
                i += 1
            session_dir = out_dir = self._setup_outdir_combined(["all_sessions"])
            make_pymol_session_on_top(decoy_list, load_as, out_dir, session_dir, score.get_outname(), None, native_path)

    def run_clustal_omega(self, processors, output_format="fasta", extra_options="", native_path=None):

        def run_clustal_omega_on_strategies():

            if not self.out_dir_name.get():
                print "No root name set!"
                return

            scores = self._setup_scores()
            top_n = self.top_n.get()

            for strategy in self.strategies:
                print "\nRunning Clustal on: " + strategy
                score_zero = scores[0]
                if isinstance(score_zero, DecoyData): pass

                # Output per strategy ALL
                decoys = score_zero.get_strategy_data(strategy)
                decoy_header_dict = defaultdict()

                for decoy in decoys:
                    decoy_header_dict[decoy] = os.path.basename(decoy)

                fasta_dir = self._setup_outdir_individual(["sequences"])
                fasta_path = fasta_dir + "/" + strategy + "_all.fasta"
                fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native",
                                                       self.is_camelid.get())

                aln_dir = self._setup_outdir_individual(["clustal"])
                aln_name = strategy + "_all.clus"
                clustal_runner = ClustalRunner(fasta_path)
                clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
                clustal_runner.set_threads(processors)
                clustal_runner.set_extra_options(extra_options)
                clustal_runner.set_output_format(output_format)
                clustal_runner.output_alignment(aln_dir, aln_name)

                # Output on Top Scoring:

                for score in scores:
                    print score
                    if isinstance(score, DecoyData): pass

                    basename = strategy + "_" + score.get_outname() + "_top_" + str(top_n)
                    fasta_path = fasta_dir + "/" + basename + ".fasta"
                    aln_name = basename + ".clus"

                    decoy_header_dict = defaultdict()
                    decoys = score.get_top_strategy_data(strategy, top_n)
                    for decoy in decoys:
                        header = "v" + get_str(decoys[decoy].score) + "::" + os.path.basename(decoy)
                        decoy_header_dict[decoy] = header
                    fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native",
                                                           self.is_camelid.get())
                    clustal_runner.set_fasta_path(fasta_path)
                    clustal_runner.set_threads(processors)
                    clustal_runner.set_extra_options(extra_options)
                    clustal_runner.set_output_format(output_format)
                    clustal_runner.output_alignment(aln_dir, aln_name)

            print "Complete"

        def run_clustal_omega_on_top_combined():

            if not self.out_dir_name.get():
                print "No root name set!"
                return

            scores = self._setup_scores()
            top_n = self.top_n.get()

            for score in scores:
                print "Running Clustal Omega on " + score.get_outname()
                if isinstance(score, DecoyData): pass
                decoys = score.get_top_all_data(top_n, False)

                decoy_header_dict = defaultdict()
                for decoy in decoys:
                    header = "v" + get_str(decoys[decoy].score) + "::" + os.path.basename(decoy)
                    decoy_header_dict[decoy] = header

                root_name = self.out_dir_name.get()
                basename = root_name + "_" + score.get_outname() + "_top_" + repr(top_n)
                fasta_dir = self._setup_outdir_combined(["sequences"])
                fasta_path = fasta_dir + "/" + basename + ".fasta"

                clustal_dir = self._setup_outdir_combined(["clustal"])
                clustal_name = basename + ".aln"

                fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "Native",
                                                       self.is_camelid.get())
                clustal_runner = ClustalRunner(fasta_path)
                clustal_runner.set_threads(processors)
                clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
                clustal_runner.set_extra_options(extra_options)
                clustal_runner.set_output_format(output_format)
                clustal_runner.output_alignment(clustal_dir, clustal_name)

        if self.individual_analysis.get():
            run_clustal_omega_on_strategies()
        if self.combined_analysis.get():
            run_clustal_omega_on_top_combined()

        print "Complete"

    def run_clustal_omega_on_all_combined(self, processors, output_format, extra_options="", native_path=None):

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        print "Running Clustal on All Combined"
        scores = self._setup_scores()
        score = scores[0]
        if isinstance(score, DecoyData): pass

        root_name = self.out_dir_name.get()
        fasta_dir = self._setup_outdir_combined(["sequences"])
        fasta_path = fasta_dir + "/" + root_name + "_all.fasta"

        clustal_dir = self._setup_outdir_combined(["clustal"])
        clustal_name = root_name + "_all.aln"

        all_data = score.get_concatonated_map(False)

        all_data_array = []
        for s in scores:
            if s.name == "combined_str_score":
                continue
            all_data_array.append(s.get_concatonated_map(False))

        decoy_header_dict = defaultdict()
        for decoy in all_data:
            a = all_data_array[0]
            header = "v" + get_str(a[decoy].score)

            for a in all_data_array[1:]:
                header = header + ":" + get_str(a[decoy].score)
            header = header + ":" + os.path.basename(decoy)
            decoy_header_dict[decoy] = header

        fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, native_path, "native",
                                               self.is_camelid.get())
        clustal_runner = ClustalRunner(fasta_path)
        clustal_runner.set_threads(processors)
        clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
        clustal_runner.set_extra_options(extra_options)
        clustal_runner.set_output_format(output_format)
        clustal_runner.output_alignment(clustal_dir, clustal_name)

        print "Complete"

    def output_len_or_clus_alignment(self, alignment_type, features_type='antibody',
                                     native_path=None):

        is_camelid = self.is_camelid.get()

        top_n = self.top_n.get()
        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type, native_path)
        len_class = self._setup_cdr_types("length", self.is_camelid.get(), features_type, native_path)

        def _output_alignment(self, outdir, top_decoys, score, type_data, strategy=None, extra_name="top"):

            if isinstance(type_data, CDRData): pass
            if isinstance(score, DecoyData): pass

            if not strategy:
                outname = "cdr_type_alignments_" + alignment_type + "_" + score.get_outname() + "_" + extra_name + ".txt"
            else:
                outname = "cdr_type_alignments_" + alignment_type + "_" + score.get_outname() + "_" + strategy + "_" + extra_name + "_.txt"

            print "Outputting cdr type alignment: " + outdir + "/" + outname

            cdr_names = [cdr for cdr in type_data.cdrs if cdr in self._setup_cdrs()]

            OUTFILE = open(outdir + "/" + outname, 'w')
            if native_path:
                header = "#score\t\tmatches"
            else:
                header = "#score\t"
            for cdr in cdr_names:
                header = header + "\t\t" + cdr

            header += "\tdecoy"
            OUTFILE.write(header + "\n")

            all_data = score.get_concatonated_map()
            all_type_data = type_data.get_concatonated_map()

            if native_path:
                line = "..." + "\t\tNA"

                if all_type_data.has_key((cdr, "native")):
                    native_info = type_data[(cdr, "native")]
                else:
                    native_info = type_data.get_native_data()

                if isinstance(native_info, CDRDataInfo): pass

                for cdr in cdr_names:
                    info = str(native_info.get_value_for_cdr(cdr))
                    if info == "NA":
                        info = cdr + "-" + str(len_class.get_native_data().get_value_for_cdr(cdr)) + "-NA"
                    line = line + "\t\t" + str(info)
                line = line + "\tnative"
                OUTFILE.write(line + "\n")

                for decoy in top_decoys:
                    decoy_info = all_type_data[decoy]
                    score_info = all_data[decoy]
                    counts = count_native_matches(decoy_info, native_info, type_data.cdrs)
                    line = get_str(score_info.score) + "\t\t" + repr(counts)
                    for cdr in cdr_names:
                        line = line + "\t\t" + str(get_star_if_native(decoy_info, native_info, cdr))

                    line = line + "\t" + os.path.basename(decoy)
                    OUTFILE.write(line + "\n")
            else:
                for decoy in top_decoys:
                    print decoy
                    decoy_info = all_type_data[decoy]
                    score_info = all_data[decoy]
                    print repr(score_info)
                    print repr(score_info.score)

                    line = get_str(score_info.score)
                    for cdr in cdr_names:

                        info = str(decoy_info.get_value_for_cdr(cdr))
                        if info == "NA":
                            info = cdr + "-" + str(
                                len_class.get_concatonated_map()[decoy].get_value_for_cdr(cdr)) + "-NA"
                        line = line + "\t\t" + str(info)

                    line = line + "\t" + os.path.basename(decoy)
                    OUTFILE.write(line + "\n")

            OUTFILE.close()

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        scores = self._setup_scores(features_type)
        if isinstance(data_class, CDRData): pass

        ### Top/All  each Strategy
        if self.individual_analysis.get():
            outdir = self._setup_outdir_individual(['cdr_alignments'])
            for strategy in self.strategies:

                for score in scores:
                    if isinstance(score, DecoyData): pass
                    top_decoys = score.get_ordered_decoy_list(strategy, top_n)
                    all_decoys = score.get_ordered_decoy_list(strategy)
                    # print "Top: "+repr(top_decoys)

                    print "Working on : " + strategy + " " + score.name
                    _output_alignment(self, outdir, top_decoys, score, data_class, strategy, "top")
                    _output_alignment(self, outdir, all_decoys, score, data_class, strategy, "all")


        ### Top/All each combined
        if self.combined_analysis.get():
            outdir = self._setup_outdir_combined(['cdr_alignments'])
            for score in scores:
                top_decoys = score.get_ordered_decoy_list_all(top_n)
                all_decoys = score.get_ordered_decoy_list_all()

                print "Top: " + repr(len(top_decoys)) + score.name
                _output_alignment(self, outdir, top_decoys, score, data_class, None, "top")
                _output_alignment(self, outdir, all_decoys, score, data_class, None, "all")

        print "Complete"

    def output_len_or_clus_recovery(self, alignment_type, features_type='antibody', native_path=None):
        len_class = self._setup_cdr_types("length", self.is_camelid.get(), features_type, native_path)
        if not native_path:
            print "Must pass select native path to calculate recoveries"
            return

        def _get_header(self, type_data):
            header = "#avg"
            cdr_names = [cdr for cdr in type_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                header = header + "\t\t" + cdr
            header += "\tname"
            return header

        def _add_native_line(self, native_info, type_data, OUTFILE):
            line = "NA"

            if isinstance(native_info, CDRDataInfo): pass

            cdr_names = [cdr for cdr in type_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                info = str(native_info.get_value_for_cdr(cdr))
                if info == "NA":
                    info = cdr + "-" + str(len_class.get_native_data().get_value_for_cdr(cdr)) + "-NA"
                line = line + "\t\t" + str(info)
            line += "\tnative"
            OUTFILE.write(line + "\n")

        def _add_recovery_line(self, label, decoys, type_data, OUTFILE, strategy=None):

            if isinstance(type_data, CDRData): pass

            total = 0
            cdr_names = [cdr for cdr in type_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                enrichment_data = calculate_recovery(type_data.get_native_data(), type_data, cdr, decoys)
                total = total + enrichment_data.get_perc_decimal()


            native_data = type_data.get_native_data()

            # If native is NA, do not count the recovery against it for averages.  Only average over set CDRs.
            t=0
            for cdr in cdr_names:
                if not native_data.get_value_for_cdr(cdr) =="NA":
                    t+=1

            avg = total / t

            line = "%.3f" % avg
            for cdr in cdr_names:
                enrichment_data = calculate_recovery(type_data.get_native_data(), type_data, cdr, decoys)

                line = line + "\t\t%.3f" % enrichment_data.get_perc_decimal()

            line += "\t"+label
            OUTFILE.write(line + "\n")

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        top_n = self.top_n.get()

        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type, native_path)
        scores = self._setup_scores(features_type)
        if isinstance(data_class, CDRData): pass

        # Each strategy + Top Values
        if self.individual_analysis.get():
            outdir = self._setup_outdir_individual(['enrichment'])
            OUTFILE = open(outdir + "/" + "cdr_type_recoveries_" + alignment_type + "_.txt", 'w')
            OUTFILE.write(_get_header(self, data_class) + "\n")
            _add_native_line(self, data_class.get_native_data(), data_class, OUTFILE)
            for strategy in self.strategies:
                outname = strategy
                decoys = scores[0].get_strategy_data(strategy).keys()
                _add_recovery_line(self, outname, decoys, data_class, OUTFILE, strategy)
                for score in scores:
                    if isinstance(score, DecoyData): pass
                    if score.name == "combined_str_score": continue

                    top_decoys = score.get_ordered_decoy_list(strategy, top_n)
                    outname = strategy + "_" + score.get_outname() + "_top_" + repr(top_n)
                    _add_recovery_line(self, outname, top_decoys, data_class, OUTFILE, strategy)
            OUTFILE.close()


        # Combined + Top Values
        if self.combined_analysis.get():
            outdir = self._setup_outdir_combined(['enrichment'])
            OUTFILE = open(outdir + "/" + "cdr_type_recoveries_all_" + alignment_type + "_.txt", 'w')
            OUTFILE.write(_get_header(self, data_class) + "\n")
            _add_native_line(self, data_class.get_native_data(), data_class, OUTFILE)
            outname = "combined_all"
            _add_recovery_line(self, outname, scores[0].get_concatonated_map().keys(), data_class, OUTFILE)
            for score in scores:
                if score.name == "combined_str_score": continue
                top_decoys = score.get_ordered_decoy_list_all(top_n)
                outname = "combined_" + score.get_outname() + "_top_" + repr(top_n)
                # print "Top: "+repr(top_decoys)
                _add_recovery_line(self, outname, top_decoys, data_class, OUTFILE)
            OUTFILE.close()
        print "Complete"

    def output_len_or_clus_enrichment(self, alignment_type, features_type='antibody', native_path = None):


        def _add_enrichments(self, label, decoys, type_data):
            if isinstance(data_class, CDRData): pass
            #OUTFILE = open(outdir + "/" + "cdr_type_enrichment_" + alignment_type + "_.txt", 'w')

            cdr_names = [cdr for cdr in type_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                enrichments = calculate_enrichments(data_class, cdr, decoys)
                OUTFILE = open(outdir+"/" + "cdr_type_enrichment_"+alignment_type+"_"+label+"_"+cdr+".txt", 'w')
                OUTFILE.write("group\tcount\tperc\n")


                for c in sorted([[enrichments[c].count, c] for c in enrichments], reverse=True):
                    t = c[1]; #c[0] is the actual counts we are sorting on.
                    perc = enrichments[t]
                    if isinstance(perc, Perc): pass

                    p = perc.get_perc_decimal()
                    OUTFILE.write(str(t)+"\t\t"+str(perc.count)+"\t\t"+perc.get_formated_perc(p)+"\n")

                OUTFILE.close()


        ##Individual Enrichments

        top_n = self.top_n.get()
        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type, native_path)

        scores = self._setup_scores(features_type)

        if self.individual_analysis.get():
            outdir = self._setup_outdir_individual(["enrichment"])
            if isinstance(data_class, CDRData): pass



            for strategy in self.strategies:
                outname = strategy
                decoys = scores[0].get_strategy_data(strategy).keys()
                _add_enrichments(self, outname, decoys, data_class)

                for score in scores:
                    if isinstance(score, DecoyData): pass
                    if score.name == "combined_str_score": continue

                    top_decoys = score.get_ordered_decoy_list(strategy, top_n)
                    outname = strategy + "_" + score.get_outname() + "_top_" + repr(top_n)
                    _add_enrichments(self, outname, top_decoys, data_class)


        #Combined Enrichments
        if self.combined_analysis.get():
            outdir = self._setup_outdir_combined(["enrichment"])
            outname = "combined_all"
            _add_enrichments(self, outname, scores[0].get_concatonated_map().keys(), data_class)
            for score in scores:
                if score.name == "combined_str_score": continue
                top_decoys = score.get_ordered_decoy_list_all(top_n)
                outname = "combined_" + score.get_outname() + "_top_" + repr(top_n)
                # print "Top: "+repr(top_decoys)
                _add_enrichments(self, outname, top_decoys, data_class)

        print "Complete"


    def create_score_subset_database(self, score_name, prefix, features_type='antibody'):
        self._setup_scores()
        score = self._get_score(score_name)
        if not isinstance(score, DecoyData):
            print "Score type not found!"
            return

        if not self.out_dir_name.get():
            print "No root name set!"
            return

        top_n = self.top_n.get()
        fdir = os.path.split(os.path.abspath(__file__))[0] + "/xml_scripts"

        for strategy in self.strategies:
            temp_name = "temp_PDBLIST.txt"
            OUTFILE = open(temp_name, 'w')
            decoys = score.get_ordered_decoy_list(strategy, top_n)
            for decoy in decoys:
                OUTFILE.write(decoy + "\n")
            OUTFILE.close()

            out_db_name = prefix + "_" + strategy
            out_db_batch = "Subset"

            # This should be done manually by getting struct_id and copying all the data in a new database via python sqlite3
            # I do not have time to figure that out and get it working right now, so this will have to do.

            analyze_strat.create_features_db(temp_name, fdir, features_type + "_features", self.rosetta_extension.get(),
                                             self.scorefxn, out_db_name, out_db_batch, self.db_dir.get(), False)

            os.remove(temp_name)


class Perc:
    """
    Simple class for holding enrichment/recovery information
    """

    def __init__(self, count, total):
        self.count = count
        self.total = total
        self._calc_per()

    def _calc_per(self):
        self.perc = self.count / float(self.total)

    def get_count(self):
        return self.count

    def get_total(self):
        return self.total

    def get_perc_decimal(self):
        return self.perc

    def get_perc_whole(self):
        return self.perc * 100

    def get_formated_perc(self, perc):
        return "%.3f" % perc


########################################################################################################################
### Helper Functions
########################################################################################################################
def get_str(value):
    if type(value) == str:
        return value
    else:
        return "%.3f" % value


def calculate_recovery(native_data, all_decoy_data, cdr, decoy_list=None):
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
            count += 1

    enrich_info = Perc(count, len(decoy_list))
    return enrich_info

def calculate_enrichments(all_decoy_data, cdr, decoy_list = None):
    """
    Returns defaultdict of [count_type] : Perc
    """

    if isinstance(all_decoy_data, CDRData): pass

    raw_decoy_map = all_decoy_data.get_concatonated_map()
    if not decoy_list:
        decoy_list = raw_decoy_map.keys()

    raw_counts = defaultdict(int)
    final_counts = defaultdict()

    for decoy in decoy_list:

        cdr_info = raw_decoy_map[decoy]
        if isinstance(cdr_info, CDRDataInfo): pass
        raw_counts[cdr_info.get_value_for_cdr(cdr)]+=1

    for count_type in raw_counts:
        final_counts[count_type] = Perc(raw_counts[count_type], len(decoy_list))

    return final_counts

def calculate_observed_value(value, all_decoy_data, cdr, decoy_list=None):
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
            count += 1

    enrich_info = Perc(count, len(decoy_list))
    return enrich_info


def count_native_matches(decoy_data, native_data, cdrs):
    if isinstance(decoy_data, CDRDataInfo): pass
    if isinstance(native_data, CDRDataInfo): pass

    count = 0
    for cdr in cdrs:
        if native_data.get_value_for_cdr(cdr) == decoy_data.get_value_for_cdr(cdr):
            count += 1
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
