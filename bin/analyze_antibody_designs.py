#!/usr/bin/env python

# Yes, the imports are all over.  Basically everything I've coded over the past few years.
import sys
import os
import sqlite3
import pandas

from argparse import ArgumentParser


# TkInter
from Tkinter import *
import tkMessageBox

p = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(p);  # Allows all modules to use all other modules, without needing to update pythonpath

# PyIgD
from antibody.cdr_data.CDRDataTypes import *
from tools.Threader import *
from RAbD.AnalyzeAntibodyDesigns import CompareAntibodyDesignStrategies

## Window Frames
from RAbD.window_main.menu import AntibodyDesignAnalysisMenu
from RAbD.window_main.CompareStrategiesFrame import CompareStrategiesFrame
from RAbD.window_main.FeaturesFrame import FeaturesFrame
from RAbD.window_main.AnalysisFrame import AnalysisFrame

## Window Modules
from RAbD.window_modules.FilterSettingsWindow import *

from tools.ImageFrame import ImageFrame
from tools import path

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

    parser.add_argument("--pyigclassify_dir",
                        help = "DBOUT directory of PyIgClassify",
                        default = "")

    options = parser.parse_args()

    if options.root_dir != os.getcwd():
        print "Changing to root."
        os.chdir(options.root_dir)

    GUI = CompareAntibodyDesignStrategies_GUI(Tk(), options.db_dir, options.analysis_name)

    # Set any values
    if options.native:
        GUI.compare_designs.native_path = options.native

    GUI.compare_designs.set_cdrs_from_list(options.cdrs)
    GUI.compare_designs.pyigclassify_dir.set(options.pyigclassify_dir)

    GUI.run()



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
        self.clustal_output_formats = ['fasta', 'clustal', 'msf', 'phylip', 'selex', 'stockholm', 'vienna',
                                       'a2m', 'fa', 'clu', 'phy', 'st', 'vie']
        self.base_clustal_options = " -v --auto --force"

        ###Sub Windows ####
        self.filter_settings_window = FilterSettingsWindow(self.compare_designs.filter_settings)

        #atexit.register(self.exit)
        self._tk_.protocol("WM_DELETE_WINDOW", self.exit)


    def set_tk(self):


        self.compare_frame = CompareStrategiesFrame(self._tk_, self.compare_designs, self, bd=2, relief=GROOVE)
        self.analysis_frame = AnalysisFrame(self._tk_, self.compare_designs, self, bd=2, relief= GROOVE)
        self.features_frame = FeaturesFrame(self._tk_, self.compare_designs, self, bd=2, relief=GROOVE)

        self.menu_class = AntibodyDesignAnalysisMenu(self._tk_, self.compare_designs, self)

        self.image_10e9 = ImageFrame(self._tk_, path.get_database_path()+"/assets/10e9_image_55.gif", bd=2, relief=RIDGE)#

    def sho_tk(self):

        # self.root_dir_label.grid(row = r+0, column = c+0, columnspan = 2, sticky = W+E, pady = 7)
        # self.db_dir_entry.grid( row = r+1, column = c+0, columnspan = 2, sticky = W+E, padx = 5)

        self.menu_class.set_tk()
        self.menu_class.sho_tk()

        self.compare_frame.grid(row=1, column=0, rowspan=2, padx=15, pady=15 )


        self.features_frame.grid(row=3, column=0, rowspan=1, padx=15, pady=15)

        self.analysis_frame.grid(row=1, column=1, rowspan=1, padx=10, pady=10)

        self.image_10e9.grid(row=2, column=1, rowspan=2, padx=10, pady=10, sticky=N+S)

        self.compare_frame.populate_all_strategies()

    def run(self):
        self.set_tk()
        self.sho_tk()
        self._tk_.mainloop()

    def exit(self):
        if threads.n_alive() > 1:
            if tkMessageBox.askyesno(title = "Exit?", message = repr(threads.n_alive())+" Processes still running.  Exit?"):
                threads.kill_all()
                self._tk_.destroy()
            else:
                self._tk_.mainloop()
        else:
            self._tk_.destroy()

    def location(self):
        """
        Allows the script to be self-aware of it's path.
        So that it can be imported/ran from anywhere.
        """

        p = os.path.abspath(__file__)
        pathSP = os.path.split(p)
        return pathSP

    def get_full_strategy_list(self):
        strategies = self.compare_frame.current_strategies_listbox.get(0, END)
        return strategies

if __name__ == "__main__":
    # 1) Main StrategyAnalysis DIR
    # 2) Output DIR
    # 3) Any Names of strategies to run.  Else will attempt load from Main StrategyAnalysisDIR

    main()
