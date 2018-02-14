#!/usr/bin/env python
# Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)

# Yes, the imports are all over.  Basically everything I've coded over the past few years.
import os, random

from argparse import ArgumentParser


# TkInter
from Tkinter import *
import tkMessageBox

# PyIgD
from jade.antibody.cdr_data.CDRDataTypes import *
from jade.basic.threading.Threader import *
from jade.RAbD_BM.AnalysisInfo import *
from jade.RAbD.AnalyzeAntibodyDesigns import CompareAntibodyDesignStrategies

## Window Frames
from jade.RAbD.window_main.menu import AntibodyDesignAnalysisMenu
from jade.RAbD.window_main.CompareStrategiesFrame import CompareStrategiesFrame
from jade.RAbD.window_main.FeaturesFrame import FeaturesFrame
from jade.RAbD.window_main.AnalysisFrame import AnalysisFrame

## Window Modules
from jade.RAbD.window_modules.FilterSettingsWindow import *

from jade.basic.TKinter.ImageFrame import ImageFrame
from jade.basic import path


#DESCRIPTION
# Application which can analyze features databases and create PyMol sessions of top models sorted by physical attributes
#  such as Interface Energy, SASA, etc.

# The application can rename decoys into their sorted info and output all sorted data in tables with full cluster, sequence, and physical data.
#  Very useful for antibody design analysis.

#REQUIREMENTS:
#clustal_omega must be installed and the application should be in your $PATH environment.
#PyMol must be installed and able to be called from the command line.
# ETC: Requires Biopython, numpy, etc.  Run the Jade setup script to install these.

#EXAMPLE COMMAND:
#
#    RAbD_Jade.py --native input_pdbs/pareto_2j88_renum_0002.pdb --cdrs L1
#      --pyigclassify_dir /Users/jadolfbr/Documents/projects/PyIgClassify --analysis_name testing

def get_parser():
    parser = ArgumentParser(description="GUI application to analyze designs output by RosettaAntibodyDesign.  "
                            "Designs should first be analyzed by both the AntibodyFeatures and CDRClusterFeatures reporters "
                            "into sqlite3 databases.")

    parser.add_argument("--db_dir",
                        help="Directory with databases to compare. DEFAULT = databases",
                        default="databases")

    parser.add_argument("--analysis_name",
                        help="Main directory to complete analysis. DEFAULT = prelim_analysis",
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
                        help = "Optional PyIgClassify Root Directory with DBOUT. Used for debugging.",
                        default = "")

    parser.add_argument("--jsons","-j",
                        help = "Analysis JSONs to use.  See RAbD_MB.AnalysisInfo for more on what is in the JSON."
                               "The JSON allows us to specify the final name, decoy directory, and features db associated with the benchmark as well as all options that went into it.",
                        nargs = "*")

    return parser

def main():

    parser = get_parser()
    options = parser.parse_args()

    if options.root_dir != os.getcwd():
        print "Changing to root."
        os.chdir(options.root_dir)

    jsons = []
    if options.jsons:
        jsons = [AnalysisInfo(json_file) for json_file in options.jsons]

    GUI = CompareAntibodyDesignStrategies_GUI(Tk(),
                                              db_dir=options.db_dir,
                                              analysis_dir=options.analysis_name,
                                              jsons=jsons)

    # Set any values
    if options.native:
        GUI.compare_designs.native_path = options.native

    GUI.compare_designs.set_cdrs_from_list(options.cdrs)
    GUI.compare_designs.pyigclassify_dir.set(options.pyigclassify_dir)

    GUI.run()



class CompareAntibodyDesignStrategies_GUI:
    def __init__(self, main, db_dir="", analysis_dir="", jsons = [], strategies=[]):
        """
        JSON file is used for analysis or renaming things:
        Requires these settings.

        { "short_name": "min.remove_antigen-T",
          "decoy_dir":"decoy_dir",
          "features_db":"ab_features_path"

        }

        :param main:
        :param db_dir:
        :param analysis_dir:
        :param jsons: [AnalysisInfo]
        :param strategies:
        """
        self._tk_ = main
        self._tk_.title("PyIgDesign Compare")
        self.compare_designs = CompareAntibodyDesignStrategies(db_dir, analysis_dir, strategies, jsons=jsons)
        if self.compare_designs.db_dir and not self.compare_designs.strategies and not jsons:
            self.compare_designs.set_strategies_from_databases()
        elif self.compare_designs.db_dir and not self.compare_designs.strategies and jsons:
            self.compare_designs.set_strategies_from_json_infos()

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

        #self.image_10e9 = ImageFrame(self._tk_, path.get_database_path()+"/assets/10e9_image_55.gif", bd=2, relief=RIDGE)#

        images = [path.get_database_path()+ "/assets/"+r for r in ["smiling_2J88_small.gif", "10e9_image_55.gif", "2j88_native_full_sticks_lines_small.gif", "ch103_native_full_small.gif"]]

        image_to_use  = random.choice(images)
        self.image_10e9 = ImageFrame(self._tk_, image_to_use, bd=2, relief=RIDGE)

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
