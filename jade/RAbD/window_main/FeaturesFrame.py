import os
from Tkinter import *
from jade.RAbD.AnalyzeAntibodyDesigns import *


class FeaturesFrame(Frame):
    def __init__(self, main, compare_designs, main_gui, **options):
        """

        :type main: Tk
        :type compare_designs: CompareAntibodyDesignStrategies
        :type main_gui:
        """

        Frame.__init__(self, main, **options)
        self._tk_ = main
        self.compare_designs = compare_designs
        self.main_gui = main_gui

        self.set_tk()
        self.sho_tk()

    def set_tk(self):

        self.label_main = Label(self, text="Features Plots", justify=CENTER)

        self.ab_features_button = Button(self, text="Run Antibody Features",
                                         command=lambda: self.run_features_reporter("antibody"), justify=CENTER)
        self.clus_features_button = Button(self, text="Run Cluster Features",
                                           command=lambda: self.run_features_reporter("cluster"), justify=CENTER)

        self.ab_features_options_label = Label(self, text="Antibody Features Options", justify=CENTER)

        self.normal_hbond_radio = Radiobutton(self, text="All Hbond R Scripts",
                                              variable=self.compare_designs.features_hbond_set, value=0)
        self.min_hbond_radio = Radiobutton(self, text="Minimal Hbond R Scripts",
                                           variable=self.compare_designs.features_hbond_set, value=1)
        self.no_hbond_radio = Radiobutton(self, text="No Hbond R Scripts",

                                          variable=self.compare_designs.features_hbond_set, value=2)
    def sho_tk(self):
        r=0
        c=0
        self.label_main.grid(row=r+5, column=c+0, columnspan=6, sticky=W+E)
        self.ab_features_button.grid(row=r + 6, column=c + 0, columnspan=3, pady=3, sticky=W + E)
        self.clus_features_button.grid(row=r + 6, column=c + 3, columnspan=3, pady=3, sticky=W + E)

        # self.ab_features_options_label.grid(row = r+7, column = c+0, columnspan = 2, pady = 3, sticky = W+E)

        self.normal_hbond_radio.grid(row=r + 8, column=c + 0, columnspan=3, pady=1, sticky=W)
        self.min_hbond_radio.grid(row=r + 9, column=c + 0, columnspan=3, pady=1, sticky=W)
        self.no_hbond_radio.grid(row=r + 10, column=c + 0, columnspan=3, pady=1, sticky=W)

    def run_features_reporter(self, type):
        strategies = self.main_gui.get_full_strategy_list()
        if len(strategies) == 0:
            print "No strategies selected..."
            return

        self.compare_designs.set_strategies(strategies)
        self.compare_designs.run_features(type)