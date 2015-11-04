import os
from Tkinter import *
import tkFileDialog
import tkSimpleDialog

from RAbD.AnalyzeAntibodyDesigns import *


class AnalysisFrame(Frame):
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

        self.cluster = StringVar()
        self.decoy = StringVar()


        self.clusterOPT = StringVar()
        self.clusterOPT.set("Open Sequence Logo")

        self.decoyOPT = StringVar()
        self.decoyOPT.set("Print Score Data")

        self.clusterFUNCTIONS = {
            "Open Sequence Logo":lambda:self.open_seq_logo(),
            "Open MSA":lambda:self.open_msa()
        }

        self.decoyFUNCTIONS = {
            "Print Score Data":lambda:self.print_decoy_info()
        }


        #############
        self.set_tk()
        self.sho_tk()

    def set_tk(self):
        self.main_label = Label(self, text="Analysis", justify=CENTER)
        self.options_cluster = OptionMenu(self, self.clusterOPT, *(sorted(self.clusterFUNCTIONS.keys())))
        self.entry_cluster = Entry(self, textvariable=self.cluster, justify = CENTER)
        self.button_cluster = Button(self, text="Run Cluster Analysis", command=lambda: self.clusterFUNCTIONS[self.clusterOPT.get()]())

        self.options_decoy = OptionMenu(self, self.decoyOPT, *(sorted(self.decoyFUNCTIONS.keys())))
        self.entry_decoy = Entry(self, textvariable=self.decoy, justify=CENTER)
        self.button_decoy = Button(self, text="Run Decoy Analysis", command = lambda: self.decoyFUNCTIONS[self.decoyOPT.get()]())

        #self.entry_sequence = Entry(self, textvariable=self.sequence, justify=CENTER)
        #self.button_sequence = Button(self, text = "NA")

    def sho_tk(self):
        self.main_label.grid(row=0, column=0, columnspan=3, sticky=W+E, padx=3, pady=3)

        self.entry_cluster.grid(row=1, column=0, sticky=W+E, padx=2, pady=2)
        self.options_cluster.grid(row=1, column=1, sticky=W+E, padx=2, pady=2)
        self.button_cluster.grid(row=1, column=2, sticky=W+E, padx=2, pady=2)


        self.entry_decoy.grid(row=2, column=0, sticky=W+E, padx=2, pady=5)
        self.options_decoy.grid(row=2, column=1, sticky=W+E, padx=2, pady=2)
        self.button_decoy.grid(row=2, column=2, sticky=W+E, padx=2, pady=2)

    def open_seq_logo(self):

        self.check_set_pyigclassify()
        seqlogo_dir = self.compare_designs.pyigclassify_dir.get()+"/"+self.compare_designs.weblogo_rel_path
        p = seqlogo_dir+"/"+self.cluster.get()+"_weblogo.png"
        os.system("open "+p+" &")

    def open_msa(self):
        self.check_set_pyigclassify()
        seqlogo_dir = self.compare_designs.pyigclassify_dir.get()+"/"+self.compare_designs.weblogo_rel_path
        p = seqlogo_dir+"/"+self.cluster.get()+"_msa.txt"
        os.system("open "+p+" &")

    def print_decoy_info(self):
        pass

    def check_set_pyigclassify(self):
        if not os.path.exists(self.compare_designs.pyigclassify_dir.get()):
            print "PyIgClassify Dir does not exist"
            pyigclassify_dir = tkFileDialog.askdirectory(initialdir=self.main_gui.current_dir, title="PyIgClassify DIR")
            if not pyigclassify_dir:return
            else:
                self.main_gui.current_dir = os.path.dirname(pyigclassify_dir)
                self.compare_designs.pyigclassify_dir.set(pyigclassify_dir)
