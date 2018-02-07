import os
from Tkinter import *
import tkFileDialog
import tkSimpleDialog

from jade.RAbD.AnalyzeAntibodyDesigns import *


class AntibodyDesignAnalysisMenu(object):
    def __init__(self, main, compare_designs, main_gui):
        """

        :type main: Tk
        :type compare_designs: rabd.CompareAntibodyDesignStrategies
        :type main_gui:
        """

        self._tk_ = main
        self.compare_designs = compare_designs

        self.main_gui = main_gui

    def sho_tk(self):
        self._tk_.config(menu=self.main_menu)

        ###Tracers###
        self.compare_designs.is_camelid.trace_variable('w', self.camelid_tracer)

    def set_tk(self):

        self.main_menu = Menu(self._tk_)

        ## File Menu ##
        self.file_menu = Menu(self.main_menu, tearoff=0)

        self.file_menu.add_checkbutton(label="Camelid Antibody", variable=self.compare_designs.is_camelid)
        # self.dtypes_menu = Menu(self.main_menu, tearoff = 0)
        # self.dtypes_menu.add_checkbutton(label = "Group by dG", variable = self.compare_designs.group_dG)

        self.file_menu.add_command(label="Filter Models",
                                   command=lambda: self.main_gui.filter_settings_window.setup_sho_gui(Toplevel(self._tk_)))
        self.file_menu.add_command(label="Read Strategies from DB DIR",
                                   command=lambda: self.read_from_db_dir_set_strategies())
        self.file_menu.add_command(label="Add Strategy", command=lambda: self.main_gui.compare_frame.add_main_strategy())
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Set Reference Native", command=lambda: self.set_native_path())
        self.file_menu.add_command(label="Set Reference Database", command=lambda: self.set_reference_db())
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Set PyIgClassify Directory", command = lambda: self.set_pyigclassify_dir)
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

        self.main_pymol_menu = Menu(self.main_menu, tearoff=0)
        self.pymol_menu = Menu(self.main_menu, tearoff=0)
        self.pymol_menu.add_command(label="Top Models",
                                    command=lambda: self.compare_designs.copy_top())
        self.pymol_menu.add_command(label="All Models",
                                    command=lambda: self.compare_designs.copy_all_models())
        self.main_pymol_menu.add_cascade(label="Create PyMol Sessions", menu=self.pymol_menu)
        self.main_pymol_menu.add_separator()
        self.main_pymol_menu.add_checkbutton(label="Align Origin CDRs", variable=self.compare_designs.load_origin_pdbs)

        self.main_menu.add_cascade(label="PyMol", menu=self.main_pymol_menu)

        ## Clustal Menu ##
        self.clustal_menu = Menu(self.main_menu, tearoff=0)
        #self.clustal_menu.add_command(label="Set Max Processors", command=lambda: self.set_max_clustal_procs()) - Requires OpenMP support
        self.clustal_menu.add_command(label="Set Output format", command=lambda: self.set_clustal_output_format())
        self.clustal_menu.add_command(label="Set Soft Wrap", command=lambda: self.set_clustal_soft_wrap())
        self.clustal_menu.add_separator()
        self.clustal_menu.add_command(label="Run Clustal Omega on Top Decoys",
                                      command=lambda: self.run_clustal_omega())
        self.clustal_menu.add_command(label="Run Clustal Omega on ALL Combined Decoys",
                                      command=lambda: self.run_clustal_on_all_combined())
        #self.main_menu.add_cascade(label="Clustal", menu=self.clustal_menu)

        ## Alignment ##
        self.per_model_menu = Menu(self.main_menu, tearoff=0)

        '''
        self.per_model_menu.add_command(label="Output Length Alignments",
                                        command=lambda: self.compare_designs.output_len_or_clus_alignment('length', 'antibody'))
        self.per_model_menu.add_command(label="Output Cluster Alignments",
                                        command=lambda: self.compare_designs.output_len_or_clus_alignment('cluster', 'antibody'))
        self.per_model_menu.add_command(label="Output CDR Sequence Alignments",
                                     command=lambda: self.compare_designs.output_len_or_clus_alignment('aligned_sequence', 'antibody'))
        self.per_model_menu.add_separator()
        '''

        self.per_model_menu.add_command(label="Output to CSV (Top)", command = lambda: self.compare_designs.output_csv_data(top = True))
        self.per_model_menu.add_command(label="Output to CSV (ALL)", command = lambda: self.compare_designs.output_csv_data(top = False))
        self.per_model_menu.add_separator()
        self.per_model_menu.add_cascade(label="Clustal", menu = self.clustal_menu)

        self.main_menu.add_cascade(label="Model Data", menu=self.per_model_menu)


        ## Recovery ##
        self.per_strategy_menu = Menu(self.main_menu, tearoff = 0)
        self.recovery_menu = Menu(self.main_menu, tearoff=0)
        self.recovery_menu.add_command(label="Output Length Recovery",
                                       command=lambda: self.compare_designs.output_len_or_clus_recovery('length', 'antibody'))
        self.recovery_menu.add_command(label="Output Cluster Recovery",
                                       command=lambda: self.compare_designs.output_len_or_clus_recovery('cluster', 'antibody'))

        ## Enrichment ##
        self.enrichment_menu = Menu(self.main_menu, tearoff=0)
        self.enrichment_menu.add_command(label="Output Length Enrichments",
                                         command = lambda: self.compare_designs.output_len_or_clus_enrichment("length"))
        self.enrichment_menu.add_command(label="Output Cluster Enrichments",
                                         command = lambda: self.compare_designs.output_len_or_clus_enrichment("cluster"))

        self.per_strategy_menu.add_command(label="Output Score Summaries to CSV (Top)", command=lambda: self.compare_designs.output_csv_data(top=True, summary=True))
        self.per_strategy_menu.add_command(label="Output Score Summaries to CSV (All)", command=lambda: self.compare_designs.output_csv_data(top=False, summary=True))
        self.per_strategy_menu.add_cascade(label = "Recovery", menu = self.recovery_menu)
        self.per_strategy_menu.add_cascade(label = "Enrichment", menu = self.enrichment_menu)

        self.main_menu.add_cascade(label="Strategy Data", menu=self.per_strategy_menu)

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



        ######### Tracers ##############

    def camelid_tracer(self, name, index, mode):
        varValue = self.compare_designs.is_camelid.get()
        if varValue == 1:
            self.compare_designs.cdrs["L1"].set(0)
            self.compare_designs.cdrs["L2"].set(0)
            self.compare_designs.cdrs["L3"].set(0)

        return


        ######### Auxilliary Functions ###########
    def open_sequence_logo(self):
        pass

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
        self.main_gui.compare_frame.populate_all_strategies()

        strategies = self.get_full_strategy_list()
        self.compare_designs.set_strategies(strategies)

    def get_full_strategy_list(self):
        strategies = self.main_gui.compare_frame.current_strategies_listbox.get(0, END)
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
                                        initialvalue=self.main_gui.clustal_procs.get())
        if not max: return
        self.main_gui.clustal_procs.set(max)

    def set_clustal_output_format(self):
        f = tkSimpleDialog.askstring(title="Clustal output format", initialvalue=self.main_gui.clustal_output_format.get())
        if not f:
            return
        if not f in self.main_gui.clustal_output_formats:
            print "Format " + f + " not recognized.  Available formats are: \n" + repr(self.main_gui.clustal_output_formats)
            return

        self.main_gui.clustal_output_format.set(f)

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
            self.compare_designs.native_path = native_path

    def set_pyigclassify_dir(self):
        origin_path = tkFileDialog.askopenfilename(title = "PyIgClassify Directory", inititaldir = self.current_dir)
        if not origin_path:
            return
        else:
            self.current_dir = os.path.dirname(origin_path)
            self.compare_designs.pyigclassify_dir.set(origin_path)


            ######## Main Analysis ############

    def run_copy_all(self):
        self.compare_designs.copy_all_models()


    def run_clustal_on_all_combined(self):

        extra_options = tkSimpleDialog.askstring(title="Extra Options", prompt="Clustal Extra Options",
                                                 initialvalue=self.main_gui.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega_on_all_combined(self.main_gui.clustal_procs.get(),
                                                               self.main_gui.clustal_output_format.get(),
                                                               extra_options=extra_options)

    def run_clustal_omega(self):

        extra_options = tkSimpleDialog.askstring(title="Extra Options", prompt="Clustal Extra Options",
                                                 initialvalue=self.main_gui.base_clustal_options)
        if not extra_options:
            return

        self.compare_designs.run_clustal_omega(self.main_gui.clustal_procs.get(),
                                                               self.main_gui.clustal_output_format.get(),
                                                               extra_options=extra_options)

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
