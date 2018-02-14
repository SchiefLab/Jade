#!/usr/bin/env python

import glob
import sqlite3
import tkFileDialog
import tkSimpleDialog
from Tkinter import *
from tkFont import *
from argparse import ArgumentParser

import jade.rosetta_jade.FeaturesJsonCreator as json_creator

import jade.RAbD_BM.tools as tools
from jade.rosetta_jade.BenchmarkInfo import *



#Runs all analysis of a particular benchmarking experiment.

def get_parser():
    parser = ArgumentParser(description="This program is a GUI used for benchmarking Rosetta Antibody Design."
                            "Before running this application, you will probably want to run 'run_rabd_features_for_benchmarks.py to create the databases required.")

    parser.add_argument("--main_dir",
                        help = "Main working directory. Not Required.  Default = PWD",
                        default = os.getcwd())

    parser.add_argument("--out_dir",
                        help = "Output data directory.  Not Required.  Default = pooled_data",
                        default = "pooled_data")

    parser.add_argument("--jsons","-j",
                        help = "Analysis JSONs to use.  See RAbD_MB.AnalysisInfo for more on what is in the JSON."
                               "The JSON allows us to specify the final name, decoy directory, and features db associated with the benchmark as well as all options that went into it.",
                        nargs = "*",
                        required = True)

    return parser

def main():

    parser = get_parser()
    options = parser.parse_args()

    GUI = CompareBenchmarks_GUI(Tk(), options.main_dir, options.out_dir, options.benchmarks)
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

class AnalyzeBenchmarksGUI:
    """
    Deprecated
    """
    def __init__(self, analyzer, nstruct_check, add_rec_data, override_check, get_ensemble_data, out_name):
        self.base_dir = os.path.split(os.path.abspath(__file__))[0]
        self.current_dir = self.base_dir

        self._tk_ = Tk()

        self.benchmark_list = StringVar()
        self.benchmark_list.set("benchmark_list.txt")
        self.benchmarks = tools.parse_benchmark_list(self.benchmark_list.get())

        self.analyser = analyzer
        if isinstance(self.analyser, AnalyzeBenchmarks):pass

        self.native_db = StringVar(value = analyzer.native_db);
        self.lambda_db = StringVar(value = analyzer.lambda_db);
        self.kappa_db = StringVar(value = analyzer.kappa_db);

        self.db_out_name = StringVar(value = out_name)
        self.add_rec_data = IntVar(value = add_rec_data)
        self.override_check = IntVar(value = override_check)
        self.get_ensemble_data = IntVar(value = get_ensemble_data)

        self.nstruct_check = StringVar(); self.nstruct_check.set("5")

    def setTk(self):

        self.entry_db = Entry(self._tk_, textvariable = self.db_out_name)
        #self.entry_native_db = Entry(self._tk_, textvariable = self.native_db)
        #self.entry_kappa_db = Entry(self._tk_, textvariable = self.kappa_db)
        #self.entry_lambda_db = Entry(self._tk_, textvariable = self.lambda_db)
        self.label_db = Label(self._tk_, text = "Result Database")
        #self.label_native = Label(self._tk_, text = "Native DB")
        #self.label_lambda = Label(self._tk_, text = "Lambda DB")
        #self.label_kappa = Label(self._tk_, text = "Kappa DB")

        self.listbox_all = Listbox(self._tk_)
        self.listbox_current = Listbox(self._tk_)
        self.button_add_data = Button(self._tk_, text = "Add Data to DB", command = lambda: self.analyze_data())

        self.label_all = Label(self._tk_, text = "Benchmark List", justify = CENTER)
        self.label_current = Label(self._tk_, text = "Benchmarks To Add", justify = CENTER)

        self.checkbox_add_rec = Checkbutton(self._tk_, variable = self.add_rec_data, text="Add Recovery data to Features DB?")
        self.checkbox_override_check = Checkbutton(self._tk_, variable = self.override_check, text = "Override nstruct check?")
        self.checkbox_ensemble = Checkbutton(self._tk_, variable = self.get_ensemble_data, text = "Calculate data for Ensemble?")

        self.entry_nstruct = Entry(self._tk_, textvariable=self.nstruct_check)
        self.label_nstruct = Label(self._tk_, text = "nstruct check")
        self.label_all2 = Label(self._tk_, text = "DB Benchmarks", justify = CENTER)
        self.label_current2 = Label(self._tk_, text = "Benchmarks To Plot", justify = CENTER)

        self.listbox_all2 = Listbox(self._tk_)
        self.listbox_current2 = Listbox(self._tk_)

        self.listbox_groups = Listbox(self._tk_)
        self.listbox_groups_current = Listbox(self._tk_)

        self.listbox_names = Listbox(self._tk_)
        self.listbox_names_current = Listbox(self._tk_)

        #self.frame = Frame(self._tk_, relief=SUNKEN, height=5,  width = 300, bg="black")
        #self.separator = Separator(self._tk_, orient = HORIZONTAL)

        self.label_groups = Label(self._tk_, text = "Groups")
        self.label_individuals = Label(self._tk_, text = "MinType")

    def set_Menu(self):

        ### File  ####
        self.main_menu=Menu(self._tk_)
        self.file_menu=Menu(self.main_menu, tearoff=0)

        self.file_menu.add_command(label="Load Benchmark List", command=lambda: self.select_bm_list_reload())
        self.file_menu.add_command(label="Set Output plot DIR")
        self.file_menu.add_command(label="Reload DB")
        self.main_menu.add_cascade(label="File", menu=self.file_menu)

        ### DB ###
        self.db_menu=Menu(self.main_menu, tearoff=0)
        self.db_menu.add_command(label="Native DB", command = lambda: self.set_db("Native DB", self.native_db))
        self.db_menu.add_command(label="Lambda DB", command = lambda: self.set_db("Lambda DB", self.lambda_db))
        self.db_menu.add_command(label="Kappa  DB", command = lambda: self.set_db("Kappa  DB", self.kappa_db))

        self.main_menu.add_cascade(label="Input Databases", menu=self.db_menu)
        self._tk_.config(menu=self.main_menu)

    def shoTk(self):
        r = 0; c = 0;
        self.entry_db.grid(row=r, column = c, sticky = W+E+S+N, padx=6, pady=6)
        self.label_db.grid(row=r, column=c+1)

        self.label_all.grid(row=r+2, column = r)
        self.label_current.grid(row=r+2, column=r+1)

        #self.entry_native_db.grid(row = r+3, column = r)
        #self.label_native.grid(row = r+4, column = r)

        #self.entry_lambda_db.grid(row = r+5, column = r)
        #self.label_lambda.grid(row = r+6, column = r)

        #self.entry_kappa_db.grid(row = r+7, column = r)
        #self.label_kappa.grid(row = r+8, column = r)


        self.listbox_all.grid(row=r+3, column=0, rowspan=6)
        self.listbox_current.grid(row=r+3, column = c+1, rowspan = 6)
        self.button_add_data.grid(row=r+5, column = c+2, padx=6)

        self.checkbox_add_rec.grid(row=r+3, column = c+3, sticky=W)
        self.checkbox_ensemble.grid(row = r+4, column = c+3, sticky = W)
        self.checkbox_override_check.grid(row = r+5, column = c+3, sticky = W)
        self.entry_nstruct.grid(row = r+6, column = c+3)
        self.label_nstruct.grid(row = r+7, column = c+3)

        #self.frame.grid(row = r+9, column = c, columnspan = 3, sticky = W)

        #self.separator.grid(row = r+9, column = c, columnspan = 4, sticky = W+E, pady = 15)
        self.label_all2.grid(row = r+10, column = c)
        self.label_current2.grid(row = r+10, column = c+1)

        self.listbox_all2.grid(row=r+12, column = c, rowspan = 6)
        self.listbox_current2.grid(row = r+12, column = c+1, rowspan = 6)

        self.label_groups.grid(row = r+18, column = c, columnspan = 2)
        self.listbox_groups.grid(row = r+19, column = c, rowspan = 5)
        self.listbox_groups_current.grid(row = r+19, column = c+1, rowspan = 5)

        self.label_individuals.grid(row = r+25, column = c, columnspan = 2)
        self.listbox_names.grid(row = r+26, column = c, rowspan = 5)
        self.listbox_names_current.grid(row = r+26, column = c+1, rowspan = 5)


        self.listbox_all.bind("<Double-Button-1>", lambda event: self.copy_item(self.listbox_all, self.listbox_current))
        self.listbox_current.bind("<Double-Button-1>", lambda event: self.delete_item(self.listbox_current))

        #self.listbox_current.bind("<ButtonRelease-1>", lambda event: self.show_benchmark_info(self.listbox_current))
        #self.listbox_all.bind("<ButtonRelease-1>", lambda event: self.show_benchmark_info(self.listbox_all))

        self.populate_all_listbox()

    def run(self):
        self._tk_.title("Analyze Benchmarks")
        self.setTk()
        self.set_Menu()
        self.shoTk()
        self._tk_.mainloop()

    ######## Functions ############

    def analyze_data(self):
        analyzer.out_db_name = self.db_out_name.get()
        analyzer.kappa_db = self.kappa_db.get()
        analyzer.lambda_db = self.lambda_db.get()
        analyzer.native_db_path = self.native_db.get()


        names =  self.listbox_current.get(0, END)
        new_benchmarks = []
        for name in names:
            print name
            for benchmark in self.benchmarks:
                if benchmark.final_name == name:
                    new_benchmarks.append(benchmark)

        print new_benchmarks
        analyzer.analyze_benchmarks(new_benchmarks, self.nstruct_check.get(), self.add_rec_data.get(), self.get_ensemble_data.get(), self.override_check.get())

    def select_bm_list_reload(self):
        infilename = tkFileDialog.askopenfilename(initialdir=self.current_dir, title="Select Benchmark List file")
        self.benchmark_list.set(infilename)
        self.current_dir = os.path.dirname(infilename)
        self._update_bm_list_and_analyzer()


    def _update_bm_list_and_analyzer(self):
        self.benchmarks = tools.parse_benchmark_list(self.benchmark_list.get())
        self.analyzer = AnalyzeBenchmarks(self.benchmark_list.get(), self.native_db.get(), self.lambda_db.get(), self.kappa_db.get())
        self.populate_all_listbox()

    def set_db(self, title, variable):
        tkSimpleDialog.askstring(title=title, variable=variable, prompt = "Set DB")

    ########## Events #############
    def show_benchmark_info(self, listbox):
        for benchmark in self.benchmarks:
            if benchmark.final_name == listbox.get(listbox.curselection()):
                print repr(benchmark)

    def copy_item(self, from_listbox, to_listbox):
        item = from_listbox.get(from_listbox.curselection())
        to_listbox.insert(END, item)

    def delete_item(self, listbox):
        listbox.delete(listbox.curselection())

    def populate_all_listbox(self):
        for benchmark in self.benchmarks:
            self.listbox_all.insert(END, benchmark.final_name)



class CompareBenchmarks_GUI:
    def __init__(self, main, main_dir, out_dir_name, benchmarks):
        self._tk_ = main
        self._tk_.title("Analyze Antibody Design Benchmarks")
        self.compare_benchmarks = CompareBenchmarks(main_dir, out_dir_name, benchmarks)
        if self.compare_benchmarks.main_dir and not self.compare_benchmarks.benchmarks:
            self.compare_benchmarks.set_benchmarks_from_databases()

        self.current_dir = os.path.split(os.path.abspath(__file__))[0]

        ###Sub Windows ####


    def run(self):
        self.set_tk()
        self.set_menu()
        self.sho_tk()
        self._tk_.mainloop()

    def set_tk(self):
        self.main_dir_entry = Entry(self._tk_, textvariable = self.compare_benchmarks.main_dir, justify = CENTER)
        self.out_dir_entry = Entry(self._tk_, textvariable = self.compare_benchmarks.out_dir_name, justify = CENTER)

        self.main_dir_label = Label(self._tk_, text = "Main Analysis Directory", justify = CENTER)
        self.out_dir_label = Label(self._tk_, text = "Combined Output Directory name", justify = CENTER)

        self.all_benchmarks_listbox = Listbox(self._tk_)
        self.current_benchmarks_listbox = Listbox(self._tk_)

        self.ab_features_button = Button(self._tk_, text = "Run Antibody Features", command = lambda: self.compare_benchmarks.run_features("antibody"), justify = CENTER)
        self.clus_features_button = Button(self._tk_, text = "Run Cluster Features", command = lambda: self.compare_benchmarks.run_features("cluster"), justify = CENTER)


    def sho_tk(self, r = 0, c = 0):

        self.main_dir_label.grid(row = r+0, column = c+0, columnspan = 2, sticky = W+E, pady = 7)
        self.main_dir_entry.grid( row = r+1, column = c+0, columnspan = 2, sticky = W+E, padx = 5)

        self.all_benchmarks_listbox.grid(row = r+2, column = c+0, padx = 6, pady = 10)
        self.current_benchmarks_listbox.grid(row = r+2, column = c+1, padx = 6, pady = 10)



        #self.separator.grid(row = r+3, column = c, columnspan = 2, sticky = W+E, pady = 15)

        self.out_dir_label.grid(row = r+5, column = c+0, columnspan = 1, pady = 5)
        self.out_dir_entry.grid(row = r+5, column = c+1, columnspan = 1, padx = 5, pady = 5)

        self.all_benchmarks_listbox.bind("<Double-Button-1>", lambda event: self.add_to_current(self.all_benchmarks_listbox, self.current_benchmarks_listbox))
        #self.all_benchmarks_listbox.bind("<Button-3>", lambda event: self.show_strat_items())

        self.ab_features_button.grid(row = r+6, column = c+0, columnspan = 1, pady = 3, sticky = W+E)
        self.clus_features_button.grid(row = r+6, column = c+1, columnspan = 1, pady = 3, sticky = W+E)


        self.current_benchmarks_listbox.bind("<Double-Button-1>", lambda event: self.delete_current(self.current_benchmarks_listbox))


        self.populate_all_benchmarks()

    def set_menu(self):
        self.main_menu = Menu(self._tk_)

        ## File Menu ##
        self.file_menu = Menu(self.main_menu, tearoff=0)
        self.file_menu.add_command(label = "Read Benchmarks from Main", command = lambda: self.read_from_main_set_benchmarks())
        self.file_menu.add_command(label = "Add Benchmark", command = lambda: self.add_main_benchmark())
        self.file_menu.add_separator()
        self.file_menu.add_command(label = "Set Reference Database", command = lambda: self.set_reference_db())
        #self.file_menu.add_checkbutton(label = "Use Recovery database where possible", variable = self.compare_benchmarks.use_recovery_dbs)
        #self.file_menu.add_separator()
        #self.file_menu.add_command(label = "Set Scorefunction", command = self.set_scorefunction())
        self.main_menu.add_cascade(label = "File", menu = self.file_menu)

        ## Plot Menu ##
        self.plot_menu = Menu(self.main_menu, tearoff = 0)
        #self.plot_menu.add_command(label = "Plot Individual RR + Recoveries", command = lambda: self.plot_individuals())
        self.plot_menu.add_command(label = "Plot Grouped RR + Recoveries", command = lambda: self.plot_groups())

        self.main_menu.add_cascade(label = "Plot", menu = self.plot_menu)

        self.plot_cdr_data_menu = Menu(self.main_menu, tearoff = 0)
        self.plot_cdr_data_menu.add_command(label = "Set Reference Database", command = lambda : self.set_reference_db())
        self.plot_cdr_data_menu.add_command(label = "ag_ab_dSASA", command = lambda: self.compare_benchmarks.plot_cdr_data("dSASA"))
        self.plot_cdr_data_menu.add_command(label = "ab_ab_contacts_total", command = lambda: self.compare_benchmarks.plot_cdr_data("contacts"))
        self.plot_cdr_data_menu.add_command(label = "custom", command = lambda: self.setup_custom_cdr_data_plot())
        self.plot_menu.add_cascade(label = "Plot by CDR Data", menu = self.plot_cdr_data_menu)

        ## Util Menu ##
        self.util_menu = Menu(self.main_menu, tearoff = 0)
        self.util_menu.add_command(label = "Spit Recovery Database", command = lambda:self.split_recovery_db())
        self.main_menu.add_cascade(label = "Util", menu = self.util_menu)

        self._tk_.config(menu = self.main_menu)


    ########### Callbacks ############

    def show_strat_items(self):
        pass

    def add_to_current(self, from_listbox, to_listbox):
        item = from_listbox.get(from_listbox.curselection())
        to_listbox.insert(END, item)
        benchmarks = self.get_full_benchmark_list()
        self.compare_benchmarks.set_benchmarks(benchmarks)

    def delete_current(self, listbox):
        listbox.delete(listbox.curselection())
        benchmarks = self.get_full_benchmark_list()
        self.compare_benchmarks.set_benchmarks(benchmarks)

    def populate_all_benchmarks(self):
        self.all_benchmarks_listbox.delete(0, END)
        for benchmark in self.compare_benchmarks.benchmarks:
            self.all_benchmarks_listbox.insert(END, benchmark)

        self.all_benchmarks_listbox.autowidth(100)
        self.current_benchmarks_listbox.autowidth(100, self.compare_benchmarks.benchmarks)


    ######### Auxilliary Functions ###########

    def add_main_benchmark(self):
        benchmark_name = tkSimpleDialog.askstring(title = "benchmark", prompt = "benchmark Name")
        if not benchmark_name:
            return

        self.all_benchmarks_listbox.insert(END, benchmark_name)

    def read_from_main_set_benchmarks(self):
        if not self.compare_benchmarks.main_dir.get():

            strat_dir = tkFileDialog.askdirectory(initialdir = self.current_dir, title = "benchmark Analysis Directory")
            if not strat_dir:
                return
            self.current_dir = strat_dir
            self.main_dir.set(strat_dir)
        self.compare_benchmarks.set_benchmarks_from_databases()
        self.populate_all_benchmarks()

    def get_full_benchmark_list(self):
        benchmarks = self.current_benchmarks_listbox.get(0, END)
        return benchmarks

    def set_reference_db(self):
        d = tkFileDialog.askopenfilename(title = "Reference DB", initialdir = self.current_dir)
        if not d: return
        self.compare_benchmarks.reference_db.set(d)

    def set_scorefunction(self):
        score = tkSimpleDialog.askstring(title="Score", prompt = "Set Scorefunction", initialvalue = self.compare_benchmarks.scorefunction.get())
        if not score:
            return
        self.compare_benchmarks.scorefunction.set(score)

    def split_recovery_db(self):
        db = tkFileDialog.askopenfilename(initialdir = self.current_dir, title = "Combined recovery database")
        if not db: return

        self.current_dir = os.path.dirname(db)
        rec = RecoveryData()
        rec.split_into_databases(db)

    ######### Main Functions ################
    def plot_groups(self):
        self.compare_benchmarks.plot_groups()

    def plot_individuals(self):
        self.compare_benchmarks.plot_individuals()

    def setup_custom_cdr_data_plot(self):
        type = tkSimpleDialog.askstring(title = "Data Type", prompt = "Please enter a custom CDR data to plot against")
        if not type: return
        self.compare_benchmarks.plot_cdr_data(type)

class RecoveryData:
    def __init__(self):
        self.tables = ["all_cluster", "all_length", "cdr_cluster", "cdr_length", "exp_cluster", "exp_length"]

    def split_into_databases(self, db):
        con = sqlite3.connect(db)


class CompareBenchmarks:
    def __init__(self, main_analysis_dir, out_dir_name, benchmarks = defaultdict(), scorefunction = "talaris2013"):

        self.main_dir = StringVar(value = main_analysis_dir)
        self.out_dir_name = StringVar(value = out_dir_name)
        self.scorefunction = StringVar(value = scorefunction)
        self.reference_db = StringVar()
        self.use_recovery_dbs = IntVar(value = 1)
        self.benchmarks = benchmarks

    def _setup_outdir(self, subdirs = [], use_out_dir_name = True):
        """
        Sets up the main output dir in the main_analys_dir, and any subdirectories such as 'decoys' or decoys/combined_3
        Returns the final output directory
        """

        #filters = self._setup_filters()


        if self.out_dir_name.get() and use_out_dir_name:
            outdir = self.main_dir.get()+"/"+self.out_dir_name.get()
        elif not self.main_dir.get():
            sys.exit("Main BenchmarkAnalysis DIR must exist.  Please use analyze_antibody_design_benchmark to populate")
        else:
            outdir = self.main_dir.get()
        if not os.path.exists(outdir): os.mkdir(outdir)

        for subdir in subdirs:
            if not subdir: continue
            outdir = outdir+"/"+subdir
            if not os.path.exists(outdir): os.mkdir(outdir)

        return outdir

    def _setup_cdr_types(self, data_type, is_camelid, features_type, native_path = None):
        if data_type == "length":
            data_class = CDRLengthData(is_camelid)
        else:
            data_class = CDRClusterData(is_camelid)
        for benchmark in self.benchmarks:
            db_path = self.get_db_path(benchmark, features_type)
            if not os.path.exists(db_path):
                sys.exit("DB path does not exist.  Please Run Features reporter with StructureScores and ScoreTypes for this benchmark\n"+db_path)
            print db_path
            con = sqlite3.connect(db_path)
            data_class.add_data(benchmark, con)

        if native_path:
            data_class.set_native_data_from_rosetta(native_path)

        return data_class

    def _setup_db_paths_recurse(self, order = ["cluster", "recoveries", "antibody"]):
        db_paths = []
        for benchmark in self.get_benchmarks():
            db_path = self.get_db_path_recurse(benchmark)
            if not db_path:
                print "Cannot get db_path to benchmark: "+benchmark
                print "No db (cluster, recovery, or antibody)"
                continue
            db_paths.append(db_path)
        return db_paths

    def set_benchmarks_from_databases(self):
        """
        Set the benchmarks from the main_dir/databases directory
        :return:
        """
        if not os.path.exists(self.main_dir.get()+"/databases"):
            print "Could not find any databases to use for analysis.  Please make sure databases are in "+self.main_dir.get()+"/databases"
            print "Databases should end with .db or .db3 extension.  Naming format of databases: benchmark.features_type.scorefunction.db"
            return

        dbs = glob.glob(self.main_dir.get()+"/databases/*.db*")
        if len(dbs) == 0:
            print "Could not find any databases to use for analysis.  Please make sure databases are in "+self.main_dir.get()+"/databases"
            print "Databases should end with .db or .db3 extension.  Naming format of databases: benchmark.features_type.scorefunction.db"
            return

        nr_dbs = defaultdict()
        for db in sorted(dbs):
            print db

            if db.split('.')[-4] == "top":
                benchmark = os.path.basename(".".join(db.split('.')[:-3]))
            else:
                benchmark = os.path.basename(".".join(db.split('.')[:-2]))

            benchmark = benchmark+"."+db.split('.')[-2] #Add scorefunction

            if not nr_dbs.has_key(benchmark):
                self.benchmarks.append(benchmark)
                nr_dbs[benchmark] = " "

    def set_benchmarks(self, benchmarks):
        self.benchmarks = benchmarks

    def get_benchmarks(self):
        return self.benchmarks

    def get_db_path(self, benchmark, features_type = 'cluster'):
        db_dir = self.main_dir.get()+"/databases"

        ext_type = ".db"
        if features_type != "recoveries":
            features_type = features_type+"_features"
            ext_type = ".db3"
        db_path = db_dir+"/"+".".join(benchmark.split('.')[:-1])+"."+features_type+"."+benchmark.split('.')[-1]+ext_type;

        return db_path

    def get_db_path_recurse(self, benchmark, order = ["cluster", "recoveries", "antibody"]):

        for features_type in order:
            db_path = self.get_db_path(benchmark, features_type)
            if os.path.exists(db_path):
                return db_path
            else:
                print "No db: "+db_path
                continue

        return None


    ################ Main Functions #######################
    def run_features(self, type):
        outdir = self._setup_outdir(["plots", "features", self.out_dir_name.get()])
        db_dir = self.main_dir.get()+"/databases"

        if len(self.benchmarks) == 0:
            return

        if not os.path.exists(db_dir): sys.exit("Please run Features reporter and copy databases to main Benchmark Analysis DIR /databases.")

        creator = json_creator.JsonCreator(outdir, type)

        if self.reference_db.get() and os.path.exists(self.reference_db.get()):
            creator.add_sample_source_info(self.reference_db.get(), "ref", True)

        for benchmark in self.benchmarks:
            id = benchmark
            id = id.replace("talaris2013", "talaris")

            #Fix up the name so it is not too long.  Will add options to do this manually later:
            """
            id = id.replace("cluster", "clus")
            id = id.replace("exclude", "excl")
            id = id.replace("include", "incl")

            """

            db_path = self.get_db_path(benchmark, type)
            if not os.path.exists(db_path):
                print " does not exist!"
                return

            creator.add_sample_source_info(db_path, id)

        creator.save_json(outdir+"/"+type+"_"+self.out_dir_name.get()+".json")
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
        os.rmdir("build")
        print "Plots ignore any set filters.  To plot with filters, create new databases through query..."
        print "Complete..."

    def plot_groups(self):

        plot_dir = self._setup_outdir(["plots", self.out_dir_name.get()], False)
        cmd = "./pooled_data/plot_groups.R "+plot_dir

        db_paths = self._setup_db_paths_recurse()

        if not db_paths:
            print "No db paths found. Returning"
            return

        for db_path in db_paths:
            print db_path
            cmd = cmd+" "+db_path

        print cmd
        os.system(cmd)
        print "Complete. Plots are in: "+plot_dir

    def plot_cdr_data(self, data_type):
        plot_dir = self._setup_outdir(["plots", self.out_dir_name.get()], False)
        cmd = "./pooled_data/plot_cdr_data.R "+plot_dir

        if not self.reference_db.get():
            print "Reference database needs to be set"
            return

        db_paths = self._setup_db_paths_recurse()

        if not db_paths:
            print "No db paths found. Returning"
            return

        cmd = cmd +" "+self.reference_db.get()+" "+data_type
        for db_path in db_paths:
            print db_path
            cmd = cmd+" "+db_path

        print cmd
        os.system(cmd)
        print "Complete. Plots are in: "+plot_dir

    def plot_individuals(self):
        plot_dir = self._setup_outdir(["plots", self.out_dir_name.get()], False)
        cmd = "./pooled_data/plot_individuals.R "+plot_dir

        db_paths = self._setup_db_paths_recurse()

        if not db_paths:
            print "No db paths found. Returning"
            return

        for db_path in db_paths:
            cmd = cmd+" "+db_path

        print cmd
        os.system(cmd)
        print "Complete."
        print "Plots are in: "+plot_dir

########################################################################################################################
###   CDRData
########################################################################################################################
class CDRDataInfo:
    """
    Simple class for holding and accessing Cluster and Length data for a particular Decoy.
    """
    def __init__(self, name, benchmark, decoy):
        self.name = name
        self.benchmark = benchmark
        self.decoy = decoy

        self.data = defaultdict()

    def get_data(self):
        return self.data

    def get_value_for_cdr(self, cdr):
        return self.data[cdr]

    def set_data(self, data):
        """
        Dictionary for each CDR: L1, L2, L3, H1, H2, H3
        """
        self.data = data

    def set_value(self, cdr, value):
        self.data[cdr] = value

    def get_data_tuple(self):
        """
        Return tuple with data at each position: [L1, L2, L3, H1, H2, H3]
        If Camelid, will only return H data.
        """

    def has_data(self, cdr):
        if self.data.has_key(cdr):
            return True
        else:
            return False

    def is_camelid(self):
        """
        Return True if missing light chain data
        """
        if (not self.has_data('L1')) and (not self.has_data('L2')) and (not self.has_data('L3')):
            return True
        else:
            return False

class CDRData:
    """
    Class holding cluster and length data from cluster or antibody features database.
    """
    def __init__(self, name, column_name, is_camelid = False):
        self.is_camelid = is_camelid
        self.name = name
        self.column_name = column_name

        self.all_data = defaultdict()
        if is_camelid:
            self.cdrs = ["H1", "H2", "H3"]
        else:
            self.cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]


        self.native_data = None


    def _get_stmt(self, column_name):
        if not self.is_camelid:
            stmt = "SELECT "+ \
                        "structures.input_tag as decoy,"+ \
                        "H1."+column_name+" as H1," +\
                        "H2."+column_name+" as H2," +\
                        "H3."+column_name+" as H3," +\
                        "L1."+column_name+" as L1," +\
                        "L2."+column_name+" as L2," +\
                        "L3."+column_name+" as L3"

            stmt = stmt + """
                    FROM
                        cdr_clusters as L1,
                        cdr_clusters as L2,
                        cdr_clusters as L3,
                        cdr_clusters as H1,
                        cdr_clusters as H2,
                        cdr_clusters as H3,
                        structures
                    WHERE
                        structures.struct_id = L1.struct_id AND
                        structures.struct_id = L2.struct_id AND
                        structures.struct_id = L3.struct_id AND
                        structures.struct_id = H1.struct_id AND
                        structures.struct_id = H2.struct_id AND
                        structures.struct_id = H3.struct_id AND
                        L1.CDR = 'L1' AND
                        L2.CDR = 'L2' AND
                        L3.CDR = 'L3' AND
                        H1.CDR = 'H1' AND
                        H2.CDR = 'H2' AND
                        H3.CDR = 'H3'
                    """
        else:
            stmt = "SELECT "+ \
                        "structures.input_tag as decoy,"+ \
                        "H1."+column_name+" as H1," +\
                        "H2."+column_name+" as H2," +\
                        "H3."+column_name+" as H3"

            stmt = stmt + """
                    FROM
                        cdr_clusters as H1,
                        cdr_clusters as H2,
                        cdr_clusters as H3,
                        structures
                    WHERE
                        structures.struct_id = H1.struct_id AND
                        structures.struct_id = H2.struct_id AND
                        structures.struct_id = H3.struct_id AND
                        H1.CDR = 'H1' AND
                        H2.CDR = 'H2' AND
                        H3.CDR = 'H3'
                    """
        return stmt

    def add_data(self, benchmark, con):
        pass

    def get_benchmark_data(self, benchmark):
        return self.all_data[benchmark]

    def get_benchmark_data_for_decoy(self, benchmark, decoy):
        return self.all_data[benchmark][decoy]

    def _get_add_data(self, benchmark, con, column_name):
        data = defaultdict()
        cur = con.cursor()
        for row in cur.execute(self._get_stmt(column_name)):
            d = CDRDataInfo(self.name, benchmark, row[0])
            d.set_value('H1', row[1])
            d.set_value('H2', row[2])
            d.set_value('H3', row[3])

            if not self.is_camelid:
                d.set_value('L1', row[4])
                d.set_value('L2', row[5])
                d.set_value('L3', row[6])
            data[row[0]] = d

        self._add_data(benchmark, data)

    def _add_data(self, benchmark, decoy_data_map):
        """
        Add data in the form of a dict of decoy:DataTriple
        """
        if not self.all_data.has_key(benchmark):
            self.all_data[benchmark] = defaultdict()
        self.all_data[benchmark] = decoy_data_map

    def get_native_data(self):
        return self.native_data

    def set_native_data_from_rosetta(self, pdb_path):
        pass

    def set_native_data_input_tag(self, con, input_tag):
        self._set_native_data_input_tag(con, input_tag, self.column_name)

    def set_native_data(self, data):
        self.native_data = data


    def _set_native_data_input_tag(self, con, input_tag, column_name):
        pass

    def get_concatonated_map(self, cdr = None):
        """
        Returns a defaultDic:
        Default:
            decoy: CDRDataInfo

        CDR to get back cdr_value, decoy for sorting on cdr_value
        """

        result_data = defaultdict()
        for benchmark in self.all_data:
            for decoy in self.all_data[benchmark]:
                triple = self.all_data[benchmark][decoy]
                if isinstance(triple, CDRDataInfo): pass

                if cdr:
                    result_data[(triple.get_value_for_cdr(cdr), decoy)] = triple
                else:
                    result_data[decoy] = triple
        return result_data

class CDRLengthData(CDRData):
    def __init__(self, is_camelid = False):
        CDRData.__init__(self, "length", is_camelid)

    def add_data(self, benchmark, con):
        self._get_add_data(benchmark, con, "length")

    def set_native_data_from_rosetta(self, pdb_path):
        """
        p = pose_from_pdb(pdb_path)
        ab_info = AntibodyInfo(p)

        native_data = CDRDataInfo(self.name, "native", pdb_path)
        for cdr in self.cdrs:
            cdr_enum = ab_info.get_CDR_name_enum(cdr)
            value = ab_info.get_CDR_length(cdr_enum)
            native_data.set_value(cdr, value)
        """
        native_data = CDRDataInfo(self.name, "native", pdb_path)
        pose = pose_from_pdb(pdb_path)
        clusterer = CDRClusterer(pose)

        data = defaultdict()
        for cdr in self.cdrs:
            length = clusterer.get_length(cdr)
            native_data.set_value(cdr, length)
        self.native_data = native_data

class CDRClusterData(CDRData):
    def __init__(self, is_camelid = False):
        CDRData.__init__(self, "cluster", is_camelid)

    def add_data(self, benchmark, con):
        self._get_add_data(benchmark, con, "fullcluster")

    def set_native_data_from_rosetta(self, pdb_path):
        """
        p = pose_from_pdb(pdb_path)
        ab_info = AntibodyInfo(p)

        native_data = CDRDataInfo(self.name, "native", pdb_path)
        for cdr in self.cdrs:
            cdr_enum = ab_info.get_CDR_name_enum(cdr)
            value = ab_info.get_cluster_name(ab_info.get_CDR_cluster(cdr_enum).cluster())
            native_data.set_value(cdr, value)
        """
        native_data = CDRDataInfo(self.name, "native", pdb_path)
        pose = pose_from_pdb(pdb_path)
        clusterer = CDRClusterer(pose)
        cdrs = ["L1", "L2", "L3", "H1", "H2", "H3"]

        data = defaultdict()
        for cdr in self.cdrs:
            clusterer.dihedrals = []
            cluster = clusterer.get_fullcluster(cdr)[0]

            native_data.set_value(cdr, cluster)
        self.native_data = native_data

###################
## Helper Functions
###################





if __name__ == "__main__":
    main()

