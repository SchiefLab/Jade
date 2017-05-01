import sqlite3

import pandas


# PyIgD
from jade.pymol_jade.PyMolScriptWriter import *
from jade.antibody.cdr_data.CDRDataTypes import *
from jade.antibody.decoy_data.DecoyDataTypes import *
from jade.basic.sequence import fasta
from jade.basic.filters.DataFilters import *
from jade.basic.filters.FilterSettings import *
from jade.basic.threading.Threader import *
from jade.basic.pandas import PandasDataFrame
from jade.RAbD_BM.AnalysisInfo import *

# Rosetta Tools
import jade.rosetta_jade.FeaturesJsonCreator as json_creator


class CompareAntibodyDesignStrategies:
    """
    Class mainly for comparing different Antibody Design strategies using our Features Databases.
    """

    def __init__(self, db_dir, out_dir_name, strategies=[], jsons = []):
        """

        :param db_dir:
        :param out_dir_name:
        :param strategies:
        :param jsons:
        """

        #Init construction options
        self.db_dir = StringVar(value=db_dir)
        self.out_dir_name = StringVar(value=out_dir_name)
        self.native_path = None
        self.strategies = strategies
        self.jsons = jsons

        #Init Classes and data
        self.db_paths = defaultdict()
        self.filter_settings = FilterSettings()

        #Init Components
        self._init_default_options()
        self._init_paths()
        self._init_scores()
        self._init_default_scores()
        self._init_cdrs()



    def _init_default_options(self):

        self.scorefxn = "talaris2014"

        self.main_dir = StringVar()

        self.reference_db = StringVar()
        self.clustal_soft_wrap = IntVar(value=100)
        self.reload_scores = IntVar(value=1)
        self.top_n = IntVar(value=15)
        self.top_n_combined = IntVar(value=15)

        self.features_hbond_set = IntVar()
        self.features_hbond_sets = ["", "_min_hbond_analysis", "_no_hbond_analysis"]
        self.features_hbond_set.set(1)
        self.query_hbonds = IntVar(value=0)

        self.is_camelid = IntVar()
        self.is_camelid.set(0)
        self.top_total_percent = IntVar()
        self.top_total_percent.set(10)
        self.backround_features = IntVar()
        self.backround_features.set(1)

        self.rosetta_extension = StringVar()
        self.rosetta_extension.set("linuxclangrelease")

        self.individual_analysis = IntVar(value = 1)
        self.combined_analysis = IntVar(value = 0)

        self.load_origin_pdbs = IntVar(value = 1)
        self.pyigclassify_dir = StringVar()


    def _init_paths(self):
        self.cdr_rel_path = "DBOUT/cdr_pdbs_redun_by_cdr_overhang_3"
        self.weblogo_rel_path = "DBOUT/weblogos"
        self.ab_db_rel_path = "DBOUT/website"

        self.redun_db_name = "antibody_database_redundant.db"
        self.nr_db_name = "antibody_database_rosetta_design.db"

    def _init_scores(self):
        total_scores = TotalDecoyData()
        dg_scores = dGDecoyData()
        dsasa_scores = dSASADecoyData()
        top10_by_10 = dGTotalScoreSubset()

        #hbonds_int = IntHbondDecoyData()
        sc_value = SCValueDecoyData()
        unsats = DeltaUnsatsPerAreaDecoyData()

        self.scores = [total_scores, dg_scores, dsasa_scores, top10_by_10, sc_value, unsats]

    def _init_default_scores(self):
        #Setup the scores that are on
        self.score_names = self._get_score_names()
        self.scores_on = defaultdict()
        for score_name in self.score_names:
            self.scores_on[score_name] = IntVar(value=0)

        self.scores_on["dG"].set(1)
        self.scores_on["dG_top_Ptotal"].set(1)
        self.scores_on["total"].set(1)
        self.scores_on["delta_unsats_per_1000_dSASA"].set(1)

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

    def _get_score_names_on(self):
        names = []
        for name in self._get_score_names():
            if self.scores_on[name]:
                names.append(name)
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

    def _setup_outdir_individual(self, subdirs=[], use_outdir_name=False):

        return self._setup_outdir(["analysis_individual", self.out_dir_name.get()] + subdirs, use_outdir_name)

    def _setup_outdir_combined(self, subdirs=[]):
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
        :rtype: list of DecoyData
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

        hb_loader = InterfaceHBondDecoyDataLoader()

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
        """
        Get a particular score class.  May or may not be initialized.  See setup_scores.
        :param score_name:
        :rtype: DecoyData
        """
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

    def _setup_cdr_types(self, data_type, is_camelid, features_type='antibody'):
        if data_type == "length":
            data_class = CDRLengthData(self.native_path, is_camelid)

        elif data_type == "cluster":
            data_class = CDRClusterData(self.native_path, is_camelid)

        elif data_type == "sequence":
            data_class = CDRSequenceData(self.native_path, is_camelid)

        elif data_type == "aligned_sequence":
            data_class = CDRAlignedSequenceData(self._setup_outdir_individual(['clustal']),
                                                self._setup_outdir_combined(['clustal']), self.native_path, is_camelid)

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

    def set_strategies_from_json_infos(self):
        """
        Uses self.json, which are AnalysisInfo classes, to populate.

        :return:
        """

        if not self.jsons:
            print "JSONS NOT SET!"
            return

        nr_dbs = defaultdict()
        for info in self.jsons:
            if not isinstance(info, AnalysisInfo): sys.exit()

            db = info.get_features_db()
            print db

            ##Example naming convention: 'ch103_5_CDR_prelim.norm_ab_features.db'
            strategy = info.get_exp()
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
        print "Creating: "+outdir + "/" + type + "_" + plot_name + ".json"
        creator.run_json(self.backround_features.get())

        # pwd = os.getcwd()
        # os.chdir("build")
        # os.system("cp -r *../"+outdir)
        # os.system("mv build build_old")
        # os.chdir(pwd)

        print "Plots ignore any set filters.  To plot with filters, create new databases through query..."
        print "Complete..."

    def get_pandas_dataframe(self):
        """
        Gets a pandas Dataframe for all
        :rtype: pandas.DataFrame
        """
        dfs = []
        output_names = ["strategy"] #Controls the order of the output names.

        for score in self._setup_scores(use_all=True):
            if score.name == "dG_top_Ptotal":continue

            if isinstance(score, DecoyData): pass

            output_names.append(score.name)

            df = score.get_pandas_dataframe()
            dfs.append(df)

        cdr_types = ["length", "cluster", "sequence", "aligned_sequence"]

        for t in cdr_types:
            cdr_data = self._setup_cdr_types(t, self.is_camelid.get())

            cdr_names = [cdr for cdr in cdr_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                output_names.append("_".join([cdr, cdr_data.name]))
            df = cdr_data.get_pandas_dataframe(cdr_names)
            dfs.append(df)

        df = PandasDataFrame.drop_duplicate_columns(pandas.concat(dfs, axis=1, join="outer"))
        return df

    def get_top_from_dataframe(self, score_name):
        """
        Gets a pandas Dataframe for top
        :rtype: pandas.DataFrame
        """
        df = self.get_pandas_dataframe()
        dfs = []
        for strategy in self.get_strategies():
            dfs.append(df["strategy" == strategy].sort_values(score_name)[0:self.top_n.get()-1])

        df = PandasDataFrame.drop_duplicate_columns(pandas.concat(dfs))
        return df

    def get_top_dataframe_by_all_scores(self):
        """
        Get a pandas DataFrame for top, grouped by the type of score that is on.
        :rtype: pandas.DataFrame
        """
        dfs = []
        score_names = self._get_score_names_on()
        for score_name in score_names:
            if score_name == "dG_top_Ptotal":continue
            df = self.get_top_from_dataframe(score_name)
            df["by_score_group"] = score_name
            dfs.append(df)
        df = PandasDataFrame.drop_duplicate_columns(pandas.concat(dfs))
        return df

    def output_all_data_as_excel_file(self, top = True):
        final_dfs = []
        final_tab_names = []


        #All Data:
        dfs, names = self.get_csv_data(top, summary = False)
        final_dfs.extend(dfs)
        final_tab_names.extend(names)

        #Summary Data:
        dfs, names = self.get_csv_data(top, summary = True)
        final_dfs.extend(dfs)
        final_tab_names.extend(dfs)





    def output_csv_data(self, top = False, summary = False):
        """
        Output a CSV file of combined or individual data.

        """
        output_dfs, output_names = self.get_csv_data(top, summary)
        for index, df in enumerate(output_dfs):
            name = output_names[index]
            df.to_csv(name+".csv")

            print "Wrote: "+name+".csv"

    def get_csv_data(self, top = False, summary = False):
        """
        Get data by converting everything to a pandas dataframe first.
        For now, one function pretty much does everything.

        :rtype: [pandas.Dataframe],[str]

        """
        final_dfs = []
        final_names = []

        dfs = []
        venn_dfs = [] #Best decoys/data seen in all score classes of the top n together.
        output_names = ["strategy"] #Controls the order of the output names.

        venn2_cat = ['dG', 'total']
        venn2_dfs = [] #Venn on dG and Total score.
        for score in self._setup_scores(use_all=True):
            if score.name == "dG_top_Ptotal":continue

            if isinstance(score, DecoyData): pass


            output_names.append(score.name)


            venn_dfs.append(score.get_pandas_dataframe(top_n=self.top_n.get()))
            df = score.get_pandas_dataframe()
            print df.tail()
            dfs.append(df)
            if score.name in venn2_cat:
                venn2_dfs.append(df)

        cdr_types = ["length", "cluster", "sequence", "aligned_sequence"]

        for t in cdr_types:
            cdr_data = self._setup_cdr_types(t, self.is_camelid.get())

            cdr_names = [cdr for cdr in cdr_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                output_names.append("_".join([cdr, cdr_data.name]))
            df = cdr_data.get_pandas_dataframe(cdr_names)
            #df.to_csv(self._setup_outdir_individual()+"/test_cdr_"+cdr_data.name+".csv")
            dfs.append(df)
            venn_dfs.append(df)
            venn2_dfs.append(df)


        combined_scores = PandasDataFrame.drop_duplicate_columns(pandas.concat(dfs, axis=1, join="outer")) #Drop Duplicates
        combined_scores = combined_scores[output_names]
        combined_scores = combined_scores.apply(pandas.to_numeric, errors='ignore')
        combined_scores.index.name = "decoy"

        if top:
            venn_df = PandasDataFrame.drop_duplicate_columns(pandas.concat(venn_dfs, axis=1, join="inner"))
            venn_df = venn_df[output_names]
            venn_df.sort_values(['strategy', 'dG'])
            #venn_df.to_csv(open(self._setup_outdir_individual()+"/ind_per_model_ven_top_"+str(self.top_n.get())+".csv", "w"))

            venn2_df = PandasDataFrame.drop_duplicate_columns(pandas.concat(venn_dfs, axis=1, join="inner"))
            venn2_df = venn2_df[output_names]
            venn2_df.sort_values(['strategy', 'dG'])
            #venn2_df.to_csv(open(self._setup_outdir_individual()+"/ind_per_model_ven_dG_total_top_"+str(self.top_n.get())+".csv", "w"))

            score_dfs=[]
            score_names = self._get_score_names_on()
            for score_name in score_names:
                if score_name == "dG_top_Ptotal":
                    sort_name = 'dG'
                else:
                    sort_name = score_name

                strat_dfs=[]
                for strategy in self.get_strategies():

                    """
                    if score_name == "dG_top_Ptotal":
                        score = self._get_score(score_name)
                        decoy_list = score.get_ordered_decoy_list(strategy, self.top_n.get())
                        df = combined_scores[combined_scores.index.isin(decoy_list)]
                        df.sort(columns=[sort_name]) #This sort is not working, have no idea why.
                    else:
                        df = combined_scores[combined_scores['strategy'] == strategy].sort(score_name)[0:self.top_n.get()] #Best N
                        df.sort(columns=[sort_name])
                    """
                    #Generally, will use order to get top scores.
                    score = self._get_score(score_name)
                    decoy_list = score.get_ordered_decoy_list(strategy, self.top_n.get())
                    df = combined_scores[combined_scores.index.isin(decoy_list)]
                    df.sort_values([sort_name])

                    strat_dfs.append(df)
                top_df = pandas.concat(strat_dfs)

                top_df["by_score_group"] = score_name
                score_dfs.append(top_df)

            top_df = PandasDataFrame.drop_duplicate_columns(pandas.concat(score_dfs))
            name_order=["by_score_group"]
            name_order.extend(output_names)
            top_df = top_df[name_order]

            if self.individual_analysis.get():

                if summary:
                    top_df = top_df.apply(pandas.to_numeric, errors='ignore')

                    final_dfs.append(top_df.groupby(by=["strategy", "by_score_group"]).describe(exclude=['object']))
                    final_names.append(self._setup_outdir_individual()+"/per_strategy_summary_top")
                else:
                    final_dfs.append(top_df)
                    final_names.append(self._setup_outdir_individual()+"/ind_per_model_top")


            if self.combined_analysis.get():
                dfs = []
                for score_name in self._get_score_names_on():

                    score = self._get_score(score_name)
                    decoy_list = score.get_ordered_decoy_list_all(self.top_n.get())
                    df = combined_scores[combined_scores.index.isin(decoy_list)]
                    df["by_score_group"] = score_name

                    dfs.append(df)
                df = PandasDataFrame.drop_duplicate_columns(pandas.concat(dfs))
                name_order=["by_score_group"]
                name_order.extend(output_names)
                df = df[name_order]
                if summary:
                    #df.groupby(by="strategy").describe().to_csv(self._setup_outdir_combined()+"/com_summary_top_by_"+score.name+".csv")
                    df = df.apply(pandas.to_numeric, errors='ignore')

                    final_dfs.append(df.groupby(by=["by_score_group"]).describe(exclude=['object']))
                    final_names.append(self._setup_outdir_combined()+"/com_summary_top.csv")
                else:
                    final_dfs.append(df)
                    final_names.append(self._setup_outdir_combined()+"/com_per_model_top")


        else:

            if self.individual_analysis.get():
                combined_scores.sort_values(['strategy', 'dG'])
                if summary:

                    final_dfs.append(combined_scores.groupby(by="strategy").describe(exclude=['object']))
                    final_names.append(self._setup_outdir_individual()+"/per_strategy_summary_all.csv")
                else:
                    final_dfs.append(combined_scores)
                    final_names.append(self._setup_outdir_individual()+"/ind_per_model_all")

            if self.combined_analysis.get():
                combined_scores.sort_values(['dG'])
                if summary:
                    #combined_scores.groupby(by="strategy").describe().to_csv(self._setup_outdir_combined()+"/")
                    final_dfs.append(combined_scores.describe(exclude=['object']))
                    final_names.append(self._setup_outdir_combined()+"/com_summary_all")
                else:
                    final_dfs.append(combined_scores)
                    final_names.append(self._setup_outdir_combined()+"/com_per_model_all")

        #How to add Native Line?
        #Maybe a 'print native info' function...

        return final_dfs, final_names

    def output_stats(self):
        """
        Depracated in favor of dataframe summaries.
        """

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
                            data = score.get_data_for_decoy(strategy, decoy)

                            if not score.name == "dG_top_Ptotal":
                                line = os.path.basename(decoy) + "\t" + get_str(data.score)
                            else:
                                line = os.path.basename(decoy)

                            for a_score in all_scores:
                                a_data = a_score.get_data_for_decoy(strategy, decoy)
                                if a_score.name == "combined_str_score": continue
                                if a_score.name == score.name: continue

                                line = line + "\t" + get_str(a_data.score)
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
                        data = score.get_data_for_decoy(strategy, decoy)
                        if not score.name == "dG_top_Ptotal":

                            line = os.path.basename(decoy) + "\t" + get_str(data.score)
                        else:
                            line = os.path.basename(decoy)
                        for a_score in all_scores:
                            a_data = a_score.get_data_for_decoy(strategy, decoy)
                            if a_score.name == "combined_str_score": continue
                            if a_score.name == score.name: continue

                            line = line + "\t" + get_str(a_data.score)
                        OUT.write(line + "\n")
                    OUT.close()

        if self.individual_analysis.get():
            print "Outputting individual stats"
            output_score_extra_stats()
        if self.combined_analysis.get():
            print "Outputting combined stats"
            output_all_stats()
        print "Complete"

    def copy_top(self):
        def copy_top_strategy():

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

                    if self.load_origin_pdbs:
                        if not self.pyigclassify_dir.get() or not os.path.exists(self.pyigclassify_dir.get()):
                            print "Origin PDB not set or does not exist. Disable this feature or set a correct directory."
                            return

                        make_pymol_session_on_top_ab_include_native_cdrs(decoy_list, load_as, out_dir, out_dir, score.get_outname(),
                                                                         self.pyigclassify_dir.get()+"/"+self.cdr_rel_path,
                                                                         top_num=top_n, native_path=self.native_path)
                    else:
                        make_pymol_session_on_top(decoy_list, load_as, out_dir, out_dir, score.get_outname(),
                                                                         top_num=top_n, native_path=self.native_path)
                    SCORELIST.close()

        def copy_top_combined():
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

                if self.load_origin_pdbs:
                    if not self.pyigclassify_dir.get() or not os.path.exists(self.pyigclassify_dir.get()):
                        print "Origin PDB not set or does not exist. Disable this feature or set a correct directory."
                        return

                    make_pymol_session_on_top_ab_include_native_cdrs(decoy_list, load_as, outdir_top_pdbs, outdir_top_sessions,
                                                                     score.get_outname(), self.pyigclassify_dir.get()+"/"+self.cdr_rel_path,
                                                                     top_num=top_n, native_path=self.native_path)
                else:
                    make_pymol_session_on_top(decoy_list, load_as, outdir_top_pdbs, outdir_top_sessions, score.get_outname(),
                                                                     top_num=top_n, native_path=self.native_path)

                SCORELIST.close()

        if self.individual_analysis.get():
            print "Outputting individual sessions"
            copy_top_strategy()
        if self.combined_analysis.get():
            print "Outputting combined sessions"
            copy_top_combined()

    def copy_all_models(self):

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

                if self.load_origin_pdbs:
                    if not self.pyigclassify_dir.get() or not os.path.exists(self.pyigclassify_dir.get()):
                        print "Origin PDB not set or does not exist. Disable this feature or set a correct directory."
                        return

                    make_pymol_session_on_top_ab_include_native_cdrs(decoy_list, load_as, out_dir, out_dir, score.get_outname(),
                                                                     self.pyigclassify_dir.get()+"/"+self.cdr_rel_path,
                                                                     top_num=None, native_path=self.native_path)
                else:
                    make_pymol_session_on_top(decoy_list, load_as, out_dir, out_dir, score.get_outname(),
                                                                     top_num=None, native_path=self.native_path)

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
            out_dir = self._setup_outdir_combined(["all_sessions"])

            if self.load_origin_pdbs:
                if not self.pyigclassify_dir.get() or not os.path.exists(self.pyigclassify_dir.get()):
                    print "Origin PDB not set or does not exist. Disable this feature or set a correct directory."
                    return

                make_pymol_session_on_top_ab_include_native_cdrs(decoy_list, load_as, out_dir, out_dir, score.get_outname(),
                                                                 self.pyigclassify_dir.get()+"/"+self.cdr_rel_path,
                                                                 top_num=None, native_path=self.native_path)
            else:
                make_pymol_session_on_top(decoy_list, load_as, out_dir, out_dir, score.get_outname(),
                                                                 top_num=None, native_path=self.native_path)

    def run_clustal_omega(self, processors, output_format="fasta", extra_options=""):

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
                fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, self.native_path, "Native",
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
                    fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, self.native_path, "Native",
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

                fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, self.native_path, "Native",
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

    def run_clustal_omega_on_all_combined(self, processors, output_format, extra_options=""):

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

        fasta.output_fasta_from_pdbs_biopython(decoy_header_dict, fasta_path, self.native_path, "native",
                                               self.is_camelid.get())
        clustal_runner = ClustalRunner(fasta_path)
        #clustal_runner.set_threads(processors)
        clustal_runner.set_hard_wrap(self.clustal_soft_wrap.get())
        clustal_runner.set_extra_options(extra_options)
        clustal_runner.set_output_format(output_format)
        clustal_runner.output_alignment(clustal_dir, clustal_name)

        print "Complete"

    def output_len_or_clus_alignment(self, alignment_type, features_type='antibody'):

        is_camelid = self.is_camelid.get()

        top_n = self.top_n.get()
        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type)
        len_class = self._setup_cdr_types("length", self.is_camelid.get(), features_type)

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
            if self.native_path:
                header = "#score\t\tmatches"
            else:
                header = "#score\t"
            for cdr in cdr_names:
                header = header + "\t\t" + cdr

            header += "\tdecoy"
            OUTFILE.write(header + "\n")

            all_data = score.get_concatonated_map()
            all_type_data = type_data.get_concatonated_map()

            if self.native_path:
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

    def output_len_or_clus_recovery(self, alignment_type, features_type='antibody'):
        len_class = self._setup_cdr_types("length", self.is_camelid.get(), features_type)
        if not self.native_path:
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

        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type)
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

    def output_len_or_clus_enrichment(self, alignment_type, features_type='antibody'):


        def _add_enrichments(self, label, decoys, type_data):
            if isinstance(data_class, CDRData): pass
            #OUTFILE = open(outdir + "/" + "cdr_type_enrichment_" + alignment_type + "_.txt", 'w')

            cdr_names = [cdr for cdr in type_data.cdrs if cdr in self._setup_cdrs()]
            for cdr in cdr_names:
                enrichments = calculate_enrichments(data_class, cdr, decoys)
                OUTFILE = open(outdir+"/" + "cdr_type_enrichment_"+alignment_type+"_"+label+"_"+cdr+".txt", 'w')
                OUTFILE.write("#group\tcount\tperc\n")


                for c in sorted([[enrichments[c].count, c] for c in enrichments], reverse=True):
                    t = c[1]; #c[0] is the actual counts we are sorting on.
                    perc = enrichments[t]
                    if isinstance(perc, Perc): pass

                    p = perc.get_perc_decimal()
                    OUTFILE.write(str(t)+"\t\t"+str(perc.count)+"\t\t"+perc.get_formated_perc(p)+"\n")

                OUTFILE.close()


        ##Individual Enrichments

        top_n = self.top_n.get()
        data_class = self._setup_cdr_types(alignment_type, self.is_camelid.get(), features_type)

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