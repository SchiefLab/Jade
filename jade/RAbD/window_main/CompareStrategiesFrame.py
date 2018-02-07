import tkFileDialog
import tkSimpleDialog

from jade.basic.TKinter.Listbox import AutoListbox
from jade.RAbD.AnalyzeAntibodyDesigns import *


class CompareStrategiesFrame( Frame ):
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
        self.all_strategies_listbox = AutoListbox(self)
        self.current_strategies_listbox = AutoListbox(self)

        self.main_label = Label(self, text="Strategies", justify=CENTER)

        # Setup CDRs
        self.L_chain_buttons = []
        self.H_chain_buttons = []
        for cdr_name in ["L1", "L2", "L3", "H1", "H2", "H3"]:
            button = Checkbutton(self, text=cdr_name, variable=self.compare_designs.cdrs[cdr_name])
            if re.search("L", cdr_name):
                self.L_chain_buttons.append(button)
            else:
                self.H_chain_buttons.append(button)

                # self.separator = Separator(self, orient = HORIZONTAL)

        self.individual_analysis = Checkbutton(self, text = "Individual Analysis", variable=self.compare_designs.individual_analysis)
        self.combined_analysis = Checkbutton(self, text = "Combined Analysis", variable = self.compare_designs.combined_analysis)

        # self.db_dir_entry = Entry(self._tk_, textvariable = self.compare_designs.out_dir_name, justify = CENTER)
        self.out_dir_entry = Entry(self, textvariable=self.compare_designs.out_dir_name, justify=CENTER)

        # self.root_dir_label = Label(self._tk_, text = "Root Directory", justify = CENTER)
        self.out_dir_label = Label(self, text="Analysis Name", justify=CENTER)


    def sho_tk(self):
        r=0
        c=0
        self.main_label.grid(row=0,column=0, columnspan=6, padx=5, pady=5)
        self.all_strategies_listbox.grid(row=r + 1, column=c + 0, columnspan=3, padx=7, pady=5)
        self.current_strategies_listbox.grid(row=r + 1, column=c + 3, columnspan=3, padx=7, pady=5)

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


        self.out_dir_label.grid(row=r+5, column=c+0, columnspan=3, padx=5, pady=5, sticky=W+E)
        self.out_dir_entry.grid(row=r+5, column=c+3, columnspan=3, padx=5, pady=5, sticky=W+E)


        self.all_strategies_listbox.bind("<Double-Button-1>",
                                         lambda event: self.add_to_current(self.all_strategies_listbox,
                                                                           self.current_strategies_listbox))
        self.all_strategies_listbox.bind("<Button-2>", lambda event: self.show_strat_items())

        self.current_strategies_listbox.bind("<Double-Button-1>",

                                             lambda event: self.delete_current(self.current_strategies_listbox))





    ##################################
    def populate_all_strategies(self):
        self.all_strategies_listbox.delete(0, END)
        for strategy in self.compare_designs.strategies:
            self.all_strategies_listbox.insert(END, strategy)

        self.all_strategies_listbox.autowidth(100)
        self.current_strategies_listbox.autowidth(100, self.compare_designs.strategies)

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

    def add_main_strategy(self):
        strategy_name = tkSimpleDialog.askstring(title="Strategy", prompt="Strategy Name")
        if not strategy_name:
            return

        strategy_path = tkFileDialog.askdirectory(initialdir=self.main_gui.current_dir, title="Strategy Path")
        self.compare_designs.strategies.append(strategy_name)
        self.compare_designs.db_paths[strategy_name] = strategy_path

        self.all_strategies_listbox.insert(END, strategy_name)

    def get_full_strategy_list(self):
        strategies = self.current_strategies_listbox.get(0, END)
        return strategies