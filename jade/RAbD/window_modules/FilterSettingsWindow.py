from Tkinter import *
from jade.basic.filters.FilterSettings import FilterSettings



class FilterSettingsWindow:
    def __init__(self, settings):
        self.settings = settings
        #self.compare_strategies =  compare_strategies
        if isinstance(self.settings, FilterSettings):pass
        #if isinstance(self.compare_strategies, CompareAntibodyDesignStrategies):pass

        self.extra_required_tables_var = StringVar()
        self.extra_required_where_var = StringVar()

    def reset_tables(self):
        self.extra_required_tables_var.set("")
    def setup_sho_gui(self, main, r = 0, c = 0):
        self.set_tk(main)
        self.sho_tk(r, c)

    def set_tk(self, main):
        main.title("Setup Filters")
        self.members_check = []
        self.members_entry = []

        for energy in self.settings.energy_types:
            member_check = Checkbutton(main, text = energy, variable = self.settings.energies_enabled[energy], justify = LEFT)
            member_entry = Entry(main, textvariable = self.settings.energy_cutoffs[energy], justify = LEFT)

            self.members_check.append(member_check)
            self.members_entry.append(member_entry)

        self.check_h3_filter = Checkbutton(main, text="Filter out extended H3", variable = self.settings.h3_filter, justify = LEFT)
        #self.check_apply_as_group = Checkbutton(main, text="Apply filters as new group", variable = self.settings.apply_as_group, justify = LEFT)

        self.label_name = Label(main, text = "New Group/DIR Name ")
        self.entry_name = Entry(main, textvariable = self.settings.name, justify = LEFT)

        self.label_custom_filters = Label(main, text = "Custom Filters")
        self.label_required_tables = Label(main, text = "Required Tables (,)")

        #self.tables = Listbox(main)
        self.label_where = Label(main, text = "WHERE")
        self.entry_required_tables = Entry(main, textvariable = self.extra_required_tables_var)
        self.entry_required_where = Entry(main, textvariable = self.extra_required_where_var)

        self.buttom_add_cutom_filter = Button(main, text = "Add Custom Filter", command = lambda: self.add_custom_filter_settings_from_entries(), justify = CENTER)


        self.done = Button(main, text = "Done", command = lambda: main.destroy())

    def sho_tk(self, r, c):

        self.label_name.grid(row=r, column = c+0, sticky = W+E,padx=3, pady = 10); self.entry_name.grid(row = r, column = c+1, padx = 10, pady = 10)
        r+=1
        for i in range(0, len(self.members_check)):
            self.members_check[i].grid(row=r+i, column = c+0, sticky = W, padx = 3, pady = 5)
            self.members_entry[i].grid(row=r+i, column = c+1, sticky = W+E, padx = 10, pady = 5)

        current_rows = len(self.members_check)
        self.check_h3_filter.grid(row = r+current_rows+1, column = c+1, sticky = W, columnspan = 1, padx = 3, pady = 3)
        current_rows+=1

        #self.check_apply_as_group.grid(row = r+current_rows+1, column = c+1, sticky = W, columnspan = 1, padx = 3, pady = 5)

        self.label_custom_filters.grid(row=r+current_rows+1, column = c+0, columnspan = 2, sticky = W+E, padx=3, pady = 5)
        self.label_required_tables.grid(row = r+current_rows+2, column = c+0, padx = 3, pady = 2)
        self.entry_required_tables.grid(row = r+current_rows+2, column = c+1, padx = 3, pady = 2)
        self.label_where.grid(row = r+current_rows+3, column = c+0, padx = 3, pady = 2)
        self.entry_required_where.grid(row = r+current_rows+3, column = c+1, padx = 3, pady = 2)
        self.buttom_add_cutom_filter.grid(row = r+current_rows+4, column = c+0, columnspan = 2, padx = 5, pady = 2)

        current_rows = current_rows+4

        self.done.grid(row = r+current_rows+2, column = c+0, columnspan = 2, sticky = W+E, padx = 8, pady=5)

    def get_settings(self):
        return self.settings

    def populate_listbox(self):
        pass

    def add_table(self):
        pass

    def add_custom_filter_settings_from_entries(self):
        tables = self.extra_required_tables_var.get()

        for table in tables.split(','):
            self.settings.extra_required_tables.append(table)
            self.settings.extra_required_where.append(table+".struct_id = structures.struct_id")
        self.settings.extra_required_where.append(self.extra_required_where_var.get())

        self.extra_required_tables_var.set("")
        self.extra_required_where_var.set("")
