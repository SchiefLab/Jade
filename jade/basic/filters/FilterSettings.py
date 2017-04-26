from Tkinter import *
from collections import defaultdict


class FilterSettings:
    """
    Simple class for accessing energy cutoff settings for custom lists made using cutoffs.
    """
    def __init__(self):
        self.energy_types = ['dG', 'total', 'dSASA']
        #self.boolean_filters = [H3ExtendedFilter()]

        self.energies_enabled = defaultdict()
        self.energy_cutoffs = defaultdict()
        for energy in self.energy_types:
            self.energies_enabled[energy] = IntVar(value=0)
            self.energy_cutoffs[energy] = DoubleVar(value=0)

        self.name = StringVar(value="filtered_H3")
        #self.apply_as_group = IntVar(value = 1)
        #self.apply_to_all = IntVar(value = 0)
        self.h3_filter = IntVar(value = 0)

        self.extra_required_tables = []
        self.extra_required_where = []

    def get_energy_enabled(self, energy_type):
        return self.energies_enabled[energy_type].get()

    def get_energy_cutoff(self, energy_type):
        return self.energy_cutoffs[energy_type].get()

    def set_energy_enabled(self, energy_type, setting):
        self.energies_enabled[energy_type].set(int(setting))

    def set_energy_cutoff(self, energy_type, setting):
        self.energy_cutoffs[energy_type].set(float(setting))
