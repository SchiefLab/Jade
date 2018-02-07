
from collections import defaultdict

from jade.basic.filters.DataFilter import DataFilter


class StatementCreator:
    """
    Simple class for constructing a statement - allows on-the-fly addition of filters, cutoffs, etc.
    """
    def __init__(self):
        self.component_types = ["SELECT", "FROM", "WHERE", "ORDER BY"]
        self.components = defaultdict()
        for component in self.component_types:
            self.components[component] = []

        self.data_filters = []

    def _add_string_or_strings(self, component_type, string_or_strings):
        if not component_type in self.component_types:
            print "Unrecognized component type: "+component_type
            print "Recognized statement component types: "+repr(self.component_types)
            return

        t = type(string_or_strings)
        if t == str:
            self.components[component_type].append(string_or_strings)
        elif t == list:
            self.components[component_type].extend(string_or_strings)
        elif t == tuple:
            self.components[component_type].extend(list(string_or_strings))
        else:
            print "Unknown type - cannot add component: "+repr(t)
            return

    def add_SELECT_string_or_strings(self, string_or_strings):
        self._add_string_or_strings("SELECT", string_or_strings)

    def add_FROM_string_or_strings(self, strings_or_strings):
        self._add_string_or_strings("FROM", strings_or_strings)

    def add_WHERE_string_or_strings(self, string_or_strings):
        self._add_string_or_strings("WHERE", string_or_strings)

    def add_ORDER_BY_string_or_strings(self, string_or_strings):
        self._add_string_or_strings("ORDER BY", string_or_strings)

    def add_data_filter(self, data_filter):
        if isinstance(data_filter, DataFilter): pass
        self.data_filters.append(data_filter)

    def _add_data_filters_to_components(self):
        for data_filter in self.data_filters:
            if isinstance(data_filter, DataFilter):pass
            required_tables = data_filter.get_required_tables()
            required_wheres = data_filter.get_required_wheres()

            for table in required_tables:
                if not table in self.components["FROM"]:
                    self.components["FROM"].append(table)
            for select in required_wheres:
                if not select in self.components["WHERE"]:
                    self.components["WHERE"].append(select)

    def create_statement(self):
        """
        Create the statment string from its components
        """
        self._add_data_filters_to_components()
        stmt = ""
        for component_type in self.component_types:
            if not len(self.components[component_type]): continue

            stmt = stmt+component_type+"\n"
            for i in range(0, len(self.components[component_type])):
                separator = ","
                if component_type == "WHERE":
                    separator = " AND"

                item = self.components[component_type][i]
                if i != len(self.components[component_type]) - 1:
                    item = item+separator

                stmt = stmt+item+"\n"

        #print stmt
        return stmt

