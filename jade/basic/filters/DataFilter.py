



class DataFilter:
    """
    Class for storing filter information - such as not equal to values or cutoffs.  Used to create custom queries for data.
    """
    def __init__(self, name, type = "boolean"):
        self.name = name
        self.type = type
        self.required_tables = []
        self.required_wheres = []

    def get_required_tables(self):
        return self.required_tables

    def get_required_wheres(self):
        return self.required_wheres
