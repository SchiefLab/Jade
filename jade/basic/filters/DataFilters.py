from jade.basic.filters.DataFilter import DataFilter

class TotalScoreCutoffFilter(DataFilter):
    """
    Filter to remove structures with total_score greater than a particular value
    """
    def __init__(self, value):
        DataFilter.__init__(self, "total", "score")
        self.set_value(value)

    def set_value(self, value):
        self.value = value

        self.required_tables = [
            "score_types",
            "structure_scores"
        ]
        self.required_wheres = [
            "score_types.score_type_name='total_score'",
            "structure_scores.score_type_id = score_types.score_type_id",
            "structures.struct_id = structure_scores.struct_id",
            "structure_scores.score_value <= "+repr(self.value)]

class dGCutoffFilter(DataFilter):
    """
    Filter to remove structures with LH_A dG greater than a particular value
    """
    def __init__(self, value):
        DataFilter.__init__(self, "dG", "score")
        self.set_value(value)

    def set_value(self, value):
        self.value = value
        self.required_tables = [
            "interfaces"
        ]
        self.required_wheres = [
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface = 'LH_A'",
            "interfaces.dG <= "+repr(self.value)
        ]

class dSASACutoffFilter(DataFilter):
    """
    Filter to remove dSASA greater than some value
    """
    def __init__(self, value):
        DataFilter.__init__(self, "dSASA", "score")
        self.set_value(value)

    def set_value(self, value):
        self.value = value
        self.required_tables = [
            "interfaces"
        ]
        self.required_wheres = [
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface = 'LH_A'",
            "interfaces.dSASA >= "+repr(self.value)
        ]

class H3ExtendedFilter(DataFilter):
    """
    Filter to remove kinked H3 structures
    """
    def __init__(self):
        DataFilter.__init__(self, "h3_extended")
        self.required_tables = ["ab_h3_kink_metrics"]
        self.required_wheres = [
            "structures.struct_id = ab_h3_kink_metrics.struct_id",
            "ab_h3_kink_metrics.kink_type != 'EXTENDED'"]