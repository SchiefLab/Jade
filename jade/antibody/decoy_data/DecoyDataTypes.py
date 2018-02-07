import numpy

from jade.antibody.decoy_data.DecoyData import *
from jade.basic.sql.StatementCreator import *



class TotalDecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "total")

    def add_data(self, strategy, con):

        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "structure_scores.score_value as total_score"])

        stmt_creator.add_FROM_string_or_strings([
            "structure_scores",
            "score_types",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "score_types.score_type_name='total_score'",
            "structure_scores.score_type_id = score_types.score_type_id",
            "structures.struct_id = structure_scores.struct_id"])

        stmt_creator.add_ORDER_BY_string_or_strings("score_value")

        self._get_add_data(strategy, stmt_creator, con)

class dGDecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "dG")

    def add_data(self, strategy, con):

        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "interfaces.dG as dG"])

        stmt_creator.add_FROM_string_or_strings([
            "interfaces",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface = "+repr(self.interface)])

        stmt_creator.add_ORDER_BY_string_or_strings([
            "dG"])

        self._get_add_data(strategy, stmt_creator, con)

class dSASADecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "dSASA", reverse_top=True)

    def add_data(self, strategy, con):
        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "interfaces.dSASA as dSASA"])

        stmt_creator.add_FROM_string_or_strings([
            "interfaces",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface =  "+ repr(self.interface)])

        stmt_creator.add_ORDER_BY_string_or_strings([
            "dSASA DESC"])

        self._get_add_data(strategy, stmt_creator, con)

class SCValueDecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "sc_value", reverse_top=True)

    def add_data(self, strategy, con):
        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "interfaces.sc_value as sc_value"])

        stmt_creator.add_FROM_string_or_strings([
            "interfaces",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface =  "+ repr(self.interface)])

        stmt_creator.add_ORDER_BY_string_or_strings([
            "sc_value DESC"])

        self._get_add_data(strategy, stmt_creator, con)

class DeltaUnsatsPerAreaDecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "delta_unsats_per_1000_dSASA")

    def add_data(self, strategy, con):
        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "interfaces.delta_unsatHbonds*1000/interfaces.dSASA as delta_unsats_per_1000_dSASA"])

        stmt_creator.add_FROM_string_or_strings([
            "interfaces",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface =  "+ repr(self.interface)])

        stmt_creator.add_ORDER_BY_string_or_strings([
            "delta_unsats_per_1000_dSASA"])

        self._get_add_data(strategy, stmt_creator, con)

class IntHbondDecoyData(DecoyData):
    """
    New way for int hbonds - added directly from IAM.
    """
    def __init__(self):
        DecoyData.__init__(self, "interface_hbonds", reverse_top=True)

    def add_data(self, strategy, con):
        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "interfaces.hbonds_int as hbonds_int"])

        stmt_creator.add_FROM_string_or_strings([
            "interfaces",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface =  "+ repr(self.interface)])

        stmt_creator.add_ORDER_BY_string_or_strings([
            "hbonds_int DESC"])

        self._get_add_data(strategy, stmt_creator, con)

class dGTotalScoreSubset(DecoyData):
    """
    dG of the top x percent of total score (for each strategy)
    """
    def __init__(self):

        self.total_scores = TotalDecoyData()
        DecoyData.__init__(self, "dG_top_Ptotal")

    def add_data(self, strategy, con, top_total_percent):
        self.top_total_percent = top_total_percent

        self.total_scores.add_data(strategy, con)
        total_cutoff = self.total_scores.get_top_x_percent_cutoff_value(strategy, self.top_total_percent)

        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "interfaces.dG as dG"])

        stmt_creator.add_FROM_string_or_strings([
            "interfaces",
            "structure_scores",
            "score_types",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "score_types.score_type_name='total_score'",
            "structure_scores.score_type_id = score_types.score_type_id",
            "structures.struct_id = structure_scores.struct_id",
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface = "+repr(self.interface),
            "structure_scores.score_value <= "+repr(total_cutoff)])

        stmt_creator.add_ORDER_BY_string_or_strings([
            "dG"])

        self._get_add_data(strategy, stmt_creator, con)


#It takes way to long to load the hbond data.  This data needs to be added to InterfaceFeatures reporter itself....

class InterfaceHBondDecoyDataLoader(DecoyData):
    """
    DecoyData class that holds the number of LH_A or L_H interface Hbonds, and energies.
    Very Slow to get this information.
     - SO - Subsequent Hbond classes accept this on construction and then parse its information
    """

    def __init__(self):
        DecoyData.__init__(self, "interface_hbonds")
        self.energy_data = defaultdict()
        self.count_data = defaultdict()



    def add_data(self, strategy, con):
        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "don_res.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "hb.energy as energy"])

        stmt_creator.add_FROM_string_or_strings([
            "interface_residues AS don_res",
            "interface_residues AS acc_res",
            "hbond_sites AS don",
            "hbond_sites AS acc",
            "hbonds AS hb",
            "structures"])

        stmt_creator.add_WHERE_string_or_strings([
            "((don_res.side== 'side1' ",
            "acc_res.side == 'side2') OR (don_res.side=='side2' ",
            "acc_res.side=='side1'))",
            "acc_res.interface == "+repr(self.interface),
            "acc_res.interface == don_res.interface",
            "don.resNum == don_res.resNum",
            "acc.resNum == acc_res.resNum",
            "hb.don_id == don.site_id",
            "hb.acc_id == acc.site_id",
            "don_res.struct_id == acc_res.struct_id",
            "acc_res.struct_id == acc.struct_id",
            "acc.struct_id == don.struct_id",
            "don.struct_id == hb.struct_id",
            "structures.struct_id == hb.struct_id"
        ])

        self._get_add_data(strategy, stmt_creator, con)

    def _get_add_data(self, strategy, stmt_creator, con):

        if isinstance(stmt_creator, StatementCreator): pass

        if self.filters:
            for filter in self.filters:
                stmt_creator.add_data_filter(filter)


        stmt = stmt_creator.create_statement()

        energies = defaultdict()
        decoy_to_struct_id = defaultdict()

        cur = con.cursor()

        print "Executing Hbond statement...."
        for row in cur.execute(stmt):
            decoy = row[1]
            struct_id = row[0]
            energy = row[2]

            if not decoy_to_struct_id.has_key(strategy):
                decoy_to_struct_id[strategy] = defaultdict()

            decoy_to_struct_id[strategy][decoy] = struct_id

            if not energies.has_key(decoy):
                energies[decoy] = []
                energies[decoy].append(float(energy))
            else:
                energies[decoy].append(float(energy))

        print "Done executing Hbond statement...."

        for decoy in energies:
            energy_mean = numpy.mean(energies[decoy])
            hbond_n = len(energies)
            if not self.count_data.has_key(strategy):
                self.count_data[strategy] = defaultdict()
                self.energy_data[strategy] = defaultdict()

            self.count_data[strategy][decoy] = hbond_n
            self.energy_data[strategy][decoy] = energy_mean

        self.decoy_to_struct_id = decoy_to_struct_id

class InterfaceHbondCountDecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "hbond_count", True, True)

    def setup_from_loader(self, hbond_loader):
        if not isinstance(hbond_loader, InterfaceHBondDecoyDataLoader): sys.exit()

        self.hbond_loader = hbond_loader
        self._init_from_hb_loader()

    def _init_from_hb_loader(self, strategy = None):

        def _organize_data(self, local_strategy):
            for decoy in self.hbond_loader.count_data[local_strategy]:
                counts = self.hbond_loader.count_data[local_strategy][decoy]
                struct_id = self.hbond_loader.decoy_to_struct_id[local_strategy][decoy]

                triple = DecoyDataTriple(local_strategy, struct_id, decoy, counts, self.get_outname(), self.name)

                if not self.all_data.has_key(local_strategy):
                    self.all_data[local_strategy] = defaultdict()

                self.all_data[local_strategy][decoy] = triple


        #################
        if strategy:
            _organize_data(self, strategy)

        else:
            for strategy in self.hbond_loader.count_data.keys():
                _organize_data(self, strategy)

    def add_data(self, strategy, con):
        self.hbond_loader.add_data(strategy, con)
        self._init_from_hb_loader(strategy)


class InterfaceHbondEnergyDecoyData(DecoyData):
    def __init__(self):
        DecoyData.__init__(self, "hbond_energy", True, True)
        self.hbond_loader = InterfaceHBondDecoyDataLoader()

    def setup_from_loader(self, hbond_loader):
        if not isinstance(hbond_loader, InterfaceHBondDecoyDataLoader): sys.exit()

        self.hbond_loader = hbond_loader
        self._init_from_hb_loader()

    def _init_from_hb_loader(self, strategy = None):

        def _organize_data(self, local_strategy):
            for decoy in self.hbond_loader.energy_data[local_strategy]:
                energy = self.hbond_loader.energy_data[local_strategy][decoy]
                struct_id = self.hbond_loader.decoy_to_struct_id[local_strategy][decoy]

                triple = DecoyDataTriple(local_strategy, struct_id, decoy, energy, self.get_outname(), self.name)

                if not self.all_data.has_key(local_strategy):
                    self.all_data[local_strategy] = defaultdict()

                self.all_data[local_strategy][decoy] = triple


        #################
        if strategy:
            _organize_data(self, strategy)

        else:
            for strategy in self.hbond_loader.energy_data.keys():
                _organize_data(self, strategy)

    def add_data(self, strategy, con):
        self.hbond_loader.add_data(strategy, con)
        self._init_from_hb_loader(strategy)

    """
  sele = "
    SELECT
    hb.energy as energy,
    don_res.interface as interface,
    don_res.struct_id as struct_id,
    hb_geom.AHdist as dis
  FROM
    interface_residues AS don_res,
    interface_residues AS acc_res,
    hbond_sites AS don,
    hbond_sites AS acc,
    hbonds AS hb,
    hbond_geom_coords as hb_geom
  WHERE
    ((don_res.side== 'side1' AND
    acc_res.side == 'side2') OR
    (don_res.side=='side2' AND
    acc_res.side=='side1')) AND
    acc_res.interface == don_res.interface AND
    don.resNum == don_res.resNum AND
    acc.resNum == acc_res.resNum AND
    hb.don_id == don.site_id AND
    hb.acc_id == acc.site_id AND
    hb.hbond_id == hb_geom.hbond_id AND
    don_res.struct_id == acc_res.struct_id AND
    acc_res.struct_id == acc.struct_id AND
    acc.struct_id == don.struct_id AND
    don.struct_id == hb.struct_id AND
    hb.struct_id == hb_geom.struct_id

"""

class CombinedStrDecoyData(DecoyData):
    """
    DecoyData class that has value as a string of the 3 main scores.
    Value held in DecoyDataTriple is a string: dG::total::dSASA for reference.
    """
    def __init__(self, filters, filt_name):
        DecoyData.__init__(self, "combined_score_summary", has_real_values=False)
        self.add_filters(filters, filt_name)

    def add_data(self, strategy, con):
        stmt_creator = StatementCreator()
        stmt_creator.add_SELECT_string_or_strings([
            "structures.struct_id as struct_id",
            "structures.input_tag as input_tag",
            "structure_scores.score_value as total_score",
            "interfaces.dG as dG",
            "interfaces.dSASA as dSASA"])

        stmt_creator.add_FROM_string_or_strings([
            "structure_scores",
            "score_types",
            "structures",
            "interfaces"])

        stmt_creator.add_WHERE_string_or_strings([
            "score_types.score_type_name='total_score'",
            "structure_scores.score_type_id = score_types.score_type_id",
            "structures.struct_id = structure_scores.struct_id",
            "structures.struct_id = interfaces.struct_id",
            "interfaces.interface = " + repr(self.interface)])

        self._get_add_data(strategy, stmt_creator, con)

    def _get_add_data(self, strategy, stmt_creator, con):

        if isinstance(stmt_creator, StatementCreator): pass

        if self.filters:
            for filter in self.filters:
                stmt_creator.add_data_filter(filter)
        stmt = stmt_creator.create_statement()

        data = defaultdict()
        cur = con.cursor()
        for row in cur.execute(stmt):
            score = "%.3f"%row[1]+"::"+"%.3f"%row[2]+"::"+"%.3f"%row[3]
            triple = DecoyDataTriple(strategy, row[0], row[1], score, self.get_outname, self.name)
            data[row[1]] = triple
        self._add_data(strategy, data)



