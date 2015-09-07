import tools
import sys
import os


class BenchmarkInfo:
    """
    Simple Class for holding info for a particular benchmark.  This will then be used during analysis
    """
    def __init__(self):
        self._initialize_data()

    def _initialize_data(self):
        self.date = tools.get_today()
        self.full_name = None
        self.final_name = None
        self.scorefunction = "talaris2013"
        self.expected_nstruct = 200
        self.found_nstruct = 200
        self.exp_batch = "20_lam_kap"
        self.decoy_dir = None
        self.rel_log_path = None

    def get_full_name(self):
        return self.full_name

    def get_final_name(self):
        return self.final_name

    def get_exp_batch(self):
        return self.exp_batch

    def get_date(self):
        return self.date

    def get_scorefunction_name(self):
        return self.scorefunction

    def get_decoy_path(self):
        return self.decoy_path

    def has_log_path(self):
        if self.rel_log_path:
            return True
        else:
            return False

    def has_required_nstruct(self, tolerance):
        if (int(self.expected_nstruct) - int(self.found_nstruct)) <= int(tolerance):
            return True
        else:
            return False

    ####################################################################################################################
    ### Setters
    ####################################################################################################################

    def set_data(self, full_name, final_name, scorefunction, decoy_path, exp_batch, rel_log_path = None, date = None, expected_nstruct=200, found_nstruct = 200 ):
        self.set_full_name(full_name)
        self.final_name = final_name
        self.scorefunction = scorefunction
        self.decoy_dir = decoy_path
        self.exp_batch = exp_batch
        self.rel_log_path = rel_log_path

        if date:
            self.date = date

        self.expected_nstruct = expected_nstruct
        self.found_nstruct = found_nstruct

    def set_full_name(self, full_name):
        #print full_name
        self.full_name = full_name
        fnSP = full_name.split('.')
        if len(fnSP) == 4:
            print "Setting exp, pdb_type, min, etc values from full name"
            self.exp = fnSP[0]
            self.pdb_type = fnSP[1]
            self.min_type = fnSP[2]
            self.etc = fnSP[3]
        elif len(fnSP) == 5:
            self.antigen_type = fnSP[0]
            self.exp = fnSP[1]
            self.pdb_type = fnSP[2]
            self.min_type = fnSP[3]
            self.etc = fnSP[4]

    def initialize_from_list_file_line(self, lineSP):
        self.date = lineSP[0]
        self.set_full_name(lineSP[1])
        self.final_name = lineSP[2]
        self.scorefunction = lineSP[3]
        self.expected_nstruct = lineSP[4]
        self.found_nstruct = lineSP[5]
        self.exp_batch = lineSP[6]
        self.decoy_dir = lineSP[7]

        if len(lineSP) > 8:
            self.rel_log_path = lineSP[8]

def parse_benchmark_list(benchmark_list):
    if not os.path.exists(benchmark_list):
        sys.exit("Benchmark list file not found.  Could not parse")


    benchmark_array = []
    FILE = open(benchmark_list, 'r')
    for line in FILE:
        line = line.strip()
        if not line: continue
        if line.startswith("#"): continue


        lineSP = line.split()

        #print repr(lineSP)
        benchmark = BenchmarkInfo()
        #print repr(lineSP)
        benchmark.initialize_from_list_file_line(lineSP)

        benchmark_array.append(benchmark)
    FILE.close()

    return benchmark_array