#!/usr/bin/env python

import sys, os
from jade.rosetta_jade.RunRosettaBenchmarks import RunRosettaBenchmarks
from jade.rosetta_jade.RunRosetta import RunRosetta

#from overrides import overrides


class RunBenchmarksRAbD( RunRosettaBenchmarks ):
    """
    Benchmark class specifically for RAbD


    Details:

        ALL INPUT PDBs should go into

            project_root/datasets

        Typically, you will have multiple directories - native, relaxed, etc.

          This is specified as a benchmark using 'input_pdb_type' in your json file.

        ALL PDBLISTs for benchmarking should go into

            project_root/datasets/pdblists



    """
    def __init__(self):
        RunRosettaBenchmarks.__init__(self, program = "antibody_designer")

        self._current_settings["CDR"] = "ALL"
        self._current_settings_ordered_keys.append("CDR")


        self.dataset_root_dir = "datasets"
        self.pdblist_dir = self.dataset_root_dir+"/pdblists"
        self.instructions_dir = "instructions"

        if self.options.l or self.options.s:
            sys.exit("PDBLIST should be created in datasets/pdblists.  See antibody_design repo for an example.")


        if not os.path.exists(self.instructions_dir):
            os.mkdir(self.instructions_dir)

    #@overrides
    def run_benchmark(self, benchmark_names, benchmark_options):
        """
        Run a single benchmark with options.

        :param benchmark_names: List of benchmark names
        :param benchmark_options: List of benchmark options
        :return:
        """

        separate_cdrs = benchmark_options[benchmark_names.index("separate_cdrs")]

        #Special case for mintype - We now do this within code itself!
        '''
        if benchmark_options[benchmark_names.index("mintype")] != "relax":
            if benchmark_options[benchmark_names.index("inner_cycles")] == 1:
                benchmark_options[benchmark_names.index("inner_cycles")] = 2
        '''


        if separate_cdrs:
            for cdr in self._get_designable_cdrs():
                self._current_settings["CDR"] = cdr
                RunRosettaBenchmarks.run_benchmark(self, benchmark_names, benchmark_options)
        else:
            self._current_settings["CDR"] = "ALL"
            RunRosettaBenchmarks.run_benchmark(self, benchmark_names, benchmark_options)


    #@overrides
    def _get_output_string(self):

        if not self.options.separate_job_per_pdb:
            self.options.l = self._get_pdb_list_fname()

        s = RunRosettaBenchmarks._get_output_string(self)

        #Decoys
        s = s + (" -in:path "+self.dataset_root_dir+"/"+self._current_settings["input_pdb_type"])

        #Instructions
        s = s + " -cdr_instructions " + self._create_instructions(self.instructions_dir+"/"+os.path.basename(self._get_make_out_path())+".instruct")
        return s

    #@overrides
    def _get_pdb_list_fname(self):
        return ".".join([self.pdblist_dir+"/"+self._current_settings["dataset"],
                        self._current_settings["l_chain"]+".PDBLIST.txt"])

    #@overrides
    def _get_job_name(self):
        return self.extra_options.get_exp()+"."+os.path.basename(self._get_make_out_path())+"."+self._current_settings["l_chain"]

    ### Helper Functions ###


    def _create_instructions(self, output_path):
        extra_lines=[]
        extra_lines.append("\n".join( str(line) for line in self.extra_options.json_dict["base_cdr_instruction_lines"]))

        print output_path
        #print repr(self._current_settings)

        extra_lines.append("ALL MinProtocol MinType "+self._current_settings["mintype"])
        extra_lines.append("ALL FIX")

        seq_design_cdrs = self.extra_options.get_benchmarks_of_key("seq_design_cdrs")
        graft_design_cdrs = self.extra_options.get_benchmarks_of_key("graft_design_cdrs")

        current_cdr = self._current_settings["CDR"]
        if current_cdr != "ALL" and current_cdr in graft_design_cdrs:
            print "Adding " +current_cdr + " to graftdesign. "
            extra_lines.append(self._current_settings["CDR"]+" GraftDesign Allow")
        elif current_cdr == "ALL":
            for cdr in self.extra_options.json_dict["graft_design_cdrs"]:
                extra_lines.append(cdr+" GraftDesign Allow")

        if current_cdr != "ALL" and current_cdr in seq_design_cdrs:
            print "Adding " + current_cdr + " to seqdesign. "
            extra_lines.append(self._current_settings["CDR"]+" SeqDesign Allow")
        elif current_cdr == "ALL":
            for cdr in self.extra_options.json_dict["seq_design_cdrs"]:
                extra_lines.append(cdr+" SeqDesign Allow")

        FILE = open(output_path, "w")
        line = "\n".join(str(s) for s in extra_lines)


        #print line
        FILE.write(line)
        FILE.close()
        return output_path

    def _get_designable_cdrs(self):

        designable_cdrs = []
        seq_design_cdrs = self.extra_options.get_benchmarks_of_key("seq_design_cdrs")
        graft_design_cdrs = self.extra_options.get_benchmarks_of_key("graft_design_cdrs")

        for cdr in ["L1","L2","L3","H1","H2","H3"]:
            if cdr in seq_design_cdrs or cdr in graft_design_cdrs:
                designable_cdrs.append(cdr)

        return designable_cdrs