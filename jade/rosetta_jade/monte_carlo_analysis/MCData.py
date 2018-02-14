from collections import defaultdict

import os
import sys
from collections import defaultdict
from jade.RAbD_BM import tools as bm_tools
from jade.RAbD_BM.AnalysisInfo import AnalysisInfo

from jade.basic import path as path_tools
import re

class MCData:
    def __init__(self, dataset_dir, analysis_info):
        """

        :param dataset_dir: str
        :param analysis_info: AnalysisInfo
        :return:
        """

        if isinstance(analysis_info, AnalysisInfo): pass
        self.dataset_dir = dataset_dir
        self.data = defaultdict(lambda: [])
        self.analysis_info = analysis_info

        self.data = defaultdict()

        self.native_pdblist = self.dataset_dir+"/"+self.analysis_info.bm_info.get_dataset()+".PDBLIST.txt"


        if not os.path.exists(self.native_pdblist):
            sys.exit("Native PDBList path, "+self.native_pdblist+" does not exist.  Make sure to concatonate both lambda and kappa abs")

        FILE = open(self.native_pdblist)
        self.pdb_ids = [line.strip() for line in FILE]
        FILE.close()

        self._read_pdb_files()
        self._calculate_cum_accepts()
        self._calculate_delta_e()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.get(item)

    def __setitem__(self, key, value):
        self.data[key] = value

    def get(self, pdb_id):
        return self.data[pdb_id]

    def get_pdb_id_of_path(self, file_path):
        for pdb_id in self.pdb_ids:
            pdb_id = pdb_id.split(".")[0]
            if re.search(pdb_id, file_path):
                return pdb_id

        sys.exit("unknown PDB ID of resultant benchmark PDB")

    def _calculate_cum_accepts(self):
        for pdb_id in self.data:
            print pdb_id
            for cy_data in self.data[pdb_id]:
                cum = 0
                for n in cy_data.values:
                    print repr(n)
                    cum = cum + cy_data.get_acceptance(n)
                    cy_data.set_cum_acceptance(n, cum)

    def _calculate_delta_e(self):
        for pdb_id in self.data:
            for cy_data in self.data[pdb_id]:
                for n in cy_data.values:
                    if n == 1:
                        cy_data.set_delta_e(n, 0)
                    else:
                        cy_data.set_delta_e(n, cy_data.get_final_e(n) - cy_data.get_final_e(n-1))

    def _read_pdb_files(self):
        file_paths = bm_tools.get_pdb_paths("decoys/"+self.analysis_info.get_dir(), self.analysis_info.get_exp(), self.analysis_info.get_name_match())
        for file_path in file_paths:
            self._parse_file_data(file_path)

    def _parse_file_data(self, file_path):
        print "Reading "+file_path
        INFILE = path_tools.open_file(file_path)
        pdb_id = self.get_pdb_id_of_path(file_path)
        cy_data = CycleData(file_path, pdb_id)

        for line in INFILE:
            line = line.strip()
            lineSP = line.split()

            if not lineSP:                continue
            if not lineSP[0] == "ACCEPT": continue
            if not lineSP[1] == "LOG":    continue

            if lineSP[3] == "FINAL":
                if lineSP[4] == "END":
                    cy_data.set_protocol_final_e(lineSP[5])
                else:
                    cy_data.set_final_energy(lineSP[4], lineSP[5])
            else:
                if lineSP[3] == "0":
                    cy_data.set_protocol_start_e(lineSP[4])
                elif lineSP[3] == "-1":
                    cy_data.set_protocol_native_e(lineSP[4])
                else:
                    cy_data.set_mc_data(lineSP[3], lineSP[4], lineSP[5])

        INFILE.close()
        if not self.data.has_key(pdb_id):
            self.data[pdb_id] = []

        #print cy_data.values.keys()
        #for i in cy_data.values:
        #    print repr(i)+" "+repr(cy_data[i])
        #    print repr(i)+" "+repr(cy_data.get_final_e(i))
        self.data[pdb_id].append(cy_data)