import copy

#PyIgClassify
from jade.antibody.CDRClusterer import *
from jade.basic.sequence import fasta

#Modules
from jade.antibody.cdr_data.CDRData import *
from jade.basic.sequence.ClustalRunner import *
from jade.basic.structure.BioPose import *
from jade.basic.structure.Structure import AntibodyStructure

class CDRLengthData(CDRData):
    def __init__(self, native_path, is_camelid = False):
        CDRData.__init__(self, "length", native_path, is_camelid)

    def add_data(self, strategy, con):
        self._get_add_data(strategy, con, "length")

    def _setup_native_data(self, pdb_path):
        if not pdb_path: return None
        else:
            self._set_native_data_from_biopose(pdb_path)

    def _set_native_data_from_biopose(self, pdb_path):
        """
        p = pose_from_pdb(pdb_path)
        ab_info = AntibodyInfo(p)

        native_data = CDRDataInfo(self.name, "native", pdb_path)
        for cdr in self.cdrs:
            cdr_enum = ab_info.get_CDR_name_enum(cdr)
            value = ab_info.get_CDR_length(cdr_enum)
            native_data.set_value(cdr, value)
        """
        native_data = CDRDataInfo(self.name, "native", pdb_path)
        pose = BioPose(pdb_path)
        clusterer = CDRClusterer(pose)

        data = defaultdict()
        for cdr in self.cdrs:
            length = int(clusterer.get_length(cdr))
            native_data.set_value(cdr, length)
        self.native_data = native_data

class CDRClusterData(CDRData):
    def __init__(self, native_path, is_camelid = False):
        CDRData.__init__(self, "cluster", native_path, is_camelid)

    def add_data(self, strategy, con):
        self._get_add_data(strategy, con, "fullcluster")

    def _setup_native_data(self, pdb_path):
        if not pdb_path: return None
        else:
            self._set_native_data_from_biopose(pdb_path)

    def _set_native_data_from_biopose(self, pdb_path):
        """
        p = pose_from_pdb(pdb_path)
        ab_info = AntibodyInfo(p)

        native_data = CDRDataInfo(self.name, "native", pdb_path)
        for cdr in self.cdrs:
            cdr_enum = ab_info.get_CDR_name_enum(cdr)
            value = ab_info.get_cluster_name(ab_info.get_CDR_cluster(cdr_enum).cluster())
            native_data.set_value(cdr, value)
        """
        native_data = CDRDataInfo(self.name, "native", pdb_path)
        pose = BioPose(pdb_path)
        clusterer = CDRClusterer(pose)

        for cdr in self.cdrs:
            clusterer.dihedrals = []
            cluster = clusterer.get_fullcluster(cdr)[0]

            native_data.set_value(cdr, cluster)
        self.native_data = native_data

class CDRSequenceData(CDRData):
    def __init__(self, native_path, is_camelid = False):
        CDRData.__init__(self, "sequence", native_path, is_camelid)

    def add_data(self, strategy, con):
        self._get_add_data(strategy, con, "sequence")

    def _setup_native_data(self, pdb_path):
        if not pdb_path: return None
        else:
            self._set_native_data_from_biopose(pdb_path)

    def _set_native_data_from_biopose(self, pdb_path):
        native_data = CDRDataInfo(self.name, "native", pdb_path)
        pose = BioPose(pdb_path)
        ab_structure = AntibodyStructure()
        clusterer = CDRClusterer(pose)
        for cdr in self.cdrs:
            seq = ab_structure.get_cdr_seq(pose, cdr)


            native_data.set_value(cdr, seq)

        self.native_data = native_data

class CDRdSASAData(CDRData):
    def __init__(self, native_path, is_camelid = False):
        CDRData.__init__(self, "dSASA", native_path, is_camelid)

    def add_data(self, strategy, con):
        self._get_add_data(strategy, con, "ab_ag_dSASA")

class CDRAlignedSequenceData(CDRSequenceData):
    """
    Uses Clustal Omega to align sequences from each database.
    """
    def __init__(self, individual_clustal_outdir, combined_clustal_outdir, native_path, is_camelid = False):
        CDRSequenceData.__init__(self, native_path, is_camelid)
        self.name = "aligned_sequence"

        self.ind_clustal_outdir = individual_clustal_outdir
        self.com_clustal_outdir = combined_clustal_outdir

        self.original_seq_data = defaultdict()
        self.concatonated_data = None

        if self.native_data:
            self.old_native_data = copy.deepcopy(self.native_data)

    def add_data(self, strategy, con):
        CDRSequenceData.add_data(self, strategy, con)
        self.original_seq_data[strategy] = copy.deepcopy(self.all_data[strategy])

        self._run_clustal_set_data(strategy)
        self.concatonated_data = None

    def get_concatonated_map(self, only_cdr = None, decoy_list = None, use_saved_data = True):

        if self.concatonated_data and use_saved_data:
            return self.concatonated_data

        final_result_data = defaultdict()

        result_data = defaultdict()
        for strategy in self.original_seq_data:
            for decoy in self.original_seq_data[strategy]:
                if decoy_list and decoy not in decoy_list:
                    continue
                triple = self.original_seq_data[strategy][decoy]
                if isinstance(triple, CDRDataInfo): pass

                if only_cdr:
                    result_data[(triple.get_value_for_cdr(only_cdr), decoy)] = triple
                else:
                    result_data[decoy] = triple

        for cdr in self.cdrs:
            if only_cdr and cdr != only_cdr: continue

            fasta_path = self._make_fasta_for_concatonated(result_data, cdr, "concatonated")

            clustal_outname = ".".join(os.path.basename(fasta_path).split(".")[0:-1])+".aln"

            runner = ClustalRunner(fasta_path)
            runner.set_extra_options("--force")
            runner.set_hard_wrap(500)
            runner.output_alignment(self.com_clustal_outdir, clustal_outname)
            clustal_path = self.com_clustal_outdir+"/"+clustal_outname

            new_data = self._parse_clustal_for_concatonated(result_data, cdr, clustal_path)
            print "Created Clustal file for reference: "+clustal_path

            for decoy in new_data:
                if only_cdr:
                    cdr_decoy = (new_data[decoy].get_value_for_cdr(cdr), decoy)
                else:
                    cdr_decoy = decoy

                if not final_result_data.has_key(cdr_decoy):
                    final_result_data[cdr_decoy] = new_data[decoy]
                else:
                    final_result_data[cdr_decoy].set_value_for_cdr(cdr, new_data[decoy].get_value_for_cdr(cdr))

        self.concatonated_data = final_result_data
        return final_result_data

    def _run_clustal_set_data(self, strategy):

        if not os.path.exists(self.ind_clustal_outdir):
            os.mkdir(self.ind_clustal_outdir)

        #Make Fasta, run clasta


        for cdr in self.cdrs:
            print "working on: "+strategy+" "+cdr
            fasta_path = self._make_fasta(cdr, strategy)
            clustal_outname = ".".join(os.path.basename(fasta_path).split(".")[0:-1])+".aln"

            runner = ClustalRunner(fasta_path)
            runner.set_extra_options("--force")
            runner.set_hard_wrap(500)
            runner.output_alignment(self.ind_clustal_outdir, clustal_outname)

            clustal_path = self.ind_clustal_outdir+"/"+clustal_outname

            self._parse_clustal_set_data(cdr, strategy, clustal_path)
            print "Created Clustal file for reference: "+clustal_path

    def _make_fasta_for_concatonated(self, concatonated_data, cdr, name):

        outpath = self.com_clustal_outdir+"/list_"+name+"_"+cdr+".fasta"


        OUTFILE = open(outpath, 'w')

        if self.old_native_data:
            seq = self.old_native_data.get_value_for_cdr(cdr)
            fasta.write_fasta(seq, "native", OUTFILE)

        for decoy in concatonated_data.keys():

            #print repr(decoy)
            data = self.original_seq_data[ concatonated_data[decoy].strategy][decoy]
            #print repr(data)
            if not isinstance(data, CDRDataInfo): sys.exit()

            seq = data.get_value_for_cdr(cdr)
            #print seq
            fasta.write_fasta(seq, decoy, OUTFILE)

        OUTFILE.close()
        return outpath


    def _make_fasta(self, cdr, strategy):
        """
        Makes Fasta for clustal, returns path
        """

        outpath = self.ind_clustal_outdir+"/all_"+strategy+"_"+cdr+".fasta"


        OUTFILE = open(outpath, 'w')

        if self.old_native_data:
            seq = self.old_native_data.get_value_for_cdr(cdr)
            fasta.write_fasta(seq, "native", OUTFILE)

        for decoy in self.original_seq_data[ strategy ]:

            data = self.original_seq_data[ strategy][decoy]
            #print repr(data)
            if not isinstance(data, CDRDataInfo): sys.exit()

            seq = data.get_value_for_cdr(cdr)
            #print seq
            fasta.write_fasta(seq, decoy, OUTFILE)

        OUTFILE.close()
        return outpath

    def _parse_clustal_for_concatonated(self, concatonated_data, cdr, clustal_path):

        new_data = defaultdict()
        if not os.path.exists(clustal_path):
            sys.exit("clustal path not good")

        INFILE = open(clustal_path, 'r')
        INFILE.readline()
        INFILE.readline()

        for line in INFILE:
            line = line.strip()
            lineSP  = line.split()
            if len(lineSP) < 2:
                continue

            decoy = lineSP[0]
            aligned_seq = lineSP[1]

            #A way to skip the ending *** alignment stuff hopefully.  Lets see if this works.
            if decoy == "native":
                if not new_data.has_key(decoy):
                    new_info = CDRDataInfo(self.name, "native", decoy)
                    new_info.set_value(cdr, aligned_seq)

                    new_data[decoy] = new_info
                else:
                    new_data[decoy].set_value(cdr, aligned_seq)
                continue

            if not concatonated_data.has_key(decoy):
                continue

            if not new_data.has_key(decoy):
                new_info = CDRDataInfo(self.name, concatonated_data[decoy].strategy, decoy)
                new_info.set_value(cdr, aligned_seq)

                new_data[decoy] = new_info
            else:
                new_data[decoy].set_value(cdr, aligned_seq)

        INFILE.close()
        #print "Setting new data for "+strategy+" - "+cdr

        return new_data

    def _parse_clustal_set_data(self, cdr, strategy, clustal_path):

        new_data = defaultdict()
        if not os.path.exists(clustal_path):
            sys.exit("clustal path not good")

        INFILE = open(clustal_path, 'r')
        INFILE.readline()
        INFILE.readline()

        for line in INFILE:
            line = line.strip()
            lineSP  = line.split()
            if len(lineSP) < 2:
                continue

            decoy = lineSP[0]
            aligned_seq = lineSP[1]

            #A way to skip the ending *** alignment stuff hopefully.  Lets see if this works.
            if decoy == "native":
                self.native_data.set_value(cdr, aligned_seq)

            if not decoy in self.original_seq_data[strategy]:
                continue

            self.all_data[strategy][decoy].set_value(cdr, aligned_seq)

            #if not new_data.has_key(decoy):
            #    new_info = CDRDataInfo(self.name, strategy, decoy)
            #    new_info.set_value(cdr, aligned_seq)

            #    new_data[decoy] = new_info
            #else:
            #    new_data[decoy].set_value(cdr, aligned_seq)

        INFILE.close()
        #print "Setting new data for "+strategy+" - "+cdr
        #self.all_data[strategy] = new_data
