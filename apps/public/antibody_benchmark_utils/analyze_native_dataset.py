#!/usr/bin/env python

import os
import re
import sqlite3
import sys

from collections import defaultdict
from optparse import OptionParser
import numpy

#from python_modules import AbDbFunctions as ab_db

#Gets data from the Native DB for organized output tables into my thesis and benchmarking paper.

#####################################
##
## Experimenting with a basic python enum.  Its not quite an enum, but it should work instead of strings.
##  Add new values here.

length_e = 1
cluster_e = 2
dis_e = 3
sequence_e = 4

ag_contacts_e = 5
ag_contacts_nres_e = 6
dSASA_e = 7
dG_e = 8
charge_e = 9
aromatics_e = 10

delta_unsats_e = 11
sc_value_e = 12
packstat_e = 13
nres_int_e = 14
#####################################


class MetaData:
    def __init__(self, struct_id):

        self.struct_id = struct_id

        self.cluster_data = defaultdict()
        self.cdr_data = defaultdict()
        self.interface_data = None

    def get_struct_id(self):
        return self.struct_id


    def get_cdr_cluster_data(self, cdr):
        return self.cluster_data[cdr]

    def get_cdr_physical_data(self, cdr):
        return self.cdr_data[ cdr ]

    def get_interface_data(self):
        return self.interface_data


    def set_cdr_cluster_data(self, cdr, cluster_data):
        assert(isinstance(cluster_data, CDRClusterData))
        self.cluster_data[cdr] = cluster_data

    def set_cdr_physical_data(self, cdr, cdr_data):
        assert(isinstance(cdr_data, CDRPhysicalData))
        self.cdr_data[ cdr ] = cdr_data

    def set_interface_data(self, interface_data):
        self.interface_data = interface_data


class CDRClusterData:
    """
    Class for holding cluster information.
    """
    def __init__(self, CDR, length, cluster, dis, sequence):
        self.CDR = CDR
        self.data = defaultdict()
        self.data[length_e] = length
        self.data[cluster_e] = cluster
        self.data[dis_e] = dis
        self.data[sequence_e] = sequence

    def get_cdr(self):
        return self.CDR

    def get_data(self):
        return self.data

    def has_value(self, enum):
        return self.data.has_key(enum)

    def get_value(self, enum):
        return self.data[enum]

class CDRPhysicalData:
    """
    Class for holding cdr physical information.
    """
    def __init__(self, CDR, ag_contacts, ag_contact_nres, dSASA, dG, charge, aromatics):
        self.CDR = CDR
        self.data = defaultdict()
        self.data[ag_contacts_e] = ag_contacts
        self.data[ag_contacts_nres_e] = ag_contact_nres
        self.data[dSASA_e] = dSASA
        self.data[dG_e] = dG
        self.data[charge_e] = charge
        self.data[aromatics_e] = aromatics

    def get_cdr(self):
        return self.CDR

    def has_value(self, enum):
        return self.data.has_key(enum)

    def get_value(self, enum):
        return self.data[ enum ]

    def get_data(self):
        return self.data

class InterfaceData:
    """
    Class to hold interface data
    """
    def __init__(self, dSASA, dG, delta_unsats, sc_value, packstat, nres_int):
        self.data = defaultdict()
        self.data[ dSASA_e ] = dSASA
        self.data[ dG_e ] = dG
        self.data[ delta_unsats_e ] = delta_unsats
        self.data[ sc_value_e ] = sc_value
        self.data[ packstat_e ] = packstat
        self.data[ nres_int_e ] = nres_int

    def has_value(self, enum):
        return self.data.has_key(enum)

    def get_value(self, enum):
        return self.data[ enum ]

    def get_data(self):
        return self.data


def get_cdr_cluster_mean_sd(data, cdr, enum):
    total = 0
    n = 0
    all = []
    for i in data:
        d = data[ i ]
        if isinstance(d, MetaData): pass
        cluster_data = d.get_cdr_cluster_data(cdr)
        assert(isinstance(cluster_data, CDRClusterData))

        total = total + cluster_data.get_value(enum)
        all.append(cluster_data.get_value(enum))
        n+=1

    return (total/float(n)), numpy.std(numpy.array(all)), all

def get_cdr_physical_mean_sd(data, cdr, enum):
    total = 0
    n = 0
    all = []
    for i in data:
        d = data[ i ]
        assert(isinstance(d, MetaData))
        physical_data = d.get_cdr_physical_data(cdr)
        assert(isinstance(physical_data, CDRPhysicalData))

        total = total + physical_data.get_value(enum)
        all.append(physical_data.get_value(enum))
        n+=1

    return (total/float(n)), numpy.std(numpy.array(all)), all

def get_interface_mean_sd(data, enum):
    total = 0
    n = 0
    all = []
    for i in data:
        d = data[ i ]
        assert(isinstance(d, MetaData))
        interface_data = d.get_interface_data()
        assert(isinstance(interface_data, InterfaceData))

        total = total + interface_data.get_value(enum)
        all.append(interface_data.get_value(enum))
        n+=1

    return (total/float(n)), numpy.std(numpy.array(all)), all


def get_data_for_natives(db, struct_id_map):
    """
    Gets data for all the natives from the database.
    Returns a dictionary of struct_id: MetaData
        The MetaData can be queried for cluster, physical, and interface data.
    """
    cdrs = ["L1", "L2", "L3", "H1", "H2"]

    data = defaultdict()

    for struct_id in sorted(struct_id_map.keys()):
        struct_data = MetaData(struct_id)
        for cdr in cdrs:

            ## CDR Cluster Data
            c = db.cursor()
            stmt = """
                    SELECT
                        length,
                        fullcluster,
                        normDis_deg,
                        sequence
                    FROM
                        cdr_clusters
                    WHERE
                        CDR=? AND struct_id=?
                    """
            for row in c.execute(stmt, [cdr, struct_id]):
                print repr(row)
                cluster_data = CDRClusterData(cdr, int(row[0]), row[1], float(row[2]), row[3])
                struct_data.set_cdr_cluster_data(cdr, cluster_data)

            c.close()

            ## CDR Physical Data
            c = db.cursor()
            stmt = """
                    SELECT
                        ag_ab_contacts_total,
                        ag_ab_contacts_nres,
                        ag_ab_dSASA,
                        ag_ab_dG,
                        charge,
                        aromatic_nres
                    FROM
                        cdr_metrics
                    WHERE
                        CDR=? AND struct_id=?
                    """
            for row in c.execute(stmt, [cdr, struct_id]):
                physical_data = CDRPhysicalData(cdr, int(row[0]), int(row[1]), float(row[2]), float(row[3]), int(row[4]), int(row[5]))
                struct_data.set_cdr_physical_data(cdr, physical_data)

            c.close()

        ## Overall Data
        c = db.cursor()
        stmt = """
                SELECT
                    dSASA,
                    dG,
                    delta_unsatHbonds,
                    sc_value,
                    packstat,
                    nres_int
                FROM
                    interfaces
                WHERE
                    interface='LH_A' AND
                    struct_id=?
                """
        for row in c.execute(stmt, [struct_id]):
            interface_data = InterfaceData(float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), int(row[5]))
            struct_data.set_interface_data(interface_data)

        c.close()

        data[struct_id] = struct_data

    return data

def get_outlier_definition_string(outlier_definition, rmsd_cutoff = 1.5, dihdist_cutoff = 40):
    """
    Returns a string for adding to a database query which removes outliers.  Need to add AND manually to the string.
    """
    if outlier_definition == "conservative":
        s = " (bb_rmsd_cdr_align < "+repr(rmsd_cutoff)+" AND DistDegree < "+repr(dihdist_cutoff)+") AND DistDegree != -1 AND bb_rmsd_cdr_align != -1"
        return s
    elif outlier_definition == "liberal":
        s = " (bb_rmsd_cdr_align < "+repr(rmsd_cutoff)+" OR DistDegree < "+repr(dihdist_cutoff)+")  AND DistDegree != -1 AND bb_rmsd_cdr_align != -1"
        return s
    else:
        sys.exit("Could not understand outlier definition:"+outlier_definition)


def get_clusters_or_lengths(data, struct_id_map, cdr, enum):

    clusters = defaultdict()
    for i in data:
        d = data[ i ]
        cluster_data = d.get_cdr_cluster_data(cdr)
        l_or_c = cluster_data.get_value(enum)
        if not clusters.has_key(l_or_c):
            clusters[l_or_c] = []
        clusters[l_or_c].append(struct_id_map[ i ])

    #print repr(clusters)
    return clusters

def get_total_cluster_or_length_in_db(db, c_or_l, cdr, enum, include_outliers = False, length_max = 30):

    total = 0
    c = db.cursor()

    if enum == length_e:
        entry = "length"
    elif enum == cluster_e:
        entry = "fullcluster"

    if not include_outliers:
        for row in c.execute("SELECT length FROM cdr_data WHERE datatag!='loopKeyNotInPaper' AND CDR=? AND "+entry+"=? AND length <=? AND "+get_outlier_definition_string("liberal"), [cdr, c_or_l, length_max]):
            total+=1
    else:
        for row in c.execute("SELECT length FROM cdr_data WHERE datatag!='loopKeyNotInPaper' AND CDR=? AND "+entry+"=? AND length <=?", [cdr, c_or_l, length_max]):
            total+=1

    return total

def get_total_cdrs(db, cdr, include_outliers = False, length_max = 30):

    total = 0
    c = db.cursor()

    if not include_outliers:
        for row in c.execute("SELECT length FROM cdr_data WHERE datatag!='loopKeyNotInPaper' AND CDR=? AND length <=? AND "+get_outlier_definition_string("liberal"), [cdr, length_max]):
            total+=1
    else:
        for row in c.execute("SELECT length FROM cdr_data WHERE datatag!='loopKeyNotInPaper' AND CDR=? AND length <=?", [cdr, length_max]):
            total+=1

    return total


###############################################
def add_m_line_cdr(data, line, cdr, enums):

    for enum in enums:
        m, sd, all = get_cdr_physical_mean_sd(data, cdr, enum)
        line = line+"\t%.2f"%m

    return line

def add_sd_line_cdr(data, line, cdr, enums):
    for enum in enums:
        m, sd, all = get_cdr_physical_mean_sd(data, cdr, enum)
        line = line+"\t%.2f"%sd

    return line

def add_min_line_cdr(data, line, cdr, enums):
    for enum in enums:
        m, sd, all = get_cdr_physical_mean_sd(data, cdr, enum)
        line = line+"\t%.2f"%min(all)

    return line

def add_max_line_cdr(data, line, cdr, enums):
    for enum in enums:
        m, sd, all = get_cdr_physical_mean_sd(data, cdr, enum)
        line = line+"\t%.2f"%max(all)

    return line

def add_m_line_int(data, line, enums):
    for enum in enums:
        m, sd, all = get_interface_mean_sd(data, enum)
        line = line+"\t%.2f"%m
    return line

def add_sd_line_int(data, line, enums):
    for enum in enums:
        m, sd, all = get_interface_mean_sd(data, enum)
        line = line+"\t%.2f"%sd
    return line

def add_min_line_int(data, line, enums):
    for enum in enums:
        m, sd, all = get_interface_mean_sd(data, enum)
        line = line+"\t%.2f"%min(all)
    return line

def add_max_line_int(data, line, enums):
    for enum in enums:
        m, sd, all = get_interface_mean_sd(data, enum)
        line = line+"\t%.2f"%max(all)
    return line
##########################################


def get_struct_id_transform(db):
    """
    Get map of struct_id to PDB
    """
    data = defaultdict()

    c = db.cursor()

    for row in c.execute("SELECT struct_id, tag FROM structures"):
        data[int(row[0])] = str(row[1]).split("_")[0].lower();

    c.close()

    return data

if __name__ == "__main__":

    parser = OptionParser()
    args = sys.argv

    ########## DB Options ##########
    parser.add_option("--db", "-i",
                      help = "rel DB path for the native features for cluster recovery",
                      default ="databases/natives_baseline2.native.talaris2013.db3")

    parser.add_option("--ab_db", "-a",
                      help = "Ab Db used for benchmarking",
                      default = os.environ[ "ROSETTA3_DB" ]+"/sampling/antibodies/antibody_database_rosetta_design.db")

    ####### Output Options #########
    parser.add_option("--outpath", "-o",
                      help = "Output directory",
                      default = "datasets/analysis")

    parser.add_option("--outname", "-n",
                      help = "Name to use for the output analysis files",
                      default = "benchmark2_20")





    ####################################################################################################################
    (options, args) = parser.parse_args(args=args[1:])

    if not os.path.exists(options.db):
        sys.exit("Native DB path does not exist")

    if not os.path.exists(options.outpath):
        os.mkdir(options.outpath)

    print options.ab_db
    if not os.path.exists(options.ab_db):
        sys.exit("AbDb path does not exist")

    db = sqlite3.connect(options.db)

    cdrs = ["L1", "L2", "L3", "H1", "H2"]

    struct_id_map = get_struct_id_transform(db)
    print repr(struct_id_map)
    data = get_data_for_natives(db, struct_id_map)
    #print repr(data)

    ####################  Output Cluster Data ############################################################
    print "Outputting Cluster Data"
    OUTFILE = open(options.outpath+"/"+options.outname+"_cluster_data.txt", 'w')
    OUTFILE.write("#PDB\tL1\tL2\tL3\tH1\tH2\n")

    line = "AVG_DIS"
    for cdr in cdrs:

        m, sd, all = get_cdr_cluster_mean_sd(data, cdr, dis_e)
        line = line+"\t%.2f"%m
    OUTFILE.write(line+"\n")

    line = "SD_DIS"
    for cdr in cdrs:

        m, sd, all = get_cdr_cluster_mean_sd(data, cdr, dis_e)
        line = line+"\t%.2f"%sd

    OUTFILE.write(line+"\n")

    for i in sorted(struct_id_map.keys()):
        pdb = struct_id_map[ i ]
        meta_data = data[ i ]
        if isinstance(meta_data, MetaData):pass



        line = pdb
        for cdr in cdrs:
            line = line + "\t"+meta_data.get_cdr_cluster_data(cdr).get_value(cluster_e)
        OUTFILE.write(line +"\n")

        #line = pdb+"_dis"
        #for cdr in cdrs:
        #    line = line + "\t%.2f"%meta_data.get_cdr_cluster_data(cdr).get_value(dis_e)
        #OUTFILE.write(line+"\n")

    OUTFILE.close()

    ################### Output Each length and cluster with frequencies and PDBs ###########################
    rosetta_db = sqlite3.connect(options.ab_db)
    print "Outputing Consensus cluster and lengths"

    OUTFILE = open(options.outpath+"/"+options.outname+"_consensus_cluster.txt", 'w')
    OUTFILE.write("#cdr\tcluster\tpercent\tfreq\ttotal\tPDBs\n")

    for cdr in cdrs:

        if cdr == "H3":
            include_outliers = True
        else:
            include_outliers = False

        unique_clusters = get_clusters_or_lengths(data, struct_id_map, cdr, cluster_e)
        total_cdrs = get_total_cdrs(rosetta_db, cdr, include_outliers)
        print "Total " +cdr+" CDRs: "+repr(total_cdrs)
        for cluster in sorted(unique_clusters):

            line = cdr
            total_clus = get_total_cluster_or_length_in_db(rosetta_db, cluster, cdr, cluster_e, include_outliers)


            print "Total clus: "+repr(total_clus)

            line = line+"\t"+cluster
            line = line+"\t%.2f"%((float(total_clus)/float(total_cdrs))*100)
            line = line+"\t"+repr(total_clus)
            line = line+"\t"+repr(total_cdrs)
            line = line+"\t"+",".join(unique_clusters[cluster])

            OUTFILE.write(line+"\n")
    OUTFILE.close()

    OUTFILE = open(options.outpath+"/"+options.outname+"_consensus_length.txt", 'w')
    OUTFILE.write("#cdr\tlength\tpercent\tfreq\ttotal\tPDBs\n")

    for cdr in cdrs:

        if cdr == "H3":
            include_outliers = True
        else:
            include_outliers = False

        unique_lengths= get_clusters_or_lengths(data, struct_id_map, cdr, length_e)
        for length in sorted(unique_lengths):
            line = cdr
            total_cdrs = get_total_cdrs(rosetta_db, cdr, include_outliers)
            total_clus = get_total_cluster_or_length_in_db(rosetta_db, length, cdr, length_e, include_outliers)

            line = line+"\t"+repr(length)
            line = line+"\t%.2f"%((float(total_clus)/float(total_cdrs))*100)
            line = line+"\t"+repr(total_clus)
            line = line+"\t"+repr(total_cdrs)
            line = line+"\t"+",".join(unique_lengths[length])

            OUTFILE.write(line+"\n")
    OUTFILE.close()

    ################### Output physical stats of each CDR ##################################################
    print "Outputing Physical stats"
    OUTFILE = open(options.outpath+"/"+options.outname+"_cdr_phys_data.txt", 'w')
    OUTFILE.write("PDB_CDR\tdSASA\tdG\tcontacts\tcontacts_nres\tcharge\taromatics\n")

    #AVGS
    enums = [dSASA_e, dG_e, ag_contacts_e, ag_contacts_nres_e, charge_e, aromatics_e]

    for cdr in cdrs:
        line = "AVG_"+cdr
        line = add_m_line_cdr(data, line, cdr, enums)
        OUTFILE.write(line+"\n")

        line = "_SD_"+cdr
        line = add_sd_line_cdr(data, line, cdr, enums)
        OUTFILE.write(line+"\n")

        line = "MIN_"+cdr
        line = add_min_line_cdr(data, line, cdr, enums)
        OUTFILE.write(line+"\n")

        line = "MAX_"+cdr
        line = add_max_line_cdr(data, line, cdr, enums)
        OUTFILE.write(line+"\n\n")

    for i in sorted(struct_id_map.keys()):
        pdb = struct_id_map[ i ]
        meta_data = data[ i ]

        for cdr in cdrs:
            line = pdb+"_"+cdr
            phys_data = meta_data.get_cdr_physical_data(cdr)
            if isinstance(phys_data, CDRPhysicalData): pass

            line = line+"\t%.2f"%phys_data.get_value(dSASA_e)
            line = line+"\t%.2f"%phys_data.get_value(dG_e)
            line = line+"\t%.2f"%phys_data.get_value(ag_contacts_e)
            line = line+"\t%.2f"%phys_data.get_value(ag_contacts_nres_e)
            line = line+"\t%.2f"%phys_data.get_value(charge_e)
            line = line+"\t%.2f"%phys_data.get_value(aromatics_e)

            OUTFILE.write(line+"\n")
        OUTFILE.write("\n")
    OUTFILE.close()



    #################### Output physical stats of interface ##################################################
    print "Outputting interface stats"
    OUTFILE = open(options.outpath+"/"+options.outname+"_interface_data.txt", 'w')
    OUTFILE.write("PDB\tdSASA\tdG\tdelta_unsats\tnres_int\tsc_value\tpackstat\n")

    enums = [dSASA_e, dG_e, delta_unsats_e, nres_int_e, sc_value_e, packstat_e]
    line = "AVG"
    line = add_m_line_int(data, line, enums)
    OUTFILE.write(line+"\n")

    line = "_SD"
    line = add_sd_line_int(data, line, enums)
    OUTFILE.write(line+"\n")

    line = "MIN"
    line = add_min_line_int(data, line, enums)
    OUTFILE.write(line+"\n")

    line = "MAX"
    line = add_max_line_int(data, line, enums)
    OUTFILE.write(line+"\n\n")

    for i in sorted(struct_id_map.keys()):
        pdb = struct_id_map[ i ]
        meta_data = data[ i ]
        int_data = meta_data.get_interface_data()
        if isinstance(int_data, InterfaceData):pass

        line = pdb
        line = line+"\t%.2f"%int_data.get_value(dSASA_e)
        line = line+"\t%.2f"%int_data.get_value(dG_e)
        line = line+"\t%.2f"%int_data.get_value(delta_unsats_e)
        line = line+"\t%.2f"%int_data.get_value(nres_int_e)
        line = line+"\t%.2f"%int_data.get_value(sc_value_e)
        line = line+"\t%.2f"%int_data.get_value(packstat_e)
        OUTFILE.write(line+"\n")
    OUTFILE.close()

    print "Complete."

