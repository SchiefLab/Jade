import sys
import os
import sqlite3
import pandas
from collections import defaultdict

from RAbD_BM.AnalysisInfo import AnalysisInfo
from RAbD_BM.AnalysisInfo import NativeInfo

#Tools for the features database.


def get_cdr_cluster_df(db_path):
    """
    Get a dataframe with typical cluster info in it
    :param db_con: sqlite3.con
    :rtype: pandas.DataFrame
    """
    query="SELECT "\
            "structures.input_tag," \
            "cdr_clusters.CDR," \
            "cdr_clusters.fullcluster," \
            "cdr_clusters.length," \
            "cdr_clusters.normDis_deg," \
            "cdr_clusters.sequence "\
        "FROM " \
            "cdr_clusters," \
            "structures "\
        "WHERE " \
            "cdr_clusters.struct_id = structures.struct_id "

    con = sqlite3.connect(db_path)
    df = pandas.read_sql_query(query, con)
    con.close()

    return df


def get_all_entries(df, pdbid, cdr):
    """
    Get all entries of a given PDBID and CDR.
    :param df: pandas.DataFrame
    :rtype: pandas.DataFrame
    """
    return df[(df['input_tag'].str.contains(pdbid)) & (df['CDR'] == cdr)]

def get_total_entries(df, pdbid, cdr):
    """
    Get the total number of entries of the particular CDR and PDBID in the database
    :param df: pandas.DataFrame
    :rtype: int
    """
    return len(get_all_entries(df, pdbid, cdr))


############### Data ###############

def get_cluster(df, pdbid, cdr):
    """
    Get the fullcluster from the dataframe for native or experimental data

    :param df: pandas.DataFrame
    :rtype: str
    """
    d = df[(df['input_tag'].str.contains(pdbid) ) & (df['CDR'] == cdr)]
    return str(d.iloc[0]['fullcluster'])

def get_length(df, pdbid, cdr):
    """
    Get the length from the dataframe for native or experimental data

    :param df: pandas.DataFrame
    :rtype: int
    """
    d = df[(df['input_tag'].str.contains(pdbid)) & (df['CDR'] == cdr)]
    return d.iloc[0]['length']


############### Matching ###############

def get_length_matches(df, pdbid, cdr, length):
    """
    Get a dataframe of the matching ("Recovered") rows (DataFrame).

    :param df: pandas.DataFrame
    :param length: int
    :rtype: pandas.DataFrame
    """
    return df[(df['input_tag'].str.contains(pdbid)) & (df['CDR'] == cdr) & (df['length'] == length )]

def get_cluster_matches(df, pdbid, cdr, cluster):
    """
    Get a dataframe of the matching ("Recovered") rows (DataFrame).

    :param df: pandas.DataFrame
    :rtype: pandas.DataFrame:
    """
    return df[(df['input_tag'].str.contains(pdbid)) & (df['CDR'] == cdr) & (df['fullcluster'] == cluster )]


############### Recovery ###############

def get_length_recovery(df, pdbid, cdr, length):
    """
    Get the number of matches in the df and pdbid to the cdr and length

    :param df: pandas.DataFrame
    :param length: int
    :rtype: int
    """
    return len(get_length_matches(df, pdbid, cdr, length)['length'])

def get_cluster_recovery(df, pdbid, cdr, cluster):
    """
    Get the number of matches in the df and pdbid to the cdr and cluster
    :param df: pandas.DataFrame
    :rtype: int
    """
    return len(get_cluster_matches(df, pdbid, cdr, cluster)['fullcluster'])
