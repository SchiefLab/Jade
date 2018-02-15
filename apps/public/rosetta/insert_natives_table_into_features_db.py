#!/usr/bin/env python

from __future__ import print_function
import sqlite3, re
from argparse import ArgumentParser


def get_parser():
    parser = ArgumentParser(description="This script takes a PDBLIST of natives and then adds a new table to the database with "
                            "struct_id as proper foreign primary key and the native structure based solely on a search of the name tag. ")

    parser.add_argument("--pdblist", help = "PDBLIST of native structures used.")
    parser.add_argument("--db", help = "The database we are working on.")
    return parser

if __name__=="__main__":
    parser = get_parser()
    options = parser.parse_args()



    pdbs=[]
    pdblist = open(options.pdblist, 'r')
    for line in pdblist:
        line = line.strip()
        if not line: continue
        if line.startswith('#'): continue
        pdbs.append(line.split('.')[0])
    pdblist.close()

    testing = sqlite3.connect(options.db)

    cur = testing.cursor()
    cur.execute("select struct_id, tag from structures")
    data = cur.fetchall()

    cur.execute("DROP TABLE IF EXISTS natives")
    cur.execute("CREATE TABLE natives(struct_id INTEGER NOT NULL, "
                "tag TEXT, "
                "native TEXT, "
                "FOREIGN KEY (struct_id) REFERENCES structures(struct_id) DEFERRABLE INITIALLY DEFERRED,"
                "PRIMARY KEY (struct_id) )")

    final_data = []
    for item in data:
        for pdb in pdbs:
            if re.search(pdb, str(item[1])):
                new_item = []
                new_item.append(item[0])
                new_item.append(str(item[1]))
                new_item.append(pdb)
                final_data.append(new_item)
                continue

    for new_data in final_data:
        cur.execute("INSERT INTO natives VALUES(?,?,?)", (new_data[0], new_data[1], new_data[2]))
        testing.commit()

    cur.close()
    testing.close()

    print("Done")