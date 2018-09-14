#!/usr/bin/env python
import os,sys

def insert (source_str, insert_str, pos):
    return source_str[:pos]+insert_str+source_str[pos:]

if __name__ == "__main__":

    if len(sys.argv) == 1:
        sys.exit("Fixes symmetry perl script issue with chains that are spaces.  First argument is the pdb file")
    input_pdb_file = sys.argv[1]

    lines = open(input_pdb_file, 'r').readlines()
    new_lines = []
    for line in lines:
        #print line
        new_line = line
        if len(line) > 4:
            print line[0:6]
            if line[0:6] == "HETATM":
                chain_column = 22-1
                chain = line[chain_column]
                print chain
                if chain == ' ' or len(line.strip()) == 77:
                    new_line = insert(new_line, " ", chain_column)


        new_lines.append(new_line)


    OUTFILE = open("test_pdb.pdb", 'w')
    for line in new_lines:
        OUTFILE.write(line)
    OUTFILE.close()
