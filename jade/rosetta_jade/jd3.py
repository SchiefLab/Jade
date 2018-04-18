from __future__ import print_function
import os,re
from collections import defaultdict

def create_substituted_jd_string(input_jd, substitutions):
    """
    Create a substituted Job Definition for a JD3 job.  Creates file in /jd3 as _substituted file 
    
    :param input_jd: file path
    :param substitutions: list of defaultdicts for each job. Don't forget %%name%% if you are using it that way.
    :return: 
    """

    new_lines = []
    unparsed = open(input_jd, 'r').readlines()

    new_lines.append("<JobDefinitionFile>\n")

    for sub_set in substitutions:
        #print(sub_set)
        for line in unparsed:

            if re.search("JobDefinitionFile", line): continue
            #print(line)
            new_line = line.format(**sub_set)
            new_lines.append(new_line)

    new_lines.append("</JobDefinitionFile>\n")
    print("Using input:",input_jd)
    x = os.path.basename(input_jd).split('.')
    print(x)
    outname = ".".join(x[0:-1]) + "_substituted.xml"

    if not os.path.exists('job_definitions'):
        os.mkdir('job_definitions')

    OUT = open('job_definitions/'+outname, 'w')
    for line in new_lines:
        OUT.write(line)
    OUT.close()

    print("Done.\nSubstituted file writted to: " + "job_definitions/" + outname)

    return new_lines

