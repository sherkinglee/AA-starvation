#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2021-01-07 14:14:09
LastEditors: Li Fajin
LastEditTime: 2021-01-07 19:31:39
Description: file content
'''

import os
import sys
# RNA_seq_summary_statistics.txt

def mergelogs(inputDir,output):
    files=os.listdir(inputDir)
    cutadaptLogs=[os.path.join(inputDir,filename) for filename in files if "Cutadapt" in filename]
    filteringLogs=[os.path.join(inputDir,filename) for filename in files if "Filtering" in filename]
    rRNAContamLogs=[os.path.join(inputDir,filename) for filename in files if "RRNA_contam" in filename]
    starMappingLogs=[os.path.join(inputDir,filename) for filename in files if "Star_mapping" in  filename]
    statisticsLogs=[os.path.join(inputDir,filename) for filename in files if "Statistics" in  filename]
    ## cutadapt
    os.system("echo -e '# cutadapt\nsample\tTotal\tTrimmed(Percent)\tshortNum(Percentage)\tLeftNum(Percentage)' >> "+output)
    for filename in cutadaptLogs:
        os.system("cat "+filename+" >>" +output)
        os.system("rm "+filename)

    ## filtering
    os.system("echo -e '# filtering\nsample\tinputNum\tRemained\tdiscarded(Percent)' >> "+output)
    for filename in filteringLogs:
        os.system("cat "+filename+" >>" +output)
        os.system("rm "+filename)

    ## rRNA contam
    os.system("echo -e '# remove RNA contamination\nsample\tProcessedNum\trRNA(Percent)\tnoContamRNA(Percent)' >> "+output)
    for filename in rRNAContamLogs:
        os.system("cat "+filename+" >>" +output)
        os.system("rm "+filename)

    ## mapping
    os.system("echo -e '# Star mapping\nsample\tinput\tUniquelyMapped(Percent)\tMutipulMapped(Percent)' >> "+output)
    for filename in starMappingLogs:
        os.system("cat "+filename+" >>" +output)
        os.system("rm "+filename)

    ## statistics
    os.system("echo -e '# DNA contamination\nsample\tExon\tDNA\tIntron\tambiguous_RNA' >> "+output)
    for filename in statisticsLogs:
        os.system("cat "+filename+" >>" +output)
        os.system("rm "+filename)


def main():
    inputdir=sys.argv[1]
    output=sys.argv[2]
    print("Start...")
    mergelogs(inputdir,output)
    print("Stop!")

if __name__=="__main__":
    main()