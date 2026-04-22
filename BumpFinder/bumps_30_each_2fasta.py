#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Description: bed2fasta
Author: Li Fajin
Email: lifajin@westlake.edu.cn
Date: 2024-03-27 13:36:29
LastEditors: Li Fajin
LastEditTime: 2024-03-27 13:45:46
FilePath: \\ORFreferenced:\\Projects\\02. AA starvation\\GSE142822\\13.metagenePlot\\BumpFinder\\Trp_48h\\minus\\bumps_30_each_2fasta.py
'''

import sys

bed=sys.argv[1]
fasta=sys.argv[2]
def bed2fasta(bed,fasta):
    with open(bed,'r') as fin, open(fasta,'w') as fout:
        for line in fin.readlines():
            id=line.strip().split("\t")[0]
            bumpstart=line.strip().split("\t")[1]
            bumpemd=line.strip().split("\t")[2]
            bumpSeq=line.strip().split("\t")[3]
            fout.write("%s%s\n" %(">",id+":"+bumpstart+"-"+bumpemd))
            fout.write("%s\n" % (bumpSeq))
    print("Finish!",file=sys.stderr)
    
bed2fasta(bed,fasta)