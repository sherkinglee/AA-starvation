#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Description: 
Author: Li Fajin
Email: lifajin@westlake.edu.cn
Date: 2024-03-27 11:10:54
LastEditors: Li Fajin
LastEditTime: 2024-03-27 11:19:13
FilePath: \\ORFreferenced:\\Projects\\02. AA starvation\\GSE142822\\13.metagenePlot\\BumpFinder\\Trp_48h\\minus\\sequence_extractor.py
'''
import subprocess
import sys 
with open("bumps_30each.bed", "r") as fh,open("bumps_30_each_tmp.txt","w") as fout:
    file_content = fh.readlines()
    for line in file_content:
        ENST = ""
        start = ""
        width = 58
        end = ""
        if "\t" in line:
            ENST, start, end = line.strip().split("\t")
            out = subprocess.getoutput(f"grep -A 1 {ENST} all_amino_acid_sequences.fa")
            line_out = out.split("\n")

            aa_1 = line_out[1][int(start):int(start)+width]
    
            fout.write(ENST + "\t" + start + "\t" + end  + "\t" + aa_1)
            fout.write("\n")
