#!/usr/bin/env python
### merge bowtie2 log
### tianliu@genomics.cn
### 20191125 v1.0
### 20200707 v1.1 add PE mapping_ratio

import argparse
import sys
import os

parser = argparse.ArgumentParser( 
    prog = "merge_bowtie2_log.py",
    description = "merge bowite2 logs (PE reads)",
    usage = "merge_bowite2_log.py *.bowti2.log > merge.bowtie2.log")

parser.add_argument( 
    "infiles", 
    metavar = "input.txt", 
    nargs = "+", 
    help = "One or more bowtie2 log to join" )

def merge(infile_lst):
    print("sample\treads\tPE_mapping_ratio\tmapping_ratio")
    for f in infile_lst:
        sample_id = os.path.split(os.path.basename(f))[1].replace('.map2scaftigs.log','')
        with open(f) as f_i:
            item = f_i.readlines()
            reads = float(item[0].strip().split(" ")[0])*2
            pe_mapped = float(item[3].strip().split(" ")[0])*2 #pairs aligned concordantly exactly 1 time
            pe_mapped += float(item[4].strip().split(" ")[0])*2 #pairs aligned concordantly >1 times
            pe_mapped += float(item[7].strip().split(" ")[0])*2 # pairs aligned discordantly 1 time
            mapped = pe_mapped + float(item[12].strip().split(" ")[0]) #aligned exactly 1 time
            mapped += float(item[13].strip().split(" ")[0]) #aligned >1 times
            mapping_ratio = mapped / reads
            pe_mapping_ratio = pe_mapped / reads 
        print("%s\t%s\t%.4f\t%.4f" %(sample_id, int(reads), pe_mapping_ratio, mapping_ratio))

def _main():
    args = parser.parse_args()
    merge(args.infiles)

if __name__ == "__main__":
    _main()
