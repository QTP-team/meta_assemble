### megahit_summary
### tianliu@genomics.cn
### 2020-07-13

import argparse
import sys
import os
import re

parser = argparse.ArgumentParser(
    prog = "megahit_summary.py",
    description = "merge megahit logs",
    usage = "megahit_summary.py */log > megahit.merge.txt")

parser.add_argument(
    "infiles",
    metavar = "input.txt",
    nargs = "+",
    help = "One or more megahit log to join" )

def merge(infile_lst):
    print("sample_ID\tcontigs\ttotal\tmin\tmax\tavg\tN50")
    for f in infile_lst:
        sample_id = os.path.split(os.path.basename(f))[1].replace('.megahit.log','')
        with open(f) as f_i:
            item = f_i.readlines()
            res = re.findall('\d+', item[-2].split(" - ")[-1])
            print("%s\t%s\t%s" %(sample_id, "\t".join(res[:-2]), res[-1]))

def _main():
    args = parser.parse_args()
    merge(args.infiles)

if __name__ == "__main__":
    _main()
