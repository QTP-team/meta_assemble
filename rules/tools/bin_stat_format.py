### tianliu@genomics.cn
### 2020-07-12
### format checkm bins_stats.analyze.tsv
import sys
import json
if len(sys.argv) != 2:
    print("Usage: python bin_stat_format.py bin_stats.analyze.tsv > demo.bins.stat.txt")
else:
    infile = sys.argv[1]

dset = open(infile).readlines()
title_dict = eval(dset[0].strip().split("\t")[1])

print("bins_id\t%s" %("\t".join(title_dict.keys())))

for rec in dset:
    bins_id, stat = rec.strip().split("\t")
    stat_lst = [str(x) for x in eval(stat).values()]
    print("%s\t%s" %(bins_id,'\t'.join(stat_lst)))

