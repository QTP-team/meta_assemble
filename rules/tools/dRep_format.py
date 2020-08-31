### to dRep's genomeInfo.csv
import os
import re
import sys

if len(sys.argv) != 3:
    print("usage: python dRep_format.py infile.tsv outfile.csv")
else:
    infile, outfile = sys.argv[1:]

dset = open(infile).readlines()[2:]

f_o = open(outfile, "w")
f_o.write("genome,completeness,contamination,length,N50\n")

for i in dset:
    item = i.strip().split("\t")
    genome = item[0].split("/")[-1]
    comp = item[2]
    cont = item[3]
    length = item[7]
    n50 = item[13]
    f_o.write("%s,%s,%s,%s,%s\n" %(genome, comp, cont, length, n50))
f_o.close()
