### to dRep's genomeInfo.csv
import os
import re
import sys

if len(sys.argv) != 3:
    print("usage: python dRep_format.py infile.tsv outfile.csv")
else:
    infile, outfile = sys.argv[1:]

dset = open(infile).readlines()[1:]

f_o = open(outfile, "w")
f_o.write("genome,completeness,contamination,length,N50\n")

for i in dset:
    item = i.strip().split("\t")
    genome = item[1].split("/")[-1]
    comp = item[3]
    cont = item[4]
    length = item[8]
    n50 = item[14]
    f_o.write("%s,%s,%s,%s,%s\n" %(genome, comp, cont, length, n50))
f_o.close()
