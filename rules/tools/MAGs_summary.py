# tianliu@genomics.cn
# 2020-06-12
# MAGs_summary

import os
import argparse
import re

def summary_report(infiles, MAGs_per_sample, MAGs_quality):
    f_o1 = open(MAGs_per_sample, "w")
    f_o1.write("sample_id\tHQ_MAGs_n\tMQ_MAGs_n\tLQ_MAGs_n\n")

    f_o2 = open(MAGs_quality, "w")
    f_o2.write("sample_id\tbin_path\tquality\tcompleteness\tcontamination\tstrain_heterogeneity\n")

    for f in infiles:
        sample_id = f.split("/")[-1].split(".")[0]
        dset = open(f).readlines()
        item = re.findall(r'\d+',dset[0])
        all_MAGs = item[0]
        HQ_MAGs = item[-2]
        MQ_MAGs = item[-1]
        LQ_MAGs = int(all_MAGs) - int(HQ_MAGs) - int(MQ_MAGs)
        f_o1.write("%s\t%s\t%s\t%s\n" %(sample_id, HQ_MAGs, MQ_MAGs, LQ_MAGs))

        for line in dset[2:]:
            f_o2.write("%s\t%s" %(sample_id, line))
    f_o1.close()
    f_o2.close()

def main():
    parser = argparse.ArgumentParser(
            prog = 'MAGs_summary.py',
            usage = 'python MAGs_summary.py infile1 ... infileN -o MAGs_per_samples.txt -O picked_MAGs_quality.txt ',
            description = "MAGs summary for Liu's meta_assemble pipline"
            )
    parser.add_argument("infiles", metavar = "*.txt or 1.txt 2.txt 3.txt etc.", nargs = "+", help = "One or more files to join")
    parser.add_argument("-o", help = "report of MAGs per samples")
    parser.add_argument("-O", help = "report of picked MAGs's quality")

    args = parser.parse_args()
    summary_report(args.infiles, args.o, args.O)

main()

