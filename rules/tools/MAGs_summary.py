# tianliu@genomics.cn
# 2020-06-12
# MAGs_summary

import os
import argparse
import re

def summary_report(infiles, all_MAGs_stat, MAGs_per_sample, MAGs_quality):
    MAG_stat_dict = {}
    dset = open(all_MAGs_stat).readlines()
    for rec in dset[1:]:
        item = rec.strip().split("\t")
        MAG_id = item[0]+'.fa'
        MAG_stat = '\t'.join(item[1:])
        MAG_stat_dict[MAG_id] = MAG_stat

    f_o1 = open(MAGs_per_sample, "w")
    f_o1.write("Sample_ID\tHQ_MAGs_n\tMQ_MAGs_n\tLQ_MAGs_n\n")

    f_o2 = open(MAGs_quality, "w")
    dset_title = '\t'.join(dset[0].strip().split("\t")[1:])
    f_o2.write("Sample_ID\tbin_path\tquality\tcompleteness\tcontamination\tstrain_heterogeneity\t%s\n" %(dset_title))

    for f in infiles:
        sample_id = f.split("/")[-1].split(".")[0]
        picked_dset = open(f).readlines()
        item = re.findall(r'\d+',picked_dset[0])
        all_MAGs = item[0]
        HQ_MAGs = item[-2]
        MQ_MAGs = item[-1]
        LQ_MAGs = int(all_MAGs) - int(HQ_MAGs) - int(MQ_MAGs)
        f_o1.write("%s\t%s\t%s\t%s\n" %(sample_id, HQ_MAGs, MQ_MAGs, LQ_MAGs))

        for line in picked_dset[2:]:
            MAG_id = line.strip().split("\t")[0].split("/")[-1]
            MAG_stat_res = MAG_stat_dict[MAG_id]
            f_o2.write("%s\t%s\t%s\n" %(sample_id, line.strip(), MAG_stat_res))
    f_o1.close()
    f_o2.close()

def main():
    parser = argparse.ArgumentParser(
            prog = 'MAGs_summary.py',
            usage = 'python MAGs_summary.py infile1 ... infileN -a All_bins_stat.txt -o MAGs_per_samples.txt -O picked_MAGs_quality.txt ',
            description = "MAGs summary for Liu's meta_assemble pipline"
            )
    parser.add_argument("infiles", metavar = "*.txt or 1.txt 2.txt 3.txt etc.", nargs = "+", help = "One or more files to join")
    parser.add_argument("-a", help = "MAGs' assembly summary.")
    parser.add_argument("-o", help = "report of MAGs per samples")
    parser.add_argument("-O", help = "report of picked MAGs's quality")

    args = parser.parse_args()
    summary_report(args.infiles, args.a, args.o, args.O)

main()
