#########################################
### tianliu@genomics.cn
### 2020-06-05
##########################################

import re
import os
import argparse

def filter_checkm(checkm_file, high_standard, medium_standard):
    h_cp, h_ct, h_sh = [float(x) for x in high_standard.split(",")]
    m_cp, m_ct, m_sh = [float(x) for x in medium_standard.split(",")]

    if h_sh < 0:
        h_sh = float('inf')
    if m_sh < 0:
        m_sh = float('inf')

    HQ_set = set()
    tmp_set = set()

    checkm_res = open(checkm_file).readlines()
    MAG_num = len(checkm_res) - 4

    for i in checkm_res[3:-1]:
        item = i.strip().split()
        bin_id, cp, ct, sh = item[0], float(item[-3]), float(item[-2]), float(item[-1])
        record = "\t".join([item[0], item[-3], item[-2], item[-1]])

        if cp > h_cp and ct < h_ct and sh < h_sh:
            HQ_set.add(record)
        if cp > m_cp and ct < m_ct and sh < m_sh:
            tmp_set.add(record)

    MQ_set = tmp_set - HQ_set
    print("#%s MAGs in %s, HQ_MAG.n = %s, MQ_MAG.n = %s" %(MAG_num, checkm_file, len(HQ_set), len(MQ_set)))
    return HQ_set, MQ_set

def summary(MAGs_dir, HQ_set, MQ_set):
    wd = os.getcwd()
    print("bin_id\tquality\tcompleteness\tcontamination\tstrain_heterogeneity")
    for record in HQ_set:
        MAG = record.split("\t")
        MAG_path = os.path.join(wd+"/"+MAGs_dir, MAG[0] + ".fa")
        print("%s\tHigh\t%s" %(MAG_path, "\t".join(MAG[1:])))

    for record in MQ_set:
        MAG = record.split("\t")
        MAG_path = os.path.join(wd+"/"+MAGs_dir, MAG[0] + ".fa")
        print("%s\tMedium\t%s" %(MAG_path, "\t".join(MAG[1:])))

def main():
    parser = argparse.ArgumentParser(
        usage = 'python pick_MAGs.py --high 90,5,-1 --medium 50,10,-1 MAGs_dir checkm.txt > MAG.quality.summary.txt',
        description = 'Select high quality(HQ) and medium quality(MQ) MAGs.\
        The three numbers after "--high/medium" are completeness, contamination, strain heterogeneity from checkm result(< 0 will be ignored).'
        )
    parser.add_argument("MAGs_dir", help = "directory of MAG.fa")
    parser.add_argument("checkm_res", help = "checkm_result.txt")
    parser.add_argument("--high", help = "standard of high quality MAGs", default="90,5,-1")
    parser.add_argument("--medium", help = "standard of medium quality MAGs", default="50,10,-1")

    args = parser.parse_args()

    HQ_set, MQ_set = filter_checkm(args.checkm_res, args.high, args.medium)
    
    summary(args.MAGs_dir, HQ_set, MQ_set)

main()