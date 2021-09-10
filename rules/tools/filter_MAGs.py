import argparse
import pandas as pd

def load_checkm(checkm_file):
    checkm_res = open(checkm_file).readlines()
    MAG_num = len(checkm_res) - 4
    for i in checkm_res[3:-1]:
        item = i.strip().split()
        bin_id, cp, ct, sh = item[0], float(item[-3]), float(item[-2]), float(item[-1])
        record = "\t".join([item[0], item[-3], item[-2], item[-1]])
    #od
    checkm_df
    print("#%s MAGs in %s" %(MAG_num, checkm_file)
    return checkm_df

def load_stat():



def main()
    parser = argparse.ArgumentParser(
        usage = 'python rules/tools/filter_MAGs.py bins_dir checkm_summary.txt bin_stats.analyze.tsv --min_comp 50 --min_cont 10 --min_qs 50 --out_all MAGs_all.summary.txt --out_filtered MAGs_filtered.summary.txt')

    parser.add_argument("MAGs_dir", help = "bins dir")
    parser.add_argument("checkm_res", help = "checkm_summary.txt from checkm")
    parser.add_argument("bin_stat", help = "bin_stats.analyze.tsv from checkm")
    parser.add_argument("--min_comp", help = "Minimal completeness from checkm", default = 50, type = float)
    parser.add_argument("--min_cont", help = "Minimal contamination from checkm", default = 10, type = float)
    parser.add_argument("--min_strh", help = "Minimal strain heterogeneity from checkm", default = 0, type = float)
    parser.add_argument("--min_len", help = "Minimal length of genome", default=0, type = float)
    parser.add_argument("--min_N50", help = "Minimal N50 of genome", default=0, type = float)
    parser.add_argument("--out_all", help = "summary for all bins")
    parser.add_argument("--out_filtered", help = "summary for filtered bins")

    args = parser.parse_args()
    checkm_df
