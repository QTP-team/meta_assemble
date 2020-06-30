###################################
# filter summary
# tianliu@genomics.cn 2019-11-26
###################################

import os
import argparse
import json
import pandas

def trim_stat(fastp_json):
    fastp_sample = os.path.basename(fastp_json).replace(".fastp.json","")
    with open(fastp_json, 'r') as f:
        OA_json = json.load(f)
    reads_before = OA_json.get('summary').get('before_filtering').get('total_reads')
    bases_before = OA_json.get('summary').get('before_filtering').get('total_bases')
    read_len = OA_json.get('summary').get('before_filtering').get('read1_mean_length') #PE only
    q30_before = float(OA_json.get('summary').get('before_filtering').get('q30_rate'))*100
    duplication = float(OA_json.get('duplication').get('rate'))*100
    res_before = "%s\t%s\t%s\t%s\t%.2f%%\t%.2f%%\t" %(fastp_sample, reads_before, bases_before, read_len, q30_before, duplication)

    reads_pass = OA_json.get('summary').get('after_filtering').get('total_reads')
    bases_pass = OA_json.get('summary').get('after_filtering').get('total_bases')
    q30_after = float(OA_json.get('summary').get('after_filtering').get('q30_rate')) * 100
    reads_pass_ratio = reads_pass/reads_before * 100
    bases_pass_ratio = bases_pass/bases_before * 100
    res_after = "%s\t%.2f%%\t%.2f%%\t%s\t%.2f%%" %(reads_pass, reads_pass_ratio, q30_after, bases_pass, bases_pass_ratio)

    fastp_res = res_before + res_after
    print(fastp_res)

def rmhost_stat(seqkit_file):
    dset = open(seqkit_file).readlines()
    r1 = dset[1].strip().split()
    r2 = dset[2].strip().split()
    rmhost_sample = os.path.basename(r1[0]).split(".rmhost.1.fq.gz")[0]
    num_seqs = int(r1[3].replace(",","")) + int(r2[3].replace(",",""))
    sum_len = int(r1[4].replace(",","")) + int(r2[4].replace(",",""))
    print("%s\t%s\t%s" %(rmhost_sample, num_seqs, sum_len))

def main():
    parser = argparse.ArgumentParser(
            prog = 'filter_summary.py',
            usage = 'python trim_summary.py [options] infile1 [infile2 ...] > summary.txt',
            description = "trimming/rmhost summary. Only used in liu's meta pipeline"
            )
    parser.add_argument("-t", "--trim", help = "trimming summary", action='store_true')
    parser.add_argument("-r", "--rmhost", help = "rmhost summary", action='store_true')
    parser.add_argument("infiles", 
        metavar = "*.json or 1.json 2.json 3.json etc.", 
        nargs = "+", 
        help = "One or more fastp.json to join")

    args = parser.parse_args()

    if args.trim:
        print("Library_ID\tRaw_reads_count\tRaw_base_count\tRead_length\tRaw_Q30_ratio\tDuplication\tClean_reads_count\tClean_reads_ratio\tClean_Q30_ratio\tClean_bases_count\tClean_bases_ratio")
        for f in args.infiles:
            trim_stat(f)

    if args.rmhost:
        print("Library_ID\tNohomo_reads_count\tNohomo_bases_count")
        for f in args.infiles:
            rmhost_stat(f)

main()
