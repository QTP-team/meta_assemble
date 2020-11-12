### s10 metabat2
### tianliu@genomics.cn
### 2020-09-20

import os
import random
import sys
import re
import yaml

if len(sys.argv) != 3:
  print("Usage: python mulit_samples_metabat2.py sample.txt config.yaml")
  exit()
else:
  samplefile, configfile = sys.argv[1:]

config = yaml.safe_load(open(configfile))
all_samples_lst = [x.split("\t")[0] for x in open(samplefile).readlines()[1:]]

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir

def map2scaftigs(sample_id, index_prefix, random_sample_lst, rmhost_dir, run_dir, run_log_dir, metabat2_cpu):
    for i in random_sample_lst:
        r1 = os.path.join(rmhost_dir, i+".rmhost.1.fq.gz")
        r2 = os.path.join(rmhost_dir, i+".rmhost.2.fq.gz")
        map_log = os.path.join(run_log_dir, i+".map2scaftigs.log")
        tmp_prefix = os.path.join(run_dir, "tmp" + i)
        tmp_bam = os.path.join(run_dir, i+".sorted.bam")
        depth = os.path.join(run_dir, i+".depth.txt")
        print("bowtie2 -p %s -x %s -1 %s -2 %s 2> %s | samtools sort -T %s -@ %s -O BAM -o %s - && jgi_summarize_bam_contig_depths %s --outputDepth %s && rm %s" %(metabat2_cpu, index_prefix, r1, r2, map_log, tmp_prefix, metabat2_cpu, tmp_bam, tmp_bam, depth, tmp_bam))

def metabat2(sample_id, sample_n, run_dir, run_log_dir, contig_path, metabat2_cpu):
    merge_depth = os.path.join(run_dir, sample_id + "_merge.depths")
    merge_log = os.path.join(run_log_dir, "merge_depths.log")
    bin_basename = os.path.join(run_dir, sample_id+"_binning", sample_id + "_bin")
    bin_log = os.path.join(run_log_dir, sample_id + ".metabat2.log")
    work_sh = "if [ `ls %s/*.depth.txt | wc -l` -eq %s ];then \
                python rules/tools/merge_depth.py %s/*.depth.txt -o %s; else echo lack depth files;fi > %s\n" %(run_dir, sample_n, run_dir, merge_depth, merge_log)
    work_sh += "metabat2 -i %s -a %s -o %s -m 1500 -t %s > %s" %(contig_path, merge_depth, bin_basename, metabat2_cpu, bin_log)
    print(work_sh)

def checkm(sample_id, checkm_run_dir,picked_dir, run_dir, logs_checkm_dir, checkm_cpu, high_quality, medium_quality):
    bins_dir = os.path.join(run_dir, sample_id+"_binning")
    bin_stat = os.path.join(run_dir, sample_id+".bins.stat.txt")
    checkm_log = os.path.join(logs_checkm_dir, sample_id+"checkm.log")
    pick_log = os.path.join(picked_dir, sample_id+'.picked.summary.txt')
    work_sh = "checkm lineage_wf -t %s -x fa %s %s 2> %s | grep -v INFO > %s/checkm_summary.txt\n" %(checkm_cpu, bins_dir, checkm_run_dir, checkm_log, checkm_run_dir)
    work_sh += 'python rules/tools/bin_stat_format.py %s/storage/bin_stats.analyze.tsv > %s\n' %(checkm_run_dir, bin_stat)
    work_sh += 'python rules/tools/pick_MAGs.py --high %s --medium %s %s %s/checkm_summary.txt > %s && rm -r %s/bins %s/storage' %(high_quality, medium_quality, bins_dir, checkm_run_dir, pick_log, checkm_run_dir, checkm_run_dir)
    print(work_sh)

def main():
    sample_n = config['params']['metabat2']['multi_samples']
    metabat2_cpu = config['params']['metabat2']['threads']
    checkm_cpu = config['params']['checkm']['threads']
    high_quality = config['params']['pick']['HQ']
    medium_quality = config['params']['pick']['MQ']

    rmhost_dir = config['assay']['rmhost']
    megahit_dir = config['assay']['megahit']
    metabat2_dir = mkdir(config['assay']['metabat2']+"_s"+str(sample_n))
    picked_dir = mkdir(config['assay']['picked']+"_s"+str(sample_n))
    logs_dir = mkdir(config['logs']['metabat2']+"_s"+str(sample_n))
    logs_checkm_dir = mkdir(config['logs']['checkm']+"_s"+str(sample_n))
    result_dir = mkdir(os.path.dirname(config['results']['binning_s1'])+'/binning_s'+str(sample_n))
    #multi_samples binning and checkm
    for sample_id in all_samples_lst:
        wd = os.getcwd()
        run_dir = mkdir(os.path.join(metabat2_dir, sample_id))
        run_log_dir = mkdir(os.path.join(logs_dir, sample_id))
        tmp_lst = [] + all_samples_lst
        tmp_lst.remove(sample_id)
        random_sample_lst = random.sample(tmp_lst, int(sample_n) - 1)
        contig_path = os.path.join(megahit_dir, sample_id+".megahit_out", sample_id + ".megahit.contigs.fa.gz")
        index_prefix = os.path.join(megahit_dir, sample_id+".megahit_out", sample_id+"_index", "g_index")

        print("ln -s %s/%s/%s/%s.depth.txt %s" %(wd, config['assay']['metabat2'], sample_id, sample_id, run_dir))
        map2scaftigs(sample_id, index_prefix, random_sample_lst, rmhost_dir, run_dir, run_log_dir, metabat2_cpu)
        metabat2(sample_id, sample_n, run_dir, run_log_dir, contig_path, metabat2_cpu)

        checkm_run_dir = mkdir(os.path.join(config['assay']['checkm']+"_s"+str(sample_n), sample_id))
        checkm(sample_id, checkm_run_dir, picked_dir, run_dir, logs_checkm_dir, checkm_cpu, high_quality, medium_quality)

    # summary
    all_bins_stat = "%s/*/*.bins.stat.txt" %(metabat2_dir)
    bins_summary = "%s/All_bins_stat.txt" %(result_dir)
    checkm_stat = "%s/*.picked.summary.txt" %(picked_dir)
    print('''cat %s | awk -F'\t' 'NR==1 || $2!="GC"' > %s''' %(all_bins_stat, bins_summary))
    print('python rules/tools/MAGs_summary.py %s -a %s -o %s/MAGs_per_sample.txt -O %s/picked_MAGs_quality.txt' %(checkm_stat, bins_summary, result_dir, result_dir))

main()
