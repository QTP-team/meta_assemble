import os
import sys
import random

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir

done_dict = {}
for i in open('00.work/qsub_map.sh'):
    item = i.strip().split(" ")
    index_id = item[4].split('/')[1]
    rmhost_id = item[-1].split('/')[-1].replace('.sorted.bam', '')
    if index_id not in done_dict.keys():
        done_dict[index_id] = {rmhost_id}
    else:
        done_dict[index_id].add(rmhost_id)

all_samples_lst = done_dict.keys()
wk_dir = '/ldfssz1/ST_META/P19Z10200N0314_LXP/tianliu/3.co-binning/Horse'
rmhost_dir = wk_dir + '/02.rmhost'

for sample_id in all_samples_lst:
    tmp_lst = all_samples_lst - done_dict[sample_id]
    random_sample_lst = random.sample(tmp_lst, 10)
    run_dir = mkdir(os.path.join(wk_dir, '06.s20_binning', sample_id))
    run_log_dir = mkdir(os.path.join(wk_dir, '06.s20_binning', sample_id, 'logs'))
    index_prefix = '/ldfssz1/ST_META/P19Z10200N0314_LXP/tianliu/3.co-binning/Horse/04.s10_binning/%s/1.index/%s_bw2_index' %(sample_id, sample_id)

    for i in random_sample_lst:
        r1 = os.path.join(rmhost_dir, i + ".rmhost.1.fq.gz")
        r2 = os.path.join(rmhost_dir, i + ".rmhost.2.fq.gz")
        map_log = os.path.join(run_log_dir, i+".map2scaftigs.log")
        tmp_prefix = os.path.join(run_dir, "tmp" + i)
        tmp_bam = os.path.join(run_dir, i+".sorted.bam")
        depth = os.path.join(run_dir, i+".depth.txt")
        print("bowtie2 -p 4 -x %s -1 %s -2 %s 2> %s | samtools sort -T %s -@ 4 -O BAM -o %s - && jgi_summarize_bam_contig_depths %s --outputDepth %s && rm %s" %(index_prefix, r1, r2, map_log, tmp_prefix, tmp_bam, tmp_bam, depth, tmp_bam))
