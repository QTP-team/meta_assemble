#############################################
### tianliu@genomics.cn
### 2020-06-13
### pipline for stool metagenomics assemble
#############################################

### First, modify rules/assemble.config.yaml and check the bowtie index of host genome (default: hg38)
### step1.a run at test node
snakemake \
--snakefile rules/assemble.smk \
--configfile rules/assemble.config.yaml \
--core 32 2> smk.log

### step1.b qsub to clusters
### Modify rules/assemble.cluster.yaml before qsub, especially for queue and project_ID.
if [ ! -d 1.assay/cluster_logs ];then mkdir -p 1.assay/cluster_logs;fi

snakemake \
--snakefile rules/assemble.smk \
--configfile rules/assemble.config.yaml \
--cluster-config rules/assemble.cluster.yaml \
--jobs 100 \
--keep-going \
--rerun-incomplete \
--latency-wait 360 \
--cluster "qsub -S /bin/bash -cwd -q {cluster.queue} -P {cluster.project} -l vf={cluster.mem},p={cluster.cores} -binding linear:{cluster.cores} -o {cluster.output} -e {cluster.error}"

#### step2 select picked MAGs and rename them.
mkdir 2.result/all_picked_MAGs
echo -e "id\tfna" > MAGs.txt
for bins in `awk 'NR>1{print $2}' 2.result/step1_assembly/picked_MAGs_quality.txt`;do 
tmp=`dirname $bins`
id=`basename $tmp _binning`
bin_id=`basename $bins .fa`;
awk 'BEGIN{count=1}{if(/>/){$0 = ">""'$id'""_""'$bin_id'""_contig"count; count++ ; print $0}else{print $0}}' $bins > 2.result/all_picked_MAGs/${id}_${bin_id}.fa
echo -e "${id}_${bin_id}\t2.result/all_picked_MAGs/${id}_${bin_id}.fa" >> MAGs.txt
done
