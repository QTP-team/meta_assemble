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

