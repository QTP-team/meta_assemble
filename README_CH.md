## 1. 简介
适用于粪便样本的大规模宏基因组组装流程。整个流程包括：fastp和Bowtie 2过滤低质量的reads和宿主的reads，MEGAHIT组装，MetaBAT2分箱，CheckM评估基因组质量。

本流程使用单样本组装，多样本对单个样本的contigs分箱的策略获得原核生物基因组。

该策略有如下优点
1. 单个样本组装可以避免多个样本中菌株的混杂，造成contig的碎裂以及基因组的菌株异质性高，组装出质量更高的基因组(图1)，同时降低计算资源的消耗。
![图1](pic/fig1.png)
<center>图1 菌株混杂会导致组装质量下降。</center>
<center>来源: https://drep.readthedocs.io/en/latest/overview.html#genome-de-replication </center>

2. 相较于常用的单个样本分箱的策略，使用多个样本对单个样本的contigs进行分箱得到基因组，可以充分利用contig在多个样本中的丰度以及相关性的信息，从而得到更多更好的基因组。

我们使用了马(n = 79)的数据，分别使用了1，4，6，8，10，15，20个样本(简称为s1-s20组)评估了样本数量对分箱结果的影响。随着使用的样本数越多，分箱的效果越好，并且在s20时仍未达到饱和。相较于s1组，s20组获得的总的bins数，中高质量的MAGs数，高质量的MAGs数分别增加了89%，115%，134% (Fig2a)，且绝大多数的样本都能受益于多样本分箱(Fig2b)。对所有基因组使用dRep聚类得到物种的代表基因组，使用的样本数越多，获得的物种也越多。相较于s1组，s20组的物种数增加了60%，有高质量基因组的物种数也增加了83%(Fig2c)。并且使用更多的样本能在保留大部分s1的物种的情况下，获得更多的新物种(Fig2d)。

![fig2](pic/Fig2.png)

<center>图2 更多的样本可以更好的改善分箱的质量</center>

metaSPAdes(不包含在本流程中)相较于MEGAHIT组装效果更好，但会消耗大量的计算资源。若去完宿主后的reads测序量小于10G，可以考虑选择metaSPAdes进行组装。MetaBAT2也可以更换为VAMB。

本流程不支持更改软件，若想换用其他软件进行组装，推荐使用[metapi](https://github.com/ohmeta/metapi)。

若使用了本流程，请引用文献：
[Todo]

### 2 Install
使用conda安装：
```shell
conda env create -n meta_assembly -f rules/environment.yaml
conda active meta_assembly
```

下载宿主的基因组，并构建bowtie2索引：
```shell
# GRCh38.p13 (human)
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39
mkdir -p 0.data/host_index && cd 0.data/host_index
curl -O https://ftp.cngb.org/pub/Assembly/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
bowtie2-build --threads 4 GCA_000001405.28_GRCh38.p13_genomic.fna.gz hg38.p13
```

## 3. Usage

### 3.1 测试数据集
使用三个单菌的基因组生成模拟的测试集。
```shell
sh demo/work.iss.sh
```

### 3.2 输入

样本输入格式见sample.txt，以制表符分割。仅支持双端测序样本。

```shell
$ cat sample.txt
id	fq1	fq2
test1	demo/test1_data_R1.fastq.gz	demo/test1_data_R2.fastq.gz
test2	demo/test2_data_R1.fastq.gz	demo/test2_data_R2.fastq.gz
test3	demo/test3_data_R1.fastq.gz	demo/test3_data_R2.fastq.gz
```

### 3.3 参数设置

软件运行参数和流程的文件结构都在```rules/assemble.config.yaml```中更改。
更新bowtie2_index的路径为"/your_path/0.data/host_index/hg38.p13"。

```yaml
params:
    ### step1 assembly
    fastp:
      min_len : 70 # Recommended value >= 30
    #rmhost
    rmhost:
      bowtie2_index : "/path/to/host_bowtie2_index"
      threads: 8
    #megahit
    megahit:
      threads: 16
      min_contigs_len: 200
    #binning
    metabat2:
      multi_samples: 3 # Recommended 10
      threads: 16
      minContig: 1500 #should be >=1500
    #quality
    checkm:
      threads: 8
    pick:
    	# "completeness, contamination, strain_heterogeneity" from CheckM
    	# Default : MIMAG standard.
    	# <0 : no filter
      HQ: "90,5,-1"
      MQ: "50,10,-1"

...
```

### 3.4a. 本地运行

先运行单个样本分箱。

```shell
snakemake \
--snakefile rules/assemble.smk \
--configfile rules/assemble.config.yaml \
--core 1 \
--dry-run
```

若测试正常，则使用--core指定CPU资源在本地运行。

```shell
snakemake \
--snakefile rules/assemble.smk \
--configfile rules/assemble.config.yaml \
--core 32 2> smk.log
```

### 3.4b. 投递任务(qsub环境)

投递任务相关参数在```rules/cluster.yaml```，可在该文件中指定每一步所需的资源。投任务前，需**修改项目编号，集群队列**。

```yaml
localrules: all

__default__:
  queue: "st.q"
  project: "your_project_ID"
  workdir: "./"
  mem: "1G"
  cores: 1

filter:
  mem: "4G"
  cores: 4
  output: "1.assay/cluster_logs/{rule}.{wildcards}.o"
  error: "1.assay/cluster_logs/{rule}.{wildcards}.e"

...
```

修改好上述文件后，使用下面的命令投递任务。

```shell
if [ ! -d 1.assay/cluster_logs ];then mkdir -p 1.assay/cluster_logs;fi

snakemake \
--snakefile rules/assemble.smk \
--configfile rules/config.yaml \
--cluster-config rules/cluster.yaml \
--jobs 80 \
--keep-going \
--rerun-incomplete \
--latency-wait 360 \
--cluster "qsub -S /bin/bash -cwd -q {cluster.queue} -P {cluster.project} -l vf={cluster.mem},p={cluster.cores} -binding linear:{cluster.cores} -o {cluster.output} -e {cluster.error}"
```

### 3.5. 多样本分箱

运行完上述的snakemake流程，即可得到单装单bin的结果。

在rules/config.yaml中更改多样本分箱使用的样本数量。测试数据默认为2，实际应用中推荐使用10个以上样本。

运行多样本分箱时，由于存在随机抽样的过程，为了使计算可重复，使用mulit_samples_binning.py生成所有的计算脚本。

```shell
python rules/tools/multi_samples_binning.py sample.txt rules/config.yaml > work_multi_binning.sh
```

work_multi_binning.sh中包含binning, checkm步骤，然后将work_multi_binning.sh投上计算任务即可。

### 4. 统计结果

过滤与去宿主的统计 : 01.assembly/filter_summary.txt

MEGAHIT组装统计: 01.assembly/contigs_stat.txt

Reads比对到contigs的回比率：02.binning_s1/map2scaftigs_summary.txt

MetaBAT2得到的所有的Bins: 02.binning_s1/All_bins_stat.txt

中高质量的MAGs : 02.binning_s1/picked_MAGs_quality.txt

每个样本得到的MAGs: 02.binning_s1/MAGs_per_sample.txt

多样本分箱统计结果文件命名同单样本分箱。


## 5. FAQ

### 1. 无需过滤宿主基因组？是单端测序？想换流程中的软件？

本流程不支持，请使用[metapi](https://github.com/ohmeta/metapi)。

### 2. 后续分析流程
dRep聚类: https://github.com/MrOlm/drep
基于唯一比对的profile：TBD
