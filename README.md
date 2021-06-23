**注意：本工作尚未发表，仅限华大内部使用**

## 1. Metagenome assembly
大规模宏基因组样本的组装流程。整个流程包括：fastp和Bowtie 2过滤低质量的reads和宿主的reads，MEGAHIT组装，MetaBAT 2分箱，CheckM评估基因组质量。

本流程使用单样本组装，多样本分箱的策略。单个样本组装可以避免多样本组装，组装出质量更高的基因组(图1)。而在分箱过程中给予更多的样本可以使metabat2更为准确:

>When the number of samples increases, abundance (ABD) becomes more reliable, so *w* effectively decreases the relative weight of tetra-nucleotide frequency (TNF) and increases the weight of ABD. Whenever there are three or more samples available, an ABD correlation score (COR) is also calculated using the Pearson correlation coefficient and then rank-normalized using ABD. In this case, a composite score (S) is calculated as the geometric mean of TNF, ABD, and COR.

![图1](pic/fig1.png)

<center>图1 菌株混杂会导致组装质量下降。</center>
<center>来源: https://drep.readthedocs.io/en/latest/overview.html#genome-de-replication </center>

我们分别使用了1，4，6，8，10个样本评估了样本数量对分箱结果的影响。4个样本即可大大增加中高质量的MAGs的数量，随着使用样本数的增加，得到的高中质量的MAG数量也一并增加(图2)。推荐使用10个样本进行多样本分箱。

![fig2](pic/Fig2.png)

<center>图2 更多的样本可以更好的改善分箱的质量</center>

metaSPAdes(不包含在本流程中)相较于MEGAHIT组装效果更好，但会消耗大量的计算资源。若去完宿主后的reads测序量小于10G，可以考虑选择metaSPAdes进行组装。本流程不支持更改软件，若想换用其他软件进行组装，推荐使用[metapi](https://github.com/ohmeta/metapi)。

若使用了本流程，请引用文献：

If you used this pipeline please cite the following paper:



### 1.1 Install

使用conda安装所有依赖：

```shell
# TBD
#conda env create -n meta_assembly -f rules/environment.yaml
```



若在华大集群上，可使用已安装好的环境变量。
```
export PATH="/ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin:$PATH"
```

### 1.2 Usage

#### 1.2.1. 测试数据集

使用三个单菌的基因组生成模拟的测试集。

```shell
sh demo/work.iss.sh
```

#### 1.2.2. 样本

样本输入格式见sample.txt，以制表符分割。仅支持双端测序样本。

```shell
$ cat sample.txt 
id	fq1	fq2
test1	demo/test1_data_R1.fastq.gz	demo/test1_data_R2.fastq.gz
test2	demo/test2_data_R1.fastq.gz	demo/test2_data_R2.fastq.gz
test3	demo/test3_data_R1.fastq.gz	demo/test3_data_R2.fastq.gz
```

#### 1.2.3 参数设置

软件运行参数和流程的文件结构都在```rules/assemble.config.yaml```。

默认宿主为人的hg38。若宿主为其他物种，需下载该物种的基因组建Bowtie 2的索引，并在config文件中替换。

```yaml
params:
    ### step1 assembly
    fastp:
      min_len : 30 # Recommended value >= 30
    #rmhost
    rmhost:
      bowtie2_index : "/ldfssz1/ST_META/share/User/tianliu/database/bowtie2_index/hg38/hg38"
      threads: 8
    #megahit
    megahit:
      threads: 16
      min_contigs_len: 200
    #binning
    metabat2:
      multi_samples: 2 # Recommended 10
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

#### 1.2.4a. 本地运行

在运行前，先使用-n进行测试。

```shell
snakemake \
--snakefile rules/assemble.smk \
--configfile rules/assemble.config.yaml \
-n
```

若测试正常，则使用--core指定CPU资源在本地运行。

```shell
snakemake \
--snakefile rules/assemble.smk \
--configfile rules/assemble.config.yaml \
--core 32 2> smk.log
```

#### 1.2.4b. 投递任务(qsub环境)

投递任务相关参数在```rules/assemble.cluster.yaml```，可在该文件中指定每一步所需的资源。投任务前，需**修改项目编号，集群队列**。

```yaml
localrules: all

__default__:
  queue: "st.q"
  project: "P19Z10200N0314"
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
--configfile rules/assemble.config.yaml \
--cluster-config rules/assemble.cluster.yaml \
--jobs 80 \
--keep-going \
--rerun-incomplete \
--latency-wait 360 \
--cluster "qsub -S /bin/bash -cwd -q {cluster.queue} -P {cluster.project} -l vf={cluster.mem},p={cluster.cores} -binding linear:{cluster.cores} -o {cluster.output} -e {cluster.error}"
```

### 4. 多样本分箱

运行完上述的snakemake流程，即可得到单装单bin的结果。

当使用10个样本对样本X进行分箱时，将样本X的contigs结果，以及从所有样本中随机抽取的10个样本(包括X自己)作为metabat2的输入。由于存在随机抽样的过程，为了使计算可重复，使用mulit_samples_binning.py而非snakemake生成所有的计算脚本。

```shell
python rules/tools/multi_samples_binning.py sample.txt rules/assemble.config.yaml > work_multi_binning.sh
```

然后将work_multi_binning.sh拆分投上任务即可。

### 5. 统计结果

过滤与去宿主的统计 : 01.assembly/filter_summary.txt

MEGAHIT组装统计: 01.assembly/contigs_stat.txt

Reads比对到contigs的回比率：02.binning_s1/map2scaftigs_summary.txt

MetaBAT 2得到的所有的Bins: 02.binning_s1/All_bins_stat.txt

中高质量的MAGs : 02.binning_s1/picked_MAGs_quality.txt

每个样本得到的MAGs: 02.binning_s1/MAGs_per_sample.txt

多样本分箱统计结果文件命名同单样本分箱。

## 2. Genome de-replication
## 3. Genome annation
## 4. Figure

## FAQ

### 1. WDL版本

WDL版本源码见rules/WDL文件夹，可以在自动化平台(http://10.225.5.11:8080/Bioinfo/pages/automation/taskAdd.jsp)调用。

### 2. 无需过滤宿主基因组？是单端测序？想换流程中的软件？

本流程不支持，请出门左转metapi。

