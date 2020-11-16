大规模宏基因组样本的组装流程。整个流程包括：fastp和bowtie2过滤低质量的reads和宿主的reads，megahit组装，metabat2分箱，checkm质量评估。(中文文档)

本流程使用单样本组装，多样本分箱的策略。单个样本组装可以避免多样本组装，组装出质量更高的基因组(图1)。

If there is only one sample available, TetraNucleotides Frequency(TNF) and abundance for each contig were used to calculated a composite score(S) in Metabat2 . When the number of samples increases, abundance becomes more reliable. When there are three or more samples available, an abundance correlation score(COR) is calculated and the more accurate S is calculated as the geometric mean of TNF, abundance, and COR.

我们使用了1，4，6，8，10个样本进行了样本数量对组装结果影响的测试。4个样本即可大大增加中高质量的MAGs的数量，随着使用样本数的增加，得到的高中质量的MAG数量也一并增加(图2)。推荐使用10个样本进行多样本分箱。

MetaSpades(不包含在本流程中)相较于megahit组装效果更好，但会消耗大量的计算资源。若去完宿主后的reads测序量小于10G，可以考虑选择Metaspades进行组装。推荐Metapi，支持更全面的软件选择。

If you used this pipeline please cite the following paper:



## Install

使用conda安装所有依赖：



若在华大集群上，可使用已安装好的环境变量。
```
export PATH="$PATH:/ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin"
```

## Usage

### 0.demo

使用xxx，xxx,xxx的基因组生成测试数据集。



### 1. 样本

样本输入格式见sample.txt，以制表符分割。

### 2. 参数设置

### 3a. 本地运行

### 3b. 投递任务

### 4. 多样本分箱

运行完上述的snakemake流程，即可得到单装单bin的结果。

当使用10个样本对样本X进行分箱时，将样本X的contigs结果，以及从所有样本中随机抽取的10个样本(包括X自己)作为metabat2的输入。

为了该过程可重复，使用mulit_samples_metabat2.py而非snakemake生成所有的计算脚本。

## FAQ

### 1. WDL版本

WDL版本源码见rules/WDL文件夹，可以在自动化平台(http://10.225.5.11:8080/Bioinfo/pages/automation/taskAdd.jsp)调用。

### 2. 没有宿主基因组？是单端测序？想换流程中的软件？

请使用Metapi。

