适用于肠道，舌苔等高宿主率样本组装。15G数据量megahit需要40G内存，30G数据量需要80G。同等数据量下，metaspades需要的内存以及计算量都是megahit的3-4倍。

![pipline](dag.svg)

输入raw data，输出高中质量的bins。输入文件格式见sample.txt。内有demo数据可以直接运行测试。

使用方法见work.sh。

WDL版本源码见rules/WDL文件夹，可以在自动化平台()调用。

