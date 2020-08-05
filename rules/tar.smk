import os
import sys
import pandas

shell.executable("bash")
#configfile: "config.yaml"

def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]


_samples = parse_samples("sample.txt")

rule all:
    input:
      expand("{megahit}/{sample}.megahit_out.tar.gz", megahit = config["assay"]["megahit"], sample = _samples.index),
      expand("{metabat2}/{sample}.tar.gz", metabat2 = config["assay"]["metabat2"], sample = _samples.index),
      expand("{checkm}/{sample}.tar.gz", checkm = config["assay"]["checkm"], sample = _samples.index),
      os.path.join(config["assay"]["all"], "01.trimming.tar.gz"),
      os.path.join(config["assay"]["all"], "logs.tar.gz"),
      os.path.join(config["assay"]["all"], "picked.tar.gz")


### step6 tar #####
rule tar_sample:
    input:
        assembly_dir = os.path.join(config["assay"]["megahit"], "{sample}.megahit_out"),
        binning_dir = os.path.join(config["assay"]["metabat2"], "{sample}"),
        quality_dir = os.path.join(config["assay"]["checkm"], "{sample}")
    output:
        assembly_tar = protected(os.path.join(config["assay"]["megahit"], "{sample}.megahit_out.tar.gz")),
        binning_tar = protected(os.path.join(config["assay"]["metabat2"], "{sample}.tar.gz")),
        quality_tar = protected(os.path.join(config["assay"]["checkm"], "{sample}.tar.gz"))
    shell:
        '''
        tar -zcf {output.assembly_tar} {input.assembly_dir} && rm -rf {input.assembly_dir}
        tar -zcf {output.binning_tar} {input.binning_dir} && rm -rf {input.binning_dir}
        tar -zcf {output.quality_tar} {input.quality_dir} && rm -rf {input.quality_dir}
        '''

rule tar_log:
    input:
        trim_dir = config["assay"]["trimming"],
        logs_dir = config["logs"]["all"],
        picked_log = config["assay"]["picked"],
    output:
        trim_tar = protected(os.path.join(config["assay"]["all"], "01.trimming.tar.gz")),
        logs_tar = protected(os.path.join(config["assay"]["all"], "logs.tar.gz")),
        picked_tar = protected(os.path.join(config["assay"]["all"], "picked.tar.gz"))
    shell:
        '''
        tar -zcf {output.trim_tar} {input.trim_dir} && rm -rf {input.trim_dir}
        tar -zcf {output.logs_tar} {input.logs_dir} && rm -rf {input.logs_dir}
        tar -zcf {output.picked_tar} {input.picked_log} && rm -rf {input.picked_log}
        '''
