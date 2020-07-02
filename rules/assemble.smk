##############################################################
### tianliu@genomics.cn
### 2020/05/06
### metagenomics assembly
##############################################################

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
        os.path.join(config["results"]["assembly"], "filter_summary.txt"),
        os.path.join(config["results"]["assembly"], "All_bins_stat.txt"),
        os.path.join(config["results"]["assembly"], "MAGs_per_sample.txt"),
        os.path.join(config["results"]["assembly"], "picked_MAGs_quality.txt")

### step1 : trimming & remove host reads
### To reduce disk storage usage, merge trimming and remove host together.

rule filter:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        trim_r1 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.1.fq.gz")),
        trim_r2 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.2.fq.gz")),
        html = os.path.join(config["assay"]["trimming"], "html/{sample}.fastp.html"),
        json = os.path.join(config["assay"]["trimming"], "json/{sample}.fastp.json"),
        rmhost_r1 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
        rmhost_r2 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
    params:
        min_len = config["params"]["fastp"]["min_len"],
        index = config["params"]["rmhost"]["bowtie2_index"]
    threads:
        config["params"]["rmhost"]["threads"]
    log:
        fastp_log = os.path.join(config["logs"]["trimming"], "{sample}.fastp.log"),
        bowtie2_log = os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.log")
    run:
        shell(
        '''
        fastp -i {input.r1} -I {input.r2} -o {output.trim_r1} -O {output.trim_r2} -w 6 --length_required {params.min_len} --disable_adapter_trimming -j {output.json} -h {output.html} 2> {log.fastp_log}

        bowtie2 --very-sensitive -p {threads} -x {params.index} -1 {output.trim_r1} -2 {output.trim_r2} 2> {log.bowtie2_log} | samtools fastq -N -c 5 -f 12 -F 256 -1 {output.rmhost_r1} -2 {output.rmhost_r2} -
        ''')

rule seqkit_stat:
    input:
        expand("{rmhost_dir}/{{sample}}.rmhost.{reads}.fq.gz", rmhost_dir = config["assay"]["rmhost"], reads = ["1","2"])
    output:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.reads.summary")
    shell:
        "seqkit stat {input} > {output}"

### step2 : assemble by megahit
rule megahit:
    input:
        r1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        r2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    output:
        protected(os.path.join(config["assay"]["megahit"], "{sample}.megahit_out/{sample}.megahit.contigs.fa.gz"))
    params:
        outdir = directory(os.path.join(config["assay"]["megahit"], "{sample}.megahit_out")),
        min_c_len = config["params"]["megahit"]["min_contigs_len"]
    log:
        os.path.join(config["logs"]["megahit"], "{sample}.megahit.log")
    threads: 
        config["params"]["megahit"]["threads"]
    priority: 10
    shell:
        '''
        if [ -d {params.outdir} ];then rm -r {params.outdir};fi
        megahit -1 {input.r1} -2 {input.r2} -t {threads} --min-contig-len {params.min_c_len} -o {params.outdir} 2> {log}
        pigz -p {threads} -c {params.outdir}/final.contigs.fa > {output}
        rm -r {params.outdir}/intermediate_contigs
        '''

### step3 : binning
rule metabat2:
    input:
        r1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        r2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"),
        scaftigs = os.path.join(config["assay"]["megahit"], "{sample}.megahit_out/{sample}.megahit.contigs.fa.gz")
    output:
        bam = temp(os.path.join(config["assay"]["megahit"], "{sample}.megahit_out/{sample}.sorted.bam")),
        depth = protected(os.path.join(config["assay"]["metabat2"], "{sample}/{sample}.depth.txt")),
        bins_dir = directory(os.path.join(config["assay"]["metabat2"], "{sample}/{sample}_binning")),
    params:
        megahit_dir = os.path.join(config["assay"]["megahit"], "{sample}.megahit_out"),
        index_dir = directory(os.path.join(config["assay"]["megahit"], "{sample}.megahit_out/{sample}_index")),
        minContig = config["params"]["metabat2"]["minContig"],
        megahit_log = os.path.join(config["assay"]["megahit"], "{sample}.megahit_out/log")
    log:
        index_log = os.path.join(config["logs"]["metabat2"], "index/megahit/{sample}.index.log"),
        map_log = os.path.join(config["logs"]["metabat2"], "map2scaftigs/megahit/{sample}.map2scaftigs.log"),
        bin_log = os.path.join(config["logs"]["metabat2"], "metabat2/megahit/{sample}.metabat2.log")
    threads: 
        config["params"]["metabat2"]["threads"]
    priority: 20
    shell:
        '''
        ### prepare
        if [ ! -d {params.index_dir} ]
        then 
            mkdir {params.index_dir}
        else
            rm {params.index_dir}/*
        fi
        rm {params.megahit_dir}/*.bam

        ### fixed large-index bug
        config_size=`tail -n 2 {params.megahit_log} | head -n 1 | grep -Po '(?<=total )\d+(?= bp,)'`
        limit_size=40000000000
        if [ $config_size -gt $limit_size ]
        then
            bowtie2-build --large-index --threads {threads} {input.scaftigs} {params.index_dir}/g_index 1> {log.index_log} 2>&1
        else
            bowtie2-build --threads {threads} {input.scaftigs} {params.index_dir}/g_index 1> {log.index_log} 2>&1
        fi

        bowtie2 -p {threads} -x {params.index_dir}/g_index -1 {input.r1} -2 {input.r2} 2> {log.map_log} | samtools sort -@ {threads} -o {output.bam} -
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}
        metabat2 -i {input.scaftigs} -a {output.depth} -o {output.bins_dir}/bin -m {params.minContig} -t {threads} > {log.bin_log}
        rm -rf {params.index_dir}
        '''

rule seqkit_bins:
    input:
        os.path.join(config["assay"]["metabat2"], "{sample}/{sample}_binning")
    output:
        os.path.join(config["assay"]["metabat2"], "{sample}/{sample}.MAGs.stat")
    shell:
        "seqkit stat -a {input}/*.fa > {output}"

### step4: MAG quality
rule checkm:
    input:
        os.path.join(config["assay"]["metabat2"], "{sample}/{sample}_binning")
    output:
        directory(os.path.join(config["assay"]["checkm"], "{sample}"))
    log:
        os.path.join(config["logs"]["checkm"], "{sample}.checkm.log")
    threads: 
        config["params"]["checkm"]["threads"]
    shell:
        '''
        mkdir -p {output}
        checkm lineage_wf -t {threads} -x fa {input} {output} 2> {log} | grep -v "INFO" > {output}/checkm_summary.txt
        rm -r {output}/bins {output}/storage
        '''

rule pick_MAGs:
    input:
        bins_dir = os.path.join(config["assay"]["metabat2"], "{sample}/{sample}_binning"),
        checkm_dir = os.path.join(config["assay"]["checkm"], "{sample}")
    output:
        os.path.join(config["assay"]["picked"], "{sample}.picked.summary.txt")
    params:
        high = config["params"]["pick"]["HQ"],
        medium = config["params"]["pick"]["MQ"]
    shell:
        "python rules/tools/pick_MAGs.py --high {params.high} --medium {params.medium} {input.bins_dir} {input.checkm_dir}/checkm_summary.txt > {output}"

### step5 summary
rule filter_summary:
    input:
        trim = expand("{trim_dir}/json/{sample}.fastp.json", trim_dir = config["assay"]["trimming"], sample = _samples.index),
        rmhost = expand("{rmhost_dir}/{sample}.rmhost.reads.summary", rmhost_dir = config["logs"]["rmhost"], sample = _samples.index)
    output:
        protected(os.path.join(config["results"]["assembly"], "filter_summary.txt"))
    params:
        trim_summary = os.path.join(config["results"]["assembly"], "trim_summary.txt"),
        rmhost_summary = os.path.join(config["results"]["assembly"], "rmhost_summary.txt")
    run:
        shell(
        '''
        python rules/tools/filter_summary.py -t {input.trim} > {params.trim_summary}
        python rules/tools/filter_summary.py -r {input.rmhost} > {params.rmhost_summary}
        python rules/tools/merge_summary.py {params.trim_summary} {params.rmhost_summary} {output}
        rm {params.trim_summary} {params.rmhost_summary}
        ''')

rule MAGs_summary:
    input:
        MAGs_stat = expand("{picked_log}/{sample}.picked.summary.txt", picked_log = config["assay"]["picked"], sample = _samples.index),
        bins_stat = expand("{bin_dir}/{sample}/{sample}.MAGs.stat", bin_dir = config["assay"]["metabat2"], sample = _samples.index)
    output:
        bins_stat = protected(os.path.join(config["results"]["assembly"], "All_bins_stat.txt")),
        MAGs_per_sample = protected(os.path.join(config["results"]["assembly"], "MAGs_per_sample.txt")),
        MAGs_quality = protected(os.path.join(config["results"]["assembly"], "picked_MAGs_quality.txt"))
    shell:
        '''
        cat {input.bins_stat} > {output.bins_stat}
        python rules/tools/MAGs_summary.py {input.MAGs_stat} -o {output.MAGs_per_sample} -O {output.MAGs_quality}
        '''
