##############################################################
### tianliu@genomics.cn
### 2020/07/16
### v1.1
### add summary_report
### fixed bug of bowtie2_index
### remove metaphlan3
###
### 2020/05/13
### metagenomics assembly v1.0
### input: raw fastq
### WDL PATH : /hwfssz5/ST_BIGDATA/USER/st_bigdata/test/Meta_Assemble_default_1.1.wdl
##############################################################
version 1.0
workflow meta_assembly_1 {
    input {
        String SampleID
        String ProjectID
        String Outdir
        String bowtie2_index_with_prefix
        String fastp_min_len = "30"
        String megahit_minContig_len = "200"
        String metabat2_minContig_len = "1500"
        String HQ_MAGs = "90,5,-1"
        String MQ_MAGs = "50,10,-1"
        Array[Array[String]] Data
        #[library,fq1,fq2]
    }

    String sample_ID=SampleID
    String work_dir=Outdir
    scatter(fq in Data){
        String sample_r1=fq[1]
        String sample_r2=fq[2]
        call filter{
            input:
                sample_r1 = sample_r1,
                sample_r2 = sample_r2,
                sample_ID = sample_ID,
                index = bowtie2_index_with_prefix,
                min_len = fastp_min_len,
                work_dir = work_dir
                }
    }
    call megahit{
        input:
            rmhost_r1 = filter.rmhost_r1,
            rmhost_r2 = filter.rmhost_r2,
            min_contigs_len = megahit_minContig_len,
            sample_ID = sample_ID,
            work_dir = work_dir,
            PID = bigCompute
            }

    call metabat2{
        input:
            rmhost_r1 = filter.rmhost_r1,
            rmhost_r2 = filter.rmhost_r2,
            sample_ID = sample_ID,
            contigs = megahit.contigs,
            minContig = metabat2_minContig_len,
            work_dir = work_dir
            }

    call checkm{
        input:
            bins_dir = metabat2.bins_dir,
            HQ = HQ_MAGs,
            MQ = MQ_MAGs,
            sample_ID = sample_ID,
            work_dir = work_dir,
            PID=bigCompute
            }
    call checkmStat{
        input:
            bins_dir = metabat2.bins_dir,
            HQ = HQ_MAGs,
            MQ = MQ_MAGs,
            sample_ID = sample_ID,
            work_dir = work_dir
            }
}

### task1 trimming & remove host reads
### To reduce disk storage usage, merge trimming and remove host together.
task filter {
    input {
        String sample_r1
        String sample_r2
        String sample_ID
        String index
        String min_len
        String work_dir
        String trim_dir = "1.assay/01.trimming"
        String log_trim_dir = "1.assay/logs/01.trimming"
        String rmhost_dir = "1.assay/02.rmhost"
        String log_rmhost_dir = "1.assay/logs/02.rmhost"
        Int cpu = 4
        Int mem = 4
    }

    command {
        set -e
        mkdir -p ${work_dir}/${trim_dir} ${work_dir}/${log_trim_dir} ${work_dir}/${rmhost_dir} ${work_dir}/${log_rmhost_dir}
        mkdir -p ${work_dir}/${trim_dir}/json ${work_dir}/${trim_dir}/html
        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/fastp \
        -i ${sample_r1} \
        -I ${sample_r2} \
        -o ${work_dir}/${trim_dir}/${sample_ID}.trimmed.1.fq.gz \
        -O ${work_dir}/${trim_dir}/${sample_ID}.trimmed.2.fq.gz \
        -w ${cpu} \
        --length_required ${min_len} \
        --disable_adapter_trimming \
        -j ${work_dir}/${trim_dir}/json/${sample_ID}.fastp.json \
        -h ${work_dir}/${trim_dir}/html/${sample_ID}.fastp.html \
        2> ${work_dir}/${log_trim_dir}/${sample_ID}.fastp.log

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/bowtie2 \
        --very-sensitive \
        -p ${cpu} \
        -x ${index} \
        -1 ${work_dir}/${trim_dir}/${sample_ID}.trimmed.1.fq.gz \
        -2 ${work_dir}/${trim_dir}/${sample_ID}.trimmed.2.fq.gz \
        2> ${work_dir}/${log_rmhost_dir}/${sample_ID}.rmhost.log | \
        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/samtools fastq \
        -N -c 5 -f 12 -F 256 \
        -1 ${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.1.fq.gz \
        -2 ${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.2.fq.gz - && \
        rm ${work_dir}/${trim_dir}/${sample_ID}.trimmed.*.fq.gz

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/seqkit stat ${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.*.fq.gz -j ${cpu} > ${work_dir}/${log_rmhost_dir}/${sample_ID}.rmhost.reads.summary
    }
    runtime{
        backend:"SGE"
        cpu:cpu
        memory:"${mem} GB"
    }
    output {
        String rmhost_r1 = "${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.1.fq.gz"
        String rmhost_r2 = "${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.2.fq.gz"
        File rmhost_r1a = "${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.1.fq.gz"
        File rmhost_r2a = "${work_dir}/${rmhost_dir}/${sample_ID}.rmhost.2.fq.gz"
    }
}

### task2 assemble
task megahit {
    input {
        Array[String] rmhost_r1
        Array[String] rmhost_r2
        String sample_ID
        String work_dir
        String min_contigs_len
        String assemble_dir = "1.assay/03.assembly/megahit"
        String log_assemble = "1.assay/logs/03.assembly/megahit"
        String PID
        Int cpu = 36
        Int mem = 50
    }

    command{
        set -e
        if [ -d ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out ];then rm -r ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out;fi
        mkdir -p ${work_dir}/${assemble_dir} ${work_dir}/${log_assemble}

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/megahit \
        -1 ${sep=',' rmhost_r1} \
        -2 ${sep=',' rmhost_r2} \
        -t ${cpu} \
        --min-contig-len ${min_contigs_len} \
        -o ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out \
        2> ${work_dir}/${log_assemble}/${sample_ID}.megahit.log

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/pigz -p ${cpu} -c ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/final.contigs.fa > ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}.megahit.contigs.fa.gz && rm ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/final.contigs.fa

        rm -r ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/intermediate_contigs
    }
    runtime{
        backend:"SGE"
        cpu:cpu
        memory:"${mem} GB"
        sge_queue:PID
    }
    output {
        String contigs = "${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}.megahit.contigs.fa.gz"
        File contigs1 = "${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}.megahit.contigs.fa.gz"
    }
}

### task3 binning
task metabat2 {
    input {
        Array[String] rmhost_r1
        Array[String] rmhost_r2
        String contigs
        String minContig
        String sample_ID
        String work_dir
        String assemble_dir = "1.assay/03.assembly/megahit"
        String binning_dir = "1.assay/04.binning/metabat2"
        String log_binning = "1.assay/logs/04.binning/metabat2"
        Int cpu = 8
        Int mem = 8
    }

    command{
        set -e
        ### prepare
        mkdir -p ${work_dir}/${binning_dir}/${sample_ID} ${work_dir}/${log_binning}
        if [ -d ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index ];then \
        rm -rf ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index;fi

        find ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out -name "*.bam" | xargs rm -f
        mkdir ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index

        config_size=`tail -n 2 ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/log|head -n 1|grep -Po '(?<=total )\d+(?= bp,)'`
        limit_size=4000000000

        if [ $config_size -gt $limit_size ];then \
          /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/bowtie2-build \
          --large-index \
          --threads ${cpu} \
          ${contigs} \
          ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index/g_index \
          > ${work_dir}/${log_binning}/${sample_ID}.index.log
        else
          /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/bowtie2-build \
          --threads ${cpu} \
          ${contigs} \
          ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index/g_index \
          > ${work_dir}/${log_binning}/${sample_ID}.index.log
        fi

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/bowtie2 -p ${cpu} \
        -x ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index/g_index \
        -1 ${sep=',' rmhost_r1} \
        -2 ${sep=',' rmhost_r2} \
        2> ${work_dir}/${log_binning}/${sample_ID}.map2scaftigs.log | \
        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/samtools sort \
        -@ ${cpu} \
        -o ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}.sorted.bam -

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/jgi_summarize_bam_contig_depths \
        --outputDepth ${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}.depth.txt \
        ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}.sorted.bam

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/metabat2 -i ${contigs} \
        -a ${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}.depth.txt \
        -o ${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}_binning/bin \
        -m ${minContig} \
        -t ${cpu} > ${work_dir}/${log_binning}/${sample_ID}.metabat2.log

        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/pigz -p ${cpu} ${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}.depth.txt


        rm -f ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}.sorted.bam
        rm -rf ${work_dir}/${assemble_dir}/${sample_ID}.megahit_out/${sample_ID}_index
    }
    runtime{
        backend:"SGE"
        cpu:cpu
        memory:"${mem} GB"
    }
    output {
        String bins_dir = "${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}_binning"
        File depth = "${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}.depth.txt.gz"
    }
}

### task4 quality
task checkm {
    input {
        String bins_dir
        String log_checkm = "1.assay/logs/05.MAG_quality/checkm"
        String quality_dir = "1.assay/05.quality/checkm"
        String binning_dir = "1.assay/04.binning/metabat2"
        String HQ
        String MQ
        String sample_ID
        String work_dir
        Int cpu = 4
        Int mem = 40
        String PID
    }

    command{
        set -e
        mkdir -p ${work_dir}/${quality_dir}/${sample_ID} ${work_dir}/${log_checkm}
        export TMPDIR="/tmp"
        export PATH=/ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin:$PATH
        /ldfssz1/ST_META/share/User/tianliu/bioenv/conda/envs/meta_assemble_wdl/bin/checkm \
        lineage_wf -t ${cpu} -x fa \
        ${bins_dir} ${work_dir}/${quality_dir}/${sample_ID} \
        2> ${work_dir}/${log_checkm}/${sample_ID}.checkm.log | \
        grep -v "INFO" > ${work_dir}/${quality_dir}/${sample_ID}/checkm_summary.txt
    }
    runtime{
        backend:"SGE"
        cpu:cpu
        memory:"${mem} GB"
        sge_queue:PID
    }
    output {
        String checkm_res = "${work_dir}/${quality_dir}/${sample_ID}/checkm_summary.txt"
        File checkm_res1 = "${work_dir}/${quality_dir}/${sample_ID}/checkm_summary.txt"
    }
}
task checkmStat {
    input {
        String bins_dir
        String log_checkm = "1.assay/logs/05.MAG_quality/checkm"
        String quality_dir = "1.assay/05.quality/checkm"
        String binning_dir = "1.assay/04.binning/metabat2"
        String HQ
        String MQ
        String sample_ID
        String work_dir
        Int cpu = 1
        Int mem = 4
    }
    command{
        python /ldfssz1/ST_META/share/User/tianliu/pipline/stable/meta_assemble_v0.2/rules/tools/bin_stat_format.py ${work_dir}/${quality_dir}/${sample_ID}/storage/bin_stats.analyze.tsv > ${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}.bins.stat.txt

        python /ldfssz1/ST_META/share/User/tianliu/pipline/stable/meta_assemble_v0.2/rules/tools/pick_MAGs.py --high ${HQ} --medium ${HQ} ${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}_binning ${work_dir}/${quality_dir}/${sample_ID}/checkm_summary.txt > ${work_dir}/${quality_dir}/${sample_ID}/${sample_ID}.picked.summary.txt

        rm -r ${work_dir}/${quality_dir}/${sample_ID}/storage ${work_dir}/${quality_dir}/${sample_ID}/bins
    }
    runtime{
        backend:"SGE"
        cpu:cpu
        memory:"${mem} GB"
    }
    output {
        File checkmStat1 = "${work_dir}/${quality_dir}/${sample_ID}/${sample_ID}.picked.summary.txt"
        File checkmStat2 = "${work_dir}/${binning_dir}/${sample_ID}/${sample_ID}.bins.stat.txt"
    }
}
