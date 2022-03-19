#quality control on BAM files: coverage, duplicates, etc..

rule coverage:
    input:
        bam = input_data_dir + "{patient}.{sample_type}.bam",
        bai = input_data_dir + "{patient}.{sample_type}.bai",
    output:
        qc_dir + "{patient}.{sample_type}.mosdepth.region.dist.txt",
        qc_dir + "{patient}.{sample_type}.regions.bed.gz",
        qc_dir + "{patient}.{sample_type}.mosdepth.global.dist.txt",
        qc_dir + "{patient}.{sample_type}.mosdepth.summary.txt"
    threads: 4
    params:
        by = 500
        #by=lambda wildcards, input: '500' if seqtype == 'WGS' else input.regions
    shell:
        """
        mosdepth --by {params.by} -t {threads} --fast-mode {qc_dir}{wildcards.patient}.{wildcards.sample_type} \
            {input.bam}
        """

rule stats:
    input:
        bam = input_data_dir + "{patient}.{sample_type}.bam",
        bai = input_data_dir + "{patient}.{sample_type}.bai",
    output:
        qc_dir + "{patient}.{sample_type}.flagstat"
    shell:
        """
        samtools flagstat {input.bam} > {output}
        """

rule fastqc:
    # FastQC gives inadequate estimate of duplicates number
    input:
        bam = input_data_dir + "{patient}.{sample_type}.bam",
        bai = input_data_dir + "{patient}.{sample_type}.bai",
    output:
        html=qc_dir + "fastqc/{patient}.{sample_type}_fastqc.html",
        zipdata=qc_dir + "fastqc/{patient}.{sample_type}_fastqc.zip"
    shell:
        """
        tmpdir={qc_dir}fastqc/{wildcards.patient}-{wildcards.sample_type}.tmp
        mkdir $tmpdir
        fastqc --outdir $tmpdir {input}
        mv $tmpdir/{wildcards.patient}.{wildcards.sample_type}_fastqc.html {output.html}
        mv $tmpdir/{wildcards.patient}.{wildcards.sample_type}_fastqc.zip {output.zipdata}
        rm -r $tmpdir
        """

rule multiqc:
    input:
        expand(qc_dir + "fastqc/{patient}.{sample_type}_fastqc.zip", patient=patients, sample_type=sample_types),
        expand(qc_dir + "{patient}.{sample_type}.mosdepth.region.dist.txt", patient=patients, sample_type=sample_types),
        expand(qc_dir + "{patient}.{sample_type}.flagstat", patient=patients, sample_type=sample_types)
    output:
        qc_dir + "multiqc_report.html"
    log:
        logs_dir + 'multiqc.log'
    wrapper:
        "0.50.4/bio/multiqc"

rule seq_depths:
    input:
        expand(qc_dir + "{patient}.{sample_type}.mosdepth.summary.txt", patient=patients, sample_type=sample_types)
    output:
        qc_dir + "depths.csv"
    log:
        logs_dir + 'seq_depths.log'
    script:
        "../scripts/gather_depths.py"

# rule plot_depths:
#     #R script doesn't work((
#     input:
#         qc_dir + "depths.csv"
#     output:
#         qc_dir + "depths.svg"
#     script:
#         "../scripts/plot_depth.R"
