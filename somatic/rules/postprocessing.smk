rule annotate_gnomAD:
    '''
    Add population allele frequency (AF) from the gnomAD database

    Needs a lightweight gnomAD VCF, with all variants in a single file,
    the INFO field having only the gnomAD AF record,
    variants not passing gnomAD filters should be removed
    '''
    input:
        vcf = filtered_dir + "{patient}.filtered.vcf.gz",
        tbi = filtered_dir + "{patient}.filtered.vcf.gz.tbi",
    output:
        vcf = temp(progress_dir + "gnomAD/{patient}.vcf.gz"),
    params:
        gnomAD_vcf = config["gnomAD_light_vcf"],
        workdir = progress_dir + "gnomAD"
    shell:
        r'''
        echo '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="Alternative allele frequency as in GNOMAD v.3.1.1">' > {params.workdir}/gnomad_header;
        bcftools annotate --threads 4 \
        -c 'ID,INFO/gnomAD_AF:=INFO/AF' \
        -h {params.workdir}/gnomad_header \
        -a {params.gnomAD_vcf} \
        {input.vcf} \
        -Oz -o {output.vcf}
        '''

rule add_tumor_VAF:
    '''
    Add tumor variant allele fraction to the INFO field
    '''
    input:
        vcf = progress_dir + "gnomAD/{patient}.vcf.gz",
    output:
        vcf = temp(progress_dir + 'with_vaf/{patient}.vcf'),
        vaf = temp(progress_dir + 'with_vaf/{patient}.vaf'),
    params:
        workdir = progress_dir + 'with_vaf/'
    shell:
        r'''
        patient={wildcards.patient}
        workdir={params.workdir}

        tumor_name=$(bcftools view -h {input.vcf}|grep  "tumor_sample"|sed 's/[^=]*=//');

        bcftools query -s "$tumor_name" -f "[%AD\t]\n" {input.vcf} |awk 'BEGIN {{FS=","}} {{if ($1+$2>0) {{VAF=$2/($1+$2)}} else {{VAF=0}};printf "AD_ref=%d;AD_alt=%d;VAF=%.2f\n",$1,$2,VAF}}' > $workdir/$patient.vaf;

        bcftools view -h {input.vcf}|sed -E -e 's/#CHROM*/##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Tumor variant allele fraction\">\n\0/' \
        -e 's/#CHROM*/##INFO=<ID=AD_ref,Number=1,Type=Integer,Description=\"Reference allele depth at the variant site\">\n\0/' \
        -e 's/#CHROM*/##INFO=<ID=AD_alt,Number=1,Type=Integer,Description=\"Alternative allele depth at the variant site\">\n\0/' > {output.vcf};
        bcftools view -H {input.vcf}|paste - $workdir/$patient.vaf | awk 'BEGIN {{OFS="\t"}} {{$8=$8";"$NF; $NF=""; print $0}}' >> {output.vcf};
        '''

rule remove_FORMAT:
    '''
    Remove FORMAT columns of tumor and normal samples
    '''
    input:
        vcf = progress_dir + 'with_vaf/{patient}.vcf',
        vaf = progress_dir + 'with_vaf/{patient}.vaf',
    output:
        vcf = temp(progress_dir + "without_format/{patient}.vcf"),
    shell:
        r'''
            bcftools view -h {input.vcf}|sed '/#CHROM/ s/INFO.*/INFO/' > {output.vcf}

            bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n" {input.vcf} >> {output.vcf}

        '''

rule mark_sample_and_project:
    '''
    Insert project name and the name of the BAM file that has the given variant into the INFO field
    This will help to identify projects and samples if we concate single-sample VCFs into a larger VCF later
    '''
    input:
        vcf = progress_dir + 'without_format/{patient}.vcf',
    output:
        vcf = progress_dir + 'per_sample/{patient}.vcf.gz',
        tbi = progress_dir + 'per_sample/{patient}.vcf.gz.tbi'
    params:
        project_name = 'MLL',
        bam_file_name = lambda wildcards:wildcards.patient+'.tumor',
        workdir = progress_dir + 'per_sample/'
    shell:
        r"""
        patient={wildcards.patient}
        workdir={params.workdir}

        bcftools view -h {input.vcf} \
        |sed -Ee '10i ##INFO=<ID=Project,Number=.,Type=String,Description=\"Project name\">' \
        -e '10i ##INFO=<ID=BAM,Number=.,Type=String,Description=\"BAM file name without extension\">' > $workdir/$patient.vcf;

        bcftools view -H --max-alleles 2 {input.vcf} \
        | awk -v bam_file_name={params.bam_file_name} -v project_name={params.project_name} 'BEGIN {{OFS="\t"}} {{$8=$8";BAM="bam_file_name";Project="project_name; print $0}}' >> $workdir/$patient.vcf;

        bgzip $workdir/$patient.vcf;
        tabix -f $workdir/$patient.vcf.gz;
        """



rule concat_vcfs:
    '''
    Concatenate calls from all samples into a single VCF
    '''
    input:
        vcf = expand(progress_dir + 'per_sample/{patient}.vcf.gz', patient=all_patients),
        tbi = expand(progress_dir + 'per_sample/{patient}.vcf.gz.tbi', patient=all_patients)
    output:
        vcf = temp(progress_dir + 'results/all.vcf.gz'),
        tbi = temp(progress_dir + 'results/all.vcf.gz.tbi')
    shell:
        r"""
        bcftools concat --no-version {input.vcf} --no-version --threads 4  -a| bcftools sort - | bgzip -c > {output.vcf};
        tabix -f {output.vcf};
        """

rule select_calls:
    '''
    Split SNPs and INDELs
    '''
    input:
        vcf = progress_dir + 'results/all.vcf.gz',
        tbi = progress_dir + 'results/all.vcf.gz.tbi'
    output:
        vcf = progress_dir + 'results/{vartype}/{vartype}.vcf.gz',
        tbi = progress_dir + 'results/{vartype}/{vartype}.vcf.gz.tbi',
    params:
        vartype = lambda w: "snps" if w.vartype == "snvs" else "indels",
    log:
        logs_dir + 'select_calls/{vartype}.log'
    shell:
        r'''
        bcftools view -v {params.vartype} {input.vcf} -Oz -o {output.vcf} > {log} 2>&1
        tabix -f {output.vcf}
        '''
