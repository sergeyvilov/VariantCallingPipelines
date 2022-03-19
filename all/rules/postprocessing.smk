rule concat_vcfs:
    '''
    Concatenate calls from different chromosomes of the same sample
    '''
    input:
        vcfs = expand(progress_dir + "called/{{bam_file_name}}/{chrom}.vcf.gz", chrom = chroms_GRCh37),
    output:
        vcf = temp(progress_dir + "called_concat/{bam_file_name}.vcf.gz"),
        tbi = temp(progress_dir + "called_concat/{bam_file_name}.vcf.gz.tbi"),
    shell:
        '''
        bcftools concat --no-version {input.vcfs}  --threads 4  | bcftools sort - |bgzip -c  > {output.vcf}
        tabix -f  {output.vcf}
        '''

rule remove_somatic:
    '''
    Remove somatic variants from each sample using somatic VCFs
    '''
    input:
        vcf = progress_dir + "called_concat/{bam_file_name}.vcf.gz",
        tbi = progress_dir + "called_concat/{bam_file_name}.vcf.gz.tbi",
        somatic_vcf = config['somatic_vcf'],
        somatic_vcf_tbi = config['somatic_vcf']+'.tbi',
    output:
        vcf = temp(progress_dir + "somatic_removed/{bam_file_name}.vcf"),
    shell:
        r'''
        if test -f "{input.somatic_vcf}"; then
            bcftools isec -C -w1 {input.vcf} {input.somatic_vcf} > {output.vcf}
        else
            bcftools view {input.vcf} > {output.vcf}
        fi
        '''

rule select_bialleic:
    '''
    Select only bialleic variants
    '''
    input:
        vcf = progress_dir + "somatic_removed/{bam_file_name}.vcf",
    output:
        vcf = temp(progress_dir + "bialleic/{bam_file_name}.vcf.gz"),
        tbi = temp(progress_dir + "bialleic/{bam_file_name}.vcf.gz.tbi"),
    shell:
        r'''
        bcftools view --max-alleles 2 {input.vcf} -Oz -o {output.vcf};
        tabix -f {output.vcf}
        '''

rule annotate_gnomAD:
    '''
    Add population allele frequency (AF) from the gnomAD database

    Needs a lightweight gnomAD VCF, with all variants in a single file,
    the INFO field having only the gnomAD AF record,
    variants not passing gnomAD filters should be removed
    '''
    input:
        vcf = progress_dir + "bialleic/{bam_file_name}.vcf.gz",
        tbi = progress_dir + "bialleic/{bam_file_name}.vcf.gz.tbi",
    output:
        vcf = progress_dir + "gnomAD/{bam_file_name}.vcf.gz",
        tbi = progress_dir + "gnomAD/{bam_file_name}.vcf.gz.tbi",
    params:
        gnomAD_vcf = config["gnomAD_light_vcf"],
        workdir = progress_dir + "gnomAD"
    shell:
        r'''
        if test -f "{params.gnomAD_vcf}"; then
            echo '##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="Alternative allele frequency as in GNOMAD v.3.1.1">' > {params.workdir}/gnomad_header;
            bcftools annotate \
            -c 'ID,INFO/gnomAD_AF:=INFO/AF' \
            -h {params.workdir}/gnomad_header \
            -a {params.gnomAD_vcf} \
            {input.vcf} \
            -Oz -o {output.vcf};
            tabix -f {output.vcf}
        else
            cp {input.vcf} {output.vcf}
            cp {input.vcf}.tbi {output.vcf}.tbi
        fi
        '''

rule mark_germline:
    '''
    Mark germline variants using germline VCFs
    '''
    input:
        vcf = progress_dir + "gnomAD/{bam_file_name}.vcf.gz",
        tbi = progress_dir + "gnomAD/{bam_file_name}.vcf.gz.tbi",
    output:
        vcf = progress_dir + "germline_annotated/{bam_file_name}.vcf.gz",
    params:
        germline_vcfs_dir = config["germline_vcfs_dir"],
        workdir = progress_dir + "germline_annotated"
    log:
        progress_dir + 'logs/mark_germline/{bam_file_name}.log'
    shell:
        r'''
        echo '##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description="Germline variant">' > {params.workdir}/germline_header;
        bam_sample_name=$(bcftools query -l {input.vcf});
        #germline_vcf=$(ls {params.germline_vcfs_dir}${{bam_sample_name}}*|grep dkfz-snv.*vcf.gz$); #NAME PATTERN FOR GERMLINE VCF
        germline_vcf={params.germline_vcfs_dir}${{bam_sample_name}}.vcf.gz; #NAME PATTERN FOR GERMLINE VCF
        if test -f $germline_vcf; then
            bcftools view -f .,PASS $germline_vcf -Oz -o {params.workdir}/{wildcards.bam_file_name}.germline.vcf.gz;
            tabix -f {params.workdir}/{wildcards.bam_file_name}.germline.vcf.gz;
            bcftools annotate \
            -a {params.workdir}/{wildcards.bam_file_name}.germline.vcf.gz \
            -m +GERMLINE \
            {input.vcf} \
            -h {params.workdir}/germline_header \
            -Oz -o {output.vcf};
            rm {params.workdir}/{wildcards.bam_file_name}.germline.vcf.gz;
            rm {params.workdir}/{wildcards.bam_file_name}.germline.vcf.gz.tbi;
        else
            echo "Germline vcf $germline_vcf not found" > {log};
            cp {input.vcf} {output.vcf}
        fi
        '''


rule add_VAF_DP:
    '''
    Calculate VAF based on the FORMAT field and add it to the INFO field
    Replace DP in the INFO field with the DP in the FORMAT field
    '''
    input:
        vcf = progress_dir + "germline_annotated/{bam_file_name}.vcf.gz",
    output:
        vcf = temp(progress_dir + 'with_dp_vaf/{bam_file_name}.vcf'),
    params:
        workdir = progress_dir + 'with_dp_vaf/'
    shell:
        r'''
        patient={wildcards.bam_file_name}
        workdir={params.workdir}

        bcftools query -f "[%AD]\n" {input.vcf} |awk 'BEGIN {{FS=","}} {{if ($1+$2>0) {{VAF=$2/($1+$2)}} else {{VAF=0}};printf "DP=%d;VAF=%.2f\n",$1+$2,VAF}}' > $workdir/$patient.dpvaf;

        bcftools view -h {input.vcf}|sed -E -e 's/#CHROM*/##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fraction\">\n\0/' > {output.vcf};
        bcftools view -H {input.vcf}|sed -E 's/DP=[^;\t]+;?//'|paste - $workdir/$patient.dpvaf | awk 'BEGIN {{OFS="\t"}} {{$8=$8";"$NF; $NF=""; print $0}}' >> {output.vcf};

        rm $workdir/$patient.dpvaf;
        '''

rule mark_sample_and_project:
    '''
    Insert project name and the name of the BAM file that has the given variant into the INFO field
    This will help to identify projects and samples if we concate single-sample VCFs into a larger VCF later
    '''
    input:
        vcf = progress_dir + "with_dp_vaf/{bam_file_name}.vcf",
    output:
        vcf =  temp(progress_dir + 'per_sample/{bam_file_name}.vcf.gz'),
        tbi =  temp(progress_dir + 'per_sample/{bam_file_name}.vcf.gz.tbi'),
    params:
        project_name = config['project_name'],
        workdir = progress_dir + 'per_sample/'
    shell:
        r"""
        bam_file_name={wildcards.bam_file_name}
        workdir={params.workdir}

        bcftools view -h {input.vcf} \
        |sed -Ee '10i ##INFO=<ID=Project,Number=.,Type=String,Description=\"Project name\">' \
        -e '10i ##INFO=<ID=BAM,Number=.,Type=String,Description=\"BAM file name without extension\">' > $workdir/$bam_file_name.vcf;

        bcftools view -H {input.vcf} \
        | awk -v bam_file_name=$bam_file_name -v project_name={params.project_name} 'BEGIN {{OFS="\t"}} {{$8=$8";BAM="bam_file_name";Project="project_name; print $0}}' >> $workdir/$bam_file_name.vcf;

        bgzip $workdir/$bam_file_name.vcf;
        tabix -f $workdir/$bam_file_name.vcf.gz;
        """

rule merge_train_vcfs:
    '''
    Merge calls of all TRAIN samples
    '''
    input:
        vcf = expand(progress_dir + 'per_sample/{bam_file_name}.vcf.gz', bam_file_name = train_bams),
        tbi = expand(progress_dir + 'per_sample/{bam_file_name}.vcf.gz.tbi', bam_file_name = train_bams)
    output:
        vcf = progress_dir + 'results/train.vcf.gz',
        tbi = progress_dir + 'results/train.vcf.gz.tbi',
    shell:
        r'''
        bcftools merge {input.vcf} -Oz -o {output.vcf}
        tabix -f {output.vcf}
        '''
#
# rule generate_test_intervals:
#     '''
#     Generate BED files defining non-overlapping intervals from TEST samples which span a single genome
#     '''
#     input:
#         train_bams = expand(input_data_dir + '{bam_file_name}.bam', bam_file_name = test_bams),
#     output:
#         temp(expand(progress_dir + 'test_samples_intervals/{bam_file_name}.bed', bam_file_name = test_bams))
#     params:
#         train_bams = test_bams,
#         test_bams = [],
#         workdir = progress_dir + 'test_samples_intervals/',
#         window_size = 500000, #should be small enough s.t. the composite genome contains parts of all the train samples and large enough to minimize the number of step_size-window_size gaps
#         step_size = 500000 + 150, #window_size + read length+1, to avoid that pieces from different train samples overlap
#         genome = 'GRCh37',
#         split_by_chromosome = False,
#     script:
#         "../scripts/make_beds.py"
#
# rule select_regions:
#     '''
#     From each TEST sample, select calls corresponding to given intervals
#     '''
#     input:
#         vcf = progress_dir + 'per_sample/{bam_file_name}.vcf.gz',
#         tbi = progress_dir + 'per_sample/{bam_file_name}.vcf.gz.tbi',
#         bed = progress_dir + 'test_samples_intervals/{bam_file_name}.bed',
#     output:
#         vcf = temp(progress_dir + 'test_selected/{bam_file_name}.vcf.gz'),
#         tbi = temp(progress_dir + 'test_selected/{bam_file_name}.vcf.gz.tbi'),
#     shell:
#         '''
#         bcftools view -R {input.bed} {input.vcf} -Oz -o {output.vcf};
#         tabix -f {output.vcf}
#         '''
#
# rule merge_test_vcfs:
#     '''
#     Merge calls of all TEST samples
#     '''
#     input:
#         vcf = expand(progress_dir + 'test_selected/{bam_file_name}.vcf.gz', bam_file_name = test_bams),
#         tbi = expand(progress_dir + 'test_selected/{bam_file_name}.vcf.gz.tbi', bam_file_name = test_bams)
#     output:
#         vcf = progress_dir + 'results/test.vcf.gz',
#         tbi = progress_dir + 'results/test.vcf.gz.tbi',
#     shell:
#         r'''
#         bcftools merge {input.vcf} -Oz -o {output.vcf}
#         tabix -f {output.vcf}
#         '''
