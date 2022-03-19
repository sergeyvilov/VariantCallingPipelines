rule HaplotypeCaller:
    '''
    Call germline variants on selected intervals
    '''
    input:
        tumor_bam = input_data_dir + "{bam_file_name}.bam",
        intervals = intervals_dir + "{bam_file_name}/{chrom}.bed",
        ref = refgen_fa,
        known = resources["dbsnp_all_vcf"], #used to annotate the ID field
    output:
        vcf = progress_dir + "called/{bam_file_name}/{chrom}.vcf.gz",
        #stats = progress_dir + "called/{bam_file_name}/{chrom}.vcf.gz.stats",
    params:
        java_options='-Xmx4g -Xms4g -Djava.io.tmpdir=/storage/groups/epigenereg01/workspace/projects/vale/tmp', #Setting a tmp dir solves AVX error!
    log:
        progress_dir + 'logs/mutect2/{bam_file_name}/{chrom}.log'
    threads: 4 #should be equal to number of calls per node
    conda:
        "../envs/gatkcondaenv.yml"
    shell:
        r'''{gatk_path}  --java-options "{params.java_options}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.tumor_bam} \
            -L {input.intervals} \
            --native-pair-hmm-threads {threads} \
            --dbsnp {input.known} \
            -O {output.vcf} > {log} 2>&1
        '''

rule concat_vcfs:
    '''
    Concatenate calls from different chromosomes of the same sample
    '''
    input:
        vcf = expand(progress_dir + "called/{{bam_file_name}}/{chrom}.vcf.gz", chrom = chroms),
    output:
        vcf = temp(progress_dir + "called_concat/{bam_file_name}.vcf.gz"),
        tbi = temp(progress_dir + "called_concat/{bam_file_name}.vcf.gz.tbi"),
    shell:
        '''
        bcftools concat --no-version {input.vcf}  --threads 4  | bcftools sort - |bgzip -c  > {output.vcf}
        tabix -f  {output.vcf}
        '''

rule merge_vcfs:
    '''
    Merge vcfs from all calls
    '''
    input:
        vcf = expand(progress_dir + "called_concat/{bam_file_name}.vcf.gz", bam_file_name = all_bams),
        tbi = expand(progress_dir + "called_concat/{bam_file_name}.vcf.gz.tbi", bam_file_name = all_bams),
    output:
        vcf = progress_dir + "merged.vcf.gz",
        tbi = progress_dir + "merged.vcf.gz.tbi",
    shell:
        '''
        bcftools merge {input.vcf}  -Oz -o {output.vcf}
        tabix -f  {output.vcf}
        '''
