rule concat_snvs_indels_and_split:
    '''
    Concatenate SNVs and INDELs
    '''
    input:
        vcf =  expand(filtered_dir + "{vartype}.recalibrated.vcf.gz", vartype=['snvs','indels']),
    output:
        vcf = filtered_dir + "all.vcf.gz"
    shell:
        r'''
        bcftools concat -a {input.vcf} -Oz -o {output.vcf};
        '''

rule split_by_sample:
    '''
    Split by BAM sample name
    '''    
    input:
        vcf = filtered_dir + "all.vcf.gz"
    output:
        placeholder = filtered_dir + "by_sample/splitted.ok" #we use an empty placeholder as output since the BAM sample names are not known a priori
    params:
        progress_dir = filtered_dir + "by_sample"
    shell:
        r'''
        echo splitting_placeholder > {output}
        for bam_sample_name in $(bcftools query -l {input.vcf}); do
            bcftools view -c 1 -s $bam_sample_name {input.vcf} -Oz -o {params.progress_dir}/$bam_sample_name.vcf.gz || rm {output} #remove placeholder if error: this will trigger snakemake error
            tabix -f {params.progress_dir}/$bam_sample_name.vcf.gz
        done
        '''
