#filter calls using hard filtering or VQSR
#choose filtering strategy in config.yaml

#for hard filtering see
#https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset,
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

#for VQSR see
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
#https://gatk.broadinstitute.org/hc/en-us/articles/360036712971-ApplyVQSR
#https://gatk.broadinstitute.org/hc/en-us/articles/360036727711-VariantRecalibrator#--trust-all-polymorphic
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

#VQSR: Minimum 30 exomes should be used (with exome size is 30M bases, this amounts to 900M bases)


rule select_calls:
    input:
        ref = refgen_fa,
        vcf = progress_dir + "merged.vcf.gz",
        tbi = progress_dir + "merged.vcf.gz.tbi",
    output:
        vcf = temp(progress_dir + "{vartype}.vcf.gz"),
    params:
        extra = lambda w: "--select-type-to-include {} --restrict-alleles-to BIALLELIC".format(
            "SNP" if w.vartype == "snvs" else "INDEL"
        ), #remove bialleic variants at this stage
    log:
        logs_dir + "gatk/selectvariants/{vartype}.log",
    conda:
        "../envs/gatkcondaenv.yml"
    shell:
        r'''{gatk_path} SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        {params.extra} \
        -O {output.vcf} > {log} 2>&1
        '''

if not config["filtering"]["type"]=="vqsr":

    rule hard_filter_calls:
        input:
            ref = refgen_fa,
            vcf = progress_dir + "{vartype}.vcf.gz",
        output:
            vcf = temp(filtered_dir + "{vartype}.hardfiltered.vcf.gz"),
        params:
            filters = " ".join([f'--filter-expression "{e}" --filter-name "{n}"' for n,e in config['filtering']['hard'].items()]),
        log:
            logs_dir + "gatk/variantfiltration/{vartype}.log",
        conda:
            "../envs/gatkcondaenv.yml"
        shell:
            r'''{gatk_path} --java-options "-Xmx3g -Xms3g" VariantFiltration \
             -R {input.ref} \
             -V {input.vcf} \
             {params.filters} \
             -O {output.vcf} > {log} 2>&1
             '''

else:

    #VQSR pipeline:
    #https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
    #Attention to different parameters for WXS and WGS samples, snps and indels

    rule excess_het:
        input:
            progress_dir + "{vartype}.vcf.gz",
        output:
            temp(filtered_dir + "{vartype}.excesshet.vcf.gz"),
        log:
            logs_dir + "gatk/excesshet/{vartype}.log"
        conda:
            "../envs/gatkcondaenv.yml"
        shell:
            r'''{gatk_path} --java-options "-Xmx3g -Xms3g" VariantFiltration \
             -V {input} \
             --filter-expression "ExcessHet > 54.69" \
             --filter-name ExcessHet \
             -O {output} > {log} 2>&1
             '''

    rule sites_only:
        input:
            filtered_dir + "{vartype}.excesshet.vcf.gz",
        output:
            temp(filtered_dir + "{vartype}.sitesonly.vcf.gz"),
        log:
            logs_dir + "gatk/sitesonly/{vartype}.log"
        conda:
            "../envs/gatkcondaenv.yml"
        shell:
            r'''{gatk_path} MakeSitesOnlyVcf \
             -I {input} \
             -O {output} > {log} 2>&1
            '''


    rule recalibrate_snv_calls:
        input:
            ref = refgen_fa,
            vcf = filtered_dir + "snvs.sitesonly.vcf.gz",
            hapmap = resources['hapmap_vcf'],
            omni = resources['omni_vcf'],
            g1k = resources['g1k_vcf'],
            dbsnp = resources['dbsnp_all_vcf'],
        output:
            vcf = filtered_dir + "snvs.recal",
            tranches = filtered_dir + "snvs.tranches",
            plots = filtered_dir + "snvs.plots.R",
        params:
            extra="--trust-all-polymorphic --max-gaussians 6 \
                   -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
                   -tranche 95.0 -tranche 90.0 \
                   --target-titv 2.15",#titv is different for WGS and WXS

            java_options='-Xmx24g -Xms24g', #this memory should also be allocated in cluster.yaml file
            mode="SNP",
            resources = {"hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
                   "omni":   {"known": False, "training": True, "truth": False, "prior": 12.0},
                   "g1k":   {"known": False, "training": True, "truth": False, "prior": 10.0},
                   "dbsnp":  {"known": True, "training": False, "truth": False, "prior": 2.0}},
            annotation = (
                    ' '.join(f"-an {an}" for an in ["QD", "FS", "ReadPosRankSum", "MQ", "MQRankSum", "SOR"])
                    if config['exp_strategy']=='WXS' #don't use DP for exome data
                    else ' '.join(f"-an {an}" for an in ["QD", "FS", "DP", "ReadPosRankSum", "MQ", "MQRankSum", "SOR"])
                )
        log:
            logs_dir + "gatk/variantrecalibrator/snvs.log",
        conda:
            "../envs/gatkcondaenv.yml"
        shell:
            r'''{gatk_path} --java-options "{params.java_options}" \
            VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
            --resource:g1k,known=false,training=true,truth=false,prior=10.0 {input.g1k} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
            {params.annotation} \
            {params.extra} \
            -mode {params.mode} \
            -O {output.vcf} \
            --tranches-file {output.tranches} \
            --rscript-file {output.plots} > {log} 2>&1
            '''


    rule recalibrate_indel_calls:
        input:
            ref = refgen_fa,
            vcf = filtered_dir + "indels.sitesonly.vcf.gz",
            dbsnp = resources['dbsnp_all_vcf'],
            mills = resources['mills_vcf'],
            axexome = resources['axexome_vcf'],
        output:
            vcf = filtered_dir + "indels.recal",
            tranches = filtered_dir + "indels.tranches",
            plots = filtered_dir + "indels.plots.R",
        params:
            extra="--trust-all-polymorphic --max-gaussians 4 \
                   -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 \
                   -tranche 95.0 -tranche 90.0 \
                   --target-titv 2.15",#titv is different for WGS and WXS
            java_options = '-Xmx24g -Xms24g', #this memory should also be allocated in cluster.yaml file
            mode = "INDEL",
            resources = {"mills": {"known": False, "training": True, "truth": True, "prior": 12.0},
                   "axexome":   {"known": False, "training": True, "truth": False, "prior": 10.0},
                   "dbsnp":  {"known": True, "training": False, "truth": False, "prior": 2.0}},
            annotation = (
                    ' '.join(f"-an {an}" for an in ["QD", "FS", "ReadPosRankSum", "MQRankSum", "SOR"])
                    if config['exp_strategy']=='WXS' #don't use DP for exome data
                    else ' '.join(f"-an {an}" for an in ["QD", "FS", "DP", "ReadPosRankSum", "MQRankSum", "SOR"])
                )
        log:
            logs_dir + "gatk/variantrecalibrator/indels.log",
        conda:
            "../envs/gatkcondaenv.yml"
        shell:
            r'''{gatk_path} --java-options "{params.java_options}" \
            VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
            --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
            --resource:axexome,known=false,training=true,truth=false,prior=10.0 {input.axexome} \
            {params.annotation} \
            {params.extra} \
            -mode {params.mode} \
            -O {output.vcf} \
            --tranches-file {output.tranches} \
            --rscript-file {output.plots} > {log} 2>&1
            '''


    rule apply_vqsr:
        input:
            ref = refgen_fa,
            vcf = filtered_dir + "{vartype}.excesshet.vcf.gz",
            recal = filtered_dir + "{vartype}.recal",
            tranches = filtered_dir + "{vartype}.tranches",
        output:
            filtered_dir + "{vartype}.recalibrated.vcf.gz",
        params:
            java_options = '-Xmx5g -Xms5g',#this memory should also be allocated in cluster.yaml file
            mode = lambda w:"SNP" if w.vartype == "snvs" else "INDEL",
            exclude_filtered = config["filtering"]["vqsr"]["exclude_filtered"],
            truth_sensitivity_filter_level = config["filtering"]["vqsr"]["truth_sensitivity_filter_level"],
        log:
            logs_dir + "gatk/applyvqsr/{vartype}.log",
        conda:
            "../envs/gatkcondaenv.yml"
        shell:
            r'''{gatk_path} --java-options "{params.java_options}" \
            ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --exclude-filtered {params.exclude_filtered} \
            --recal-file {input.recal} \
            --tranches-file {input.tranches} \
            --truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} \
            --create-output-variant-index true \
            -mode {params.mode} \
            -O {output} > {log} 2>&1
            '''
