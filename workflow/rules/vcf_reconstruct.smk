localrules:
    link_altref_to_genomes_dir,

if config["vcf_gatk"]:
    rule picard_index:
        input:
            "{reference_dir}/{reference}.fasta",
        output:
            dict="{reference_dir}/{reference}.dict",
            fai="{reference_dir}/{reference}.fasta.fai",
        log:
            std= "{reference_dir}/picard_index.{reference}.log",
            cluster_log= "{reference_dir}/picard_index.{reference}.cluster.log",
            cluster_err= "{reference_dir}/picard_index.{reference}.cluster.err",
        benchmark:
            "{reference_dir}/picard_index.{reference}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"]),
        resources:
            queue=config["vcf_processing_queue"],
            cpus=config["vcf_processing_threads"],
            time=config["vcf_processing_time"],
            mem_mb=config["vcf_processing_mem_mb"],
        shell:
            "picard CreateSequenceDictionary -R {input} && "
            "samtools faidx {input} "
            "1> {log.std} 2>&1" 


    rule gatk_vcf_index:
        input:
            vcf="{vcf_dir}/{vcf}",
        output:
            tbi="{vcf_dir}/{vcf}.tbi",
        params:
            gatk_path=config["gatk_path"],
        log:
            std_gatk= "{vcf_dir}/gatk_vcf_index.{vcf}.log",
            cluster_log="{vcf_dir}/gatk_vcf_index.{vcf}.cluster.log",
            cluster_err="{vcf_dir}/gatk_vcf_index.{vcf}.cluster.err",
        benchmark:
            "{vcf_dir}/gatk_vcf_index.{vcf}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
        resources:
            queue=config["vcf_processing_queue"],
            cpus=config["vcf_processing_threads"],
            time=config["vcf_processing_time"],
            mem_mb=config["vcf_processing_mem_mb"],
        shell:
            """
            {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m IndexFeatureFile -I {input.vcf} > {log.std_gatk} 2>&1; 
            """

    

    rule vcf_separation:
        input:
            vcf="{vcf_dir}/{vcf_prefix}.vcf.gz",
        output:
            separated_vcf="{vcf_dir}/tmp/{sample}.{vcf_prefix}.vcf.gz",
        params:
            gatk_path=config["gatk_path"],
        log:
            std="{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.log",
            std_bcf_view="{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.bcf_view.log",
            std_bcf_filter="{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.bcf_filter.log",
            std_zip="{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.zip.log",
            cluster_log="{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.cluster.log",
            cluster_err="{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.cluster.err",
        benchmark:
            "{vcf_dir}/vcf_separation.{sample}.{vcf_prefix}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
        resources:
            queue=config["vcf_separation_queue"],
            cpus=config["vcf_separation_threads"],
            time=config["vcf_separation_time"],
            mem_mb=config["vcf_separation_mem_mb"],
        shell:
            """
            bcftools view --threads {resources.cpus} --min-ac 1 -s {wildcards.sample} -Ov {input.vcf} 2> {log.std_bcf_view} |\
            bcftools filter -e 'ALT ~ "[^ATCGN]"' 2>{log.std_bcf_filter} |\
            bgzip -c > "{wildcards.vcf_dir}/tmp/{wildcards.sample}.{wildcards.vcf_prefix}.vcf.gz"  2>{log.std_zip}
            """

    rule gatk_altref:
        input:
            ref=lambda wc: "{0}/{1}.fasta".format(vcf_species_location_dict["{0}.{1}.{2}".format(wc.sample, wc.vcf_prefix, wc.reference)], wc.reference),
            refidx=lambda wc: "{0}/{1}.dict".format(vcf_species_location_dict["{0}.{1}.{2}".format(wc.sample, wc.vcf_prefix, wc.reference)], wc.reference),
            vcf=lambda wc: "{0}/tmp/{1}.{2}.vcf.gz".format(vcf_species_location_dict["{0}.{1}.{2}".format(wc.sample, wc.vcf_prefix, wc.reference)], wc.sample, wc.vcf_prefix),
            vcfidx=lambda wc: "{0}/tmp/{1}.{2}.vcf.gz.tbi".format(vcf_species_location_dict["{0}.{1}.{2}".format(wc.sample, wc.vcf_prefix, wc.reference)], wc.sample, wc.vcf_prefix),
        output:
            altref_dir_path / "{sample}.{vcf_prefix}.{reference}.fasta",
        params:
            gatk_path=config["gatk_path"],
        log:
            std=log_dir_path / "gatk_altref.{sample}.{vcf_prefix}.{reference}.log",
            cluster_log=cluster_log_dir_path / "gatk_altref.{sample}.{vcf_prefix}.{reference}.cluster.log",
            cluster_err=cluster_log_dir_path / "gatk_altref.{sample}.{vcf_prefix}.{reference}.cluster.err",
        benchmark:
            benchmark_dir_path / "gatk_altref.{sample}.{vcf_prefix}.{reference}.benchmark.txt",
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"]),
        resources:
            queue=config["vcf_altref_queue"],
            cpus=config["vcf_altref_threads"],
            time=config["vcf_altref_time"],
            mem_mb=config["vcf_altref_mem_mb"],
        shell:
            " {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m FastaAlternateReferenceMaker "
            " --output {output} --reference {input.ref} --variant {input.vcf} --showHidden true --use-iupac-sample {wildcards.sample} 1> {log.std} 2>&1; "


    rule move_altref_to_genomes_dir:
        input:
            altref_dir_path / "{species}.fasta",
        output:
            genome_dir_path / "{species}.fasta",
        log:
            std=log_dir_path / "link_altref_to_genomes_dir.{species}.log",
            cluster_log=cluster_log_dir_path / "link_altref_to_genomes_dir.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "link_altref_to_genomes_dir.{species}.cluster.err",
        resources:
            queue=config["vcf_separation_queue"],
            cpus=config["vcf_separation_threads"],
            time=config["vcf_separation_time"],
            mem_mb=config["vcf_separation_mem_mb"],
        shell:
            " cp -a {input} {output} 1> {log.std} 2>&1; "

if config["vcf2phylip"]:
    from pathlib import Path, PurePath

    vcf_dir = Path(config["vcf_reconstruct_dir"])
    vcf_list = list(vcf_dir.rglob("*.vcf.gz"))
    if len(vcf_list) != 1:
        raise ValueError(f"vcf2phylip requires exactly one VCF in {vcf_dir}, found: {vcf_list}")
    vcf_file = vcf_list[0]
    prefix   = vcf_file.name.replace(".vcf.gz", "")

    checkpoint vcf2phylip:
        input:
            vcf=str(vcf_file),
        output:
            f"{prefix}.min4.fasta",
        log:
            std=log_dir_path / f"vcf2phylip.{prefix}.log",
            cluster_log=cluster_log_dir_path / f"vcf2phylip.{prefix}.cluster.log",
            cluster_err=cluster_log_dir_path / f"vcf2phylip.{prefix}.cluster.err",
        benchmark:
            benchmark_dir_path / f"vcf2phylip.{prefix}.benchmark.txt",
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"]
            else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"]),
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "workflow/scripts/vcf2phylip.py -i {input.vcf} --fasta "

    rule consolidate_vcf2phylip:
        input:
            fasta=f"{prefix}.min4.fasta",
        output:
            altref_dir_path / "consensus.fasta",
        log:
            std=log_dir_path / "consolidate_vcf2phylip.log",
            cluster_log=cluster_log_dir_path / "consolidate_vcf2phylip.cluster.log",
            cluster_err=cluster_log_dir_path / "consolidate_vcf2phylip.cluster.err",
        benchmark:
            benchmark_dir_path / f"consolidate_vcf2phylip.benchmark.txt",
        conda:
            config["conda"]["buscoclade_main"]["name"]
            if config["use_existing_envs"]
            else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"]),
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "mv {input.fasta} {output}"
