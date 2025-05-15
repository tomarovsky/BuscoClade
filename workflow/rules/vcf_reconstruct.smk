localrules:
    link_altref_to_genomes_dir,

if config["vcf_gatk"]:
    rule picard_index:
        input:
            "{reference}.fasta",
        output:
            dict="{reference}.dict",
            fai="{reference}.fasta.fai",
        log:
            std=log_dir_path / "picard_index.{reference}.log",
            cluster_log=cluster_log_dir_path / "picard_index.{reference}.cluster.log",
            cluster_err=cluster_log_dir_path / "picard_index.{reference}.cluster.err",
        benchmark:
            benchmark_dir_path / "picard_index.{reference}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"]),
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "picard CreateSequenceDictionary -R {input} && "
            "samtools faidx {input} "
            "1> {log.std} 2>&1" 



rule gatk_vcf_index:
    input:
        vcf="{vcf_dir}/{vcf}",
    output:
        tbi="{vcf_dir}/{vcf}.tbi",
        sample_list=temp("{vcf_dir}/{vcf}.sample_list.txt"),
    params:
        gatk_path=config["gatk_path"],
    log:
        std_gatk= "{vcf_dir}/gatk_vcf_index.{vcf}.log",
        std_bcf= "{vcf_dir}/gatk_vcf_index.{vcf}.bcf.log",
        std_grep= "{vcf_dir}/gatk_vcf_index.{vcf}.grep.log",
        cluster_log="{vcf_dir}/gatk_vcf_index.{vcf}.cluster.log",
        cluster_err="{vcf_dir}/gatk_vcf_index.{vcf}.cluster.err",
    benchmark:
        "{vcf_dir}/gatk_vcf_index.{vcf}.benchmark.txt"
    conda:
        config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        """
        {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m IndexFeatureFile -I {input.vcf} > {log.std_gatk} 2>&1; 
        bcftools query -l {input.vcf} 2>{log.std_bcf}|
        grep -v '^$' > {output.sample_list} 2>{log.std_grep}; 
        """

checkpoint vcf_separation:
    input:
        vcf="{vcf_dir}/{vcf}",
        tbi="{vcf_dir}/{vcf}.tbi",
        sample_list="{vcf_dir}/{vcf}.sample_list.txt",
    output:
        separated="{vcf_dir}/{vcf}.separation.done",
    params:
        gatk_path=config["gatk_path"],
        vcf_basename=lambda wildcards: os.path.basename(wildcards.vcf).replace('.vcf.gz', ''),
        output_dir=lambda wildcards: os.path.dirname(wildcards.vcf),
        original_vcf_dir="results/original_vcf",
    log:
        std="{vcf_dir}/vcf_separation.{vcf}.log",
        std_bcf_view="{vcf_dir}/vcf_separation.{vcf}.bcf_view.log",
        std_bcf_filter="{vcf_dir}/vcf_separation.{vcf}.bcf_filter.log",
        std_zip="{vcf_dir}/vcf_separation.{vcf}.zip.log",
        std_tabix="{vcf_dir}/vcf_separation.{vcf}.tabix.log",
        std_iff="{vcf_dir}/vcf_separation.{vcf}.iff.log",
        cluster_log="{vcf_dir}/vcf_separation.{vcf}.cluster.log",
        cluster_err="{vcf_dir}/vcf_separation.{vcf}.cluster.err",
    benchmark:
        "{vcf_dir}/vcf_separation.{vcf}.benchmark.txt"
    conda:
        config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        r"""
        SAMPLES=($(grep -v '^$' {input.sample_list}))

        for SAMPLE in "${{SAMPLES[@]}}"; do

            bcftools view \
                --threads {resources.cpus} \
                --min-ac 1 \
                -s "$SAMPLE" \
                -Ov -o "{params.output_dir}/$SAMPLE.{params.vcf_basename}.unfiltered.vcf" \
                {input.vcf} > {log.std_bcf_view} 2>&1

            bcftools filter \
                -e 'ALT ~ "[^ATCGN]"' \
                "{params.output_dir}/$SAMPLE.{params.vcf_basename}.unfiltered.vcf" \
                > "{params.output_dir}/$SAMPLE.{params.vcf_basename}.filtered.vcf" 2>{log.std_bcf_filter}

            bgzip -c "{params.output_dir}/$SAMPLE.{params.vcf_basename}.filtered.vcf" \
                > "{params.output_dir}/$SAMPLE.{params.vcf_basename}.vcf.gz"  2>{log.std_zip}

            tabix -p vcf "{params.output_dir}/$SAMPLE.{params.vcf_basename}.vcf.gz" > {log.std_tabix} 2>&1

            {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m \
                IndexFeatureFile -I "{params.output_dir}/$SAMPLE.{params.vcf_basename}.vcf.gz" >{log.std_iff} 2>&1
            
            rm -f \
                "{params.output_dir}/$SAMPLE.{params.vcf_basename}.unfiltered.vcf" \
                "{params.output_dir}/$SAMPLE.{params.vcf_basename}.filtered.vcf"

        done

        touch {output.separated}
        """

if config["vcf_gatk"]:
    rule gatk_altref:
        input:
            ref=lambda wc: vcf_reconstruct_map[wc.species]["reference"],
            refidx=lambda wc: vcf_reconstruct_map[wc.species]["reference"].with_suffix(".dict"),
            vcf=lambda wc: checkpoints.vcf_separation.get(vcf=Path(vcf_reconstruct_map[wc.species]["vcf"]).name,
                                                          vcf_dir=Path(vcf_reconstruct_map[wc.species]["vcf"]).parent).output.separated,
            vcfidx=lambda wc: vcf_reconstruct_map[wc.species]["vcf"].with_name(vcf_reconstruct_map[wc.species]["vcf"].name + ".tbi"),
        output:
            altref_dir_path / "{species}.fasta",
        params:
            gatk_path=config["gatk_path"],
            sample=lambda wc: wc.species.split(".")[0],
        log:
            std=log_dir_path / "gatk_altref.{species}.log",
            cluster_log=cluster_log_dir_path / "gatk_altref.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "gatk_altref.{species}.cluster.err",
        benchmark:
            benchmark_dir_path / "gatk_altref.{species}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            " {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m FastaAlternateReferenceMaker "
            " --output {output} --reference {input.ref} --variant {input.vcf} --showHidden true --use-iupac-sample {params.sample} 1> {log.std} 2>&1; "




    checkpoint move_altref_to_genomes_dir:
        input:
            altref_dir_path / "{species}.fasta",
        output:
            genome_dir_path / "{species}.fasta",
        log:
            std=log_dir_path / "link_altref_to_genomes_dir.{species}.log",
            cluster_log=cluster_log_dir_path / "link_altref_to_genomes_dir.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "link_altref_to_genomes_dir.{species}.cluster.err",
        shell:
            " cp -a {input} {output} 1> {log.std} 2>&1; "

if config["vcf2phylip"]:
    from pathlib import Path, PurePath

    vcf_dir = Path(config["vcf_reconstruct_dir"])
    vcf_list = list(vcf_dir.rglob("*.vcf.gz"))
    if len(vcf_list) != 1:
        raise ValueError(f"vcf2phylip requires exactly one VCF in {vcf_dir}, found: {vcf_list}")
    vcf_file = vcf_list[0]
    prefix   = vcf_file.stem

    checkpoint vcf2phylip:
        input:
            vcf=str(vcf_file),
        output:
            altref_dir_path / f"{prefix}.fasta",
        params:
            prefix=prefix,
        log:
            std=log_dir_path / f"vcf2phylip.{prefix}.log",
            cluster_log=cluster_log_dir_path / f"vcf2phylip.{prefix}.cluster.log",
            cluster_err=cluster_log_dir_path / f"vcf2phylip.{prefix}.cluster.err",
        benchmark:
            benchmark_dir_path / f"vcf2phylip.{prefix}.benchmark.txt",
        conda:
            config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"]
            else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"]),
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "workflow/scripts/vcf2phylip.py -i {input.vcf} --fasta --output-prefix {params.prefix} --resolve-IUPAC"

    rule consolidate_vcf2phylip:
        input:
            fasta=altref_dir_path / f"{prefix}.fasta",
        output:
            altref_dir_path / "consensus.fasta",
        log:
            std=log_dir_path / "consolidate_vcf2phylip.log",
        run:
            shell("mv {input.fasta} {output}")
