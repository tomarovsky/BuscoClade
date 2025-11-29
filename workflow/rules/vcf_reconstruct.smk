localrules:
    link_altref_to_genomes_dir,


rule samtools_index:
    input:
        "{reference}.fasta",
    output:
        "{reference}.fasta.fai",
    log:
        std=log_dir_path / "samtools_index.{reference}.log",
        cluster_log=cluster_log_dir_path / "samtools_index.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "samtools_index.{reference}.cluster.err",
    benchmark:
        benchmark_dir_path / "samtools_index.{reference}.benchmark.txt"
    conda:
        config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"],
    shell:
        " samtools faidx {input} 1> {log.std} 2>&1; "


rule picard_index:
    input:
        "{reference}.fasta",
    output:
        "{reference}.dict",
    log:
        std=log_dir_path / "picard_index.{reference}.log",
        cluster_log=cluster_log_dir_path / "picard_index.{reference}.cluster.log",
        cluster_err=cluster_log_dir_path / "picard_index.{reference}.cluster.err",
    benchmark:
        benchmark_dir_path / "picard_index.{reference}.benchmark.txt"
    conda:
        config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"],
    shell:
        " picard CreateSequenceDictionary -R {input} 1> {log.std} 2>&1; "


rule gatk_vcf_index:
    input:
        "{vcf}",
    output:
        "{vcf}.tbi",
    params:
        gatk_path=config["gatk_path"],
    log:
        std=log_dir_path / "gatk_vcf_index.{vcf}.log",
        cluster_log=cluster_log_dir_path / "gatk_vcf_index.{vcf}.cluster.log",
        cluster_err=cluster_log_dir_path / "gatk_vcf_index.{vcf}.cluster.err",
    benchmark:
        benchmark_dir_path / "gatk_vcf_index.{vcf}.benchmark.txt"
    conda:
        config["conda"]["buscoclade_gatk"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_gatk"]["yaml"])
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"],
    shell:
        " {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m IndexFeatureFile -I {input}"


rule gatk_altref:
    input:
        ref=lambda wc: vcf_reconstruct_map[wc.species]["reference"],
        refidx_picard=lambda wc: vcf_reconstruct_map[wc.species]["reference"].with_suffix(".dict"),
        refidx_samtools=lambda wc: vcf_reconstruct_map[wc.species]["reference"].with_suffix(".fasta.fai"),
        vcf=lambda wc: vcf_reconstruct_map[wc.species]["vcf"],
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
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"],
    shell:
        " {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m FastaAlternateReferenceMaker "
        " --output {output} --reference {input.ref} --variant {input.vcf} --showHidden true --use-iupac-sample {params.sample} 1> {log.std} 2>&1; "


rule link_altref_to_genomes_dir:
    input:
        altref_dir_path / "{species}.fasta",
    output:
        genome_dir_path / "{species}.fasta",
    log:
        std=log_dir_path / "link_altref_to_genomes_dir.{species}.log",
        cluster_log=cluster_log_dir_path / "link_altref_to_genomes_dir.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "link_altref_to_genomes_dir.{species}.cluster.err",
    shell:
        " ln -sr {input} {output} 1> {log.std} 2>&1; "


if config.get("vcf2phylip"):

    rule vcf2phylip:
        input:
            vcf=str(vcf_file),
        output:
            concat_alignments_dir_path / fasta_filename,
        params:
            prefix=f"{prefix}.min4.fasta",
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
            slurm_partition=config["processing_queue"],
            runtime=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"],
        shell:
            " workflow/scripts/vcf2phylip.py -i {input.vcf} --phylip-disable --fasta 1> {log.std} 2>&1; "
            " mv {params.prefix} {output}"
