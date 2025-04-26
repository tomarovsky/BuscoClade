localrules:
    link_altref_to_genomes_dir,


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
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
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
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " {params.gatk_path}/gatk --java-options -Xmx{resources.mem_mb}m IndexFeatureFile -I {input}"


rule gatk_altref:
    input:
        ref=lambda wc: vcf_reconstruct_map[wc.species]["reference"],
        refidx=lambda wc: vcf_reconstruct_map[wc.species]["reference"].with_suffix(".dict"),
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
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
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
