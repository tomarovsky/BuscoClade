localrules:
    tabix_index,
    link_ref_to_genomes_dir,


rule tabix_index:
    input:
        "{vcf}.vcf.gz",
    output:
        "{vcf}.vcf.gz.tbi",
    log:
        std=log_dir_path / "tabix_index.{vcf}.log",
        cluster_log=cluster_log_dir_path / "tabix_index.{vcf}.cluster.log",
        cluster_err=cluster_log_dir_path / "tabix_index.{vcf}.cluster.err",
    benchmark:
        benchmark_dir_path / "tabix_index.{vcf}.benchmark.txt"
    conda:
        main_env
    shell:
        " tabix -p vcf {input} 1> {log.std} 2>&1; "


rule link_ref_to_genomes_dir:
    """Creates a symlink for the VCF reference genome in genomes/ so busco_metaeuk can find it."""
    wildcard_constraints:
        ref_prefix="|".join(altref_refs) if altref_refs else "^$"
    input:
        lambda wc: next(
            v["reference"] for v in altref_map.values()
            if v["ref_prefix"] == wc.ref_prefix
        ),
    output:
        genome_dir_path / "{ref_prefix}.fasta",
    log:
        std=log_dir_path / "link_ref_to_genomes_dir.{ref_prefix}.log",
        cluster_log=cluster_log_dir_path / "link_ref_to_genomes_dir.{ref_prefix}.cluster.log",
        cluster_err=cluster_log_dir_path / "link_ref_to_genomes_dir.{ref_prefix}.cluster.err",
    shell:
        " ln -sr {input} {output} 1> {log.std} 2>&1; "


rule apply_vcf_to_busco:
    wildcard_constraints:
        species="|".join(altref_map.keys()) if altref_map else "^$"
    input:
        single_copy_busco_sequences=lambda wc: expand(
            rules.busco_metaeuk.output.single_copy_busco_sequences,
            species=altref_map[wc.species]["ref_prefix"]
        )[0],
        metaeuk_output=lambda wc: expand(
            rules.busco_metaeuk.output.metaeuk_output,
            species=altref_map[wc.species]["ref_prefix"]
        )[0],
        busco_summary=lambda wc: expand(
            rules.busco_metaeuk.output.summary,
            species=altref_map[wc.species]["ref_prefix"]
        )[0],
        vcf=lambda wc: altref_map[wc.species]["vcf"],
        vcf_tbi=lambda wc: str(altref_map[wc.species]["vcf"]) + ".tbi",
    output:
        seqs=directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences"),
        multi_copy=directory(busco_dir_path / "{species}/busco_sequences/multi_copy_busco_sequences"),
        summary=busco_dir_path / "{species}/short_summary_{species}.txt",
    params:
        sample=lambda wc: altref_map[wc.species]["vcf_id"],
        iupac="--iupac" if config.get("apply_vcf_iupac") else "",
    log:
        std=log_dir_path / "apply_vcf_to_busco.{species}.log",
        cluster_log=cluster_log_dir_path / "apply_vcf_to_busco.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "apply_vcf_to_busco.{species}.cluster.err",
    benchmark:
        benchmark_dir_path / "apply_vcf_to_busco.{species}.benchmark.txt"
    conda:
        main_env
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/apply_vcf_to_busco.py "
        " --single_copy_busco_sequences {input.single_copy_busco_sequences} "
        " --metaeuk_output {input.metaeuk_output} "
        " --vcf {input.vcf} "
        " --sample {params.sample} "
        " --output_dir {output.seqs} "
        " {params.iupac} "
        " 1> {log.std} 2>&1; "
        " mkdir -p {output.multi_copy}; "
        " ln -r -s {input.busco_summary} {output.summary}; "
