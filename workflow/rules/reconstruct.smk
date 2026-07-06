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
    """Symlinks the reconstruct reference genome into genomes/ under its original
    basename (extension incl. .gz preserved) so busco_metaeuk finds it and treats it
    like a native genome assembly."""
    wildcard_constraints:
        ref_basename="|".join(re.escape(b) for b in ref_link_basenames) if ref_link_basenames else "^$"
    input:
        lambda wc: ref_link_source_by_basename[wc.ref_basename],
    output:
        genome_dir_path / "{ref_basename}",
    log:
        std=log_dir_path / "link_ref_to_genomes_dir.{ref_basename}.log",
        cluster_log=cluster_log_dir_path / "link_ref_to_genomes_dir.{ref_basename}.cluster.log",
        cluster_err=cluster_log_dir_path / "link_ref_to_genomes_dir.{ref_basename}.cluster.err",
    shell:
        " ln -sr {input} {output} 1> {log.std} 2>&1; "


rule prepare_reconstruct_consensus:
    """Normalizes a consensus genome to a plain, faidx-indexed FASTA in results/
    (input/ is read-only, so the .fai cannot be written next to the source)."""
    wildcard_constraints:
        species="|".join(consensus_species) if consensus_species else "^$"
    input:
        lambda wc: reconstruct_map[wc.species]["consensus_fasta"],
    output:
        fasta=reconstruct_consensus_dir_path / "{species}.fasta",
        fai=reconstruct_consensus_dir_path / "{species}.fasta.fai",
    log:
        std=log_dir_path / "prepare_reconstruct_consensus.{species}.log",
        cluster_log=cluster_log_dir_path / "prepare_reconstruct_consensus.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "prepare_reconstruct_consensus.{species}.cluster.err",
    benchmark:
        benchmark_dir_path / "prepare_reconstruct_consensus.{species}.benchmark.txt"
    conda:
        main_env
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " case \"{input}\" in "
        "   *.gz) zcat {input} > {output.fasta} 2> {log.std} ;; "
        "   *) ln -sr {input} {output.fasta} 2> {log.std} ;; "
        " esac; "
        " samtools faidx {output.fasta} 2>> {log.std}; "


rule apply_vcf_to_busco:
    wildcard_constraints:
        species="|".join(vcf_species) if vcf_species else "^$"
    input:
        single_copy_busco_sequences=lambda wc: expand(
            rules.busco_metaeuk.output.single_copy_busco_sequences,
            species=reconstruct_map[wc.species]["ref_prefix"]
        )[0],
        metaeuk_output=lambda wc: expand(
            rules.busco_metaeuk.output.metaeuk_output,
            species=reconstruct_map[wc.species]["ref_prefix"]
        )[0],
        busco_summary=lambda wc: expand(
            rules.busco_metaeuk.output.summary,
            species=reconstruct_map[wc.species]["ref_prefix"]
        )[0],
        vcf=lambda wc: reconstruct_map[wc.species]["vcf"],
        vcf_tbi=lambda wc: str(reconstruct_map[wc.species]["vcf"]) + ".tbi",
    output:
        seqs=directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences"),
        multi_copy=directory(busco_dir_path / "{species}/busco_sequences/multi_copy_busco_sequences"),
        summary=busco_dir_path / "{species}/short_summary_{species}.txt",
    params:
        sample=lambda wc: reconstruct_map[wc.species]["sample"],
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


rule apply_consensus_to_busco:
    wildcard_constraints:
        species="|".join(consensus_species) if consensus_species else "^$"
    input:
        single_copy_busco_sequences=lambda wc: expand(
            rules.busco_metaeuk.output.single_copy_busco_sequences,
            species=reconstruct_map[wc.species]["ref_prefix"]
        )[0],
        metaeuk_output=lambda wc: expand(
            rules.busco_metaeuk.output.metaeuk_output,
            species=reconstruct_map[wc.species]["ref_prefix"]
        )[0],
        busco_summary=lambda wc: expand(
            rules.busco_metaeuk.output.summary,
            species=reconstruct_map[wc.species]["ref_prefix"]
        )[0],
        consensus_fasta=reconstruct_consensus_dir_path / "{species}.fasta",
        consensus_fai=reconstruct_consensus_dir_path / "{species}.fasta.fai",
    output:
        seqs=directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences"),
        multi_copy=directory(busco_dir_path / "{species}/busco_sequences/multi_copy_busco_sequences"),
        summary=busco_dir_path / "{species}/short_summary_{species}.txt",
    log:
        std=log_dir_path / "apply_consensus_to_busco.{species}.log",
        cluster_log=cluster_log_dir_path / "apply_consensus_to_busco.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "apply_consensus_to_busco.{species}.cluster.err",
    benchmark:
        benchmark_dir_path / "apply_consensus_to_busco.{species}.benchmark.txt"
    conda:
        main_env
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/apply_consensus_to_busco.py "
        " --single_copy_busco_sequences {input.single_copy_busco_sequences} "
        " --metaeuk_output {input.metaeuk_output} "
        " --consensus {input.consensus_fasta} "
        " --output_dir {output.seqs} "
        " 1> {log.std} 2>&1; "
        " mkdir -p {output.multi_copy}; "
        " ln -r -s {input.busco_summary} {output.summary}; "
