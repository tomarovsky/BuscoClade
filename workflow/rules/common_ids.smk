localrules:
    species_single_copy_ids,
    species_multi_copy_ids,
    common_ids,
    gene_count_report,


rule species_single_copy_ids:  # get single copy ids for each species
    input:
        busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences",
    output:
        species_ids_dir_path / "single_copy/{species}.ids",
    log:
        std=log_dir_path / "species_ids.{species}.log",
        cluster_log=cluster_log_dir_path / "species_ids.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "species_ids.{species}.cluster.err",
    benchmark:
        benchmark_dir_path / "species_ids.{species}.benchmark.txt"
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " ls {input} | grep -P '.fna$' | sed 's/.fna//' > {output} 2> {log.std}; "


rule species_multi_copy_ids:  # get multi copy ids for each species
    input:
        busco_dir_path / "{species}/busco_sequences/multi_copy_busco_sequences",
    output:
        species_ids_dir_path / "multi_copy/{species}.ids",
    log:
        std=log_dir_path / "species_ids.{species}.log",
        cluster_log=cluster_log_dir_path / "species_ids.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "species_ids.{species}.cluster.err",
    benchmark:
        benchmark_dir_path / "species_ids.{species}.benchmark.txt"
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " ls {input}/*.fna 2>/dev/null | xargs -I{{}} basename {{}} .fna > {output} 2> {log.std} || touch {output}; "

rule common_ids:  # get common IDs for all species given config["gene_blacklist"] and split them into files
    input:
        expand(species_ids_dir_path / "single_copy/{species}.ids", species=config["species_list"]),
    output:
        common_ids_dir_path / "common.ids",
    params:
        nfiles=len(config["species_list"]),
    log:
        std=log_dir_path / "common_ids.log",
        cluster_log=cluster_log_dir_path / "common_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "common_ids.cluster.err",
    benchmark:
        benchmark_dir_path / "common_ids.benchmark.txt"
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " cat {input} | sort | uniq -c | awk '{{if($1=={params.nfiles}){{print $2}}}}' > {output} 2> {log.std}; "


rule gene_count_report:  # per-species single-copy counts vs the common set (surfaces gene loss)
    input:
        single_copy=expand(species_ids_dir_path / "single_copy/{species}.ids", species=config["species_list"]),
        multi_copy=expand(species_ids_dir_path / "multi_copy/{species}.ids", species=config["species_list"]),
        common=common_ids_dir_path / "common.ids",
    output:
        main_ids_dir_path / "gene_counts.tsv",
    params:
        single_copy_dir=species_ids_dir_path / "single_copy",
        multi_copy_dir=species_ids_dir_path / "multi_copy",
    log:
        std=log_dir_path / "gene_count_report.log",
        cluster_log=cluster_log_dir_path / "gene_count_report.cluster.log",
        cluster_err=cluster_log_dir_path / "gene_count_report.cluster.err",
    benchmark:
        benchmark_dir_path / "gene_count_report.benchmark.txt"
    conda:
        main_env
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/gene_count_report.py "
        " --single_copy_dir {params.single_copy_dir} "
        " --multi_copy_dir {params.multi_copy_dir} "
        " --common {input.common} "
        " --output {output} "
        " 1> {log.std} 2>&1; "


checkpoint merged_sequences:  # get merged sequences by common IDs
    input:
        rules.common_ids.output,
    output:
        directory(merged_sequences_dir_path),
    params:
        single_copy_files=expand(
            busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences",
            species=species_list_for_raw_alignment if reconstruct_gapaware and reconstruct_map else config["species_list"],
        ),
    log:
        std=log_dir_path / "merged_sequences.log",
        cluster_log=cluster_log_dir_path / "merged_sequences.cluster.log",
        cluster_err=cluster_log_dir_path / "merged_sequences.cluster.err",
    benchmark:
        benchmark_dir_path / "merged_sequences.benchmark.txt"
    conda:
        main_env
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/merge_common_ids.py -c {input} -s {params.single_copy_files} -o {output} > {log.std} 2>&1; "
