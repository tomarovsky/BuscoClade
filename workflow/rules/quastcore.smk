rule quastcore:
    input:
        lambda w: get_all_genome_files(),
    output:
        quastcore_dir_path / "assembly_stats.csv",
    params:
        config["quastcore_params"],
    log:
        std=log_dir_path / "quastcore.log",
        cluster_log=cluster_log_dir_path / "quastcore.cluster.log",
        cluster_err=cluster_log_dir_path / "quastcore.cluster.err",
    benchmark:
        benchmark_dir_path / "quastcore.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"]),
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/quastcore.py -i {input} {params} -o {output} > {log.std} 2>&1; "
