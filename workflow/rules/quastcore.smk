rule quastcore:
    input:
        expand(genome_dir_path / "{species}.fasta", species=config["species_list"]),
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
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"]),
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " quast_core.py -i {input} {params} -b 10000000 -o {output} > {log.std} 2>&1; "
