rule rapidnj:
    input:
        concat_alignments_dir_path / stockholm_filename,
    output:
        tree=rapidnj_dir_path / rapidnj_tree,
        matrix=rapidnj_dir_path / rapidnj_matrix,
    params:
        config["rapidnj_params"],
    log:
        std=log_dir_path / "rapidnj_tree.log",
        cluster_log=cluster_log_dir_path / "rapidnj_tree.cluster.log",
        cluster_err=cluster_log_dir_path / "rapidnj_tree.cluster.err",
    benchmark:
        benchmark_dir_path / "rapidnj_tree.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        slurm_partition=config["rapidnj_queue"],
        runtime=config["rapidnj_time"],
        mem_mb=config["rapidnj_mem_mb"],
    threads: config["rapidnj_threads"]
    shell:
        " rapidnj -i sth -c {threads} -o m {input} > {output.matrix} 2>{log.std}; "
        " rapidnj -i sth -c {threads} {params} {input} > {output.tree} 2>>{log.std}; "
