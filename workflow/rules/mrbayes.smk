rule mrbayes:
    input:
        concat_alignments_dir_path / nexus_filename,
    output:
        directory(mrbayes_dir_path / "fna"),
    params:
        mrbayes_path=config["mrbayes_path"],
        options=config["mrbayes_params"],
    log:
        std=log_dir_path / "mrbayes.log",
        cluster_log=cluster_log_dir_path / "mrbayes.cluster.log",
        cluster_err=cluster_log_dir_path / "mrbayes.cluster.err",
    benchmark:
        benchmark_dir_path / "mrbayes.benchmark.txt"
    resources:
        slurm_partition=config["mrbayes_queue"],
        runtime=config["mrbayes_time"],
        mem_mb=config["mrbayes_mem_mb"],
    threads: config["mrbayes_threads"]
    shell:
        " mkdir -p {output}; "
        " mpirun -np {resources.cpus} {params.mrbayes_path}/mb-mpi {input} {params.options} 1> {log.std} 2>&1; "
        " mv {input}.* {output}/; "
