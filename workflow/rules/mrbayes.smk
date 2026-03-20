rule mrbayes:
    input:
        concat_alignments_dir_path / nexus_filename,
    output:
        mrbayes_dir_path / f"{fasta_filename}.nex.con.tre",
    params:
        output_dir=mrbayes_dir_path,
        options=config["mrbayes_params"],
    log:
        std=log_dir_path / "mrbayes.log",
        cluster_log=cluster_log_dir_path / "mrbayes.cluster.log",
        cluster_err=cluster_log_dir_path / "mrbayes.cluster.err",
    benchmark:
        benchmark_dir_path / "mrbayes.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        slurm_partition=config["mrbayes_queue"],
        runtime=config["mrbayes_time"],
        mem_mb=config["mrbayes_mem_mb"],
    threads: config["mrbayes_threads"]
    shell:
        " mkdir -p {params.output_dir}; "
        " INPUT=$(readlink -f {input}); "
        " LOGFILE=$(readlink -f {log.std}); "
        " ln -sf $INPUT {params.output_dir}/; "
        " cd {params.output_dir}; "
        " mpirun -np {threads} mb-mpi $(basename {input}) {params.options} 1> $LOGFILE 2>&1; "


rule mrbayes_convert:
    input:
        mrbayes_dir_path / f"{fasta_filename}.nex.con.tre",
    output:
        mrbayes_dir_path / f"{fasta_filename}.nex.con.tre.nwk",
    log:
        std=log_dir_path / "mrbayes_convert.log",
        cluster_log=cluster_log_dir_path / "mrbayes_convert.cluster.log",
        cluster_err=cluster_log_dir_path / "mrbayes_convert.cluster.err",
    benchmark:
        benchmark_dir_path / "mrbayes_convert.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " gotree reformat newick -i {input} -f nexus -o {output} > {log.std} 2>&1; "
