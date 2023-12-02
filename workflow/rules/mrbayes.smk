if config["mrbayes_dna"]:
    rule mrbayes_dna:
        input:
            concat_alignments_dir_path / nexus_dna_filename
        output:
            directory(mrbayes_dir_path / "fna")
        params:
            mrbayes_path=config["mrbayes_path"],
            options=config["mrbayes_dna_params"]
        log:
            std=log_dir_path / "mrbayes_dna.log",
            cluster_log=cluster_log_dir_path / "mrbayes_dna.cluster.log",
            cluster_err=cluster_log_dir_path / "mrbayes_dna.cluster.err"
        benchmark:
            benchmark_dir_path / "mrbayes_dna.benchmark.txt"
        resources:
            queue=config["mrbayes_queue"],
            cpus=config["mrbayes_threads"],
            time=config["mrbayes_time"],
            mem_mb=config["mrbayes_mem_mb"]
        shell:
            " mkdir -p {output}; "
            " mpirun -np {resources.cpus} {params.mrbayes_path}/mb-mpi {input} {params.options} 1> {log.std} 2>&1; "
            " mv {input}.* {output}/; "

if config["mrbayes_protein"]:
    rule mrbayes_protein:
        input:
            concat_alignments_dir_path / nexus_protein_filename
        output:
            directory(mrbayes_dir_path / "faa")
        params:
            mrbayes_path=config["mrbayes_path"],
            options=config["mrbayes_protein_params"]
        log:
            std=log_dir_path / "mrbayes_protein.log",
            cluster_log=cluster_log_dir_path / "mrbayes_protein.cluster.log",
            cluster_err=cluster_log_dir_path / "mrbayes_protein.cluster.err"
        benchmark:
            benchmark_dir_path / "mrbayes_protein.benchmark.txt"
        resources:
            queue=config["mrbayes_queue"],
            cpus=config["mrbayes_threads"],
            time=config["mrbayes_time"],
            mem_mb=config["mrbayes_mem_mb"]
        shell:
            " mkdir -p {output}; "
            " mpirun -np {resources.cpus} {params.mrbayes_path}/mb-mpi {input} {params.options} 1> {log.std} 2>&1; "
            " mv {input}.* {output}/; "


