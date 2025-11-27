rule iqtree:
    input:
        concat_alignments_dir_path / fasta_filename,
    output:
        iqtree_dir_path / "fna" / f"{fasta_filename}.bionj",
        iqtree_dir_path / "fna" / f"{fasta_filename}.ckp.gz",
        iqtree_dir_path / "fna" / f"{fasta_filename}.iqtree",
        iqtree_dir_path / "fna" / f"{fasta_filename}.log",
        iqtree_dir_path / "fna" / f"{fasta_filename}.mldist",
        iqtree_dir_path / "fna" / f"{fasta_filename}.model.gz",
        iqtree_dir_path / "fna" / f"{fasta_filename}.treefile",
    params:
        options=config["iqtree_params"],
        outdir=iqtree_dir_path / "fna",
        prefix=fasta_filename,
    log:
        std=log_dir_path / "iqtree.log",
        cluster_log=cluster_log_dir_path / "iqtree.cluster.log",
        cluster_err=cluster_log_dir_path / "iqtree.cluster.err",
    benchmark:
        benchmark_dir_path / "iqtree.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["iqtree_queue"],
        cpus=config["iqtree_threads"],
        time=config["iqtree_time"],
        mem_mb=config["iqtree_mem_mb"],
    threads: config["iqtree_threads"]
    shell:
        " mkdir -p {params.outdir}; iqtree -nt {threads} -s {input} "
        " --prefix {params.outdir}/{params.prefix} {params.options} 1> {log.std} 2>&1; "
