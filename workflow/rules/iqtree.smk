rule iqtree_dna:
    input:
        concat_alignments_dir_path / fasta_dna_filename
    output:
        directory(iqtree_dir_path / "fna")
    params:
        iqtree_path=config["iqtree_path"],
        options=config["iqtree_dna_params"],
        prefix=fasta_dna_filename
    log:
        std=log_dir_path / "iqtree_dna.log",
        cluster_log=cluster_log_dir_path / "iqtree_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "iqtree_dna.cluster.err"
    benchmark:
        benchmark_dir_path / "iqtree_dna.benchmark.txt"
    resources:
        cpus=config["iqtree_threads"],
        time=config["iqtree_time"],
        mem_mb=config["iqtree_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "{params.iqtree_path}/iqtree -s {input} -pre {params.prefix} -nt {resources.cpus} {params.options} 1> {log.std} 2>&1; "
        "mv {params.prefix}.* {output}"


rule iqtree_protein:
    input:
        concat_alignments_dir_path / fasta_protein_filename
    output:
        directory(iqtree_dir_path / "faa")
    params:
        iqtree_path=config["iqtree_path"],
        options=config["iqtree_protein_params"],
        prefix=fasta_protein_filename
    log:
        std=log_dir_path / "iqtree_protein.log",
        cluster_log=cluster_log_dir_path / "iqtree_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "iqtree_protein.cluster.err"
    benchmark:
        benchmark_dir_path / "iqtree_protein.benchmark.txt"
    resources:
        cpus=config["iqtree_threads"],
        time=config["iqtree_time"],
        mem_mb=config["iqtree_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "{params.iqtree_path}/iqtree -s {input} -pre {params.prefix} -nt {resources.cpus} {params.options} 1> {log.std} 2>&1; "
        "mv {params.prefix}.* {output}"
