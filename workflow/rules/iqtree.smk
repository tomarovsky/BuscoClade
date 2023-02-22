if config["iqtree_dna"]:
    rule iqtree_dna:
        input:
            concat_alignments_dir_path / fasta_dna_filename
        output:
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.bionj",
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.ckp.gz",
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.iqtree",
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.log",
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.mldist",
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.model.gz",
            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.treefile"
        params:
            options=config["iqtree_dna_params"],
            outdir=iqtree_dir_path / "fna",
            prefix=fasta_dna_filename
        log:
            std=log_dir_path / "iqtree_dna.log",
            cluster_log=cluster_log_dir_path / "iqtree_dna.cluster.log",
            cluster_err=cluster_log_dir_path / "iqtree_dna.cluster.err"
        benchmark:
            benchmark_dir_path / "iqtree_dna.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["iqtree_threads"],
            time=config["iqtree_time"],
            mem_mb=config["iqtree_mem_mb"]
        threads:
            config["iqtree_threads"]
        shell:
            "mkdir -p {params.outdir}; iqtree -nt {threads} -s {input} "
            "--prefix {params.outdir}/{params.prefix} {params.options} 1> {log.std} 2>&1; "

if config["iqtree_protein"]:
    rule iqtree_protein:
        input:
            concat_alignments_dir_path / fasta_protein_filename
        output:
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.bionj",
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.ckp.gz",
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.iqtree",
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.log",
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.mldist",
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.model.gz",
            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.treefile"
        params:
            options=config["iqtree_protein_params"],
            outdir=iqtree_dir_path / "faa",
            prefix=fasta_protein_filename
        log:
            std=log_dir_path / "iqtree_protein.log",
            cluster_log=cluster_log_dir_path / "iqtree_protein.cluster.log",
            cluster_err=cluster_log_dir_path / "iqtree_protein.cluster.err"
        benchmark:
            benchmark_dir_path / "iqtree_protein.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["iqtree_threads"],
            time=config["iqtree_time"],
            mem_mb=config["iqtree_mem_mb"]
        threads:
            config["iqtree_threads"]
        shell:
            "mkdir -p {params.outdir}; iqtree -nt {threads} -s {input} "
            "--prefix {params.outdir}/{params.prefix} {params.options} 1> {log.std} 2>&1; "


