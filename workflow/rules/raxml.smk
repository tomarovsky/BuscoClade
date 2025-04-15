if config["raxml"]:

    rule raxml_tree:
        input:
            concat_alignments_dir_path / fasta_dna_filename
        output:
            raxml_dir_path / f"{fasta_dna_filename}.raxml.bestTree",
            raxml_dir_path / f"{fasta_dna_filename}.raxml.log",
            raxml_dir_path / f"{fasta_dna_filename}.raxml.rba",
            raxml_dir_path / f"{fasta_dna_filename}.raxml.support"
        params: 
            options=config["raxml_params"],
            outdir=raxml_dir_path,
            prefix=fasta_dna_filename,
        log:
            std=log_dir_path / "raxml.log",
            cluster_log=cluster_log_dir_path / "raxml.cluster.log",
            cluster_err=cluster_log_dir_path / "raxml.cluster.err",
        benchmark:
            benchmark_dir_path / "raxml.benchmark.txt",
        conda:
            "../../%s" % config["conda_config"],
        resources:
            queue=config["raxml_queue"],
            cpus=config["raxml_threads"],
            time=config["raxml_time"],
            mem_mb=config["raxml_mem_mb"],
        threads: config["raxml_threads"],
        shell:
            "mkdir -p {params.outdir}; raxml-ng --all --msa {input} --prefix {params.outdir}/{params.prefix} "
            "{params.options} 1> {log.std} 2>&1; "
            
    rule raxml_tree_rename:
        input:
            raxml_dir_path / f"{fasta_dna_filename}.raxml.support",
        output:
            raxml_dir_path / f"{fasta_dna_filename}.raxml.treefile",
        params:
            outdir=raxml_dir_path,   
        log:
            std=log_dir_path / "raxml_rename.log",
            cluster_log=cluster_log_dir_path / "raxml_rename.cluster.log",
            cluster_err=cluster_log_dir_path / "raxml_rename.cluster.err",
        benchmark:
            benchmark_dir_path / "raxml.benchmark.txt",
        conda:
            "../../%s" % config["conda_config"],
        resources:
            queue=config["raxml_queue"],
            cpus=config["raxml_threads"],
            time=config["raxml_time"],
            mem_mb=config["raxml_mem_mb"],
        threads: config["raxml_threads"],
        shell:
            "mv {input} {output} > {log.std} 2>&1; "