if config["alignment_tool"] == "mafft":
    rule mafft_dna:
        input:
            merged_sequences_dir_path / "group_{N}"
        output:
            temp(directory(alignment_dir_path / "fna_tmp" / "group_{N}"))
        params:
            options=config["mafft_dna_params"]
        log:
            std=log_dir_path / "mafft_dna.{N}.log",
            cluster_log=cluster_log_dir_path / "mafft_dna.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "mafft_dna.{N}.cluster.err"
        benchmark:
            benchmark_dir_path / "mafft_dna.{N}.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["mafft_threads"],
            time=config["mafft_time"],
            mem_mb=config["mafft_mem_mb"]
        threads:
            config["mafft_threads"]
        shell:
            "mkdir -p {output}; "
            "for FILE in `ls {input}/*.fna`; do "
            "mafft --thread {threads} {params.options} $FILE > {output}/$(basename $FILE) 2> {log.std}; "
            "done; "


rule mafft_protein:
    input:
        merged_sequences_dir_path / "group_{N}"
    output:
        temp(directory(alignment_dir_path / "faa_tmp" / "group_{N}"))
    params:
        options=config["mafft_protein_params"]
    log:
        std=log_dir_path / "mafft_protein.{N}.log",
        cluster_log=cluster_log_dir_path / "mafft_protein.{N}.cluster.log",
        cluster_err=cluster_log_dir_path / "mafft_protein.{N}.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft_protein.{N}.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem_mb=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls {input}/*.faa`; do "
        "mafft --thread {threads} {params.options} $FILE > {output}/$(basename $FILE) 2> {log.std}; "
        "done; "