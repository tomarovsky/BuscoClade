if config.get("alignment") == "prank":

    rule prank:
        input:
            merged_sequences_dir_path / "{N}.fna",
        output:
            alignments_dir_path / "fna" / "{N}.fna",
        params:
            prefix=lambda wildcards, output: output[0][:-4],
            options=config["prank_params"],
        log:
            std=log_dir_path / "prank.{N}.log",
            cluster_log=cluster_log_dir_path / "prank.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "prank.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "prank.{N}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            queue=config["alignment_queue"],
            cpus=config["prank_threads"],
            time=config["prank_time"],
            mem_mb=config["prank_mem_mb"],
        threads: config["prank_threads"]
        shell:
            " prank -d={input} -o={params.prefix} {params.options} 1> {log.std} 2>&1; "
            " mv {params.prefix}.best.fas {output} >> {log.std} 2>&1; "

if config["alignment"] == "mafft":

    rule mafft:
        input:
            merged_sequences_dir_path / "{N}.fna",
        output:
            alignments_dir_path / "fna" / "{N}.fna",
        params:
            options=config["mafft_params"],
        log:
            std=log_dir_path / "mafft.{N}.log",
            cluster_log=cluster_log_dir_path / "mafft.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "mafft.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "mafft.{N}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            queue=config["alignment_queue"],
            cpus=config["mafft_threads"],
            time=config["mafft_time"],
            mem_mb=config["mafft_mem_mb"],
        threads: config["mafft_threads"]
        shell:
            " mafft --thread {threads} {params.options} {input} > {output} 2> {log.std}; "
