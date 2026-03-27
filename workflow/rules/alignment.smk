if config.get("alignment") == "prank":

    def parse_time_to_minutes(time_str):
        return sum(
            int(v) * {"h": 60, "m": 1}.get(u, 1)
            for v, u in re.findall(r"(\d+)([hm]?)", str(time_str))
            if v
        )

    rule prank:
        input:
            merged_sequences_dir_path / "{N}.fna",
        output:
            alignments_dir_path / "fna" / "{N}.fna",
        params:
            prefix=lambda wildcards, output: output[0][:-4],
            options=config["prank_params"],
            timeout=lambda wildcards: f"{parse_time_to_minutes(config['prank_time']) - 15}m",
        log:
            std=log_dir_path / "prank.{N}.log",
            cluster_log=cluster_log_dir_path / "prank.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "prank.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "prank.{N}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            slurm_partition=config["alignment_queue"],
            runtime=config["prank_time"],
            mem_mb=config["prank_mem_mb"],
        threads: config["prank_threads"]
        shell:
            " timeout --kill-after=5m {params.timeout} "
            " prank -d={input} -o={params.prefix} {params.options} 1> {log.std} 2>&1 "
            " || {{ exit_code=$?; "
            " if [ $exit_code -ne 0 ]; then "
            "     echo 'PRANK failed or timed out (exit code: '$exit_code') for {wildcards.N}' >> {log.std}; "
            "     touch {output}; "
            " fi; "
            " }}; "
            " if [ ! -f {output} ]; then "
            "     mv {params.prefix}.best.fas {output} >> {log.std} 2>&1; "
            " fi; "


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
            slurm_partition=config["alignment_queue"],
            runtime=config["mafft_time"],
            mem_mb=config["mafft_mem_mb"],
        threads: config["mafft_threads"]
        shell:
            " mafft --thread {threads} {params.options} {input} > {output} 2> {log.std}; "


if config["alignment"] == "muscle":

    rule muscle:
        input:
            merged_sequences_dir_path / "{N}.fna",
        output:
            alignments_dir_path / "fna" / "{N}.fna",
        params:
            options=config["muscle_params"],
        log:
            std=log_dir_path / "muscle.{N}.log",
            cluster_log=cluster_log_dir_path / "muscle.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "muscle.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "muscle.{N}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            slurm_partition=config["alignment_queue"],
            runtime=config["muscle_time"],
            mem_mb=config["muscle_mem_mb"],
        threads: config["muscle_threads"]
        shell:
            " muscle -in {input} -out {output} {params.options} > {log.std} 2>&1; "
