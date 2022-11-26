if "dna_alignment" in config:
    if config["dna_alignment"] == "prank":
        rule prank_dna:
            input:
                merged_sequences_dir_path / "{N}.fna"
            output:
                alignments_dir_path / "fna" / "{N}.fna"
            params:
                prefix=lambda wildcards, output: output[0][:-4],
                options=config["prank_dna_params"]
            log:
                std=log_dir_path / "prank_dna.{N}.log",
                cluster_log=cluster_log_dir_path / "prank_dna.{N}.cluster.log",
                cluster_err=cluster_log_dir_path / "prank_dna.{N}.cluster.err"
            benchmark:
                benchmark_dir_path / "prank_dna.{N}.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                cpus=config["prank_threads"],
                time=config["prank_time"],
                mem_mb=config["prank_mem_mb"]
            threads:
                config["prank_threads"]
            shell:
                "prank -d={input} -o={params.prefix} {params.options} > {log.std} 2>&1; "
                "mv {params.prefix}.best.fas {output}; "

    if config["dna_alignment"] == "mafft":
        rule mafft_dna:
            input:
                merged_sequences_dir_path / "{N}.fna"
            output:
                alignments_dir_path / "fna" / "{N}.fna"
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
                "mafft --thread {threads} {params.options} {input} > {output} 2> {log.std}; "

if "protein_alignment" in config:
    if config["protein_alignment"] == "prank":
        rule prank_protein:
            input:
                merged_sequences_dir_path / "{N}.faa"
            output:
                alignments_dir_path / "faa" / "{N}.faa"
            params:
                prefix=lambda wildcards, output: output[0][:-4],
                options=config["prank_protein_params"]
            log:
                std=log_dir_path / "prank_protein.{N}.log",
                cluster_log=cluster_log_dir_path / "prank_protein.{N}.cluster.log",
                cluster_err=cluster_log_dir_path / "prank_protein.{N}.cluster.err"
            benchmark:
                benchmark_dir_path / "prank_protein.{N}.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                cpus=config["prank_threads"],
                time=config["prank_time"],
                mem_mb=config["prank_mem_mb"]
            threads:
                config["prank_threads"]
            shell:
                "prank -d={input} -o={params.prefix} {params.options} > {log.std} 2>&1; "
                "mv {params.prefix}.best.fas {output}; "

    if config["protein_alignment"] == "mafft":
        rule mafft_protein:
            input:
                merged_sequences_dir_path / "{N}.faa"
            output:
                alignments_dir_path / "faa" / "{N}.faa"
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
                "mafft --thread {threads} {params.options} {input} > {output} 2> {log.std}; "


