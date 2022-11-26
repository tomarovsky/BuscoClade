if "dna_filtration" in config:
    if config["dna_filtration"] == "gblocks":
        rule gblocks_dna:
            input:
                alignments_dir_path / "fna" / "{N}.fna"
            output:
                filtered_alignments_dir_path / "fna" / "{N}.fna"
            params:
                outdir=filtered_alignments_dir_path / "fna",
                options=config["gblocks_dna_params"]
            log:
                std=log_dir_path / "{N}.fna.gblocks.log",
                cluster_log=cluster_log_dir_path / "{N}.fna.gblocks.cluster.log",
                cluster_err=cluster_log_dir_path / "{N}.fna.gblocks.cluster.err"
            benchmark:
                benchmark_dir_path / "{N}.fna.gblocks.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                cpus=config["gblocks_threads"],
                time=config["gblocks_time"],
                mem=config["gblocks_mem_mb"]
            shell:
                "mkdir -p {params.outdir}; set +e; " # Gblocks always returns 1 exit code
                "Gblocks {input} {params.options} 1> {log.std} 2>&1; "
                "rm -r {input}-gb.htm; mv {input}-gb {output} " # If Gblocks exits with an error, the output files will not exist

    if config["dna_filtration"] == "trimal":
        rule trimal_dna:
            input:
                alignments_dir_path / "fna" / "{N}.fna"
            output:
                filtered_alignments_dir_path / "fna" / "{N}.fna"
            params:
                prefix=lambda wildcards, output: output[0][:-4],
                options=config["trimal_dna_params"]
            log:
                std=log_dir_path / "trimal_dna.{N}.log",
                cluster_log=cluster_log_dir_path / "trimal_dna.{N}.cluster.log",
                cluster_err=cluster_log_dir_path / "trimal_dna.{N}.cluster.err"
            benchmark:
                benchmark_dir_path / "trimal_dna.{N}.benchmark.txt"
            resources:
                cpus=config["trimal_threads"],
                time=config["trimal_time"],
                mem_mb=config["trimal_mem_mb"]
            shell:
                "trimal -in {input} -out {params.prefix} {params.options} > {log.std} 2>&1; "
                "trimal -in {params.prefix} -out {params.prefix} -nogaps > {log.std} 2>&1; "

if "protein_filtration" in config:
    if config["protein_filtration"] == "gblocks":
        rule gblocks_protein:
            input:
                alignments_dir_path / "faa" / "{N}.faa"
            output:
                filtered_alignments_dir_path / "faa" / "{N}.faa"
            params:
                outdir=filtered_alignments_dir_path / "faa",
                options=config["gblocks_protein_params"]
            log:
                std=log_dir_path / "{N}.faa.gblocks.log",
                cluster_log=cluster_log_dir_path / "{N}.faa.gblocks.cluster.log",
                cluster_err=cluster_log_dir_path / "{N}.faa.gblocks.cluster.err"
            benchmark:
                benchmark_dir_path / "{N}.faa.gblocks.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                cpus=config["gblocks_threads"],
                time=config["gblocks_time"],
                mem=config["gblocks_mem_mb"]
            shell:
                "mkdir -p {params.outdir}; set +e; " # Gblocks always returns 1 exit code
                "Gblocks {input} {params.options} 1> {log.std} 2>&1; "
                "rm -r {input}-gb.htm; mv {input}-gb {output} " # If Gblocks exits with an error, the output files will not exist

    if config["protein_filtration"] == "trimal":
        rule trimal_protein:
            input:
                alignments_dir_path / "faa" / "{N}.faa"
            output:
                filtered_alignments_dir_path / "faa" / "{N}.faa"
            params:
                prefix=lambda wildcards, output: output[0][:-4],
                options=config["trimal_protein_params"]
            log:
                std=log_dir_path / "trimal_protein.{N}.log",
                cluster_log=cluster_log_dir_path / "trimal_protein.{N}.cluster.log",
                cluster_err=cluster_log_dir_path / "trimal_protein.{N}.cluster.err"
            benchmark:
                benchmark_dir_path / "trimal_protein.{N}.benchmark.txt"
            resources:
                cpus=config["trimal_threads"],
                time=config["trimal_time"],
                mem_mb=config["trimal_mem_mb"]
            shell:
                "trimal -in {input} -out {params.prefix} {params.options} > {log.std} 2>&1; "
                "trimal -in {params.prefix} -out {params.prefix} -nogaps > {log.std} 2>&1; "