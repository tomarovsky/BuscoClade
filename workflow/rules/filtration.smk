if config.get("filtration") == "gblocks":

    rule gblocks:
        input:
            alignments_dir_path / "fna" / "{N}.fna",
        output:
            filtered_alignments_dir_path / "fna" / "{N}.fna",
        params:
            outdir=filtered_alignments_dir_path / "fna",
            options=config["gblocks_params"],
        log:
            std=log_dir_path / "gblocks.{N}.log",
            cluster_log=cluster_log_dir_path / "gblocks.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "gblocks.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "gblocks.{N}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            slurm_partition=config["filtration_queue"],
            runtime=config["gblocks_time"],
            mem_mb=config["gblocks_mem_mb"],
        threads: config["gblocks_threads"]
        shell:
            " mkdir -p {params.outdir}; set +e; "
            " Gblocks {input} {params.options} 1> {log.std} 2>&1; "
            " rm -r {input}-gb.htm; mv {input}-gb {output} "  # Gblocks always returns 1 exit code
            # If Gblocks exits with an error, the output files will not exist

if config.get("filtration") == "trimal":

    rule trimal:
        input:
            alignments_dir_path / "fna" / "{N}.fna",
        output:
            filtered_alignments_dir_path / "fna" / "{N}.fna",
        params:
            prefix=lambda wildcards, output: output[0][:-4],
            options=config["trimal_params"],
        log:
            std=log_dir_path / "trimal.{N}.log",
            cluster_log=cluster_log_dir_path / "trimal.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "trimal.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "trimal.{N}.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            slurm_partition=config["filtration_queue"],
            runtime=config["trimal_time"],
            mem_mb=config["trimal_mem_mb"],
        threads: config["trimal_threads"]
        shell:
            " trimal -in {input} -out {params.prefix} {params.options} > {log.std} 2>&1; "
            " trimal -in {params.prefix} -out {output} -nogaps >> {log.std} 2>&1; "
            " rm {params.prefix} >> {log.std} 2>&1; "
