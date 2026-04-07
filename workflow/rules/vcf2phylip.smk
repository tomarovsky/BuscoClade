if config.get("vcf2phylip"):

    rule vcf2phylip:
        input:
            vcf_file,
        output:
            concat_alignments_dir_path / fasta_filename,
        params:
            prefix=f"{prefix}.min4.fasta",
            outdir=concat_alignments_dir_path,
        log:
            std=log_dir_path / f"vcf2phylip.{prefix}.log",
            cluster_log=cluster_log_dir_path / f"vcf2phylip.{prefix}.cluster.log",
            cluster_err=cluster_log_dir_path / f"vcf2phylip.{prefix}.cluster.err",
        benchmark:
            benchmark_dir_path / f"vcf2phylip.{prefix}.benchmark.txt",
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"]
            else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"]),
        resources:
            slurm_partition=config["processing_queue"],
            runtime=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"],
        shell:
            " MYPWD=$(pwd); "
            " cd {params.outdir}; "
            " python $MYPWD/workflow/scripts/vcf2phylip.py -i {input} --phylip-disable --fasta > $MYPWD/{log.std} 2>&1; "
            " mv {params.prefix} $MYPWD/{output}; "
