if config.get("vcf2phylip") != True:

    rule concat_fasta:
        input:
            lambda w: expand_fna_from_merged_sequences(w, filtered_alignments_dir_path / "fna" / "{N}.fna",
                                                    busco_blacklist=busco_blacklist),
        output:
            concat_alignments_dir_path / fasta_filename,
        log:
            std=log_dir_path / "concat_fasta.log",
            cluster_log=cluster_log_dir_path / "concat_fasta.cluster.log",
            cluster_err=cluster_log_dir_path / "concat_fasta.cluster.err",
        benchmark:
            benchmark_dir_path / "concat_fasta.benchmark.txt"
        conda:
            config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            """
            workflow/scripts/concat_fasta.py -i {input} -o {output} 1> {log.std} 2>&1
            """


rule concat_nexus:
    input:
        concat_alignments_dir_path / fasta_filename,
    output:
        concat_alignments_dir_path / nexus_filename,
    params:
        type="DNA",
        block=config["mrbayes_block"],
    log:
        std=log_dir_path / "concat_nexus.log",
        cluster_log=cluster_log_dir_path / "concat_nexus.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_nexus.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_nexus.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_nexus.py -i {input} "
        " -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1; "


rule concat_stockholm:
    input:
        concat_alignments_dir_path / fasta_filename,
    output:
        concat_alignments_dir_path / stockholm_filename,
    log:
        std=log_dir_path / "concat_stockholm.log",
        cluster_log=cluster_log_dir_path / "concat_stockholm.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_stockholm.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_stockholm.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_stockholm.py -i {input} -o {output} 1> {log.std} 2>&1 "


rule concat_phylip:
    input:
        concat_alignments_dir_path / fasta_filename,
    output:
        phy=concat_alignments_dir_path / phylip_filename,
        map=concat_alignments_dir_path / f"{phylip_filename}.map",
    log:
        std=log_dir_path / "concat_phylip.log",
        cluster_log=cluster_log_dir_path / "concat_phylip.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_phylip.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_phylip.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_phylip.py -i {input} -o {output.phy} 1> {log.std} 2>&1; "
