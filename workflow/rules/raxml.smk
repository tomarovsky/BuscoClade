rule raxml:
    input:
        concat_alignments_dir_path / fasta_filename
    output:
        bestTree=raxml_dir_path / f"{fasta_filename}.raxml.bestTree",
        log=raxml_dir_path / f"{fasta_filename}.raxml.log",
        rba=raxml_dir_path / f"{fasta_filename}.raxml.rba",
        treefile=raxml_dir_path / f"{fasta_filename}.raxml.treefile",
    params:
        support=raxml_dir_path / f"{fasta_filename}.raxml.support",
        options=config["raxml_params"],
        outdir=raxml_dir_path,
        prefix=fasta_filename,
    log:
        std=log_dir_path / "raxml.log",
        cluster_log=cluster_log_dir_path / "raxml.cluster.log",
        cluster_err=cluster_log_dir_path / "raxml.cluster.err",
    benchmark:
        benchmark_dir_path / "raxml.benchmark.txt",
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["raxml_queue"],
        cpus=config["raxml_threads"],
        time=config["raxml_time"],
        mem_mb=config["raxml_mem_mb"],
    threads: config["raxml_threads"],
    shell:
        """
        if [ -d {params.outdir} ]; then
            rm -rf {params.outdir}/*;
        fi
        mkdir -p {params.outdir};
        raxml-ng --all --msa {input} --prefix {params.outdir}/{params.prefix} {params.options} 1> {log.std} 2>&1;
        mv {params.support} {output.treefile};
        """
