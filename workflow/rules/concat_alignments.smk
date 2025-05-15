rule concat_fasta_dna:
    input:
        lambda w: (
            rules.consolidate_vcf2phylip.output
            if config["vcf2phylip"]
            else expand_fna_from_merged_sequences(w, filtered_alignments_dir_path / "fna" / "{N}.fna")),
    output:
        concat_alignments_dir_path / fasta_dna_filename,
    params:
        vcf2phylip=config["vcf2phylip"],
    log:
        std=log_dir_path / "concat_fasta_dna.log",
        cluster_log=cluster_log_dir_path / "concat_fasta_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_fasta_dna.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_fasta_dna.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        """
        if {params.vcf2phylip}; then
            mv {input} {output} 1> {log.std} 2>&1
        else
            workflow/scripts/concat_fasta.py -i {input} -o {output} 1> {log.std} 2>&1
        fi
        """


rule concat_fasta_protein:
    input:
        lambda w: expand_fna_from_merged_sequences(w, filtered_alignments_dir_path / "faa" / "{N}.faa"),
    output:
        concat_alignments_dir_path / fasta_protein_filename,
    log:
        std=log_dir_path / "concat_fasta_protein.log",
        cluster_log=cluster_log_dir_path / "concat_fasta_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_fasta_protein.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_fasta_protein.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/concat_fasta.py -i {input} -o {output} 1> {log.std} 2>&1; "


rule concat_nexus_dna:
    input:
        rules.concat_fasta_dna.output, 
    output:
        concat_alignments_dir_path / nexus_dna_filename
    params:
        type="DNA",
        block=config["mrbayes_block"],
    log:
        std=log_dir_path / "concat_nexus_dna.log",
        cluster_log=cluster_log_dir_path / "concat_nexus_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_nexus_dna.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_nexus_dna.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_nexus.py -i {input} "
        " -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1; "


rule concat_nexus_protein:
    input:
        rules.concat_fasta_protein.output,
    output:
        concat_alignments_dir_path / nexus_protein_filename,
    params:
        type="protein",
        block=config["mrbayes_block"],
    log:
        std=log_dir_path / "concat_nexus_protein.log",
        cluster_log=cluster_log_dir_path / "concat_nexus_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_nexus_protein.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_nexus_protein.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_nexus.py -i {input} "
        " -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1; "


rule concat_stockholm_dna:
    input:
        rules.concat_fasta_dna.output,
    output:
        concat_alignments_dir_path / stockholm_dna_filename,
    log:
        std=log_dir_path / "concat_stockholm_dna.log",
        cluster_log=cluster_log_dir_path / "concat_stockholm_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_stockholm_dna.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_stockholm_dna.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_stockholm.py -i {input} -o {output} 1> {log.std} 2>&1 "


rule concat_stockholm_protein:
    input:
        rules.concat_fasta_protein.output,
    output:
        concat_alignments_dir_path / stockholm_protein_filename,
    log:
        std=log_dir_path / "concat_stockholm_protein.log",
        cluster_log=cluster_log_dir_path / "concat_stockholm_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_stockholm_protein.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_stockholm_protein.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_stockholm.py -i {input} -o {output} 1> {log.std} 2>&1; "


rule concat_phylip_dna:
    input:
        rules.concat_fasta_dna.output,
    output:
        concat_alignments_dir_path / phylip_dna_filename,
    params:
        type="DNA",
    log:
        std=log_dir_path / "concat_phylip_dna.log",
        cluster_log=cluster_log_dir_path / "concat_phylip_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_phylip_dna.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_phylip_dna.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_phylip.py -i {input} "
        " -t {params.type} -o {output} 1> {log.std} 2>&1; "


rule concat_phylip_protein:
    input:
        rules.concat_fasta_protein.output,
    output:
        concat_alignments_dir_path / phylip_protein_filename,
    params:
        type="protein",
    log:
        std=log_dir_path / "concat_phylip_protein.log",
        cluster_log=cluster_log_dir_path / "concat_phylip_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_phylip_protein.cluster.err",
    benchmark:
        benchmark_dir_path / "concat_phylip_protein.benchmark.txt"
    conda:
        config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    shell:
        " workflow/scripts/fasta_to_phylip.py -i {input} "
        " -t {params.type} -o {output} 1> {log.std} 2>&1; "
