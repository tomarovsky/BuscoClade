rule iqtree_tree_visualization:
    input:
        iqtree_dir_path / "fna" / f"{fasta_filename}.treefile",
    output:
        iqtree_dir_path / "fna" / f"{fasta_filename}.length_and_support_tree.svg",
        iqtree_dir_path / "fna" / f"{fasta_filename}.only_support_tree.svg",
        iqtree_dir_path / "fna" / f"{fasta_filename}.only_tree.svg",
    params:
        prefix=iqtree_dir_path / "fna" / f"{fasta_filename}",
        options=config["tree_visualization_params"],
    log:
        std=log_dir_path / "iqtree_tree_visualization.log",
        cluster_log=cluster_log_dir_path / "iqtree_tree_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "iqtree_tree_visualization.cluster.err",
    benchmark:
        benchmark_dir_path / "iqtree_tree_visualization.benchmark.txt"
    conda:
        config["conda"]["buscoclade_ete3"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_ete3"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " QT_QPA_PLATFORM=offscreen workflow/scripts/draw_phylotrees.py -i {input} "
        " -o {params.prefix} {params.options} 1> {log.std} 2>&1; "


rule astral_tree_visualization:
    input:
        treefile=astral_dir_path / astral_tree,
        common_ids=rules.common_ids.output,
    output:
        astral_dir_path / f"{astral_tree}.svg",
        astral_dir_path / f"{astral_tree}.pp.svg",
        astral_dir_path / f"{astral_tree}.q.svg",
        astral_dir_path / f"{astral_tree}.tsv",
    params:
        prefix=astral_dir_path / astral_tree,
        options=config["tree_visualization_params"],
    log:
        std=log_dir_path / "astral_tree_visualization.log",
        cluster_log=cluster_log_dir_path / "astral_tree_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "astral_tree_visualization.cluster.err",
    benchmark:
        benchmark_dir_path / "astral_tree_visualization.benchmark.txt"
    conda:
        config["conda"]["buscoclade_ete3"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_ete3"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " QT_QPA_PLATFORM=offscreen; "
        " workflow/scripts/draw_phylotrees_from_astral.py -i {input.treefile} -o {params.prefix} -n $(cat {input.common_ids} | wc -l) {params.options} > {log.std} 2>&1; "
        " workflow/scripts/draw_phylotrees_from_astral_pp.py -i {input.treefile} -o {params.prefix}.pp {params.options} > {log.std} 2>&1; "
        " workflow/scripts/draw_phylotrees_from_astral_q.py -i {input.treefile} -o {params.prefix}.q {params.options} > {log.std} 2>&1; "
        " workflow/scripts/astral_metrics.py -i {input.treefile} -o {params.prefix}.tsv {params.options} > {log.std} 2>&1; "


rule rapidnj_tree_visualization:
    input:
        rapidnj_dir_path / rapidnj_tree,
    output:
        rapidnj_dir_path / f"{fasta_filename}.length_and_support_tree.svg",
        rapidnj_dir_path / f"{fasta_filename}.only_support_tree.svg",
        rapidnj_dir_path / f"{fasta_filename}.only_tree.svg",
    params:
        prefix=rapidnj_dir_path / f"{fasta_filename}",
        options=config["tree_visualization_params"],
    log:
        std=log_dir_path / "rapidnj_tree_visualization.log",
        cluster_log=cluster_log_dir_path / "rapidnj_tree_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "rapidnj_tree_visualization.cluster.err",
    benchmark:
        benchmark_dir_path / "rapidnj_tree_visualization.benchmark.txt"
    conda:
        config["conda"]["buscoclade_ete3"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_ete3"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " QT_QPA_PLATFORM=offscreen workflow/scripts/draw_phylotrees.py -i {input} "
        " -o {params.prefix} {params.options} 1> {log.std} 2>&1; "


rule phylip_tree_visualization:
    input:
        phylip_dir_path / phylip_tree,
    output:
        phylip_dir_path / f"{fasta_filename}.length_and_support_tree.svg",
        phylip_dir_path / f"{fasta_filename}.only_support_tree.svg",
        phylip_dir_path / f"{fasta_filename}.only_tree.svg",
    params:
        prefix=phylip_dir_path / f"{fasta_filename}",
        options=config["tree_visualization_params"],
    log:
        std=log_dir_path / "phylip_tree_visualization.log",
        cluster_log=cluster_log_dir_path / "phylip_tree_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "phylip_tree_visualization.cluster.err",
    benchmark:
        benchmark_dir_path / "phylip_tree_visualization.benchmark.txt"
    conda:
        config["conda"]["buscoclade_ete3"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_ete3"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " QT_QPA_PLATFORM=offscreen workflow/scripts/draw_phylotrees.py -i {input} "
        " -o {params.prefix} {params.options} 1> {log.std} 2>&1; "


rule raxml_tree_visualization:
    input:
        raxml_dir_path / raxml_tree,
    output:
        raxml_dir_path / f"{fasta_filename}.length_and_support_tree.svg",
        raxml_dir_path / f"{fasta_filename}.only_support_tree.svg",
        raxml_dir_path / f"{fasta_filename}.only_tree.svg"
    params:
        prefix=raxml_dir_path / f"{fasta_filename}",
        options=config["tree_visualization_params"],
    log:
        std=log_dir_path / "raxml_tree_visualization.log",
        cluster_log=cluster_log_dir_path / "raxml_tree_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "raxml_tree_visualization.cluster.err",
    benchmark:
        benchmark_dir_path / "raxml_tree_visualization.benchmark.txt",
    conda:
        config["conda"]["buscoclade_ete3"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_ete3"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"],
    shell:
        "QT_QPA_PLATFORM=offscreen workflow/scripts/draw_phylotrees.py -i {input} "
        "-o {params.prefix} {params.options} 1> {log.std} 2>&1; "


rule species_ids_plot:
    input:
        single_copy_ids=expand(species_ids_dir_path / "single_copy/{species}.ids", species=config["species_list"]),
        multi_copy_ids=expand(species_ids_dir_path / "multi_copy/{species}.ids", species=config["species_list"]),
    output:
        species_ids_dir_path / "unique_species_ids.svg",
    log:
        std=log_dir_path / "species_ids_plot.log",
        cluster_log=cluster_log_dir_path / "species_ids_plot.cluster.log",
        cluster_err=cluster_log_dir_path / "species_ids_plot.cluster.err",
    conda:
        config["conda"]["buscoclade_ete3"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_ete3"]["yaml"])
    benchmark:
        benchmark_dir_path / "species_ids_plot.benchmark.txt"
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/unique_ids_plot.py --single_copy_ids_files {input.single_copy_ids} --multi_copy_ids_files {input.multi_copy_ids}"
        " --outplot {output} > {log.std} 2>&1 "


rule busco_summaries_to_tsv:
    input:
        expand(busco_dir_path / "{species}/short_summary_{species}.txt", species=config["species_list"]),
    output:
        busco_dir_path / "busco_summaries.tsv",
    log:
        std=log_dir_path / "busco_summaries_to_tsv.log",
        cluster_log=cluster_log_dir_path / "busco_summaries_to_tsv.cluster.log",
        cluster_err=cluster_log_dir_path / "busco_summaries_to_tsv.cluster.err",
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/busco_summaries_to_tsv.py -i {input} -o {output} > {log.std} 2>&1; "


rule busco_histogram:
    input:
        busco_dir_path / "busco_summaries.tsv",
    output:
        busco_dir_path / "busco_summaries.svg",
    params:
        colors=config["busco_histogram_colors"],
    log:
        std=log_dir_path / "busco_histogram.log",
        cluster_log=cluster_log_dir_path / "busco_histogram.cluster.log",
        cluster_err=cluster_log_dir_path / "busco_histogram.cluster.err",
    benchmark:
        benchmark_dir_path / "busco_histogram.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        queue=config["processing_queue"],
        cpus=config["processing_threads"],
        time=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/busco_histogram.py -i {input} -o {output} -c '{params.colors}' > {log.std} 2>&1; "
