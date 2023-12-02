if "astral" in config:
    rule iqtree_per_fna:
        input:
            filtered_alignments_dir_path / "fna" / "{N}.fna"
        output:
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.bionj",
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.ckp.gz",
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.iqtree",
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.log",
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.mldist",
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.model.gz",
            astral_dir_path / "iqtree_per_fna" / "{N}.fna.treefile"
        params:
            options=config["iqtree_per_fna_params"],
            outdir=astral_dir_path / "iqtree_per_fna",
            prefix=lambda wildcards, output: output[0][:-6],
        log:
            std=log_dir_path / "iqtree_per_fna.{N}.log",
            cluster_log=cluster_log_dir_path / "iqtree_per_fna.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "iqtree_per_fna.{N}.cluster.err"
        benchmark:
            benchmark_dir_path / "iqtree_per_fna.{N}.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["iqtree_per_fna_threads"],
            time=config["iqtree_per_fna_time"],
            mem_mb=config["iqtree_per_fna_mem_mb"]
        threads:
            config["iqtree_per_fna_threads"]
        shell:
            " mkdir -p {params.outdir}; "
            " iqtree -nt {threads} -s {input} --prefix {params.prefix} {params.options} 1> {log.std} 2>&1; "

    rule concat_newick_files:
        input:
            lambda w: expand_fna_from_merged_sequences(w, astral_dir_path / "iqtree_per_fna" / "{N}.fna.treefile")
        output:
            astral_dir_path / astral_input_trees
        log:
            std=log_dir_path / "concat_newick_files.log",
            cluster_log=cluster_log_dir_path / "concat_newick_files.cluster.log",
            cluster_err=cluster_log_dir_path / "concat_newick_files.cluster.err"
        benchmark:
            benchmark_dir_path / "concat_newick_files.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["concat_newick_files_threads"],
            time=config["concat_newick_files_time"],
            mem_mb=config["concat_newick_files_mem_mb"]
        threads:
            config["concat_newick_files_threads"]
        shell:
            " cat {input} > {output} 2> {log.std}; "

    rule nodes_filtrataion_by_support:
        input:
            astral_dir_path / astral_input_trees
        output:
            astral_dir_path / astral_filtered_trees
        params:
            support=config["nodes_filtrataion_by_support"]
        log:
            std=log_dir_path / "nodes_filtrataion_by_support.log",
            cluster_log=cluster_log_dir_path / "nodes_filtrataion_by_support.cluster.log",
            cluster_err=cluster_log_dir_path / "nodes_filtrataion_by_support.cluster.err"
        benchmark:
            benchmark_dir_path / "nodes_filtrataion_by_support.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["concat_newick_files_threads"],
            time=config["concat_newick_files_time"],
            mem_mb=config["concat_newick_files_mem_mb"]
        threads:
            config["concat_newick_files_threads"]
        shell:
            " nw_ed {input} 'i & b<{params.support}' o > {output} 2> {log.std}; "

    rule astral_tree:
        input:
            astral_dir_path / astral_filtered_trees
        output:
            astral_dir_path / astral_tree
        params:
            config["astral_params"]
        log:
            std=log_dir_path / "astral_tree.log",
            cluster_log=cluster_log_dir_path / "astral_tree.cluster.log",
            cluster_err=cluster_log_dir_path / "astral_tree.cluster.err"
        benchmark:
            benchmark_dir_path / "astral_tree.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["astral_threads"],
            time=config["astral_time"],
            mem_mb=config["astral_mem_mb"]
        threads:
            config["astral_threads"]
        shell:
            " java -jar $CONDA_PREFIX/share/astral-tree-5.7.8-0/astral.5.7.8.jar " # 'astral' conda bin file is broken
            " -i {input} -o {output} {params} 1> {log.std} 2>&1"


