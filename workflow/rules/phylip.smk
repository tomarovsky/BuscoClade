localrules:
    phylip_tree_namefix,


rule phylip_dnadist:
    input:
        concat_alignments_dir_path / phylip_filename,
    output:
        dnadist=phylip_dir_path / f"{phylip_filename}.dnadist",
    params:
        outdir=phylip_dir_path,
        dnadist_params=config["phylip_dnadist_params"],
    log:
        std=log_dir_path / "phylip_dnadist.log",
        cluster_log=cluster_log_dir_path / "phylip_dnadist.cluster.log",
        cluster_err=cluster_log_dir_path / "phylip_dnadist.cluster.err",
    benchmark:
        benchmark_dir_path / "phylip_dnadist.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        slurm_partition=config["phylip_queue"],
        runtime=config["phylip_time"],
        mem_mb=config["phylip_mem_mb"],
    threads: config["phylip_threads"]
    shell:
        " MYPWD=$(pwd); cd {params.outdir} 2> $MYPWD/{log.std}; "
        ' echo -e "$MYPWD/{input}\n{params.dnadist_params}Y\n" | dnadist > $MYPWD/{log.std} 2>&1; '
        " mv outfile $MYPWD/ 2> $MYPWD/{log.std}; "
        " cd $MYPWD 2> $MYPWD/{log.std}; "
        " mv outfile {output.dnadist} 2> $MYPWD/{log.std}; "


rule phylip_neighbor:
    input:
        phylip_dir_path / f"{phylip_filename}.dnadist",
    output:
        treefile=phylip_dir_path / f"{phylip_filename}.treefile",
        outfile=phylip_dir_path / f"{phylip_filename}.outfile",
    params:
        outdir=phylip_dir_path,
        neighbor_params=config["phylip_neighbor_params"],
    log:
        std=log_dir_path / "phylip_neighbor.log",
        cluster_log=cluster_log_dir_path / "phylip_neighbor.cluster.log",
        cluster_err=cluster_log_dir_path / "phylip_neighbor.cluster.err",
    benchmark:
        benchmark_dir_path / "phylip_neighbor.benchmark.txt"
    conda:
        config["conda"]["buscoclade_main"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade_main"]["yaml"])
    resources:
        slurm_partition=config["phylip_queue"],
        runtime=config["phylip_time"],
        mem_mb=config["phylip_mem_mb"],
    threads: config["phylip_threads"]
    shell:
        " MYPWD=$(pwd); cd {params.outdir} 2> $MYPWD/{log.std}; "
        ' echo -e "$MYPWD/{input}\n{params.neighbor_params}Y\n" | neighbor > $MYPWD/{log.std} 2>&1; '
        " mv outfile $MYPWD/ 2> $MYPWD/{log.std}; "
        " cd $MYPWD 2> $MYPWD/{log.std}; "
        " mv outfile {output.outfile} 2> $MYPWD/{log.std}; "
        " cd {params.outdir} 2> $MYPWD/{log.std}; "
        " mv outtree $MYPWD/ 2> $MYPWD/{log.std}; "
        " cd $MYPWD 2> $MYPWD/{log.std}; "
        " mv outtree {output.treefile} 2> $MYPWD/{log.std}; "


rule phylip_tree_namefix:
    input:
        phy=phylip_dir_path / f"{phylip_filename}.treefile",
        map=concat_alignments_dir_path / f"{phylip_filename}.map",
    output:
        phylip_dir_path / phylip_tree,
    log:
        std=log_dir_path / "phylip_tree_namefix.log",
        cluster_log=cluster_log_dir_path / "phylip_tree_namefix.cluster.log",
        cluster_err=cluster_log_dir_path / "phylip_tree_namefix.cluster.err",
    benchmark:
        benchmark_dir_path / "phylip_tree_namefix.benchmark.txt"
    resources:
        slurm_partition=config["processing_queue"],
        runtime=config["processing_time"],
        mem_mb=config["processing_mem_mb"],
    threads: config["processing_threads"]
    shell:
        " workflow/scripts/phylip_tree_namefix.py -i {input.phy} -m {input.map} -o {output} > {log.std} 2>&1; "
