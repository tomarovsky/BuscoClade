rule phylip_tree:
    input:
        concat_alignments_dir_path / phylip_dna_filename
    output:
        treefile = phylip_dir_path / phylip_tree,
        dnadist = phylip_dir_path / f"{phylip_dna_filename}.dnadist",
        outfile = phylip_dir_path / f"{phylip_dna_filename}.outfile"
    params:
        outdir = phylip_dir_path,
        dnadist_params = config["phylip_dnadist_params"],
        neighbor_params = config["phylip_neighbor_params"]
    log:
        std=log_dir_path / "phylip_tree.log",
        cluster_log=cluster_log_dir_path / "phylip_tree.cluster.log",
        cluster_err=cluster_log_dir_path / "phylip_tree.cluster.err"
    benchmark:
        benchmark_dir_path / "phylip_tree.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["phylip_threads"],
        time=config["phylip_time"],
        mem_mb=config["phylip_mem_mb"]
    threads:
        config["phylip_threads"]
    shell:
        "MYPWD=$(pwd); cd {params.outdir} 2> $MYPWD/{log.std}; "
        "echo -e \"$MYPWD/{input}\n{params.dnadist_params}Y\n\" | dnadist > $MYPWD/{log.std} 2>&1; "
        "mv outfile $MYPWD/ 2> $MYPWD/{log.std}; "
        "cd $MYPWD 2> $MYPWD/{log.std}; "
        "mv outfile {output.dnadist} 2> $MYPWD/{log.std}; "
        "cd {params.outdir} 2> $MYPWD/{log.std}; "
        "echo -e \"$MYPWD/{output.dnadist}\n{params.neighbor_params}Y\n\" | neighbor > $MYPWD/{log.std} 2>&1; "
        "mv outfile $MYPWD/ 2> $MYPWD/{log.std}; "
        "cd $MYPWD 2> $MYPWD/{log.std}; "
        "mv outfile {output.outfile} 2> $MYPWD/{log.std}; "
        "cd {params.outdir} 2> $MYPWD/{log.std}; "
        "mv outtree $MYPWD/ 2> $MYPWD/{log.std}; "
        "cd $MYPWD 2> $MYPWD/{log.std}; "
        "mv outtree {output.treefile} 2> $MYPWD/{log.std}; "


