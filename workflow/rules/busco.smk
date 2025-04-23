if config["busco_gene_prediction_tool"] == "metaeuk":

    rule busco_metaeuk:
        input:
            genome_dir_path / "{species}.fasta",
        output:
            busco_outdir=directory(busco_dir_path / "{species}"),
            single_copy_busco_sequences=directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences"),
            multi_copy_busco_sequences=directory(busco_dir_path / "{species}/busco_sequences/multi_copy_busco_sequences"),
            summary=busco_dir_path / "{species}/short_summary_{species}.txt",
        params:
            mode=config["busco_mode"],
            busco_dataset_path=config["busco_dataset_path"],
            busco_options=config["busco_options"],
            output_prefix="{species}",
        log:
            std=log_dir_path / "busco.{species}.log",
            cluster_log=cluster_log_dir_path / "busco.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "busco.{species}.cluster.err",
        benchmark:
            benchmark_dir_path / "busco.{species}.benchmark.txt"
        conda:
            config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
        resources:
            queue=config["busco_queue"],
            cpus=config["busco_threads"],
            time=config["busco_time"],
            mem_mb=config["busco_mem_mb"],
        threads: config["busco_threads"]
        shell:
            " MYPWD=$(pwd); mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; "
            " busco --metaeuk -m {params.mode} -i {input} -c {threads} -l {params.busco_dataset_path} "
            " -o {params.output_prefix} {params.busco_options} 1> $MYPWD/{log.std} 2>&1; "
            " mv {params.output_prefix}/* . 1> $MYPWD/{log.std} 2>&1; "
            " rm -r {params.output_prefix}/ 1> $MYPWD/{log.std} 2>&1; "
            " rm -r busco_sequences/ 1> $MYPWD/{log.std} 2>&1; "
            " mv run*/* . 1> $MYPWD/{log.std} 2>&1; "
            " rm -r run* 1> $MYPWD/{log.std} 2>&1; "
            " mv full_table.tsv full_table_{params.output_prefix}.tsv 1> $MYPWD/{log.std} 2>&1; "
            " mv missing_busco_list.tsv missing_busco_list_{params.output_prefix}.tsv 1> $MYPWD/{log.std} 2>&1; "
            " mv short_summary.txt short_summary_{params.output_prefix}.txt 1> $MYPWD/{log.std} 2>&1; "


elif config["busco_gene_prediction_tool"] == "augustus":

    rule busco_augustus:
        input:
            genome_dir_path / "{species}.fasta",
        output:
            busco_outdir=directory(busco_dir_path / "{species}"),
            single_copy_busco_sequences=directory(busco_dir_path / "{species}/single_copy_busco_sequences"),
            multi_copy_busco_sequences=directory(busco_dir_path / "{species}/multi_copy_busco_sequences"),
            summary=busco_dir_path / "{species}/short_summary_{species}.txt",
        params:
            mode=config["busco_mode"],
            species=config["busco_augustus_species"],
            busco_dataset_path=config["busco_dataset_path"],
            busco_options=config["busco_options"],
            output_prefix="{species}",
        log:
            std=log_dir_path / "busco.{species}.log",
            cluster_log=cluster_log_dir_path / "busco.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "busco.{species}.cluster.err",
        benchmark:
            benchmark_dir_path / "busco.{species}.benchmark.txt"
        conda:
            config["conda"]["buscoclade"]["name"] if config["use_existing_envs"] else ("../../%s" % config["conda"]["buscoclade"]["yaml"])
        resources:
            queue=config["busco_queue"],
            cpus=config["busco_threads"],
            time=config["busco_time"],
            mem_mb=config["busco_mem_mb"],
        threads: config["busco_threads"]
        shell:
            " MYPWD=$(pwd); mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; "
            " busco --augustus --augustus_species {params.species} -m {params.mode} -i {input} -c {threads} "
            " -l {params.busco_dataset_path} -o {params.output_prefix} {params.busco_options} 1> $MYPWD/{log.std} 2>&1; "
            " mv {params.output_prefix}/* ./ 1> $MYPWD/{log.std} 2>&1; "
            " rm -r {params.output_prefix}/ 1> $MYPWD/{log.std} 2>&1; "
            " rm -r augustus_output/ 1> $MYPWD/{log.std} 2>&1; "
            " mv run*/* . ; rm -r run* 1> $MYPWD/{log.std} 2>&1; "
            " mv full_table.tsv full_table_{params.output_prefix}.tsv 1> $MYPWD/{log.std} 2>&1; "
            " mv missing_busco_list.tsv missing_busco_list_{params.output_prefix}.tsv 1> $MYPWD/{log.std} 2>&1; "
            " mv short_summary.txt short_summary_{params.output_prefix}.txt 1> $MYPWD/{log.std} 2>&1; "
