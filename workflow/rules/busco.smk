if config['busco_version'] == 3:
    rule busco3:
        input:
            genome_dir_path / "{species}.fasta"
        output:
            busco_outdir=protected(directory(busco_dir_path / "{species}")),
            single_copy_files_dir=protected(directory(busco_dir_path / "{species}/busco_sequences//single_copy_busco_sequences")),
            summary=protected(busco_dir_path / "{species}/short_summary_{species}.txt")
        params:
            busco_path=config["busco_path"],
            mode=config["busco_mode"],
            species=config["augustus_species"],
            busco_dataset_path=config["busco_dataset_path"],
            output_prefix="{species}"
        log:
            std=log_dir_path / "busco.{species}.log",
            cluster_log=cluster_log_dir_path / "busco.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "busco.{species}.cluster.err"
        benchmark:
            benchmark_dir_path / "busco.{species}.benchmark.txt"
        resources:
            cpus=config["busco_threads"],
            time=config["busco_time"],
            mem_mb=config["busco_mem_mb"],
        threads:
            config["busco_threads"]
        shell:
            "mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; {params.busco_path}/run_BUSCO.py -m {params.mode} -sp {params.species}"
            " -i {input} -c {threads} -l {params.busco_dataset_path} -o {params.output_prefix} 1>../../../{log.std} 2>&1;"
            " mv run_{params.output_prefix}/* ./; rm -r run_{params.output_prefix} tmp/; "
            "mkdir busco_sequences/; mv single_copy_busco_sequences/ busco_sequences/ "

elif config['busco_version'] == 5:
    if config['gene_prediction_tool'] == "metaeuk":
        rule busco5_metaeuk:
            input:
                genome_dir_path / "{species}.fasta"
            output:
                busco_outdir=protected(directory(busco_dir_path / "{species}")),
                single_copy_busco_sequences=protected(directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences")),
                summary=protected(busco_dir_path / "{species}/short_summary_{species}.txt")
            params:
                mode=config["busco_mode"],
                busco_dataset_path=config["busco_dataset_path"],
                output_prefix="{species}"
            log:
                std=log_dir_path / "busco.{species}.log",
                cluster_log=cluster_log_dir_path / "busco.{species}.cluster.log",
                cluster_err=cluster_log_dir_path / "busco.{species}.cluster.err"
            benchmark:
                benchmark_dir_path / "busco.{species}.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                cpus=config["busco_threads"],
                time=config["busco_time"],
                mem_mb=config["busco_mem_mb"],
            threads:
                config["busco_threads"]
            shell:
                "mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; "
                "busco -m {params.mode} -i {input} -c {threads} "
                "-l {params.busco_dataset_path} -o {params.output_prefix} 1>../../../{log.std} 2>&1; "
                "mv {params.output_prefix}/* . 1>../../../{log.std} 2>&1; "
                "rm -r {params.output_prefix}/ 1>../../../{log.std} 2>&1; "
                "rm -r busco_sequences/ 1>../../../{log.std} 2>&1; " # empty directory
                "mv run*/* . 1>../../../{log.std} 2>&1; "
                "rm -r run* 1>../../../{log.std} 2>&1; "
                "mv full_table.tsv full_table_{params.output_prefix}.tsv 1>../../../{log.std} 2>&1; "
                "mv missing_busco_list.tsv missing_busco_list_{params.output_prefix}.tsv 1>../../../{log.std} 2>&1; "
                "mv short_summary.txt short_summary_{params.output_prefix}.txt 1>../../../{log.std} 2>&1; "

    elif config['gene_prediction_tool'] == "augustus":
        rule busco5_augustus:
            input:
                genome_dir_path / "{species}.fasta"
            output:
                busco_outdir=protected(directory(busco_dir_path / "{species}")),
                single_copy_busco_sequences=protected(directory(busco_dir_path / "{species}/single_copy_busco_sequences")),
                augustus_gff=protected(directory(busco_dir_path / "{species}/augustus_output/gff")),
                summary=protected(busco_dir_path / "{species}/short_summary_{species}.txt")
            params:
                mode=config["busco_mode"],
                species=config["augustus_species"],
                busco_dataset_path=config["busco_dataset_path"],
                output_prefix="{species}"
            log:
                std=log_dir_path / "busco.{species}.log",
                cluster_log=cluster_log_dir_path / "busco.{species}.cluster.log",
                cluster_err=cluster_log_dir_path / "busco.{species}.cluster.err"
            benchmark:
                benchmark_dir_path / "busco.{species}.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                cpus=config["busco_threads"],
                time=config["busco_time"],
                mem_mb=config["busco_mem_mb"],
            threads:
                config["busco_threads"]
            shell:
                "mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; "
                "busco --augustus --augustus_species {params.species} -m {params.mode} "
                "-i {input} -c {threads} -l {params.busco_dataset_path} -o {params.output_prefix} 1>../../../{log.std} 2>&1; "
                "mv {params.output_prefix}/* ./ ; rm -r {params.output_prefix}/ ; "
                "rm -r augustus_output/ ; " # empty directory
                "mv run*/* ./ ; rm -r run* ; "
                "mv full_table.tsv full_table_{params.output_prefix}.tsv ; "
                "mv missing_busco_list.tsv missing_busco_list_{params.output_prefix}.tsv ; "
                "mv short_summary.txt short_summary_{params.output_prefix}.txt ; "
    else:
        print("Specify the tool name in 'gene_prediction_tool' parameter! Use 'metaeuk' or 'augustus'")
else:
    print("Specify the version of BUSCO in 'busco_version' parameter! Use '3' or '5'")