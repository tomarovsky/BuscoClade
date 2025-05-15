if config["vcf2phylip"] != True:
    localrules:
        species_single_copy_ids,
        species_multi_copy_ids,
        common_ids,


    rule species_single_copy_ids:  # get single copy ids for each species
        input:
            busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences",
        output:
            species_ids_dir_path / "single_copy/{species}.ids",
        log:
            std=log_dir_path / "species_ids.{species}.log",
            cluster_log=cluster_log_dir_path / "species_ids.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "species_ids.{species}.cluster.err",
        benchmark:
            benchmark_dir_path / "species_ids.{species}.benchmark.txt"
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"]
        shell:
            " ls {input} | grep -P '.fna$' | sed 's/.fna//' > {output} 2> {log.std}; "


    rule species_multi_copy_ids:  # get multi copy ids for each species
        input:
            busco_dir_path / "{species}/busco_sequences/multi_copy_busco_sequences",
        output:
            species_ids_dir_path / "multi_copy/{species}.ids",
        log:
            std=log_dir_path / "species_ids.{species}.log",
            cluster_log=cluster_log_dir_path / "species_ids.{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "species_ids.{species}.cluster.err",
        benchmark:
            benchmark_dir_path / "species_ids.{species}.benchmark.txt"
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"]
        shell:
            " ls {input} | grep -P '.fna$' | sed 's/.fna//' > {output} 2> {log.std}; "


    rule common_ids:  # get common IDs for all species given config["gene_blacklist"] and split them into files
        input:
            expand(species_ids_dir_path / "single_copy/{species}.ids", species=config["species_list"]),
        output:
            common_ids_dir_path / "common.ids",
        params:
            nfiles=len(config["species_list"]),
            gene_blacklist=config["gene_blacklist"],
        log:
            std=log_dir_path / "common_ids.log",
            cluster_log=cluster_log_dir_path / "common_ids.cluster.log",
            cluster_err=cluster_log_dir_path / "common_ids.cluster.err",
        benchmark:
            benchmark_dir_path / "common_ids.benchmark.txt"
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"]
        shell:
            " cat {input} | sort | uniq -c | awk '{{if($1=={params.nfiles}){{print $2}}}}' > {output} 2> {log.std}; "
            ' for id in {params.gene_blacklist}; do sed -i "/$id$/d" {output}; done 2> {log.std}; '


    checkpoint merged_sequences:  # get merged sequences by common IDs
        input:
            rules.common_ids.output,
        output:
            directory(merged_sequences_dir_path),
        params:
            single_copy_files=expand(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences", species=config["species_list"]),
        log:
            std=log_dir_path / "merged_sequences.log",
            cluster_log=cluster_log_dir_path / "merged_sequences.cluster.log",
            cluster_err=cluster_log_dir_path / "merged_sequences.cluster.err",
        benchmark:
            benchmark_dir_path / "merged_sequences.benchmark.txt"
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"]
        shell:
            " workflow/scripts/merge_common_ids.py -c {input} -s {params.single_copy_files} -o {output} > {log.std} 2>&1; "
