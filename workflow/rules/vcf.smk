if config["vcf"]:
    
    checkpoint vcf_discovery:
        input:
            vcf_dir_path
        output:
            touch(vcf_dir_path / "vcf_files.discovered")
        shell:
            "find {input} -name '*.vcf*' > {output}"
    
    if config["vcf_masking"]:
        rule vcf_masking:
            input:
                vcf_file=lambda w: f"vcf_dir_path/{w.vcf_directory}/{w.sample}.vcf.gz",
                mask_file=lambda w: f"vcf_dir_path/{w.vcf_directory}/mask.bed"
            output:
                vcf_dir_path / "{vcf_directory}/{sample}.masked.vcf.gz"
            params:
                type="DNA"
            log:
                std=log_dir_path / "{vcf_directory}/{sample}_vcf_masking.log",
                cluster_log=cluster_log_dir_path / "{vcf_directory}/{sample}_vcf_masking.cluster.log",
                cluster_err=cluster_log_dir_path / "{vcf_directory}/{sample}_vcf_masking.cluster.err",
            benchmark:
                benchmark_dir_path / "{vcf_directory}/{sample}_vcf_masking.benchmark.txt"
            conda:
                "../../%s" % config["conda_config"]
            resources:
                queue=config["processing_queue"],
                cpus=config["processing_threads"],
                time=config["processing_time"],
                mem_mb=config["processing_mem_mb"],
            shell:
                "bedtools intersect -header -v -a {input.vcf_file} -b {input.mask_file} > {output} 2> {log.std}"
                

    rule vcf_extract_sample_names:
        input:
            vcf=lambda w: f"vcf_dir_path/{w.vcf_directory}/sample.masked.vcf.gz",
            checkpoint_file=checkpoints.vcf_discovery.get().output[0]
        output:
            f"vcf_dir_path/{{vcf_directory}}/sample_list.txt",
        params:
            type="DNA"
        log:
            std=log_dir_path / "{vcf_directory}/vcf_extract_sample_names.log",
            cluster_log=cluster_log_dir_path / "{vcf_directory}/vcf_extract_sample_names.cluster.log",
            cluster_err=cluster_log_dir_path / "{vcf_directory}/vcf_extract_sample_names.cluster.err",
        benchmark:
            benchmark_dir_path / "{vcf_directory}/vcf_extract_sample_names.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "bcftools query -l {input.vcf} > {output}"


    rule vcf_separation:
        input:
            vcf=lambda w: f"vcf_dir_path/{w.vcf_directory}/sample.masked.vcf.gz",
            sample_list="vcf_dir_path/{vcf_directory}/sample_list.txt"
        output:
            directory(f"vcf_dir_path/{{vcf_directory}}/separated"),
        params:
            type="DNA"
        log:
            std=log_dir_path / "{vcf_directory}/vcf_separation.log",
            cluster_log=cluster_log_dir_path / "{vcf_directory}/vcf_separation.cluster.log",
            cluster_err=cluster_log_dir_path / "{vcf_directory}/vcf_separation.cluster.err",
        benchmark:
            benchmark_dir_path / "{vcf_directory}/vcf_separation.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            """
            mkdir -p {output}
            while read -r SAMPLE; do
                [ -z "$SAMPLE" ] && continue
                bcftools view --threads 10 --min-ac 1 \\
                    -s "$SAMPLE" \\
                    -O z -o {output}/"$SAMPLE".separated.vcf.gz \\
                    {input.vcf}
            done < {input.sample_list}
            """
    rule vcf_filtering:
        input:
            "vcf_dir_path/{vcf_directory}/separated/{sample}.separated.vcf.gz",
        output:
            "vcf_dir_path/{vcf_directory}/filtered/{sample}.filtered.vcf.gz",
        params:
            vcf_directory="{vcf_directory}",
            common_filters=config["common_filters"],
            haploid_filters=config["haploid_filters"],
            diploid_filters=config["diploid_filters"],
        log:
            std=log_dir_path / "{vcf_directory}/{sample}_vcf_filtering.log",  # Removed {ext}
            cluster_log=cluster_log_dir_path / "{vcf_directory}/{sample}_vcf_filtering.cluster.log",
            cluster_err=cluster_log_dir_path / "{vcf_directory}/{sample}_vcf_filtering.cluster.err",
        benchmark:
            benchmark_dir_path / "{vcf_directory}/{sample}_vcf_filtering.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            """
            mkdir -p {vcf_dir_path}/{wildcards.vcf_directory}/filtered
            ploidy=$(bcftools query -f '%%GT\\n' {input} | \\
                awk '$1 ~ /\// || $1 ~ /\|/ {{print "diploid"; exit}} END {{print "haploid"}}')

            if [ "$ploidy" = "diploid" ]; then
                bcftools filter -i {params.common_filters} \\
                -e {params.diploid_filters} {input} > {output}
            else
                bcftools filter -i {params.common_filters} \\
                -e {params.haploid_filters} {input} > {output}
            fi
            """
    
    rule vcf_sequence_dictionary:
        input:
            "vcf_dir_path/{vcf_directory}/{reference}.{ext}",  # Simplified pattern
        output:
            "vcf_dir_path/{vcf_directory}/{reference}.dict"
        log:
            std=log_dir_path / "{vcf_directory}/{reference}_vcf_sequence_dictionary.log",
            cluster_log=cluster_log_dir_path / "{vcf_directory}/{reference}_vcf_sequence_dictionary.cluster.log",
            cluster_err=cluster_log_dir_path / "{vcf_directory}/{reference}_vcf_sequence_dictionary.cluster.err",
        benchmark:
            benchmark_dir_path / "{vcf_directory}/{reference}_vcf_sequence_dictionary.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "java -Xmx5g -jar picard.jar CreateSequenceDictionary R={input} O={output} 2> {log.std}"
    

    rule vcf_indexing:
        input:
            "vcf_dir_path/{vcf_directory}/filtered/{sample}.filtered.{ext}",
        output:
            "vcf_dir_path/{vcf_directory}/filtered/{sample}.filtered.{ext}.idx",
        log:
            std=log_dir_path / "{vcf_directory}/{sample}.filtered.{ext}_vcf_indexing.log",
            cluster_log=cluster_log_dir_path / "{vcf_directory}/{sample}.filtered.{ext}_vcf_indexing.cluster.log",  # Fixed typo
            cluster_err=cluster_log_dir_path / "{vcf_directory}/{sample}.filtered.{ext}_vcf_indexing.cluster.err",
        benchmark:
            benchmark_dir_path / "{vcf_directory}/{sample}.filtered.{ext}_vcf_indexing.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "java -Xmx10g -jar gatk IndexFeatureFile -I {input} 1> {log.std} 2>&1"
    
    rule reference_indexing:
        input:
            "vcf_dir_path/{vcf_directory}/{reference}.{ext}",
        output:
            "vcf_dir_path/{vcf_directory}/{reference}.{ext}.fai",
        log:
            std=log_dir_path / "{vcf_directory}/{reference}.{ext}_reference_indexing.log",
            cluster_log=cluster_log_dir_path / "{vcf_directory}/{reference}.{ext}_reference_indexing.cluster.log",
            cluster_err=cluster_log_dir_path / "{vcf_directory}/{reference}.{ext}_reference_indexing.cluster.err",
        benchmark:
            benchmark_dir_path / "{vcf_directory}/{reference}.{ext}_reference_indexing.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "samtools faidx {input} 2> {log.std}"

    rule vcf_to_fasta:
        input:
            vcf_file="vcf_dir_path/{vcf_directory}/filtered/{sample}.filtered.vcf.gz",
            reference_file="vcf_dir_path/{vcf_directory}/reference.fna"
        output:
            genome_dir_path / "{sample}.AltRef.fasta",  # Keep flat structure
        params:
            vcf_directory=lambda w: os.path.basename(os.path.dirname(os.path.dirname(input.vcf_file))),
        log:
            std=log_dir_path / "{sample}_vcf_to_fasta.log",
            cluster_log=cluster_log_dir_path / "{sample}_vcf_to_fasta.cluster.log",
            cluster_err=cluster_log_dir_path / "{sample}_vcf_to_fasta.cluster.err",
        benchmark:
            benchmark_dir_path / "{sample}_vcf_to_fasta.benchmark.txt"
        params:
            vcf_directory="{sample}"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            queue=config["processing_queue"],
            cpus=config["processing_threads"],
            time=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        shell:
            "java -Xmx10g -jar gatk FastaAlternateReferenceMaker "
            "--output {output} "
            "--reference {input.reference_file} "
            "--variant {input.vcf_file} "
            "--showHidden true "
            "--use-iupac-sample {wildcards.sample} 2> {log.std}"