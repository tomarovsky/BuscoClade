if config["alignment_tool"] == "prank":
    rule prank_dna:
        input:
            merged_sequences_dir_path / "group_{N}"
        output:
            temp(directory(alignment_dir_path / "fna_tmp" / "group_{N}"))
        params:
            options=config["prank_dna_params"]
        log:
            std=log_dir_path / "prank_dna.{N}.log",
            cluster_log=cluster_log_dir_path / "prank_dna.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "prank_dna.{N}.cluster.err"
        benchmark:
            benchmark_dir_path / "prank_dna.{N}.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["prank_threads"],
            time=config["prank_time"],
            mem_mb=config["prank_mem_mb"]
        threads:
            config["prank_threads"]
        shell:
            "mkdir -p {output}; "
            "for FILE in `ls {input}/*.fna`; do "
            "prank -d=$FILE -o={output}/$(basename $FILE) {params.options} > {log.std} 2>&1; "
            "mv {output}/$(basename $FILE).best.fas {output}/$(basename $FILE); "
            "done; "