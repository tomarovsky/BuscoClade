rule trimal_dna:
    input:
        alignment_dir_path / "fna_tmp" / "group_{N}"
    output:
        temp(directory(trimal_dir_path / "fna_tmp" /"group_{N}"))
    params:
        trimal_path=config["trimal_path"],
        options=config["trimal_dna_params"]
    log:
        std=log_dir_path / "trimal_dna.{N}.log",
        cluster_log=cluster_log_dir_path / "trimal_dna.{N}.cluster.log",
        cluster_err=cluster_log_dir_path / "trimal_dna.{N}.cluster.err"
    benchmark:
        benchmark_dir_path / "trimal_dna.{N}.benchmark.txt"
    resources:
        cpus=config["trimal_threads"],
        time=config["trimal_time"],
        mem_mb=config["trimal_mem_mb"]
    shell:
        "mkdir -p {output}; for FILE in `ls {input}/*.fna`; do "
        "{params.trimal_path}/trimal -in $FILE -out {output}/$(basename $FILE) {params.options} > {log.std} 2>&1; "
        "{params.trimal_path}/trimal -in {output}/$(basename $FILE) -out {output}/$(basename $FILE) -nogaps > {log.std} 2>&1; "
        "done"


rule trimal_protein:
    input:
        alignment_dir_path / "faa_tmp" / "group_{N}"
    output:
        temp(directory(trimal_dir_path / "faa_tmp" /"group_{N}"))
    params:
        trimal_path=config["trimal_path"],
        options=config["trimal_protein_params"]
    log:
        std=log_dir_path / "trimal_protein.{N}.log",
        cluster_log=cluster_log_dir_path / "trimal_protein.{N}.cluster.log",
        cluster_err=cluster_log_dir_path / "trimal_protein.{N}.cluster.err"
    benchmark:
        benchmark_dir_path / "trimal_protein.{N}.benchmark.txt"
    resources:
        cpus=config["trimal_threads"],
        time=config["trimal_time"],
        mem_mb=config["trimal_mem_mb"]
    shell:
        "mkdir -p {output}; for FILE in `ls {input}/*`; do "
        "{params.trimal_path}/trimal -in $FILE -out {output}/$(basename $FILE) {params.options} > {log.std} 2>&1; "
        "{params.trimal_path}/trimal -in {output}/$(basename $FILE) -out {output}/$(basename $FILE) -nogaps > {log.std} 2>&1; "
        "done"