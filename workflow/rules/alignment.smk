# ---- Gap-aware AltRef insertion setup ----
# altref_gapaware and ref_to_altrefs are expected to be defined in Snakefile:
#   altref_gapaware = config.get("altref_gapaware_insertion", False)
#   ref_to_altrefs  = {ref_prefix: [altref_sp, ...], ...}  (built from altref_map)
#
# In gap-aware mode, merged_sequences/ contains only non-AltRef species (refs + FASTA
# assemblies). AltRef sequences are inserted into the finished alignment afterwards by
# add_altref_to_alignment.py. In standard mode, merged_sequences/ contains all species.
#
# IMPORTANT: input must be a lambda so Snakemake defers file-existence checks
# until after the checkpoint has run. A static Path is evaluated at DAG-build
# time — before merged_sequences/ exists — causing MissingInputException.
def _aligner_input(wildcards):
    checkpoints.merged_sequences.get(**wildcards)
    return str(merged_sequences_dir_path / f"{wildcards.N}.fna")

_aligner_output_dir = (
    pre_altref_alignments_dir_path
    if altref_gapaware and altref_map
    else alignments_dir_path / "fna"
)

# ---- Aligners ----

if config.get("alignment") == "prank":

    def parse_time_to_minutes(time_str):
        return sum(
            int(v) * {"h": 60, "m": 1}.get(u, 1)
            for v, u in re.findall(r"(\d+)([hm]?)", str(time_str))
            if v
        )

    rule prank:
        input:
            _aligner_input,
        output:
            _aligner_output_dir / "{N}.fna",
        params:
            prefix=lambda wildcards, output: output[0][:-4],
            options=config["prank_params"],
            timeout=lambda wildcards: f"{parse_time_to_minutes(config['prank_time']) - 15}m",
        log:
            std=log_dir_path / "prank.{N}.log",
            cluster_log=cluster_log_dir_path / "prank.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "prank.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "prank.{N}.benchmark.txt"
        conda:
            main_env
        resources:
            slurm_partition=config["alignment_queue"],
            runtime=config["prank_time"],
            mem_mb=config["prank_mem_mb"],
        threads: config["prank_threads"]
        shell:
            " timeout --kill-after=5m {params.timeout} "
            " prank -d={input} -o={params.prefix} {params.options} 1> {log.std} 2>&1 "
            " || {{ exit_code=$?; "
            " if [ $exit_code -ne 0 ]; then "
            "     echo 'PRANK failed or timed out (exit code: '$exit_code') for {wildcards.N}' >> {log.std}; "
            "     touch {output}; "
            " fi; "
            " }}; "
            " if [ ! -f {output} ]; then "
            "     mv {params.prefix}.best.fas {output} >> {log.std} 2>&1; "
            " fi; "


if config["alignment"] == "mafft":

    rule mafft:
        input:
            _aligner_input,
        output:
            _aligner_output_dir / "{N}.fna",
        params:
            options=config["mafft_params"],
        log:
            std=log_dir_path / "mafft.{N}.log",
            cluster_log=cluster_log_dir_path / "mafft.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "mafft.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "mafft.{N}.benchmark.txt"
        conda:
            main_env
        resources:
            slurm_partition=config["alignment_queue"],
            runtime=config["mafft_time"],
            mem_mb=config["mafft_mem_mb"],
        threads: config["mafft_threads"]
        shell:
            " mafft --thread {threads} {params.options} {input} > {output} 2> {log.std}; "


if config["alignment"] == "muscle":

    rule muscle:
        input:
            _aligner_input,
        output:
            _aligner_output_dir / "{N}.fna",
        params:
            options=config["muscle_params"],
        log:
            std=log_dir_path / "muscle.{N}.log",
            cluster_log=cluster_log_dir_path / "muscle.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "muscle.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "muscle.{N}.benchmark.txt"
        conda:
            main_env
        resources:
            slurm_partition=config["alignment_queue"],
            runtime=config["muscle_time"],
            mem_mb=config["muscle_mem_mb"],
        threads: config["muscle_threads"]
        shell:
            " muscle -in {input} -out {output} {params.options} > {log.std} 2>&1; "


# ---- Gap-aware AltRef insertion ----
# Active only when altref_gapaware_insertion: True and there are AltRef species.
# Reads raw alignment (refs + non-AltRef species), inserts each AltRef sequence
# by copying gap positions from its corresponding ref, then optionally removes
# ref sequences if vcf_reconstruct_ref_as_species: False.

if altref_gapaware and altref_map:

    import json as _json

    rule add_altref_to_alignment:
        input:
            raw_aln=pre_altref_alignments_dir_path / "{N}.fna",
            # Declare dependency on all AltRef busco dirs so Snakemake
            # knows to wait for apply_vcf_to_busco to finish.
            altref_seqs=expand(
                str(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences"),
                species=altref_species,
            ),
        output:
            alignments_dir_path / "fna" / "{N}.fna",
        params:
            busco_dir=str(busco_dir_path),
            ref_to_altrefs=_json.dumps(ref_to_altrefs),
            keep_refs="--keep_refs" if config.get("vcf_reconstruct_ref_as_species") else "",
        log:
            std=log_dir_path / "add_altref_to_alignment.{N}.log",
            cluster_log=cluster_log_dir_path / "add_altref_to_alignment.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "add_altref_to_alignment.{N}.cluster.err",
        benchmark:
            benchmark_dir_path / "add_altref_to_alignment.{N}.benchmark.txt"
        conda:
            main_env
        resources:
            slurm_partition=config["processing_queue"],
            runtime=config["processing_time"],
            mem_mb=config["processing_mem_mb"],
        threads: config["processing_threads"]
        shell:
            " workflow/scripts/add_altref_to_alignment.py "
            " --raw_aln {input.raw_aln} "
            " --busco_dir {params.busco_dir} "
            " --ref_to_altrefs '{params.ref_to_altrefs}' "
            " --gene_id {wildcards.N} "
            " {params.keep_refs} "
            " --output {output} "
            " 1> {log.std} 2>&1; "
