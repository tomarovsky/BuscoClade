from pathlib import Path
import os

#---- setup config ----
configfile: "config/default.yaml"

#---- setup paths ----
cluster_log_dir_path = Path(config["cluster_log_dir"])
genome_dir_path = Path(config["genome_dir"])
log_dir_path = Path(config["log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])
output_dir_path = Path(config["output_dir"])

busco_dir_path = output_dir_path / config["busco_dir"]
species_ids_dir_path = output_dir_path / config["species_ids_dir"]
common_ids_dir_path = output_dir_path / config["common_ids_dir"]
merged_sequences_dir_path = output_dir_path / config["merged_sequences_dir"]
alignment_dir_path = output_dir_path / config["alignment_dir"]
trimal_dir_path = output_dir_path / config["trimal_dir"]
concat_alignments_dir_path = output_dir_path / config["concat_alignments_dir"]
iqtree_dir_path = output_dir_path / config["iqtree_dir"]
mrbayes_dir_path = output_dir_path / config["mrbayes_dir"]

#---- setup filenames ----
alignment_file_prefix = config["alignment_file_prefix"]
fasta_dna_filename = f"{alignment_file_prefix}.fna"
fasta_protein_filename = f"{alignment_file_prefix}.faa"
nexus_dna_filename = f"{fasta_dna_filename}.nex"
nexus_protein_filename = f"{fasta_protein_filename}.nex"

if "species_list" not in config:
    config["species_list"] = [f.name[:-6] for f in genome_dir_path.iterdir() if f.is_file() and f.suffix == ".fasta"]

#---- necessary functions ----
def expand_template_from_common_ids(wildcards, template):
    checkpoint_output = checkpoints.common_ids.get(**wildcards).output[0]
    N = glob_wildcards(os.path.join(checkpoint_output, "group_{N}")).N
    return expand(str(template), N=N)

############################
# ---- the "all" rule ---- #
############################
#
# By default, all steps are run
#
output_files = [
    # busco:
    expand(busco_dir_path / "{species}/short_summary_{species}.txt",species=config["species_list"]),

    # common ids and merged sequences:
    lambda w: expand_template_from_common_ids(w,merged_sequences_dir_path / "group_{N}"),

    # mafft:
    alignment_dir_path / "fna",
    alignment_dir_path / "faa",

    # trimal:
    trimal_dir_path / "fna",
    trimal_dir_path / "faa",

    # concat alignments:
    concat_alignments_dir_path / fasta_dna_filename,
    concat_alignments_dir_path / fasta_protein_filename,
    concat_alignments_dir_path / nexus_dna_filename,
    concat_alignments_dir_path / nexus_protein_filename
]

phylogenetic_methods_dict = {
    "iqtree_dna" : iqtree_dir_path / "fna",
    "iqtree_protein": iqtree_dir_path / "faa",
    "mrbayes_dna" : mrbayes_dir_path / "fna",
    "mrbayes_protein" : mrbayes_dir_path / "faa"
}

if config["iqtree_dna_method"] in ["yes", "+", "true"]:
    output_files.append(phylogenetic_methods_dict["iqtree_dna"])
if config["iqtree_protein_method"] in ["yes", "+", "true"]:
    output_files.append(phylogenetic_methods_dict["iqtree_protein"])
if config["mrbayes_dna_method"] in ["yes", "+", "true"]:
    output_files.append(phylogenetic_methods_dict["mrbayes_dna"])
if config["mrbayes_protein_method"] in ["yes", "+", "true"]:
    output_files.append(phylogenetic_methods_dict["mrbayes_protein"])


localrules: all, files_transfer

rule all:
    input:
        output_files


rule files_transfer:
    input:
        mafft_fna_dirs=lambda w: expand_template_from_common_ids(w, alignment_dir_path / "fna_tmp" / "group_{N}"),
        mafft_faa_dirs=lambda w: expand_template_from_common_ids(w, alignment_dir_path / "faa_tmp" / "group_{N}"),
        trimal_fna_dirs=lambda w: expand_template_from_common_ids(w, trimal_dir_path / "fna_tmp" / "group_{N}"),
        trimal_faa_dirs=lambda w: expand_template_from_common_ids(w, trimal_dir_path / "faa_tmp" / "group_{N}")
    output:
        mafft_fna_dir=directory(alignment_dir_path / "fna"),
        mafft_faa_dir=directory(alignment_dir_path / "faa"),
        trimal_fna_dir=directory(trimal_dir_path / "fna"),
        trimal_faa_dir=directory(trimal_dir_path / "faa"),
    shell:
        "mkdir -p {output.mafft_fna_dir}; for i in {input.mafft_fna_dirs}; do mv $i/*.fna {output.mafft_fna_dir}/; done; "
        "mkdir -p {output.mafft_faa_dir}; for i in {input.mafft_faa_dirs}; do mv $i/*.faa {output.mafft_faa_dir}/; done; "
        "mkdir -p {output.trimal_fna_dir}; for i in {input.trimal_fna_dirs}; do mv $i/*.fna {output.trimal_fna_dir}/; done; "
        "mkdir -p {output.trimal_faa_dir}; for i in {input.trimal_faa_dirs}; do mv $i/*.faa {output.trimal_faa_dir}/; done; "


#---- load rules ----
include: "workflow/rules/busco.smk"
include: "workflow/rules/common_ids.smk"
include: "workflow/rules/mafft.smk"
include: "workflow/rules/prank.smk"
include: "workflow/rules/trimal.smk"
include: "workflow/rules/concat_alignments.smk"
include: "workflow/rules/iqtree.smk"
include: "workflow/rules/mrbayes.smk"