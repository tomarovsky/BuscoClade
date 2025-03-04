from pathlib import Path
import os

import pandas as pd

# ---- Setup config ----
configfile: "config/default.yaml"


# ---- Setup paths ----
cluster_log_dir_path = Path(config["cluster_log_dir"])
genome_dir_path = Path(config["genome_dir"]).resolve()
log_dir_path = Path(config["log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])
output_dir_path = Path(config["output_dir"])

quastcore_dir_path = output_dir_path / config["quastcore_dir"]
busco_dir_path = output_dir_path / config["busco_dir"]
species_ids_dir_path = output_dir_path / config["species_ids_dir"]
common_ids_dir_path = output_dir_path / config["common_ids_dir"]
merged_sequences_dir_path = output_dir_path / config["merged_sequences_dir"]
alignments_dir_path = output_dir_path / config["alignments_dir"]
filtered_alignments_dir_path = output_dir_path / config["filtered_alignments_dir"]
concat_alignments_dir_path = output_dir_path / config["concat_alignments_dir"]
iqtree_dir_path = output_dir_path / config["iqtree_dir"]
mrbayes_dir_path = output_dir_path / config["mrbayes_dir"]
astral_dir_path = output_dir_path / config["astral_dir"]
rapidnj_dir_path = output_dir_path / config["rapidnj_dir"]
phylip_dir_path = output_dir_path / config["phylip_dir"]

if "species_list" not in config:
    config["species_list"] = [f.stem for f in genome_dir_path.iterdir() if f.is_file() and f.suffix == ".fasta"]

# ---- Setup filenames ----
fasta_dna_filename = "{}.fna".format(config["alignment_file_prefix"])
fasta_protein_filename = "{}.faa".format(config["alignment_file_prefix"])
nexus_dna_filename = "{}.fna.nex".format(config["alignment_file_prefix"])
nexus_protein_filename = "{}.faa.nex".format(config["alignment_file_prefix"])
phylip_dna_filename = "{}.fna.phy".format(config["alignment_file_prefix"])
phylip_protein_filename = "{}.faa.phy".format(config["alignment_file_prefix"])
stockholm_dna_filename = "{}.fna.sth".format(config["alignment_file_prefix"])
stockholm_protein_filename = "{}.faa.sth".format(config["alignment_file_prefix"])
astral_input_trees = "{}.iqtree_per_fna.concat.treefile".format(config["alignment_file_prefix"])
astral_filtered_trees = "{0}.iqtree_per_fna.concat.{1}.treefile".format(config["alignment_file_prefix"], config["nodes_filtrataion_by_support"])
astral_tree = "{0}.{1}.fna.astral.treefile".format(config["alignment_file_prefix"], config["nodes_filtrataion_by_support"])
rapidnj_tree = "{}.fna.rapidnj.treefile".format(config["alignment_file_prefix"])
phylip_tree = "{}.fna.phy.namefix.treefile".format(config["alignment_file_prefix"])


# ---- Necessary functions ----
def expand_fna_from_merged_sequences(wildcards, template, busco_blacklist=None):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    N = glob_wildcards(os.path.join(checkpoint_output, "{N}.fna")).N
    print("test")
    print(busco_blacklist)
    print(len(N))
    print(N[0])
    if busco_blacklist is not None:
        N = set(N) - set(list(map(lambda s: f"{s}.merged", busco_blacklist)))
    print(len(N))
    return expand(str(template), N=N)


def expand_faa_from_merged_sequences(wildcards, template, busco_blacklist=None):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    N = glob_wildcards(os.path.join(checkpoint_output, "{N}.faa")).N
    if busco_blacklist is not None:
        N = set(N) - set(list(map(lambda s: f"{s}.merged", busco_blacklist)))
    return expand(str(template), N=N)

#------------------TEMPORARY CODE!!!!!!!!!!!!! -----------------------
# blacklist is applied at the concatenation stage
busco_blacklist = None
busco_blacklist_path = Path("input/BUSCO.blacklist")

if busco_blacklist_path.exists() and (busco_dir_path.stat().st_size > 0):
    busco_blacklist = pd.read_csv("input/BUSCO.blacklist", sep="\t", header=None).squeeze()

#---------------------------------------------------------------------

# +-----------------+
# |  the "All" rule |
# +-----------------+

output_files = [
    # ---- Busco ----
    expand(busco_dir_path / "{species}/short_summary_{species}.txt", species=config["species_list"]),
    # ---- Merge sequences with common ids ----
    lambda w: expand_fna_from_merged_sequences(w, merged_sequences_dir_path / "{N}.fna"),
    lambda w: expand_faa_from_merged_sequences(w, merged_sequences_dir_path / "{N}.faa"),
    species_ids_dir_path / "unique_species_ids.png",
    busco_dir_path / "busco_summaries.svg",
]
if "quastcore" in config:
    if config["quastcore"]:
        output_files.append(quastcore_dir_path / "assembly_stats.csv")
if "dna_alignment" in config:
    if config["dna_alignment"]:
        output_files.append(lambda w: expand_fna_from_merged_sequences(w, alignments_dir_path / "fna" / "{N}.fna"))
        if "dna_filtration" in config:
            if config["dna_filtration"]:
                output_files.append(lambda w: expand_fna_from_merged_sequences(w, filtered_alignments_dir_path / "fna" / "{N}.fna"))
                output_files.append(concat_alignments_dir_path / fasta_dna_filename)
                output_files.append(concat_alignments_dir_path / nexus_dna_filename)
                if "iqtree_dna" in config:
                    if config["iqtree_dna"]:
                        output_files.append(iqtree_dir_path / "fna" / f"{fasta_dna_filename}.treefile")
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(iqtree_dir_path / "fna" / f"{fasta_dna_filename}.length_and_support_tree.svg")
                if "astral" in config:
                    if config["astral"]:
                        output_files.append(astral_dir_path / astral_tree)
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(astral_dir_path / f"{astral_tree}.svg")
                if "rapidnj" in config:
                    if config["rapidnj"]:
                        output_files.append(concat_alignments_dir_path / stockholm_dna_filename)
                        output_files.append(rapidnj_dir_path / rapidnj_tree)
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(rapidnj_dir_path / f"{fasta_dna_filename}.only_tree.svg")
                if "phylip" in config:
                    if config["phylip"]:
                        output_files.append(concat_alignments_dir_path / phylip_dna_filename)
                        output_files.append(phylip_dir_path / phylip_tree)
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(phylip_dir_path / f"{fasta_dna_filename}.only_tree.svg")
if "protein_alignment" in config:
    if config["protein_alignment"]:
        output_files.append(lambda w: expand_faa_from_merged_sequences(w, alignments_dir_path / "faa" / "{N}.faa"))
        if "protein_filtration" in config:
            if config["protein_filtration"]:
                output_files.append(lambda w: expand_fna_from_merged_sequences(w, filtered_alignments_dir_path / "faa" / "{N}.faa"))
                output_files.append(concat_alignments_dir_path / fasta_protein_filename)
                # output_files.append(concat_alignments_dir_path / nexus_protein_filename)
                # output_files.append(concat_alignments_dir_path / stockholm_protein_filename)
                # output_files.append(concat_alignments_dir_path / phylip_protein_filename)
                if "iqtree_protein" in config:
                    if config["iqtree_protein"]:
                        output_files.append(iqtree_dir_path / "faa" / f"{fasta_protein_filename}.treefile")
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(iqtree_dir_path / "faa" / f"{fasta_protein_filename}.length_and_support_tree.svg")
if "mrbayes_dna" in config:
    if config["mrbayes_dna"]:  # to-do: upgrade
        output_files.append(mrbayes_dir_path / "fna")
if "mrbayes_protein" in config:
    if config["mrbayes_protein"]:
        output_files.append(mrbayes_dir_path / "faa")


localrules:
    all,


rule all:
    input:
        output_files,


# ---- Load rules ----
include: "workflow/rules/quastcore.smk"
include: "workflow/rules/busco.smk"
include: "workflow/rules/common_ids.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/filtration.smk"
include: "workflow/rules/concat_alignments.smk"
include: "workflow/rules/iqtree.smk"
include: "workflow/rules/mrbayes.smk"
include: "workflow/rules/visualization.smk"
include: "workflow/rules/astral.smk"
include: "workflow/rules/rapidnj.smk"
include: "workflow/rules/phylip.smk"
