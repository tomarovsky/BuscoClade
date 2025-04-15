from pathlib import Path
import os

# ---- Setup config ----
configfile: "config/default.yaml"


# ---- Setup paths ----
cluster_log_dir_path = Path(config["cluster_log_dir"])
genome_dir_path = Path(config["genome_dir"]).resolve()
log_dir_path = Path(config["log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])
output_dir_path = Path(config["output_dir"])
vcf_dir_path = Path(config["vcf_dir"]).resolve()


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
raxml_dir_path = output_dir_path / config["raxml_dir"]

# if "species_list" not in config:
#    config["species_list"] = [f.stem for f in genome_dir_path.iterdir() if f.is_file() and f.suffix == ".fasta"]


# ---- Setup filenames ----
vcf_filename = "{}.{ext,vcf|vcf.gz}"
vcf_mask_filename = "{}.{ext,bed|bed.gz}"
vcf_reference_filename = "{}.{ext,fna|fasta}"
fasta_dna_filename = f"concat_alignment{config['altref_suffix']}.fna"
fasta_protein_filename = "{}.faa".format(config["alignment_file_prefix"])
nexus_dna_filename = f"{fasta_dna_filename}.nex"
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
raxml_tree = f"{fasta_dna_filename}.raxml.treefile"



# ---- Necessary functions ----
def expand_fna_from_merged_sequences(wildcards, template):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    N = glob_wildcards(os.path.join(checkpoint_output, "{N}.fna")).N
    return expand(str(template), N=N)


def expand_faa_from_merged_sequences(wildcards, template):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    N = glob_wildcards(os.path.join(checkpoint_output, "{N}.faa")).N
    return expand(str(template), N=N)

def get_samples_from_file(sample_list_path):
    with open(sample_list_path) as f:
        return [line.strip() for line in f if line.strip()]

        

#---- VCF ----
if config["vcf"]:
    include: "workflow/rules/vcf.smk"
    
    checkpoint vcf_processed:
        input:
            expand(
                genome_dir_path / "{sample}.AltRef.fasta",
                sample=glob_wildcards(os.path.join(vcf_dir_path, "*", "{sample}_PASS.vcf.gz")).sample
            )

def get_vcf_samples(vcf_dir):
    vcf_files = list(Path(vcf_dir).rglob("*.vcf.gz"))
    if not vcf_files:
        raise ValueError(f"No VCF samples found in {vcf_dir} subdirectories.")
    return sorted({f.stem.replace("_PASS", "") for f in vcf_files})

if "species_list" not in config:
    if config.get("vcf", False):
        config["species_list"] = []
    else:
        config["species_list"] = [f.stem for f in genome_dir_path.glob("*.fasta") if f.is_file()]

# +-----------------+
# |  the "All" rule |
# +-----------------+

output_files = [
    # ---- VCF ----
    checkpoints.vcf_processed.get().output[0] if config["vcf"] else [],
    *(
        expand(
            genome_dir_path / "{sample}.AltRef.fasta",
            sample=config["species_list"]
        )
        if config.get("vcf", False)
        else []
    ),
    expand(busco_dir_path / "{species}/short_summary_{species}.txt", species=config["species_list"]),
    # ---- Busco ----
    expand(busco_dir_path /
           "{species}/short_summary_{species}.txt", species=config["species_list"]),
    # ---- Merge sequences with common ids ----
    lambda w: expand_fna_from_merged_sequences(
        w, merged_sequences_dir_path / "{N}.fna"),
    lambda w: expand_faa_from_merged_sequences(
        w, merged_sequences_dir_path / "{N}.faa"),
    species_ids_dir_path / "unique_species_ids.png",
    busco_dir_path / "busco_summaries.svg",


    
]
# if "vcf" in config:
#    if config["vcf"]:
#        expand("genomes/{sample}.AltRef.fasta", sample=get_sample_names())
if "quastcore" in config:
    if config["quastcore"]:
        output_files.append(quastcore_dir_path / "assembly_stats.csv")
if "dna_alignment" in config:
    if config["dna_alignment"]:
        output_files.append(lambda w: expand_fna_from_merged_sequences(
            w, alignments_dir_path / "fna" / "{N}.fna"))
        if "dna_filtration" in config:
            if config["dna_filtration"]:
                output_files.append(lambda w: expand_fna_from_merged_sequences(
                    w, filtered_alignments_dir_path / "fna" / "{N}.fna"))
                output_files.append(
                    concat_alignments_dir_path / fasta_dna_filename)
                output_files.append(
                    concat_alignments_dir_path / nexus_dna_filename)
                if "iqtree_dna" in config:
                    if config["iqtree_dna"]:
                        output_files.append(
                            iqtree_dir_path / "fna" / f"{fasta_dna_filename}.treefile")
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(
                                    iqtree_dir_path / "fna" / f"{fasta_dna_filename}.length_and_support_tree.png")
                if "astral" in config:
                    if config["astral"]:
                        output_files.append(astral_dir_path / astral_tree)
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(
                                    astral_dir_path / f"{astral_tree}.png")
                if "rapidnj" in config:
                    if config["rapidnj"]:
                        output_files.append(
                            concat_alignments_dir_path / stockholm_dna_filename)
                        output_files.append(rapidnj_dir_path / rapidnj_tree)
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(
                                    rapidnj_dir_path / f"{fasta_dna_filename}.only_tree.png")
                if "phylip" in config:
                    if config["phylip"]:
                        output_files.append(
                            concat_alignments_dir_path / phylip_dna_filename)
                        output_files.append(phylip_dir_path / phylip_tree)
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(
                                    phylip_dir_path / f"{fasta_dna_filename}.only_tree.png")
                if "raxml" in config:
                    if config["raxml"]:
                        output_files.append(raxml_dir_path / raxml_tree)
#                        if "draw_phylotrees" in config:
#                            if config["draw_phylotrees"]:
#                                output_files.append(
#                                    raxml_dir_path / f"{fasta_dna_filename}.only_tree.png")
if "protein_alignment" in config:
    if config["protein_alignment"]:
        output_files.append(lambda w: expand_faa_from_merged_sequences(
            w, alignments_dir_path / "faa" / "{N}.faa"))
        if "protein_filtration" in config:
            if config["protein_filtration"]:
                output_files.append(lambda w: expand_fna_from_merged_sequences(
                    w, filtered_alignments_dir_path / "faa" / "{N}.faa"))
                output_files.append(
                    concat_alignments_dir_path / fasta_protein_filename)
                # output_files.append(concat_alignments_dir_path / nexus_protein_filename)
                # output_files.append(concat_alignments_dir_path / stockholm_protein_filename)
                # output_files.append(concat_alignments_dir_path / phylip_protein_filename)
                if "iqtree_protein" in config:
                    if config["iqtree_protein"]:
                        output_files.append(
                            iqtree_dir_path / "faa" / f"{fasta_protein_filename}.treefile")
                        if "draw_phylotrees" in config:
                            if config["draw_phylotrees"]:
                                output_files.append(
                                    iqtree_dir_path / "faa" / f"{fasta_protein_filename}.length_and_support_tree.png")
if "mrbayes_dna" in config:
    if config["mrbayes_dna"]:
        output_files.append(mrbayes_dir_path / "fna")
if "mrbayes_protein" in config:
    if config["mrbayes_protein"]:
        output_files.append(mrbayes_dir_path / "faa")


localrules: all

rule all:
    input:
        output_files

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
include: "workflow/rules/raxml.smk"
include: "workflow/rules/rapidnj.smk"
include: "workflow/rules/phylip.smk"
