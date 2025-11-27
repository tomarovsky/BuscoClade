from pathlib import Path
import gzip
import os

import pandas as pd


# ---- Setup paths ----
# -- Input --
genome_dir_path = Path(config["genome_dir"]).resolve()
vcf_reconstruct_dir_path = Path(config["vcf_reconstruct_dir"]).resolve()

# -- Logs and benchmarks --
cluster_log_dir_path = Path(config["cluster_log_dir"])
log_dir_path = Path(config["log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])
output_dir_path = Path(config["output_dir"])

# -- Results --
quastcore_dir_path = output_dir_path / config["quastcore_dir"]
altref_dir_path = output_dir_path / config["altref_dir"]
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

# ---- Setup filenames ----
fasta_filename = "{}.fna".format(config["alignment_file_prefix"])
nexus_filename = "{}.fna.nex".format(config["alignment_file_prefix"])
phylip_filename = "{}.fna.phy".format(config["alignment_file_prefix"])
stockholm_filename = "{}.fna.sth".format(config["alignment_file_prefix"])
astral_input_trees = "{}.iqtree_per_fna.concat.treefile".format(config["alignment_file_prefix"])
astral_filtered_trees = "{0}.iqtree_per_fna.concat.{1}.treefile".format(config["alignment_file_prefix"], config["nodes_filtrataion_by_support"])
astral_tree = "{0}.{1}.fna.astral.treefile".format(config["alignment_file_prefix"], config["nodes_filtrataion_by_support"])
rapidnj_tree = "{}.fna.rapidnj.treefile".format(config["alignment_file_prefix"])
rapidnj_matrix = "{}.fna.rapidnj.matrix".format(config["alignment_file_prefix"])
phylip_tree = "{}.fna.phy.namefix.treefile".format(config["alignment_file_prefix"])
raxml_tree = "{}.fna.raxml.treefile".format(config["alignment_file_prefix"])


# ---- Necessary functions ----
def get_vcf_reconstruct_map(vcf_dir: Path) -> dict:
    """Returns a dictionary where keys are species names and values are dictionaries of VCF and FASTA files."""
    vcf_mapping = {}

    if vcf_dir.exists():
        for vcf_subdir in vcf_dir.iterdir():
            if not vcf_subdir.is_dir():
                continue

            # Find first FASTA file
            ref_file = next(vcf_subdir.glob("*.fasta"), None)
            if ref_file is None:
                continue
            ref_prefix = ref_file.stem

            # Process VCF files
            for vcf_file in vcf_subdir.glob("*.vcf.gz"):
                vcf_id = vcf_file.stem.split('.')[0]
                alt_name = f"{vcf_id}.{ref_prefix}.AltRef"
                vcf_mapping[alt_name] = {
                    'vcf': vcf_file,
                    'reference': ref_file
                }

    return vcf_mapping


def get_species_list(vcf_species: list, genome_species: list) -> list:
    """Merges and returns final species list."""
    return sorted(set(vcf_species + genome_species))


def extract_samples_from_vcf(vcf_file: Path) -> list:
    """Extract sample names from VCF file header."""
    try:
        with gzip.open(vcf_file, 'rt') as file:
            for line in file:
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    return parts[9:]  # Sample names start at column 10
    except Exception as e:
        raise ValueError(f"Failed to read VCF file {vcf_file}: {e}")

    raise ValueError(f"No header line found in VCF file {vcf_file}")


def load_busco_blacklist(blacklist_path: Path):
    """Load BUSCO blacklist from file."""
    try:
        df = pd.read_csv(blacklist_path, sep="\t", header=None, comment="#")
        if not df.empty:
            return df.squeeze().tolist()
    except (pd.errors.EmptyDataError, Exception) as e:
        print(f"Warning: Could not load BUSCO blacklist from {blacklist_path}: {e}")
    return None


def expand_fna_from_merged_sequences(wildcards, template, busco_blacklist=None):
    """Expand fna files from merged sequences."""
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    N = glob_wildcards(os.path.join(checkpoint_output, "{N}.fna")).N
    print(f"Initial BUSCO count: {len(N)}")
    if busco_blacklist:
        blacklist_set = {f"{s}.merged" for s in busco_blacklist}
        N = [n for n in N if n not in blacklist_set]
    print(f"Final BUSCO count after filtering: {len(N)}")
    return expand(str(template), N=N)


# ---- Input data ----
genome_species = [f.stem for f in genome_dir_path.glob("*.fasta") if f.is_file()]
vcf_reconstruct_map = get_vcf_reconstruct_map(vcf_reconstruct_dir_path)
vcf_reconstruct_species = list(vcf_reconstruct_map.keys())

# Species list configuration
if "species_list" not in config:
    if not config.get("vcf2phylip"):
        config["species_list"] = get_species_list(genome_species, vcf_reconstruct_species)
    else:
        # VCF2Phylip specific processing
        vcf_files = list(vcf_reconstruct_dir_path.rglob("*.vcf.gz"))
        if len(vcf_files) != 1:
            raise ValueError(
                f"vcf2phylip requires exactly one VCF in {vcf_reconstruct_dir_path}, "
                f"found {len(vcf_files)}: {vcf_files}"
            )
        vcf_file = vcf_files[0]
        prefix = vcf_file.name.replace(".vcf.gz", "")
        config["species_list"] = extract_samples_from_vcf(vcf_file)

    print("Species list:", config["species_list"])


# Busco blacklist processing
busco_blacklist = None
if "busco_blacklist" in config:
    blacklist_path = Path(config["busco_blacklist"])
    if blacklist_path.exists():
        busco_blacklist = load_busco_blacklist(blacklist_path)


# ---- "All" rule ----

output_files = []

if config.get("vcf2phylip"):
    output_files.append(concat_alignments_dir_path / fasta_filename)
    if config.get("iqtree"):
        output_files.append(iqtree_dir_path / "fna" / f"{fasta_filename}.treefile")
        if config.get("draw_phylotrees"):
            output_files.append(iqtree_dir_path / "fna" / f"{fasta_filename}.length_and_support_tree.svg")
    if config.get("rapidnj"):
        output_files.append(concat_alignments_dir_path / stockholm_filename)
        output_files.append(rapidnj_dir_path / rapidnj_tree)
        if config.get("draw_phylotrees"):
            output_files.append(rapidnj_dir_path / f"{fasta_filename}.only_tree.svg")
    if config.get("phylip"):
        output_files.append(concat_alignments_dir_path / phylip_filename)
        output_files.append(phylip_dir_path / phylip_tree)
        if config.get("draw_phylotrees"):
            output_files.append(phylip_dir_path / f"{fasta_filename}.only_tree.svg")
    if config.get("raxml"):
        output_files.append(raxml_dir_path / raxml_tree)
        if config.get("draw_phylotrees"):
            output_files.append(raxml_dir_path / f"{fasta_filename}.only_tree.svg")
else:
    expand(busco_dir_path / "{species}/short_summary_{species}.txt", species=config["species_list"]),
    lambda w: expand_fna_from_merged_sequences(w, merged_sequences_dir_path / "{N}.fna", busco_blacklist=busco_blacklist),
    species_ids_dir_path / "unique_species_ids.svg",
    busco_dir_path / "busco_summaries.svg",

    if config.get("quastcore"):
        output_files.append(quastcore_dir_path / "assembly_stats.csv")

    if config.get("alignment"):
        output_files.append(lambda w: expand_fna_from_merged_sequences(w, alignments_dir_path / "fna" / "{N}.fna", busco_blacklist=busco_blacklist))
        if config.get("filtration"):
            output_files.append(lambda w: expand_fna_from_merged_sequences(w, filtered_alignments_dir_path / "fna" / "{N}.fna", busco_blacklist=busco_blacklist))
            output_files.append(concat_alignments_dir_path / fasta_filename)
            if config.get("iqtree"):
                output_files.append(iqtree_dir_path / "fna" / f"{fasta_filename}.treefile")
                if config.get("draw_phylotrees"):
                    output_files.append(iqtree_dir_path / "fna" / f"{fasta_filename}.length_and_support_tree.svg")
            if config.get("astral"):
                output_files.append(astral_dir_path / astral_tree)
                if config.get("draw_phylotrees"):
                    output_files.append(astral_dir_path / f"{astral_tree}.svg")
            if config.get("rapidnj"):
                output_files.append(concat_alignments_dir_path / stockholm_filename)
                output_files.append(rapidnj_dir_path / rapidnj_tree)
                if config.get("draw_phylotrees"):
                    output_files.append(rapidnj_dir_path / f"{fasta_filename}.only_tree.svg")
            if config.get("phylip"):
                output_files.append(concat_alignments_dir_path / phylip_filename)
                output_files.append(phylip_dir_path / phylip_tree)
                if config.get("draw_phylotrees"):
                    output_files.append(phylip_dir_path / f"{fasta_filename}.only_tree.svg")
            if config.get("raxml"):
                output_files.append(raxml_dir_path / raxml_tree)
                if config.get("draw_phylotrees"):
                    output_files.append(raxml_dir_path / f"{fasta_filename}.only_tree.svg")
            if config.get("mrbayes"):
                    output_files.append(concat_alignments_dir_path / nexus_filename)
                    output_files.append(mrbayes_dir_path / "fna")


localrules:
    all,


rule all:
    input:
        output_files,


# ---- Load rules ----
include: "workflow/rules/vcf_reconstruct.smk"
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
include: "workflow/rules/raxml.smk"
