from pathlib import Path
import gzip
import os
import re

import pandas as pd
from snakemake.utils import validate


configfile: "config/default.yaml"


# ---- Validate config (types + typo-catching against the schema) ----
validate(config, "workflow/schemas/config.schema.yaml")


# ---- Constants ----
FASTA_EXTENSIONS = [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]
FASTA_PATTERNS = [f"*{ext}" for ext in FASTA_EXTENSIONS]

# ---- Conda environment ----
main_env = (
    config["conda"]["buscoclade_main"]["name"]
    if config["use_existing_envs"]
    else "../../%s" % config["conda"]["buscoclade_main"]["yaml"]
)

# ---- Setup paths ----
# -- Input --
genome_dir_path = Path(config["genome_dir"]).resolve()
reconstruct_dir_path = Path(config["reconstruct_dir"]).resolve()
vcf2phylip_dir_path = Path(config["vcf2phylip_dir"]).resolve()

# -- Output --
output_dir_path = Path(config["output_dir"])

# -- Logs and benchmarks --
log_dir_path = output_dir_path / config["log_dir"]
cluster_log_dir_path = output_dir_path / config["cluster_log_dir"]
benchmark_dir_path = output_dir_path / config["benchmark_dir"]

# ---- Create directories ----
onstart:
    os.makedirs(output_dir_path, exist_ok=True)
    os.makedirs(log_dir_path, exist_ok=True)
    os.makedirs(cluster_log_dir_path, exist_ok=True)
    os.makedirs(benchmark_dir_path, exist_ok=True)

# -- Results --
quastcore_dir_path = output_dir_path / config["quastcore_dir"]
busco_dir_path = output_dir_path / config["busco_dir"]
main_ids_dir_path = output_dir_path / config["main_ids_dir"]
species_ids_dir_path = output_dir_path / config["species_ids_dir"]
common_ids_dir_path = output_dir_path / config["common_ids_dir"]
merged_sequences_dir_path = output_dir_path / config["merged_sequences_dir"]
alignments_dir_path = output_dir_path / config["alignments_dir"]
pre_altref_alignments_dir_path = output_dir_path / config["pre_altref_alignments_dir"]
reconstruct_consensus_dir_path = output_dir_path / config["reconstruct_consensus_dir"]
filtered_alignments_dir_path = output_dir_path / config["filtered_alignments_dir"]
concat_alignments_dir_path = output_dir_path / config["concat_alignments_dir"]
iqtree_dir_path = output_dir_path / config["iqtree_dir"]
mrbayes_dir_path = output_dir_path / config["mrbayes_dir"]
astral_dir_path = output_dir_path / config["astral_dir"]
rapidnj_dir_path = output_dir_path / config["rapidnj_dir"]
phylip_dir_path = output_dir_path / config["phylip_dir"]
raxml_dir_path = output_dir_path / config["raxml_dir"]


# ---- Setup filenames ----
pfx = config["alignment_file_prefix"]
sup = config["nodes_filtration_by_support"]

fasta_filename        = f"{pfx}.fna"
nexus_filename        = f"{pfx}.fna.nex"
phylip_filename       = f"{pfx}.fna.phy"
stockholm_filename    = f"{pfx}.fna.sth"
astral_input_trees    = f"{pfx}.iqtree_per_fna.concat.treefile"
astral_filtered_trees = f"{pfx}.iqtree_per_fna.concat.{sup}.treefile"
astral_tree           = f"{pfx}.{sup}.fna.astral.treefile"
rapidnj_tree          = f"{pfx}.fna.rapidnj.treefile"
rapidnj_matrix        = f"{pfx}.fna.rapidnj.matrix"
phylip_tree           = f"{pfx}.fna.phy.namefix.treefile"
raxml_tree            = f"{pfx}.fna.raxml.treefile"


# ---- Necessary functions ----
def get_fasta_stem(path: Path) -> str:
    """Returns filename stem, correctly handling multi-part extensions like .fna.gz."""
    name = path.name
    for ext in FASTA_EXTENSIONS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return path.stem


def get_fasta_files(directory: Path) -> list[Path]:
    """Returns all FASTA files in a directory, matching all known extensions."""
    seen = set()
    files = []
    for pattern in FASTA_PATTERNS:
        for f in sorted(directory.rglob(pattern)):
            if f.is_file() and f not in seen:
                seen.add(f)
                files.append(f)
    return files


def get_toplevel_fasta_files(directory: Path) -> list[Path]:
    """Returns FASTA files directly inside a directory (non-recursive), matching
    all known extensions. Used to find the reference FASTA at a reconstruct
    subdir's top level without descending into vcf/ or fasta/."""
    seen = set()
    files = []
    for pattern in FASTA_PATTERNS:
        for f in sorted(directory.glob(pattern)):
            if f.is_file() and f not in seen:
                seen.add(f)
                files.append(f)
    return files


def get_reconstruct_map(root: Path) -> dict:
    """
    Scans reconstruct/ subdirectories and returns a mapping of reconstructed
    species names to their reconstruction source, reference FASTA, and reference
    prefix. Each subdir holds one reference FASTA at its top level plus optional
    vcf/ and fasta/ subdirectories:

        <subdir>/reference.fasta          -> reference (BUSCO run once on it)
        <subdir>/vcf/<sample>.vcf.gz      -> "{sample}.{ref_prefix}.AltRef"    (source: vcf)
        <subdir>/fasta/<sample>.fasta     -> "{sample}.{ref_prefix}.Consensus" (source: consensus)
    """
    mapping = {}

    if not root.exists():
        return mapping

    for subdir in sorted(root.iterdir()):
        if not subdir.is_dir():
            continue

        ref_files = get_toplevel_fasta_files(subdir)
        if not ref_files:
            continue
        ref_file = ref_files[0]
        ref_prefix = get_fasta_stem(ref_file)

        vcf_subdir = subdir / "vcf"
        if vcf_subdir.is_dir():
            for vcf_file in sorted(vcf_subdir.glob("*.vcf.gz")):
                sample = vcf_file.stem.split(".")[0]
                mapping[f"{sample}.{ref_prefix}.AltRef"] = {
                    "source": "vcf",
                    "sample": sample,
                    "reference": ref_file,
                    "ref_prefix": ref_prefix,
                    "vcf": vcf_file,
                }

        fasta_subdir = subdir / "fasta"
        if fasta_subdir.is_dir():
            for consensus_file in get_fasta_files(fasta_subdir):
                sample = get_fasta_stem(consensus_file)
                mapping[f"{sample}.{ref_prefix}.Consensus"] = {
                    "source": "consensus",
                    "sample": sample,
                    "reference": ref_file,
                    "ref_prefix": ref_prefix,
                    "consensus_fasta": consensus_file,
                }

    return mapping


def get_species_list(reconstructed: list, genome_species: list, refs: list) -> list:
    """
    Merges genome and reconstructed species into a single sorted list.
    Reference genomes are included only if reconstruct_refs_as_species is True.
    """
    all_species = set(reconstructed + genome_species)
    if config.get("reconstruct_refs_as_species"):
        all_species.update(refs)
    return sorted(all_species)


def extract_samples_from_vcf(vcf_file: Path) -> list:
    """Extract sample names from VCF file header."""
    try:
        with gzip.open(vcf_file, "rt") as file:
            for line in file:
                if line.startswith("#CHROM"):
                    parts = line.strip().split("\t")
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


def get_genome_file(wildcards) -> Path:
    """
    Returns the genome FASTA path for a given species wildcard.
    - Reconstructed species (VCF- or consensus-derived) are not allowed here —
      BUSCO is not run on them directly.
    - Reference genomes from reconstruct/ are expected as symlinks in genomes/.
    - Regular genome assemblies are looked up directly in genomes/.
    """
    if wildcards.species in reconstruct_map:
        raise ValueError(
            f"get_genome_file called for reconstructed species '{wildcards.species}' — "
            "BUSCO should not be run directly on reconstructed species."
        )

    if wildcards.species in reconstructed_refs:
        return genome_dir_path / ref_link_name[wildcards.species]

    for f in get_fasta_files(genome_dir_path):
        if get_fasta_stem(f) == wildcards.species:
            return f

    raise ValueError(f"No genome file found for species '{wildcards.species}' in {genome_dir_path}")


def get_all_genome_files() -> list[Path]:
    """Returns genome files for all non-reconstructed species in species_list."""
    return [
        get_genome_file(type("W", (), {"species": s})())
        for s in config["species_list"]
        if s not in reconstruct_map
    ]


def get_vcf2phylip_outputs() -> list:
    """Returns output files for the vcf2phylip route."""
    files = [concat_alignments_dir_path / fasta_filename]
    if config.get("iqtree"):
        files.append(iqtree_dir_path / f"{fasta_filename}.treefile")
        if config.get("draw_phylotrees"):
            files.append(iqtree_dir_path / f"{fasta_filename}.length_and_support_tree.svg")
    if config.get("rapidnj"):
        files.append(concat_alignments_dir_path / stockholm_filename)
        files.append(rapidnj_dir_path / rapidnj_tree)
        if config.get("draw_phylotrees"):
            files.append(rapidnj_dir_path / f"{fasta_filename}.only_tree.svg")
    if config.get("phylip"):
        files.append(concat_alignments_dir_path / phylip_filename)
        files.append(phylip_dir_path / phylip_tree)
        if config.get("draw_phylotrees"):
            files.append(phylip_dir_path / f"{fasta_filename}.only_tree.svg")
    if config.get("raxml"):
        files.append(raxml_dir_path / raxml_tree)
        if config.get("draw_phylotrees"):
            files.append(raxml_dir_path / f"{fasta_filename}.only_tree.svg")
    return files


def get_busco_outputs() -> list:
    """Returns output files for the main BUSCO-based pipeline route."""
    files = []

    # Ensure BUSCO runs on reference genomes used for reconstruction
    if reconstructed_species:
        files.append(expand(busco_dir_path / "{species}/short_summary_{species}.txt", species=reconstructed_refs))

    files += [
        expand(busco_dir_path / "{species}/short_summary_{species}.txt", species=config["species_list"]),
        lambda w: expand_fna_from_merged_sequences(w, merged_sequences_dir_path / "{N}.fna", busco_blacklist=busco_blacklist),
        species_ids_dir_path / "unique_species_ids.svg",
        main_ids_dir_path / "gene_counts.tsv",
        busco_dir_path / "busco_summaries.svg",
    ]

    if config.get("quastcore") and get_all_genome_files():
        files.append(quastcore_dir_path / "assembly_stats.csv")

    if config.get("alignment"):
        files.append(lambda w: expand_fna_from_merged_sequences(
            w, alignments_dir_path / "fna" / "{N}.fna", busco_blacklist=busco_blacklist
        ))
        if config.get("filtration"):
            files += [
                lambda w: expand_fna_from_merged_sequences(
                    w, filtered_alignments_dir_path / "fna" / "{N}.fna", busco_blacklist=busco_blacklist
                ),
                concat_alignments_dir_path / fasta_filename,
            ]
            if config.get("iqtree"):
                files.append(iqtree_dir_path / f"{fasta_filename}.treefile")
                if config.get("draw_phylotrees"):
                    files.append(iqtree_dir_path / f"{fasta_filename}.length_and_support_tree.svg")
            if config.get("astral"):
                files.append(astral_dir_path / astral_tree)
                if config.get("draw_phylotrees"):
                    files.append(astral_dir_path / f"{astral_tree}.svg")
            if config.get("rapidnj"):
                files.append(concat_alignments_dir_path / stockholm_filename)
                files.append(rapidnj_dir_path / rapidnj_tree)
                if config.get("draw_phylotrees"):
                    files.append(rapidnj_dir_path / f"{fasta_filename}.only_tree.svg")
            if config.get("phylip"):
                files.append(concat_alignments_dir_path / phylip_filename)
                files.append(phylip_dir_path / phylip_tree)
                if config.get("draw_phylotrees"):
                    files.append(phylip_dir_path / f"{fasta_filename}.only_tree.svg")
            if config.get("raxml"):
                files.append(raxml_dir_path / raxml_tree)
                if config.get("draw_phylotrees"):
                    files.append(raxml_dir_path / f"{fasta_filename}.only_tree.svg")
            if config.get("mrbayes"):
                files.append(concat_alignments_dir_path / nexus_filename)
                files.append(mrbayes_dir_path / f"{fasta_filename}.nex.con.tre.nwk")

    return files


# ---- Input data ----
# reconstruct_map holds both reconstruction sources under one reference:
#   source == "vcf"       -> reconstructed via apply_vcf_to_busco (SNPs from a VCF)
#   source == "consensus" -> reconstructed via apply_consensus_to_busco (already-built FASTA)
reconstruct_map      = get_reconstruct_map(reconstruct_dir_path)
reconstructed_species = list(reconstruct_map.keys())
reconstructed_refs   = sorted({v["ref_prefix"] for v in reconstruct_map.values()})
vcf_species          = [s for s, v in reconstruct_map.items() if v["source"] == "vcf"]
consensus_species    = [s for s, v in reconstruct_map.items() if v["source"] == "consensus"]

# ---- Gap-aware reconstructed-sequence insertion ----
reconstruct_gapaware = config.get("altref_gapaware_insertion", False)

# ref_prefix → [reconstructed_sp1, ...]  (used in alignment.smk; both sources)
ref_to_reconstructed = {}
for _sp, _info in reconstruct_map.items():
    ref_to_reconstructed.setdefault(_info["ref_prefix"], []).append(_sp)

# The reference is symlinked into genomes/ under its ORIGINAL basename (keeping the
# extension, including .gz) so BUSCO treats it exactly like a native genome assembly.
# Naming the link "{ref_prefix}.fasta" for a gzipped reference would hand BUSCO gzip
# bytes under a .fasta name and crash it.
ref_link_source = {}
for _info in reconstruct_map.values():
    ref_link_source.setdefault(_info["ref_prefix"], _info["reference"])
ref_link_name = {rp: src.name for rp, src in ref_link_source.items()}
ref_link_source_by_basename = {src.name: src for src in ref_link_source.values()}
ref_link_basenames = sorted(ref_link_source_by_basename)

genome_species = sorted(
    {get_fasta_stem(f) for f in get_fasta_files(genome_dir_path)} - set(reconstructed_refs)
)

# ---- Species list ----
if "species_list" not in config:
    if not config.get("vcf2phylip"):
        config["species_list"] = get_species_list(reconstructed_species, genome_species, reconstructed_refs)
    else:
        vcf_file = list(vcf2phylip_dir_path.rglob("*.vcf.gz"))
        if len(vcf_file) != 1:
            raise ValueError(
                f"vcf2phylip requires exactly one VCF in {vcf2phylip_dir_path}, "
                f"found {len(vcf_file)}: {vcf_file}"
            )
        vcf_file = vcf_file[0]
        prefix = vcf_file.name.replace(".vcf.gz", "")
        config["species_list"] = extract_samples_from_vcf(vcf_file)

    print("Species list:", config["species_list"])


# Species without reconstructed samples for raw alignment (refs stay — gaps are
# projected onto them). In gap-aware mode reconstructed samples are inserted after
# alignment, so only refs + non-reconstructed species go into the aligner input.
species_list_for_raw_alignment = (
    sorted(
        {s for s in config["species_list"] if s not in reconstruct_map}
        | set(reconstructed_refs)  # refs are always needed as gap templates
    )
    if reconstruct_gapaware and reconstruct_map
    else config["species_list"]
)

# ---- BUSCO blacklist ----
busco_blacklist = None
if "busco_blacklist" in config:
    blacklist_path = Path(config["busco_blacklist"])
    if blacklist_path.exists():
        busco_blacklist = load_busco_blacklist(blacklist_path)


# ---- "All" rule ----
localrules:
    all,


rule all:
    input:
        get_vcf2phylip_outputs() if config.get("vcf2phylip") else get_busco_outputs(),


# ---- Load rules ----
include: "workflow/rules/reconstruct.smk"
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
include: "workflow/rules/vcf2phylip.smk"
