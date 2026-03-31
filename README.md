<p align="center">
  <h1 align="center">BuscoClade</h1>
</p>

<p align="center">
  <i>Snakemake-based workflow to construct species phylogenies using BUSCOs</i>
</p>

<p align="center" style="text-decoration: none; border: none;">
	<a href="https://snakemake.github.io" style="text-decoration: none">
		<img alt="Snakemake" src="https://img.shields.io/badge/Snakemake-9.16-dfb034?style=flat&logo=snakemake&logoColor=ffffff&labelColor=00623d"></a>
	<a href="https://opensource.org/licenses/MIT" style="text-decoration: none">
		<img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-dfb034?style=flat&logo=opensourceinitiative&logoColor=ffffff&labelColor=00623d"></a>
	<a href="https://github.com/tomarovsky/BuscoClade/stargazers" style="text-decoration: none">
		<img alt="Stars" src="https://img.shields.io/github/stars/tomarovsky/BuscoClade?style=flat&logo=starship&color=dfb034&logoColor=ffffff&labelColor=00623d"></a>
</p>

![BuscoClade logo](logo.png)

---

## Workflow

```mermaid
flowchart LR

%% ----- INPUT -----
subgraph INPUT["Input data"]
A_fa["Genome assemblies (FASTA)"]
A_vcf2["Multi-sample VCF"]
subgraph VCF_REF["VCF reconstruction"]
A_ref["Reference genome (FASTA)"]
A_vcf["Per-sample VCFs"]
end
end

%% ----- BUSCO -----
subgraph BUSCO["Ortholog extraction"]
B_busco["BUSCO"]
end

%% ----- PREPROCESSING -----
subgraph PREP["Sequence processing"]
subgraph ALN["Multiple alignment"]
C_aln["MAFFT\nMUSCLE\nPRANK"]
end
subgraph FLT["Trimming"]
C_flt["ClipKIT\nGBlocks\nTrimAl"]
end
end

%% ----- PHYLOGENY -----
subgraph PHYLO["Phylogenetic tree inference"]
subgraph CONCAT["Supermatrix approach"]
E_phy["IQTree\nMrBayes\nPHYLIP\nRAxML-NG\nRapidNJ"]
end
subgraph TREE["Multispecies coalescent"]
D_ast["Astral-IV"]
end
end

%% ----- MERGE NODE (invisible) -----
MERGE(( ))

%% ----- EDGES: MAIN -----
A_fa --> B_busco
A_ref -.-> B_busco
B_busco --> C_aln
B_busco -.-> MERGE
A_vcf --> MERGE
MERGE --> |"apply_vcf_to_busco.py"| C_aln
C_aln --> C_flt
C_flt -->|"Concat alignment"| E_phy
C_flt -->|"IQTree per gene"| D_ast

%% ----- EDGES: VCF2PHYLIP -----
A_vcf2 -. "vcf2phylip.py" .-> E_phy

%% ----- STYLE -----
classDef input fill:#e8f4ff,stroke:#2b7cd3,stroke-width:1px
classDef process fill:#eaf7ea,stroke:#2f9e44,stroke-width:1px
classDef phylo fill:#fff4e6,stroke:#e67700,stroke-width:1px
classDef optional fill:#e8f4ff,stroke:#2b7cd3,stroke-width:1px,stroke-dasharray:4 4
classDef merge fill:none,stroke:none,width:0px

class A_fa,A_ref,A_vcf input
class B_busco,C_aln,C_flt process
class D_ast,E_phy phylo
class A_vcf2 optional
class MERGE merge
```

- **Ortholog extraction:** [BUSCO](https://busco.ezlab.org/)
- **VCF-based SNP application:** [apply_vcf_to_busco.py](#vcf-based-snp-application-apply_vcf_to_buscopy), [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)
- **Alignment:** [MAFFT](https://mafft.cbrc.jp/alignment/software/), [MUSCLE](https://doi.org/10.1038/s41467-022-34630-w), [PRANK](http://wasabiapp.org/software/prank/)
- **Trimming:** [ClipKIT](https://github.com/JLSteenwyk/ClipKIT), [TrimAl](http://trimal.cgenomics.org/), [GBlocks](https://academic.oup.com/mbe/article/17/4/540/1127654)
- **Phylogenetic tree construction:** [IQTree](http://www.iqtree.org/), [MrBayes](https://nbisweden.github.io/MrBayes/), [ASTRAL-IV](https://doi.org/10.1093/molbev/msaf172), [RapidNJ](https://birc.au.dk/software/rapidnj), [PHYLIP](https://phylipweb.github.io/phylip/), [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- **Visualization:** [Etetoolkit](http://etetoolkit.org/), [Matplotlib](https://matplotlib.org/stable/)

---

## Usage

### Step 1. Deploy workflow

Clone the repository or download the [latest release](https://github.com/tomarovsky/BuscoClade/releases):

```bash
git clone https://github.com/tomarovsky/BuscoClade.git
```

### Step 2. Prepare input data

#### FASTA assemblies

Place genome assemblies into `input/genomes/`. The file prefix is used as the sample name in the output phylogeny. Supported extensions: `.fasta`, `.fna`, `.fa`, and their gzipped versions (`.fasta.gz`, `.fna.gz`, `.fa.gz`).

#### Per-sample VCFs + reference genome

If you have per-sample VCFs, the pipeline applies SNPs directly to the BUSCO sequences of a reference sample using `apply_vcf_to_busco.py`. This avoids rebuilding full pseudo-genome assemblies and re-running BUSCO for each sample — instead, BUSCO is run once on the reference, and ortholog sequences for all other samples are reconstructed by applying their SNPs to the reference BUSCO sequences exon-by-exon.

Place per-sample VCF files and the corresponding reference genome together into a subdirectory under `input/vcf_reconstruct/`. Each subdirectory is processed independently, which allows reconstructing sequences against different references in a single run:

```
input/
    genomes/
        Species1.fasta
        Species2.fasta.gz
    vcf_reconstruct/
        project_hg38/               # one reference per directory
            reference.fasta
            SampleA.vcf.gz
            SampleB.vcf.gz
        project_mm39/               # another reference
            reference.fasta
            SampleC.vcf.gz
```

The directory name is used only for organization — the VCF file prefix determines the sample name in the output phylogeny. No additional config changes are needed; the pipeline detects subdirectories automatically.

BUSCO is run once on the reference genome in each subdirectory. The resulting `single_copy_busco_sequences/` and `metaeuk_output/` are then used by `apply_vcf_to_busco.py` to produce per-sample ortholog sequences.

#### Multi-sample VCF via vcf2phylip (optional)

As an alternative to the BUSCO-based pipeline, a concatenated phylip alignment can be built directly from a multi-sample VCF using `vcf2phylip.py`, bypassing BUSCO and sequence alignment entirely. When `vcf2phylip: True` is set, only this route is executed.

Place exactly one multi-sample `.vcf.gz` file into `input/vcf2phylip/` and enable the option in the config:

```
input/
    vcf2phylip/
        all_samples.vcf.gz    # exactly one multi-sample VCF
```

```yaml
vcf2phylip: True
```

### Step 3. Configure workflow

Modify `config/default.yaml` (recommended: copy it and pass with `--configfile`). The config has four sections:

#### Pipeline configuration

Enable or disable tools and modes:

```yaml
vcf2phylip: False       # set True to use vcf2phylip instead of BUSCO-based pipeline
quastcore: True         # assembly statistics

alignment: "mafft"      # 'mafft', 'muscle' or 'prank'
filtration: "clipkit"   # 'clipkit', 'trimal' or 'gblocks'

iqtree: True
astral: True
rapidnj: True
phylip: True
raxml: True
mrbayes: False          # recommended to run GPU-compiled version separately

draw_phylotrees: True
```

#### Tool parameters

Key parameters to configure before running:

**BUSCO:**
- `busco_dataset_path`: Path to a pre-downloaded OrthoDB dataset (e.g. `"$TOOLS/busco_datasets/mammalia_odb12/"`).
- `busco_options`: Use `"--offline"` to run without internet access.
- `busco_mode`: Typically `"genome"`.
- `busco_blacklist`: Path to a file with BUSCO IDs to exclude (optional).

**apply_vcf_to_busco:**
- `apply_vcf_iupac`: Set `True` to encode heterozygous SNPs as IUPAC ambiguity codes (equivalent to GATK `--use-iupac-sample`). Default: `False` (ALT allele used for het/hom-alt calls).
- `vcf_reconstruct_ref_as_species`: Set `True` to include the reference genome itself as a sample in the phylogeny. Default: `False`.

**Alignment** (parameters passed directly to the chosen tool):
- `prank_params`, `mafft_params`, `muscle_params`

> **PRANK timeout:** PRANK can be very slow on long genes when run in codon alignment mode (`-codon`). To prevent the pipeline from stalling, PRANK is automatically terminated 15 minutes before the time limit set by `prank_time` (default: `"100h"`). Genes that hit the timeout are discarded from the analysis.

**Filtration:**
- `clipkit_params`, `gblocks_params`, `trimal_params`

**Phylogenetic inference:**
- `iqtree_params`: e.g. `"-keep-ident -m TESTNEW -bb 1000"`. Add `-o 'OUTGROUP'` to set an outgroup.
- `astral_params`: e.g. `"--support 2"`. Add `--root 'OUTGROUP'` to set an outgroup.
- `raxml_params`: e.g. `"--model GTR+G --bs-trees 100"`.
- `rapidnj_params`: e.g. `"-b 1000"`.
- `phylip_dnadist_params`: Use `"D\n"` for Kimura 2-parameter model, or `""` for F84 (default).
- `phylip_neighbor_params`: Use `"N"` for UPGMA, or `""` for NJ (default).
- `mrbayes_params`, `mrbayes_block`: MrBayes configuration block file and extra parameters.

**Visualization:**
- `tree_visualization_params`: Specify outgroup as `"--outgroup OUTGROUP"`.

#### Directory structure

Input and output paths are defined here. The defaults are:

```yaml
# Input
genome_dir:          "input/genomes/"
vcf_reconstruct_dir: "input/vcf_reconstruct/"
vcf2phylip_dir:      "input/vcf2phylip/"

# Output (all under results/)
output_dir:          "results/"
```

It is recommended to leave the directory structure unchanged.

#### Resources

Per-tool Slurm settings: partition (`*_queue`), threads (`*_threads`), memory in MB (`*_mem_mb`), and runtime (`*_time`). Adjust these to match your cluster configuration. Note that BUSCO and PRANK are the most time-consuming steps and may require generous time limits (default: `150h` and `100h` respectively).

### Step 4. Execute workflow

Install Snakemake:

```bash
mamba create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake snakemake-executor-plugin-cluster-generic
mamba activate snakemake
```

Dry run to preview all steps:

```bash
snakemake --profile profile/slurm/ --configfile config/default.yaml --dry-run
```

Remove `--dry-run` to start the actual run.

---

## Advanced usage

### Starting from completed BUSCO results

Move genome assemblies (or create empty placeholder files) into `input/genomes/`, then place BUSCO output directories under `results/busco/`. Expected structure for `Ailurus_fulgens.fasta`:

```
results/
    busco/
        Ailurus_fulgens/
            busco_sequences/
                fragmented_busco_sequences/
                multi_copy_busco_sequences/
                single_copy_busco_sequences/
            hmmer_output/
            logs/
            metaeuk_output/
            full_table_Ailurus_fulgens.tsv
            missing_busco_list_Ailurus_fulgens.tsv
            short_summary_Ailurus_fulgens.txt
            short_summary.json
            short_summary.specific.mammalia_odb10.Ailurus_fulgens.json
            short_summary.specific.mammalia_odb10.Ailurus_fulgens.txt
```

---

## Contact

Please email: <andrey.tomarovsky@gmail.com> for questions or feedback.
