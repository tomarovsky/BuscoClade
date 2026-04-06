# Usage

BuscoClade can be run either directly with Snakemake or via a pre-built
Apptainer container — see [[Apptainer]] for instructions.

---

## Step 1. Deploy workflow

### Option A: Clone repository

```bash
git clone https://github.com/tomarovsky/BuscoClade.git
```

### Option B: Apptainer container

Pull the latest image and run immediately:

```bash
apptainer pull buscoclade.sif oras://ghcr.io/tomarovsky/buscoclade:latest
```

See [[Apptainer]] for full usage instructions.

---

## Step 2. Prepare input data

### FASTA assemblies

Place genome assemblies into `input/genomes/`. The file prefix is used as the
sample name in the output phylogeny. Supported extensions: `.fasta`, `.fna`,
`.fa`, and their gzipped versions (`.fasta.gz`, `.fna.gz`, `.fa.gz`).

### Per-sample VCFs + reference genome

If you have per-sample VCFs, the pipeline applies SNPs directly to the BUSCO
sequences of a reference genome using `apply_vcf_to_busco.py`. Each VCF file
must contain exactly one sample. Only SNPs are used — indels present in the
VCF are ignored. This avoids rebuilding full pseudo-genome assemblies and
re-running BUSCO for each sample — instead, BUSCO is run once on the
reference, and ortholog sequences for all other samples are reconstructed by
applying their SNPs to the reference BUSCO sequences exon-by-exon.

Place per-sample VCF files and the corresponding reference genome together
into a subdirectory under `input/vcf_reconstruct/`. Each subdirectory is
processed independently, which allows reconstructing sequences against
different references in a single run:

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

The directory name is used only for organization — the VCF file prefix
determines the sample name in the output phylogeny. No additional config
changes are needed; the pipeline detects subdirectories automatically.

BUSCO is run once on the reference genome in each subdirectory. The resulting
`single_copy_busco_sequences/` and `metaeuk_output/` are then used by
`apply_vcf_to_busco.py` to produce per-sample ortholog sequences.

### Multi-sample VCF via vcf2phylip (optional)

As an alternative to the BUSCO-based pipeline, a concatenated phylip alignment
can be built directly from a multi-sample VCF using `vcf2phylip.py`, bypassing
BUSCO and sequence alignment entirely. When `vcf2phylip: True` is set, only
this route is executed.

Place exactly one multi-sample `.vcf.gz` file into `input/vcf2phylip/` and
enable the option in the config:

```
input/
    vcf2phylip/
        all_samples.vcf.gz    # exactly one multi-sample VCF
```

```yaml
vcf2phylip: True
```

---

## Step 3. Configure workflow

Modify `config/default.yaml` (recommended: copy it and pass with
`--configfile`). See [[Configuration]] for a full parameter reference.

### Directory structure

Input and output paths are defined in the config. The defaults are:

```yaml
# Input
genome_dir:          "input/genomes/"
vcf_reconstruct_dir: "input/vcf_reconstruct/"
vcf2phylip_dir:      "input/vcf2phylip/"

# Output (all under results/)
output_dir:          "results/"
```

It is recommended to leave the directory structure unchanged.

### Resources

Per-tool Slurm/PBS settings: partition (`*_queue`), threads (`*_threads`), memory
in MB (`*_mem_mb`), and runtime (`*_time`). Adjust these to match your cluster
configuration. Note that BUSCO and PRANK are the most time-consuming steps and
may require generous time limits (default: `150h` and `100h` respectively).

---

## Step 4. Execute workflow

Install Snakemake:

```bash
mamba create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake snakemake-executor-plugin-cluster-generic
mamba activate snakemake
```

Dry run to preview all steps:

```bash
snakemake --configfile config/default.yaml --dry-run
```

Remove `--dry-run` to start the actual run.

### Running on a cluster
 
BuscoClade includes ready-made profiles for Slurm and PBS. Pass `--profile`
instead of `--cores` to submit jobs to the scheduler:
 
```bash
# Slurm
snakemake --profile profiles/slurm/ --configfile config/default.yaml
 
# PBS
snakemake --profile profiles/pbs/ --configfile config/default.yaml
```
 
Resource allocation (partition, threads, memory, runtime) is configured
per-tool in the config file — see [[Configuration#cluster-resources]].
