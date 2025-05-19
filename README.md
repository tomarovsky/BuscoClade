# Snakemake workflow: BuscoClade

[![Snakemake](https://img.shields.io/badge/snakemake-<8.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# BuscoClade Pipeline Documentation

**BuscoClade** is a **Snakemake**-based workflow that constructs species phylogenies using *BUSCO* (Benchmarking Universal Single-Copy Orthologs). It runs multiple analysis stages—from preparing inputs and running BUSCO on each genome, through multiple-sequence alignment (MSA), trimming, and tree-building—to produce a final phylogenetic tree and related visualizations. BuscoClade is organized as modular Snakemake rules; each rule corresponds to a step (e.g. alignment, filtering, tree inference). Snakemake automatically infers dependencies between rules based on input/output files. By leveraging Snakemake and Conda, BuscoClade ensures a reproducible, scalable pipeline that can run locally or on an HPC cluster with minimal user effort.

Key components and stages of the pipeline include:

* **Input Preparation:** Genome assemblies (FASTA files) and optionally variant data (VCF files) are placed in designated input directories.
* **[BUSCO](https://busco.ezlab.org/) Analysis:** Each genome is analyzed with BUSCO to identify conserved single-copy orthologous genes. BuscoClade can use either the *MetaEuk* or *AUGUSTUS* gene predictor (configurable). The pipeline collects single-copy BUSCO gene IDs that are common across all species.
* **Sequence Extraction:** BUSCO outputs are parsed to extract the DNA or protein sequences of the common single-copy genes from each genome.
* **Alignment:** Extracted gene sequences are aligned across species. The pipeline supports [PRANK](http://wasabiapp.org/software/prank/) or [MAFFT](https://mafft.cbrc.jp/alignment/software/) for alignment; users select by configuration (e.g. `dna_alignment: 'prank'` or `'mafft'`).
* **Filtering/Trimming:** Alignments are cleaned of poorly aligned or gap-rich regions. BuscoClade supports [GBlocks](https://academic.oup.com/mbe/article/17/4/540/1127654) or [TrimAl](http://trimal.cgenomics.org/) for trimming (configurable). This reduces noise before tree inference.
* **Concatenation:** Individual gene alignments are concatenated into a supermatrix (combined alignment). A concatenated alignment file (FASTA, PHYLIP, etc.) is produced for downstream analyses.
* **Phylogenetic Inference:** The pipeline can infer phylogenies using multiple methods:

  * **[IQTree](http://www.iqtree.org/) (Maximum Likelihood):** Fast ML tree reconstruction with model testing and bootstrapping.
  * **[ASTRAL III](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) (Coalescent):** Builds a species tree from gene trees (produced by IQ-TREE) using the ASTRAL algorithm for coalescent-based inference.
  * **[RapidNJ](https://birc.au.dk/software/rapidnj) (Neighbor Joining):** Quick distance-based tree (with bootstraps) from the concatenated alignment.
  * **[PHYLIP](https://phylipweb.github.io/phylip/) (Neighbor Joining/UPGMA):** Tree building via the classic PHYLIP package.
  * **[Raxml-NG](https://github.com/amkozlov/raxml-ng) (Maximum Likelihood):** Alternative ML tree builder.
  * **[MrBayes](https://nbisweden.github.io/MrBayes/) (Bayesian MCMC):** Optional Bayesian inference (requires more time).
* **VCF-Based Alternative:** Instead of BUSCO, BuscoClade can operate on aligned variant data. Provided multi-sample VCFs can be split per sample and converted to a phylogenetic alignment (via [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)), enabling phylogeny from SNP data.
* **Visualization:** Final trees are visualized ([Etetoolkit](http://etetoolkit.org/)) and saved as image files. The pipeline also produces summary plots such as BUSCO completeness histograms and charts of ortholog presence using [Matplotlib](https://matplotlib.org/stable/) scripts.




## Quick Start Guide

Follow these steps to run BuscoClade on your data. This guide assumes basic familiarity with the command line.

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/Zubiriguitar/BuscoClade.git
   cd BuscoClade
   ```

2. **Install Requirements:**

   BuscoClade relies on Snakemake and various bioinformatics tools. It is recommended to use Conda to install dependencies. You can either install Snakemake and tools globally or let Snakemake create environments. For example:

   ```bash
   # (Optional) Create a new Conda environment with Snakemake
   conda create -n buscoclade_env snakemake python=3.8 -y
   conda activate buscoclade_env

   # Or install Snakemake and Mamba/Conda globally as needed
   ```

3. **Prepare Input Data:**

   * Create the input directory structure: put your genome FASTA files in `input/genomes/`. Each file should have a `.fasta` extension. For example:

     ```
     input/genomes/
         species_1_.fasta
         ...
         species_n_.fasta
     ```
   * (Optional) If you have variant data instead of raw genomes, place VCF files in .gz archive format with reference FASTA file if using GATK into `input/vcf_reconstruct/*any_folder*`. The pipeline will split multi-sample VCFs into per-sample VCFs and reconstruct alignments from them. Genome input is also an option combining with VCF. 
   
    ```
     input/genomes/
         species_1_.fasta
         ...
         species_n_.fasta
     input/vcf_reconstruct/*folder_1*
         reference_1.fasta
         vcf_file_1.vcf.gz
         ...
         vcf_file_m.vcf.gz
      input/vcf_reconstruct/*folder_k*
         reference_k.fasta
         vcf_file_1.vcf.gz
         ...
         vcf_file_i.vcf.gz
     ```

   * If you using vcf2phylip yo are able to place only one VCF `input/vcf_reconstruct/*any_folder*`

4. **Download BUSCO Datasets:**

   BuscoClade requires a BUSCO lineage dataset. Download the appropriate BUSCO dataset (e.g. `saccharomycetes_odb10`) from the BUSCO website and note its path. Update the `busco_dataset_path` in the config (see next step).

5. **Configure the Workflow:**

   * Open `config/default.yaml` in an editor. Adjust the parameters as needed (see next section for details). At minimum, set:

     ```yaml
     genome_dir: "input/genomes/"
     busco_dataset_path: "/path/to/busco/saccharomycetes_odb10/"
     output_dir: "results/"
     ```
   * You can leave many defaults unchanged on first run. The default config enables DNA alignment with PRANK, trimming with Gblocks, and tree inference with IQ-TREE, ASTRAL, RapidNJ, PHYLIP, and RAxML.

6. **Run the Pipeline:**

   Use Snakemake to execute the workflow. For a local run (no cluster), a simple command is:

   ```bash
   snakemake --cores 4 --configfile config/default.yaml
   ```

   This tells Snakemake to use at most 4 CPU cores and to create/use Conda environments as specified. The `--use-conda` flag enables environment management for reproducibility. Snakemake will execute the pipeline steps in order.

   If you have a SLURM cluster, you can use the provided profile:

   ```bash
   snakemake --profile profile/slurm --configfile config/default.yaml
   ```

   This submits jobs to SLURM as per the `profile/slurm/config.yaml` settings. Adjust the `profile/slurm/config.yaml` (queues, time, etc.) if needed for your cluster.

7. **Inspect Results:**

   After completion, see the `results/` directory (or the directory you set). It will contain subdirectories for alignments, phylogenetic trees, and summary outputs (see **Output Description** below).

Throughout execution, Snakemake will produce logs in the `logs/` directory (and cluster logs in `cluster_logs/`). You can resume or rerun specific steps by adjusting the command (e.g. adding `-R` to rerun certain rules).

## Contact

Please email me at: <andrey.tomarovsky@gmail.com> for any questions or feedback.

