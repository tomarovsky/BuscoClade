# Snakemake workflow: BuscoClade

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description

Pipeline to construct species phylogenies using [BUSCO](https://busco.ezlab.org/).

- Alignment: [PRANK](http://wasabiapp.org/software/prank/), [MAFFT](https://mafft.cbrc.jp/alignment/software/).
- Trimming: [GBlocks](https://academic.oup.com/mbe/article/17/4/540/1127654), [TrimAl](http://trimal.cgenomics.org/).
- Phylogenetic tree constraction: [IQTree](http://www.iqtree.org/), [MrBayes](https://nbisweden.github.io/MrBayes/), [ASTRAL III](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y), [PHYLIP](https://phylipweb.github.io/phylip/).
- Visualization: [Etetoolkit](http://etetoolkit.org/).

## Dependencies

Install Snakemake:

```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

## Usage

### Step 1. Deploy workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/tomarovsky/BuscoClade/releases). Or clone the repository:

```
git clone https://github.com/tomarovsky/BuscoClade.git
```

### Step 2. Add species genomes

Add the unpacked FASTA genome assemblies to the 'genomes/' directory. Please note that the file prefixes will be used in the output phylogeny.

### Step 3. Configure workflow

To set up the workflow, customize the config/default.yaml file according to your needs, following the explanations provided below.

- **Pipeline Configuration:**
This section outlines the workflow. By default, it includes alignments, nucleotide sequence filtering, and all tools for phylogeny reconstruction, except for MrBayes (it is recommended to run the GPU compiled version separately). To disable a tool, set its value to `False` or comment out the corresponding line.

- **Tool Parameters:**
This section contains parameters for each tool. Key considerations are:
  - `busco_dataset_path`: Download the BUSCO dataset beforehand and specify its path here.
  - `busco_params`: Use the --offline flag and the --download_path parameter, indicating the path to the busco_downloads/ directory.

- **Directory Structure:**
Define the structure of directories for output files. The default structure is as follows:

```
results/
    busco/
    ids/
        species_ids/
        merged_sequences/
        common_ids.ids
    alignments/
        raw/
        filtered/
    concat_alignments/
    phylogeny/
        astral/
        iqtree/
        phylip/
```

- **Resources:**
Specify the Slurm queue and resources for each tool that will be used during execution.

### Step 4. Execute workflow

To perform a dry run of the pipeline, execute the following command:

```
snakemake --profile profile/slurm/ --configfile config/default.yaml --dry-run
```

Snakemake will print all the rules that will be executed. To initiate the actual run, remove the `--dry-run` flag.

### FAQ

**1. How to run the workflow if I have completed BUSCOs?**

First, move the genome assemblies to the "genomes/" directory or create empty files with corresponding names. Then, create a "busco/" directory in the "results/" directory and move the BUSCO output directories into it. Note that BUSCO output must be formatted. Thus, for "Ailurus_fulgens.fasta" the output should look like this:

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
            short_summary_Ailurus_fulgens_DNAzoo.txt
            short_summary.json
            short_summary.specific.mammalia_odb10.Ailurus_fulgens.json
            short_summary.specific.mammalia_odb10.Ailurus_fulgens.txt
```

**2. Why does the tree visualization from PHYLIP return an error?**

This issue arises because PHYLIP crops the species names to the first 10 characters by default. To perform visualization, you must manually edit the output NEWICK tree and restart workflow. This will be fixed soon.

## Contact

Please email me at: <andrey.tomarovsky@gmail.com> for any questions or feedback.

