# Configuration

All parameters are set in `config/default.yaml`. It is recommended to copy
the file and pass it with `--configfile` rather than editing it directly.

---

## Pipeline configuration

| Option | Default | Choices | Description |
|---|---|---|---|
| `vcf2phylip` | `False` | `True` / `False` | Use multi-sample VCF instead of BUSCO-based pipeline |
| `quastcore` | `True` | `True` / `False` | Compute assembly statistics |
| `alignment` | `"mafft"` | `mafft`, `muscle`, `prank` | Multiple alignment tool |
| `filtration` | `"clipkit"` | `clipkit`, `trimal`, `gblocks` | Alignment trimming tool |
| `iqtree` | `True` | `True` / `False` | Run IQ-TREE |
| `astral` | `True` | `True` / `False` | Run ASTRAL-IV (multispecies coalescent) |
| `rapidnj` | `True` | `True` / `False` | Run RapidNJ |
| `phylip` | `True` | `True` / `False` | Run PHYLIP |
| `raxml` | `True` | `True` / `False` | Run RAxML-NG |
| `mrbayes` | `False` | `True` / `False` | Run MrBayes (recommended to run GPU-compiled version separately) |
| `draw_phylotrees` | `True` | `True` / `False` | Visualize output trees |

---

## Tool parameters

### BUSCO

| Parameter | Default | Description |
|---|---|---|
| `busco_dataset_path` | `"path/to/busco_datasets/ortho_odb12/"` | Path to pre-downloaded OrthoDB dataset |
| `busco_options` | `"--offline"` | Extra BUSCO flags; use `"--offline"` to run without internet |
| `busco_mode` | `"genome"` | BUSCO mode |
| `busco_blacklist` | `"input/BUSCO.blacklist"` | File with BUSCO IDs to exclude (optional) |
| `busco_histogram_colors` | `"#23b4e8,#008dbf,#fbbc04,#ea4335"` | Bar colors for S, D, F, M in the BUSCO histogram |

### VCF reconstruction

| Parameter | Default | Description |
|---|---|---|
| `apply_vcf_iupac` | `False` | Encode heterozygous SNPs as IUPAC ambiguity codes (equivalent to GATK `--use-iupac-sample`); if `False`, ALT allele is used for het/hom-alt calls |
| `vcf_reconstruct_ref_as_species` | `False` | Include the reference genome itself as a sample in the output phylogeny |
| `altref_gapaware_insertion` | `False` | Insert AltRef sequences into alignments using gap positions of the corresponding reference instead of re-aligning; recommended when working with many VCF samples; see [[Advanced-Usage#gap-aware-altref-insertion]] |

### Alignment

| Parameter | Default | Description |
|---|---|---|
| `mafft_params` | `"--reorder --auto"` | Flags passed to MAFFT |
| `muscle_params` | `""` | Flags passed to MUSCLE |
| `prank_params` | `"-codon"` | Flags passed to PRANK |
| `prank_time` | `"100h"` | PRANK is killed 15 min before this limit; genes that time out are discarded |

### Filtration

| Parameter | Default | Description |
|---|---|---|
| `clipkit_params` | `"--mode smart-gap"` | Flags passed to ClipKIT; add `--codon` for codon-aware trimming |
| `trimal_params` | `"-automated1"` | Flags passed to TrimAl |
| `gblocks_params` | `"-t=Codons"` | Flags passed to GBlocks |

### Phylogenetic inference

| Parameter | Default | Description |
|---|---|---|
| `iqtree_params` | `"-keep-ident -m TESTNEW -bb 1000"` | IQ-TREE flags for the concatenated alignment; add `-o 'OUTGROUP'` to set an outgroup |
| `iqtree_per_fna_params` | `"-keep-ident -m TESTNEW -bb 1000"` | IQ-TREE flags for per-gene trees (used as input to ASTRAL) |
| `nodes_filtration_by_support` | `70` | Minimum bootstrap support for collapsing nodes before ASTRAL |
| `astral_params` | `"--support 2"` | ASTRAL-IV flags; add `--root 'OUTGROUP'` to set an outgroup |
| `raxml_params` | `"--model GTR+G --bs-trees 100"` | RAxML-NG flags |
| `rapidnj_params` | `"-b 1000"` | RapidNJ flags |
| `phylip_dnadist_params` | `"D\n"` | `"D\n"` for Kimura 2-parameter model; `""` for F84 (default) |
| `phylip_neighbor_params` | `""` | `"N\n"` for UPGMA; `""` for NJ (default) |
| `mrbayes_params` | `""` | Extra MrBayes flags |
| `mrbayes_block` | `"resources/mrbayes.block"` | Path to MrBayes block configuration file |
| `mrbayes_path` | `"path/to/mrbayes-3.2.7/mrbayes"` | Path to MrBayes binary |

### Visualization

| Parameter | Default | Description |
|---|---|---|
| `tree_visualization_params` | `""` | Add `"--outgroup 'Species name'"` to root trees; multiple outgroups: `"--outgroup 'sp1,sp2'"` |

---
 
## Cluster resources
 
Per-tool resource allocation for Slurm and PBS cluster execution. These
parameters are passed automatically via the Snakemake profile — you do not
need to set them manually on the command line.
 
### Partitions / queues
 
| Parameter | Default | Description |
|---|---|---|
| `processing_queue` | `"main"` | Partition for lightweight processing jobs |
| `busco_queue` | `"main"` | Partition for BUSCO |
| `alignment_queue` | `"main"` | Partition for alignment |
| `filtration_queue` | `"main"` | Partition for trimming |
| `iqtree_queue` | `"main"` | Partition for IQ-TREE |
| `astral_queue` | `"main"` | Partition for ASTRAL |
| `rapidnj_queue` | `"main"` | Partition for RapidNJ |
| `phylip_queue` | `"main"` | Partition for PHYLIP |
| `raxml_queue` | `"main"` | Partition for RAxML-NG |
| `mrbayes_queue` | `"main"` | Partition for MrBayes |
 
### Threads
 
| Parameter | Default |
|---|---|
| `processing_threads` | `1` |
| `busco_threads` | `8` |
| `mafft_threads` | `1` |
| `muscle_threads` | `1` |
| `prank_threads` | `1` |
| `clipkit_threads` | `1` |
| `trimal_threads` | `1` |
| `gblocks_threads` | `1` |
| `iqtree_threads` | `8` |
| `iqtree_per_fna_threads` | `1` |
| `astral_threads` | `4` |
| `rapidnj_threads` | `4` |
| `phylip_threads` | `1` |
| `raxml_threads` | `4` |
| `mrbayes_threads` | `8` |
 
### Memory (MB)
 
| Parameter | Default |
|---|---|
| `processing_mem_mb` | `2000` |
| `busco_mem_mb` | `10000` |
| `mafft_mem_mb` | `2000` |
| `muscle_mem_mb` | `2000` |
| `prank_mem_mb` | `2000` |
| `clipkit_mem_mb` | `2000` |
| `trimal_mem_mb` | `2000` |
| `gblocks_mem_mb` | `2000` |
| `iqtree_mem_mb` | `10000` |
| `iqtree_per_fna_mem_mb` | `2000` |
| `astral_mem_mb` | `10000` |
| `rapidnj_mem_mb` | `8000` |
| `phylip_mem_mb` | `4000` |
| `raxml_mem_mb` | `10000` |
| `mrbayes_mem_mb` | `10000` |
 
### Runtime
 
| Parameter | Default |
|---|---|
| `processing_time` | `"5h"` |
| `busco_time` | `"150h"` |
| `mafft_time` | `"10h"` |
| `muscle_time` | `"10h"` |
| `prank_time` | `"100h"` |
| `clipkit_time` | `"10h"` |
| `trimal_time` | `"10h"` |
| `gblocks_time` | `"10h"` |
| `iqtree_time` | `"100h"` |
| `iqtree_per_fna_time` | `"100h"` |
| `astral_time` | `"50h"` |
| `rapidnj_time` | `"150h"` |
| `phylip_time` | `"50h"` |
| `raxml_time` | `"100h"` |
| `mrbayes_time` | `"100h"` |
