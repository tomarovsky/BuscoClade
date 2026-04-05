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
| `busco_dataset_path` | `"$TOOLS/busco_datasets/mammalia_odb12/"` | Path to pre-downloaded OrthoDB dataset |
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
| `nodes_filtrataion_by_support` | `70` | Minimum bootstrap support for collapsing nodes before ASTRAL |
| `astral_params` | `"--support 2"` | ASTRAL-IV flags; add `--root 'OUTGROUP'` to set an outgroup |
| `raxml_params` | `"--model GTR+G --bs-trees 100"` | RAxML-NG flags |
| `rapidnj_params` | `"-b 1000"` | RapidNJ flags |
| `phylip_dnadist_params` | `"D\n"` | `"D\n"` for Kimura 2-parameter model; `""` for F84 (default) |
| `phylip_neighbor_params` | `""` | `"N\n"` for UPGMA; `""` for NJ (default) |
| `mrbayes_params` | `""` | Extra MrBayes flags |
| `mrbayes_block` | `"resources/mrbayes.block"` | Path to MrBayes block configuration file |
| `mrbayes_path` | `"$TOOLS/mrbayes/mrbayes"` | Path to MrBayes binary |

### Visualization

| Parameter | Default | Description |
|---|---|---|
| `tree_visualization_params` | `""` | Add `"--outgroup 'Species name'"` to root trees; multiple outgroups: `"--outgroup 'sp1,sp2'"` |
