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
- **VCF-based SNP application:** [apply_vcf_to_busco.py](Advanced-Usage#vcf-based-snp-application), [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)
- **Alignment:** [MAFFT](https://mafft.cbrc.jp/alignment/software/), [MUSCLE](https://doi.org/10.1038/s41467-022-34630-w), [PRANK](http://wasabiapp.org/software/prank/)
- **Trimming:** [ClipKIT](https://github.com/JLSteenwyk/ClipKIT), [TrimAl](http://trimal.cgenomics.org/), [GBlocks](https://academic.oup.com/mbe/article/17/4/540/1127654)
- **Phylogenetic tree construction:** [IQTree](http://www.iqtree.org/), [MrBayes](https://nbisweden.github.io/MrBayes/), [ASTRAL-IV](https://doi.org/10.1093/molbev/msaf172), [RapidNJ](https://birc.au.dk/software/rapidnj), [PHYLIP](https://phylipweb.github.io/phylip/), [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- **Visualization:** [Etetoolkit](http://etetoolkit.org/), [Matplotlib](https://matplotlib.org/stable/)

---

## Wiki contents

- [[Usage]] — input preparation, configuration, running the pipeline
- [[Configuration]] — all config parameters with defaults and descriptions
- [[Apptainer]] — running BuscoClade via pre-built Apptainer container
- [[Advanced-Usage]] — starting from existing BUSCO results, gap-aware AltRef insertion
