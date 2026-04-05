# Advanced Usage

## Starting from completed BUSCO results

Move genome assemblies (or create empty placeholder files) into
`input/genomes/`, then place BUSCO output directories under `results/busco/`.
Expected structure for `Ailurus_fulgens.fasta`:

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

## Gap-aware AltRef insertion

When `altref_gapaware_insertion: True`, BuscoClade uses an alternative
strategy to incorporate VCF-reconstructed (AltRef) sequences into multiple
alignments. Instead of including them in the aligner input — which multiplies
wall time linearly with the number of VCF samples — the pipeline runs the
aligner only on the reference genomes and FASTA-assembled species, then inserts
each AltRef sequence into the finished alignment by copying the gap pattern of
its reference.

This is always valid in BuscoClade because `apply_vcf_to_busco.py` uses only
SNPs — indels present in the VCF are ignored. Each input VCF must contain
exactly one sample. As a result, every AltRef gene has exactly the same
ungapped length as the corresponding reference gene, and the reference gap
positions in the alignment are directly transferable to the AltRef sequence.

### Execution flow with `altref_gapaware_insertion: True`

```
merged_sequences/          ← refs + FASTA species only (no AltRef)
        │
        ▼
alignments/pre_altref/     ← aligner output (MAFFT / MUSCLE / PRANK)
        │
        │  add_altref_to_alignment.py
        ▼
alignments/raw/            ← AltRef sequences inserted using reference gap positions
```

> References removed if `vcf_reconstruct_ref_as_species: False`

The script `workflow/scripts/add_altref_to_alignment.py` handles the insertion
step. For each reference present in the alignment it:

1. Reads the aligned reference sequence and records gap positions.
2. Loads the unaligned AltRef gene from the BUSCO output directory.
3. Validates that the ungapped length matches the reference.
4. Inserts gaps at the same positions and appends the sequence to the alignment.
5. Matches the case (upper / lower) of the inserted sequence to the alignment convention.
6. Removes the reference sequence from the output if `vcf_reconstruct_ref_as_species: False`.

Within the pipeline, every gene in the alignment is guaranteed to exist for
every AltRef sample — this is ensured by `common_ids`, which is built from the
intersection of genes present across all species.

### When to use this option

| Scenario | Recommendation |
|---|---|
| Many VCF samples (>5–10) | `altref_gapaware_insertion: True` — significant speedup |
| Few VCF samples | Either setting works; `False` is simpler |
