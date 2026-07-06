#!/usr/bin/env python3
"""
busco_reconstruct_common.py

Shared helpers for reconstructing BUSCO single-copy sequences from a reference
genome's MetaEuk BUSCO output. Used by both reconstruction routes:

    - apply_vcf_to_busco.py       — substitutions come from per-sample VCF SNPs
    - apply_consensus_to_busco.py — substitutions come from an already
                                    reconstructed (consensus) genome that is
                                    coordinate-identical to the reference

Both routes rely on the same invariant: only substitutions are applied (no
indels), so every reconstructed gene keeps the reference's ungapped length and
exon structure. This module holds everything that is identical between the two:
MetaEuk GFF/CDS parsing, FNA header parsing, exon-aware coordinate mapping, and
(ambiguous) translation.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

from Bio.Data import CodonTable

# ---------------------------------------------------------------------------
# IUPAC helpers
# ---------------------------------------------------------------------------

IUPAC_MAP = {
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("AC"): "M",
    frozenset("GT"): "K",
    frozenset("AT"): "W",
    frozenset("CG"): "S",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("AGT"): "D",
    frozenset("CGT"): "B",
    frozenset("ACGT"): "N",
}
IUPAC_EXPAND = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "W": "AT",
    "S": "CG",
    "H": "ACT",
    "V": "ACG",
    "D": "AGT",
    "B": "CGT",
    "N": "ACGT",
}

STD_TABLE = CodonTable.unambiguous_dna_by_id[1]

# Complement table covering IUPAC ambiguity codes (upper and lower case).
COMP = str.maketrans(
    "ACGTacgtRYMKWSBDHVNrykmwsbdhvn",
    "TGCAtgcaYRKMWSVHDBNyrkmwsvhdbn",
)


def bases_to_iupac(alleles: set) -> str:
    """Return IUPAC code for a set of nucleotide bases."""
    alleles_upper = frozenset(a.upper() for a in alleles)
    if len(alleles_upper) == 1:
        return next(iter(alleles_upper))
    return IUPAC_MAP.get(alleles_upper, "N")


def reverse_complement(seq: str) -> str:
    """Reverse-complement a nucleotide sequence (IUPAC-aware)."""
    return seq.translate(COMP)[::-1]


def translate_ambiguous_codon(codon: str) -> str:
    """
    Translate a codon that may contain IUPAC ambiguity codes.
    Returns a single-letter amino acid if all expansions agree, else 'X'.
    """
    codon = codon.upper()
    if all(c in "ACGT" for c in codon):
        try:
            return STD_TABLE.forward_table[codon]
        except KeyError:
            return "*" if codon in STD_TABLE.stop_codons else "X"

    options = [""]
    for base in codon:
        expansions = IUPAC_EXPAND.get(base, "N")
        options = [prev + b for prev in options for b in expansions]

    aas = set()
    for opt in options:
        if opt in STD_TABLE.stop_codons:
            aas.add("*")
        else:
            aas.add(STD_TABLE.forward_table.get(opt, "X"))

    return next(iter(aas)) if len(aas) == 1 else "X"


def translate_sequence(nuc_seq: str) -> str:
    """Translate a nucleotide sequence (may contain IUPAC codes) to protein."""
    nuc_seq = nuc_seq.upper()
    protein = []
    for i in range(0, len(nuc_seq) - 2, 3):
        codon = nuc_seq[i : i + 3]
        aa = translate_ambiguous_codon(codon)
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


# ---------------------------------------------------------------------------
# MetaEuk GFF loading
# ---------------------------------------------------------------------------

def find_metaeuk_gffs(metaeuk_output: Path) -> tuple[Optional[Path], Optional[Path]]:
    """
    Locate rerun_results and initial_results GFF files inside metaeuk_output.
    Returns (rerun_gff_path, initial_gff_path), either may be None if absent.
    """
    rerun_gff   = None
    initial_gff = None

    for subdir_name, target in [("rerun_results", "rerun_gff"),
                                  ("initial_results", "initial_gff")]:
        subdir = metaeuk_output / subdir_name
        if not subdir.is_dir():
            continue
        gffs = list(subdir.glob("*.gff"))
        if not gffs:
            continue
        if len(gffs) > 1:
            print(f"  WARNING: multiple GFF files in {subdir}, using {gffs[0].name}")
        if target == "rerun_gff":
            rerun_gff = gffs[0]
        else:
            initial_gff = gffs[0]

    return rerun_gff, initial_gff


def load_cds_index(gff_path: Path) -> dict:
    """
    Parse a MetaEuk GFF file and return a CDS index:
        {(target_id, chrom): [(start1, end1), ...]}
    Coordinates are 1-based inclusive (as in GFF).
    Only 'CDS' feature rows are used.
    """
    index = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            m = re.search(r"Target_ID=([^;]+)", parts[8])
            if not m:
                continue
            tid   = m.group(1)
            chrom = parts[0]
            start = int(parts[3])
            end   = int(parts[4])
            index[(tid, chrom)].append((start, end))
    return dict(index)


def build_cds_lookup(metaeuk_output: Path) -> tuple[dict, dict]:
    """
    Load CDS indices from rerun_results and initial_results GFFs.
    Returns (rerun_index, initial_index).
    """
    rerun_gff, initial_gff = find_metaeuk_gffs(metaeuk_output)

    if rerun_gff is None and initial_gff is None:
        sys.exit(f"ERROR: No MetaEuk GFF files found in {metaeuk_output}")

    rerun_index   = load_cds_index(rerun_gff)   if rerun_gff   else {}
    initial_index = load_cds_index(initial_gff) if initial_gff else {}

    print(f"Loaded CDS index: rerun={len(rerun_index)} keys, "
          f"initial={len(initial_index)} keys")
    return rerun_index, initial_index


def get_cds_exons(tid: str, chrom: str, gene_start1: int, gene_end1: int,
                  rerun_index: dict, initial_index: dict) -> list:
    """
    Return the ordered list of CDS exon intervals (1-based inclusive) for a
    gene, using rerun_results first and initial_results as fallback.

    Exons are filtered to those that fall within [gene_start1, gene_end1] to
    handle cases where a single Target_ID has predictions on multiple loci on
    the same chromosome.

    For minus-strand genes MetaEuk lists exons in descending genomic order
    (exon_0 has the highest coordinates). We preserve that order here because
    the FNA sequence was built in the same order.

    Raises ValueError if no matching exons are found.
    """
    key = (tid, chrom)

    for source_name, index in [("rerun", rerun_index), ("initial", initial_index)]:
        if key not in index:
            continue
        exons = index[key]
        # Filter to gene's genomic range
        filtered = [(s, e) for s, e in exons
                    if s >= gene_start1 and e <= gene_end1 + 1]
        if not filtered:
            continue
        return filtered, source_name

    raise ValueError(
        f"No CDS exons found for {tid} on {chrom}:{gene_start1}-{gene_end1} "
        f"in either rerun or initial GFF."
    )


# ---------------------------------------------------------------------------
# MetaEuk header parsing
# ---------------------------------------------------------------------------

def parse_metaeuk_fna_header(header: str) -> dict:
    """
    Parse MetaEuk FNA header:
        8689at40674_1868482_0:0013eb|NC_037328.1:120896888-120904250|-
    Returns dict with keys: full_id, chrom, start0, end0, strand.
    start0/end0 are 0-based (as stored in the header).
    """
    header = header.lstrip(">").strip()
    parts = header.split("|")
    if len(parts) < 3:
        raise ValueError(f"Unexpected header format: {header}")

    full_id    = parts[0]
    coord_part = parts[1]
    strand     = parts[2]

    m = re.match(r"^(.+):(\d+)-(\d+)$", coord_part)
    if not m:
        raise ValueError(f"Cannot parse coords from: {coord_part}")

    return {
        "full_id":          full_id,
        "chrom":            m.group(1),
        "start0":           int(m.group(2)),  # 0-based
        "end0":             int(m.group(3)),  # 0-based inclusive
        "strand":           strand,
        "original_header":  header,
    }


# ---------------------------------------------------------------------------
# Exon-aware SNP application
# ---------------------------------------------------------------------------

def apply_snps_exon_aware(original_seq: str, strand: str,
                           exons: list, snps_by_region: dict) -> str:
    """
    Apply SNPs to a spliced CDS sequence using per-exon genomic coordinates.

    Parameters
    ----------
    original_seq : str
        The CDS sequence as stored in the FNA file (already on the gene's
        strand, exons concatenated in transcript order).
    strand : str
        '+' or '-'.
    exons : list of (start1, end1)
        CDS exon intervals in 1-based inclusive coordinates, in the same
        order as they appear in the FNA sequence (transcript order).
        For minus-strand genes this means descending genomic coordinates
        (highest coords first), matching MetaEuk's exon_0..exon_N order.
    snps_by_region : dict
        Pre-loaded SNPs per exon: {exon_index: {pos1: base}}.

    Returns
    -------
    str
        Modified CDS sequence with SNPs applied.
    """
    seq_list = list(original_seq.upper())
    fna_offset = 0  # running position in the FNA sequence

    for exon_idx, (exon_start1, exon_end1) in enumerate(exons):
        exon_len = exon_end1 - exon_start1 + 1
        snps = snps_by_region.get(exon_idx, {})

        if strand == "+":
            # FNA bases for this exon run left-to-right in genomic order.
            # pos1 maps directly: fna_idx = fna_offset + (pos1 - exon_start1)
            for pos1, base in snps.items():
                fna_idx = fna_offset + (pos1 - exon_start1)
                if 0 <= fna_idx < fna_offset + exon_len:
                    seq_list[fna_idx] = base

        else:  # strand == "-"
            # MetaEuk outputs the reverse-complement of each exon.
            # Within this exon's FNA segment, position fna_offset corresponds
            # to the highest genomic coordinate (exon_end1) and
            # fna_offset + exon_len - 1 corresponds to exon_start1.
            # A SNP at genomic pos1 maps to:
            #   fna_idx = fna_offset + (exon_end1 - pos1)
            # The base must also be complemented.
            for pos1, base in snps.items():
                fna_idx = fna_offset + (exon_end1 - pos1)
                if 0 <= fna_idx < fna_offset + exon_len:
                    seq_list[fna_idx] = base.translate(COMP)

        fna_offset += exon_len

    return "".join(seq_list)
