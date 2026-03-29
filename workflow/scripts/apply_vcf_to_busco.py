#!/usr/bin/env python3
"""
apply_vcf_to_busco.py

Applies SNPs from a VCF file to BUSCO single-copy sequences (FNA/FAA/GFF)
produced by MetaEuk, without rebuilding the full alternate reference.

Requirements:
    pip install biopython pysam

Usage:
    python apply_vcf_to_busco.py \\
        --busco-dir /path/to/single_copy_busco_sequences \\
        --ref reference.fa \\
        --vcf sample.vcf.gz \\
        --sample SAMPLE_NAME \\
        --output-dir /path/to/output \\
        [--iupac]

Notes:
    - Only SNPs are handled; indels are ignored (as in the original pipeline).
    - With --iupac, heterozygous SNPs are encoded as IUPAC ambiguity codes in FNA,
      and FAA uses ambiguous translation (X if ambiguous amino acid, otherwise
      the unambiguous amino acid if all codon variants agree).
    - Without --iupac, the ALT allele is used for homozygous alt and het calls
      (matching GATK FastaAlternateReferenceMaker default behaviour).
    - GFF files are copied unchanged (SNP-only: coordinates never shift).
    - The VCF must be bgzipped and tabix-indexed, or a plain .vcf (slower).
"""

import argparse
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

import pysam
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq

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


def bases_to_iupac(alleles: set) -> str:
    """Return IUPAC code for a set of nucleotide bases."""
    alleles_upper = frozenset(a.upper() for a in alleles)
    if len(alleles_upper) == 1:
        return next(iter(alleles_upper))
    return IUPAC_MAP.get(alleles_upper, "N")


def translate_ambiguous_codon(codon: str) -> str:
    """
    Translate a codon that may contain IUPAC ambiguity codes.
    Returns a single-letter amino acid if all expansions agree, else 'X'.
    """
    codon = codon.upper()
    # Fast path: no ambiguity
    if all(c in "ACGT" for c in codon):
        try:
            return STD_TABLE.forward_table[codon]
        except KeyError:
            return "*" if codon in STD_TABLE.stop_codons else "X"

    # Expand ambiguous positions
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
# VCF / SNP application
# ---------------------------------------------------------------------------


def load_vcf_snps(vcf_path: str, chrom: str, start1: int, end1: int, sample: Optional[str], iupac: bool) -> dict:
    """
    Load SNPs from VCF for a given region (1-based, inclusive).
    Returns dict {pos1 (int): base (str)}.

    With iupac=True and a heterozygous call, returns IUPAC code.
    Without iupac, returns ALT allele for het/hom-alt, REF for hom-ref.
    """
    snps = {}
    try:
        vcf = pysam.VariantFile(vcf_path)
    except Exception as e:
        sys.exit(f"ERROR: Cannot open VCF {vcf_path}: {e}")

    # Resolve sample index
    sample_idx = None
    if sample:
        samples = list(vcf.header.samples)
        if sample not in samples:
            sys.exit(f"ERROR: Sample '{sample}' not found in VCF. Available: {samples}")
        sample_idx = samples.index(sample)
    else:
        if len(vcf.header.samples) == 1:
            sample_idx = 0
        elif len(vcf.header.samples) > 1:
            sys.exit("ERROR: VCF has multiple samples. Please specify --sample.")

    try:
        region_iter = vcf.fetch(chrom, start1 - 1, end1)
    except ValueError:
        # Chromosome not in VCF index — no SNPs in this region
        return snps

    for rec in region_iter:
        # Skip indels
        if len(rec.ref) != 1:
            continue
        alts = rec.alts
        if alts is None:
            continue
        if any(len(a) != 1 for a in alts):
            continue

        pos1 = rec.pos  # pysam uses 1-based pos

        if sample_idx is not None:
            sample_name = list(vcf.header.samples)[sample_idx]
            gt = rec.samples[sample_name]["GT"]
            if gt is None or all(g is None for g in gt):
                continue
            allele_indices = [g for g in gt if g is not None]
            allele_bases = set()
            for idx in allele_indices:
                if idx == 0:
                    allele_bases.add(rec.ref.upper())
                elif idx <= len(alts):
                    allele_bases.add(alts[idx - 1].upper())

            if iupac:
                snps[pos1] = bases_to_iupac(allele_bases)
            else:
                # Use ALT if any alt allele present
                alt_alleles = [alts[i - 1].upper() for i in allele_indices if i > 0]
                if alt_alleles:
                    snps[pos1] = alt_alleles[0]
                # else hom-ref: no change needed
        else:
            # No sample info: use first ALT
            if iupac:
                snps[pos1] = bases_to_iupac({rec.ref.upper(), alts[0].upper()})
            else:
                snps[pos1] = alts[0].upper()

    vcf.close()
    return snps


def apply_snps_to_sequence(seq: str, snps: dict, region_start1: int) -> str:
    """
    Apply SNP dict {pos1: base} to a sequence string.
    region_start1: 1-based genomic start of the sequence.
    Returns modified sequence (same length, SNP-only).
    """
    seq_list = list(seq.upper())
    for pos1, base in snps.items():
        idx = pos1 - region_start1
        if 0 <= idx < len(seq_list):
            seq_list[idx] = base
    return "".join(seq_list)


# ---------------------------------------------------------------------------
# MetaEuk header parsing
# ---------------------------------------------------------------------------


def parse_metaeuk_fna_header(header: str) -> dict:
    """
    Parse MetaEuk FNA/FAA header like:
        99985at4891_1095631_0:001c89|NC_079283.1:696170-696544|-
    Returns dict with keys: full_id, chrom, start0, end0, strand
    (start0 is 0-based, end0 is 0-based inclusive, matching header).
    """
    # Strip leading >
    header = header.lstrip(">").strip()
    # Split on |
    parts = header.split("|")
    if len(parts) < 3:
        raise ValueError(f"Unexpected header format: {header}")

    full_id = parts[0]
    coord_part = parts[1]  # e.g. NC_079283.1:696170-696544
    strand = parts[2]  # + or -

    m = re.match(r"^(.+):(\d+)-(\d+)$", coord_part)
    if not m:
        raise ValueError(f"Cannot parse coords from: {coord_part}")

    chrom = m.group(1)
    start0 = int(m.group(2))  # 0-based start (BED-style)
    end0 = int(m.group(3))  # 0-based inclusive end

    return {
        "full_id": full_id,
        "chrom": chrom,
        "start0": start0,
        "end0": end0,
        "strand": strand,
        "original_header": header,
    }


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------


def process_busco_gene(busco_id: str, busco_dir: Path, output_dir: Path, vcf_path: str, ref_path: str, sample: Optional[str], iupac: bool):
    """Process one BUSCO gene: update FNA, FAA; copy GFF unchanged."""

    fna_file = busco_dir / f"{busco_id}.fna"
    faa_file = busco_dir / f"{busco_id}.faa"
    gff_file = busco_dir / f"{busco_id}.gff"

    # --- Parse FNA ---
    records = list(SeqIO.parse(fna_file, "fasta"))
    if not records:
        print(f"  WARNING: empty FNA for {busco_id}, skipping.")
        return

    new_fna_records = []
    new_faa_records = []

    for rec in records:
        info = parse_metaeuk_fna_header(rec.description)
        chrom = info["chrom"]
        start0 = info["start0"]
        end0 = info["end0"]
        strand = info["strand"]

        # 1-based coords for VCF lookup
        start1 = start0 + 1
        end1 = end0 + 1  # end0 is 0-based inclusive → end1 = end0+1

        # Load SNPs for this region
        snps = load_vcf_snps(vcf_path, chrom, start1, end1, sample, iupac)

        # The FNA sequence is already on the gene's strand (MetaEuk outputs
        # revcomp for minus-strand genes). We need to apply SNPs in genomic
        # (+strand) coordinates, so we work on the plus-strand sequence.

        original_seq = str(rec.seq).upper()

        if strand == "-":
            # Convert to plus-strand for SNP application
            plus_seq = str(Seq(original_seq).reverse_complement())
            plus_seq = apply_snps_to_sequence(plus_seq, snps, start1)
            # Convert back to minus-strand
            new_seq = str(Seq(plus_seq).reverse_complement())
        else:
            new_seq = apply_snps_to_sequence(original_seq, snps, start1)

        # Build new FNA record
        rec.seq = Seq(new_seq)
        new_fna_records.append(rec)

        # Translate to protein
        protein_seq = translate_sequence(new_seq)
        from Bio.SeqRecord import SeqRecord

        faa_rec = SeqRecord(
            Seq(protein_seq),
            id=rec.id,
            description=rec.description.split(None, 1)[1] if " " in rec.description else "",
        )
        new_faa_records.append(faa_rec)

    # Write outputs
    out_fna = output_dir / f"{busco_id}.fna"
    out_faa = output_dir / f"{busco_id}.faa"
    out_gff = output_dir / f"{busco_id}.gff"

    SeqIO.write(new_fna_records, out_fna, "fasta")
    SeqIO.write(new_faa_records, out_faa, "fasta")
    shutil.copy2(gff_file, out_gff)


def main():
    parser = argparse.ArgumentParser(description="Apply VCF SNPs to BUSCO MetaEuk sequences.")
    parser.add_argument("--busco-dir", required=True, help="Directory with .fna/.faa/.gff BUSCO files")
    parser.add_argument("--ref", required=True, help="Reference FASTA (must be faidx-indexed)")
    parser.add_argument("--vcf", required=True, help="VCF file (bgzipped+tabix preferred)")
    parser.add_argument("--sample", default=None, help="Sample name in VCF (required if multi-sample)")
    parser.add_argument("--output-dir", required=True, help="Output directory for modified BUSCO files")
    parser.add_argument("--iupac", action="store_true", help="Encode heterozygous SNPs as IUPAC ambiguity codes (like GATK --use-iupac-sample)")
    args = parser.parse_args()

    busco_dir = Path(args.busco_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect unique BUSCO IDs from .gff files
    busco_ids = sorted(p.stem for p in busco_dir.glob("*.gff"))

    if not busco_ids:
        sys.exit(f"ERROR: No .gff files found in {busco_dir}")

    print(f"Found {len(busco_ids)} BUSCO genes to process.")
    print(f"IUPAC mode: {'ON' if args.iupac else 'OFF'}")
    print(f"Output dir: {output_dir}\n")

    for i, busco_id in enumerate(busco_ids, 1):
        print(f"[{i}/{len(busco_ids)}] {busco_id}")
        try:
            process_busco_gene(
                busco_id=busco_id,
                busco_dir=busco_dir,
                output_dir=output_dir,
                vcf_path=args.vcf,
                ref_path=args.ref,
                sample=args.sample,
                iupac=args.iupac,
            )
        except Exception as e:
            print(f"  ERROR processing {busco_id}: {e}")
            continue

    print("\nDone.")


if __name__ == "__main__":
    main()
