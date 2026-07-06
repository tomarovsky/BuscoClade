#!/usr/bin/env python3
__author__ = 'tomarovsky'

"""
apply_vcf_to_busco.py

Applies SNPs from a VCF file to BUSCO single-copy sequences (FNA/FAA/GFF)
produced by MetaEuk, without rebuilding the full alternate reference.

Requirements:
    pip install biopython pysam

Usage:
    python apply_vcf_to_busco.py \\
        --single_copy_busco_sequences /path/to/single_copy_busco_sequences \\
        --metaeuk_output /path/to/run_lineage/metaeuk_output \\
        --vcf sample.vcf.gz \\
        --sample SAMPLE_NAME \\
        --output_dir /path/to/output \\
        [--iupac]

Notes:
    - Only SNPs are handled; indels are ignored.
    - With --iupac, heterozygous SNPs are encoded as IUPAC ambiguity codes in FNA,
      and FAA uses ambiguous translation (X if ambiguous amino acid, otherwise
      the unambiguous amino acid if all codon variants agree).
    - Without --iupac, the ALT allele is used for homozygous alt and het calls
      (matching GATK FastaAlternateReferenceMaker default behaviour).
    - GFF files are copied unchanged (SNP-only: coordinates never shift).
    - The VCF must be bgzipped and tabix-indexed, or a plain .vcf (slower).
    - CDS structure is resolved from MetaEuk GFF files (rerun_results preferred,
      initial_results as fallback). The correct exons are identified by matching
      (Target_ID, chrom) and filtering to the gene's genomic range from the FNA
      header, which correctly handles multi-locus Target_IDs.
"""

import argparse
import shutil
import sys
from pathlib import Path
from typing import Optional

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from busco_reconstruct_common import (
    apply_snps_exon_aware,
    bases_to_iupac,
    build_cds_lookup,
    get_cds_exons,
    parse_metaeuk_fna_header,
    translate_sequence,
)


# ---------------------------------------------------------------------------
# VCF / SNP application
# ---------------------------------------------------------------------------

def resolve_vcf_sample(vcf: pysam.VariantFile, sample: Optional[str]) -> Optional[str]:
    """
    Resolve which sample to read from an open VCF, once, up front.
    Returns the sample name, or None for a sites-only VCF (use REF/ALT[0]).
    Exits with a clear error on an unknown or ambiguous sample.
    """
    samples = list(vcf.header.samples)
    if sample:
        if sample not in samples:
            sys.exit(f"ERROR: Sample '{sample}' not found in VCF. Available: {samples}")
        return sample
    if len(samples) == 1:
        return samples[0]
    if len(samples) > 1:
        sys.exit("ERROR: VCF has multiple samples. Please specify --sample.")
    return None


def load_vcf_snps(vcf: pysam.VariantFile, chrom: str, start1: int, end1: int,
                  sample_name: Optional[str], iupac: bool) -> dict:
    """
    Load SNPs from an already-open VCF for a given region (1-based, inclusive).
    Returns dict {pos1 (int): base (str)}.
    """
    snps = {}
    try:
        region_iter = vcf.fetch(chrom, start1 - 1, end1)
    except ValueError:
        return snps

    for rec in region_iter:
        if len(rec.ref) != 1:
            continue
        alts = rec.alts
        if alts is None:
            continue
        if any(len(a) != 1 for a in alts):
            continue

        pos1 = rec.pos  # pysam uses 1-based pos

        if sample_name is not None:
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
                alt_alleles = [alts[i - 1].upper() for i in allele_indices if i > 0]
                if alt_alleles:
                    snps[pos1] = alt_alleles[0]
        else:
            if iupac:
                snps[pos1] = bases_to_iupac({rec.ref.upper(), alts[0].upper()})
            else:
                snps[pos1] = alts[0].upper()

    return snps


# ---------------------------------------------------------------------------
# Per-gene processing
# ---------------------------------------------------------------------------

def process_busco_gene(busco_id: str, single_copy_busco_sequences: Path, output_dir: Path,
                       vcf: pysam.VariantFile, sample_name: Optional[str], iupac: bool,
                       rerun_index: dict, initial_index: dict):
    """Process one BUSCO gene: update FNA and FAA; copy GFF unchanged."""

    fna_file = single_copy_busco_sequences / f"{busco_id}.fna"
    faa_file = single_copy_busco_sequences / f"{busco_id}.faa"
    gff_file = single_copy_busco_sequences / f"{busco_id}.gff"

    records = list(SeqIO.parse(fna_file, "fasta"))
    if not records:
        print(f"  WARNING: empty FNA for {busco_id}, skipping.")
        return

    new_fna_records = []
    new_faa_records = []

    for rec in records:
        info = parse_metaeuk_fna_header(rec.description)
        tid    = info["full_id"]
        chrom  = info["chrom"]
        strand = info["strand"]

        # Header stores 0-based coords; convert to 1-based for GFF/VCF.
        gene_start1 = info["start0"] + 1
        gene_end1   = info["end0"]   + 1  # end0 is 0-based inclusive

        # Resolve exon structure from MetaEuk GFF.
        try:
            exons, source = get_cds_exons(
                tid, chrom, gene_start1, gene_end1,
                rerun_index, initial_index
            )
        except ValueError as e:
            raise ValueError(str(e))

        # Load SNPs per exon (one VCF query per exon; the VCF is opened once in main).
        snps_by_region = {}
        for exon_idx, (exon_start1, exon_end1) in enumerate(exons):
            snps = load_vcf_snps(vcf, chrom, exon_start1, exon_end1,
                                  sample_name, iupac)
            if snps:
                snps_by_region[exon_idx] = snps

        original_seq = str(rec.seq).upper()

        if snps_by_region:
            new_seq = apply_snps_exon_aware(
                original_seq, strand, exons, snps_by_region
            )
        else:
            new_seq = original_seq

        rec.seq = Seq(new_seq)
        new_fna_records.append(rec)

        protein_seq = translate_sequence(new_seq)
        faa_rec = SeqRecord(
            Seq(protein_seq),
            id=rec.id,
            description=rec.description.split(None, 1)[1] if " " in rec.description else "",
        )
        new_faa_records.append(faa_rec)

    SeqIO.write(new_fna_records, output_dir / f"{busco_id}.fna", "fasta")
    SeqIO.write(new_faa_records, output_dir / f"{busco_id}.faa", "fasta")
    shutil.copy2(gff_file, output_dir / f"{busco_id}.gff")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Apply VCF SNPs to BUSCO MetaEuk sequences (exon-aware)."
    )
    parser.add_argument(
        "--single_copy_busco_sequences", required=True,
        help="Directory with .fna/.faa/.gff BUSCO files "
             "(single_copy_busco_sequences)"
    )
    parser.add_argument(
        "--metaeuk_output", required=True,
        help="MetaEuk output directory containing initial_results/ and/or "
             "rerun_results/ subdirectories (run_<lineage>/metaeuk_output)"
    )
    parser.add_argument(
        "--vcf", required=True,
        help="VCF file (bgzipped + tabix-indexed preferred)"
    )
    parser.add_argument(
        "--sample", default=None,
        help="Sample name in VCF (required if multi-sample)"
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Output directory for modified BUSCO files"
    )
    parser.add_argument(
        "--iupac", action="store_true",
        help="Encode heterozygous SNPs as IUPAC ambiguity codes "
             "(like GATK --use-iupac-sample)"
    )
    args = parser.parse_args()

    single_copy_busco_sequences   = Path(args.single_copy_busco_sequences)
    metaeuk_output = Path(args.metaeuk_output)
    output_dir  = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load CDS structure from MetaEuk GFFs once, upfront.
    rerun_index, initial_index = build_cds_lookup(metaeuk_output)

    # Open the VCF once and resolve the sample once — reused for every exon query.
    try:
        vcf = pysam.VariantFile(args.vcf)
    except Exception as e:
        sys.exit(f"ERROR: Cannot open VCF {args.vcf}: {e}")
    sample_name = resolve_vcf_sample(vcf, args.sample)

    busco_ids = sorted(p.stem for p in single_copy_busco_sequences.glob("*.gff"))
    if not busco_ids:
        sys.exit(f"ERROR: No .gff files found in {single_copy_busco_sequences}")

    print(f"Found {len(busco_ids)} BUSCO genes to process.")
    print(f"IUPAC mode: {'ON' if args.iupac else 'OFF'}")
    print(f"Output dir: {output_dir}\n")

    errors = 0
    for i, busco_id in enumerate(busco_ids, 1):
        print(f"[{i}/{len(busco_ids)}] {busco_id}")
        try:
            process_busco_gene(
                busco_id=busco_id,
                single_copy_busco_sequences=single_copy_busco_sequences,
                output_dir=output_dir,
                vcf=vcf,
                sample_name=sample_name,
                iupac=args.iupac,
                rerun_index=rerun_index,
                initial_index=initial_index,
            )
        except Exception as e:
            print(f"  ERROR: {e}")
            errors += 1
            continue

    vcf.close()
    print(f"\nDone. {len(busco_ids) - errors} succeeded, {errors} failed.")


if __name__ == "__main__":
    main()
