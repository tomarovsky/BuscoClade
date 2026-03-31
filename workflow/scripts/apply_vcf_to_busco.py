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
    - Only SNPs are handled; indels are ignored (as in the original pipeline).
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
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pysam
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
        total_len = sum(e - s + 1 for s, e in filtered)
        return filtered, source_name

    raise ValueError(
        f"No CDS exons found for {tid} on {chrom}:{gene_start1}-{gene_end1} "
        f"in either rerun or initial GFF."
    )


# ---------------------------------------------------------------------------
# VCF / SNP application
# ---------------------------------------------------------------------------

def load_vcf_snps(vcf_path: str, chrom: str, start1: int, end1: int,
                  sample: Optional[str], iupac: bool) -> dict:
    """
    Load SNPs from VCF for a given region (1-based, inclusive).
    Returns dict {pos1 (int): base (str)}.
    """
    snps = {}
    try:
        vcf = pysam.VariantFile(vcf_path)
    except Exception as e:
        sys.exit(f"ERROR: Cannot open VCF {vcf_path}: {e}")

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
                alt_alleles = [alts[i - 1].upper() for i in allele_indices if i > 0]
                if alt_alleles:
                    snps[pos1] = alt_alleles[0]
        else:
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
# Core SNP application logic
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
            COMP = str.maketrans("ACGTacgtRYMKWSBDHVNrykmwsbdhvn",
                                  "TGCAtgcaYRKMWSVHDBNyrkmwsvhdbn")
            for pos1, base in snps.items():
                fna_idx = fna_offset + (exon_end1 - pos1)
                if 0 <= fna_idx < fna_offset + exon_len:
                    seq_list[fna_idx] = base.translate(COMP)

        fna_offset += exon_len

    return "".join(seq_list)


# ---------------------------------------------------------------------------
# Per-gene processing
# ---------------------------------------------------------------------------

def process_busco_gene(busco_id: str, single_copy_busco_sequences: Path, output_dir: Path,
                       vcf_path: str, sample: Optional[str], iupac: bool,
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

        # Load SNPs per exon (one VCF query per exon).
        snps_by_region = {}
        for exon_idx, (exon_start1, exon_end1) in enumerate(exons):
            snps = load_vcf_snps(vcf_path, chrom, exon_start1, exon_end1,
                                  sample, iupac)
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
                vcf_path=args.vcf,
                sample=args.sample,
                iupac=args.iupac,
                rerun_index=rerun_index,
                initial_index=initial_index,
            )
        except Exception as e:
            print(f"  ERROR: {e}")
            errors += 1
            continue

    print(f"\nDone. {len(busco_ids) - errors} succeeded, {errors} failed.")


if __name__ == "__main__":
    main()
