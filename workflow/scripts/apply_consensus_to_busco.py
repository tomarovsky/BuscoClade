#!/usr/bin/env python3
__author__ = 'tomarovsky'

"""
apply_consensus_to_busco.py

Reconstructs BUSCO single-copy sequences (FNA/FAA/GFF) for an already
reconstructed ("consensus") genome, reusing the reference genome's MetaEuk
BUSCO output — without re-running BUSCO on the consensus genome.

This is the FASTA analogue of apply_vcf_to_busco.py. Instead of applying SNPs
from a VCF to the reference BUSCO genes, it slices the corresponding CDS exon
regions directly out of the consensus genome (reverse-complementing minus-strand
exons) and concatenates them in transcript order.

Requirements:
    pip install biopython pysam

Usage:
    python apply_consensus_to_busco.py \\
        --single_copy_busco_sequences /path/to/single_copy_busco_sequences \\
        --metaeuk_output /path/to/run_lineage/metaeuk_output \\
        --consensus consensus_genome.fasta \\
        --output_dir /path/to/output

Notes:
    - The consensus genome MUST be coordinate-identical to the reference:
      identical chromosome names, identical coordinates, substitutions only
      (no indels, no renaming). This is what guarantees every reconstructed
      gene keeps the reference's ungapped length — the same invariant the
      VCF route relies on.
    - CDS exon structure is resolved from the reference MetaEuk GFF files
      (rerun_results preferred, initial_results as fallback), exactly as in
      apply_vcf_to_busco.py.
    - GFF files are copied unchanged (SNP-only: coordinates never shift).
    - The consensus FASTA must be faidx-indexable (plain or bgzipped); a .fai
      index is created next to it if absent.
"""

import argparse
import shutil
import sys
from pathlib import Path

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from busco_reconstruct_common import (
    build_cds_lookup,
    get_cds_exons,
    parse_metaeuk_fna_header,
    reverse_complement,
    translate_sequence,
)


# ---------------------------------------------------------------------------
# Up-front coordinate-identity validation
# ---------------------------------------------------------------------------

def collect_required_contigs(single_copy_busco_sequences: Path) -> dict:
    """Scan every reference BUSCO FNA header and return {chrom: max_end1} — the
    contigs (and highest coordinate) the reconstruction will read from the
    consensus genome."""
    required = {}
    for fna in single_copy_busco_sequences.glob("*.fna"):
        with open(fna) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                try:
                    info = parse_metaeuk_fna_header(line)
                except ValueError:
                    continue
                end1 = info["end0"] + 1
                cur = required.get(info["chrom"])
                if cur is None or end1 > cur:
                    required[info["chrom"]] = end1
    return required


def validate_consensus_contigs(required: dict, fasta: pysam.FastaFile,
                               consensus_path: str) -> None:
    """Fail fast (once) if the consensus genome cannot be coordinate-identical to
    the reference: a contig used by BUSCO genes is missing or too short. This turns
    a wrong/renamed/truncated genome into one clear error instead of thousands of
    per-gene warnings."""
    refs = set(fasta.references)
    problems = []
    for chrom, max_end1 in sorted(required.items()):
        if chrom not in refs:
            problems.append(f"  missing contig '{chrom}' (BUSCO genes reference it up to pos {max_end1})")
            continue
        length = fasta.get_reference_length(chrom)
        if length < max_end1:
            problems.append(f"  contig '{chrom}' too short: length {length} < required {max_end1}")

    if problems:
        sys.exit(
            f"ERROR: consensus genome {consensus_path} is not coordinate-identical to "
            f"the reference (same contig names and coordinates, substitutions only, are "
            f"required):\n" + "\n".join(problems)
        )


# ---------------------------------------------------------------------------
# Consensus genome extraction
# ---------------------------------------------------------------------------

def build_cds_from_consensus(fasta: pysam.FastaFile, chrom: str, strand: str,
                             exons: list) -> str:
    """
    Reconstruct a spliced CDS from the consensus genome.

    exons are 1-based inclusive intervals in MetaEuk transcript order (for
    minus-strand genes, descending genomic order). Each exon region is fetched
    from the consensus genome and, on the minus strand, reverse-complemented —
    mirroring how MetaEuk builds the reference FNA — then concatenated.

    Raises ValueError if the chromosome is absent from the consensus genome.
    """
    parts = []
    for exon_start1, exon_end1 in exons:
        try:
            seg = fasta.fetch(chrom, exon_start1 - 1, exon_end1).upper()
        except (KeyError, ValueError) as e:
            raise ValueError(
                f"Cannot fetch {chrom}:{exon_start1}-{exon_end1} from consensus "
                f"genome (is it coordinate-identical to the reference?): {e}"
            )
        if strand == "-":
            seg = reverse_complement(seg)
        parts.append(seg)
    return "".join(parts)


# ---------------------------------------------------------------------------
# Per-gene processing
# ---------------------------------------------------------------------------

def process_busco_gene(busco_id: str, single_copy_busco_sequences: Path,
                       output_dir: Path, fasta: pysam.FastaFile,
                       rerun_index: dict, initial_index: dict):
    """Process one BUSCO gene: rebuild FNA and FAA from the consensus genome;
    copy GFF unchanged."""

    fna_file = single_copy_busco_sequences / f"{busco_id}.fna"
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

        # Header stores 0-based coords; convert to 1-based for GFF.
        gene_start1 = info["start0"] + 1
        gene_end1   = info["end0"]   + 1  # end0 is 0-based inclusive

        # Resolve exon structure from the reference MetaEuk GFF.
        exons, source = get_cds_exons(
            tid, chrom, gene_start1, gene_end1,
            rerun_index, initial_index
        )

        new_seq = build_cds_from_consensus(fasta, chrom, strand, exons)

        ref_len = len(str(rec.seq))
        if len(new_seq) != ref_len:
            # Coordinate-identity is assumed; a mismatch means the consensus
            # genome is not a pure-substitution copy of the reference. The
            # gap-aware insertion step would drop such a gene anyway.
            print(
                f"  WARNING: {busco_id} length mismatch vs reference "
                f"(reference {ref_len}, consensus {len(new_seq)}); "
                f"consensus genome may not be coordinate-identical."
            )

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
        description="Reconstruct BUSCO MetaEuk sequences from a consensus genome "
                    "(exon-aware, coordinate-identical to the reference)."
    )
    parser.add_argument(
        "--single_copy_busco_sequences", required=True,
        help="Reference directory with .fna/.faa/.gff BUSCO files "
             "(single_copy_busco_sequences)"
    )
    parser.add_argument(
        "--metaeuk_output", required=True,
        help="Reference MetaEuk output directory containing initial_results/ "
             "and/or rerun_results/ subdirectories (run_<lineage>/metaeuk_output)"
    )
    parser.add_argument(
        "--consensus", required=True,
        help="Consensus genome FASTA (coordinate-identical to the reference)"
    )
    parser.add_argument(
        "--output_dir", required=True,
        help="Output directory for reconstructed BUSCO files"
    )
    args = parser.parse_args()

    single_copy_busco_sequences = Path(args.single_copy_busco_sequences)
    metaeuk_output = Path(args.metaeuk_output)
    output_dir  = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Ensure the consensus FASTA has a faidx index (created alongside it).
    if not Path(str(args.consensus) + ".fai").exists():
        pysam.faidx(str(args.consensus))
    fasta = pysam.FastaFile(str(args.consensus))

    # Load CDS structure from the reference MetaEuk GFFs once, upfront.
    rerun_index, initial_index = build_cds_lookup(metaeuk_output)

    busco_ids = sorted(p.stem for p in single_copy_busco_sequences.glob("*.gff"))
    if not busco_ids:
        sys.exit(f"ERROR: No .gff files found in {single_copy_busco_sequences}")

    # Fail fast if the consensus genome is clearly not coordinate-identical.
    validate_consensus_contigs(
        collect_required_contigs(single_copy_busco_sequences), fasta, args.consensus
    )

    print(f"Found {len(busco_ids)} BUSCO genes to process.")
    print(f"Consensus genome: {args.consensus}")
    print(f"Output dir: {output_dir}\n")

    errors = 0
    for i, busco_id in enumerate(busco_ids, 1):
        print(f"[{i}/{len(busco_ids)}] {busco_id}")
        try:
            process_busco_gene(
                busco_id=busco_id,
                single_copy_busco_sequences=single_copy_busco_sequences,
                output_dir=output_dir,
                fasta=fasta,
                rerun_index=rerun_index,
                initial_index=initial_index,
            )
        except Exception as e:
            print(f"  ERROR: {e}")
            errors += 1
            continue

    fasta.close()
    print(f"\nDone. {len(busco_ids) - errors} succeeded, {errors} failed.")


if __name__ == "__main__":
    main()
