#!/usr/bin/env python3
"""
add_altref_to_alignment.py

Insert AltRef sequences into a raw multiple-sequence alignment using
gap positions from the corresponding reference sequence.

Background
----------
AltRef sequences are produced by apply_vcf_to_busco.py, which applies only
SNPs from a VCF to the reference gene sequences.  Because only substitutions
are applied (no indels), every AltRef gene has exactly the same ungapped
length as the corresponding reference gene.  This means we can place gaps
into an AltRef sequence by simply copying the gap pattern of its reference
in the finished alignment — no re-alignment is needed.

Workflow
--------
For each reference species present in the raw alignment:
  1. Read the aligned reference sequence and record the positions of every
     gap character ('-').
  2. For each AltRef species associated with that reference:
       a. Load the unaligned gene from
          <busco_dir>/<altref_species>/busco_sequences/
          single_copy_busco_sequences/<gene_id>.fna
       b. Verify that its ungapped length equals the ungapped length of the
          reference in the alignment.  Warn and skip if they differ (this
          should not happen under normal pipeline conditions).
       c. Insert gaps at the same positions as the reference → the resulting
          aligned AltRef sequence has the same total length as every other
          sequence in the alignment.
       d. Append it to the alignment.
  3. If --keep_refs is NOT given, remove all reference sequences from the
     final alignment (they were needed only as gap templates).
  4. Write the result to --output.

Usage
-----
python add_altref_to_alignment.py \
    --raw_aln   alignments/raw/1234at7742.merged.fna \
    --busco_dir results/busco \
    --ref_to_altrefs '{"ref1": ["A1.ref1.AltRef", "A2.ref1.AltRef"], \
                       "ref2": ["A1.ref2.AltRef", "A3.ref2.AltRef"]}' \
    --gene_id   1234at7742.merged \
    --keep_refs            # omit to drop refs from output
    --output    alignments/fna/1234at7742.merged.fna
"""

import argparse
import json
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> dict[str, str]:
    """
    Parse a FASTA file and return an ordered dict {header: sequence}.

    The header is everything after '>' up to the first whitespace so that
    BUSCO sequence headers (which often contain extra annotation) are matched
    by their bare ID.
    """
    sequences: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(parts)
                current_id = line[1:].split()[0]
                parts = []
            elif current_id is not None:
                parts.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(parts)

    return sequences


def write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Write sequences to a FASTA file, 60 characters per line."""
    with open(path, "w") as fh:
        for header, seq in sequences.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def insert_gaps(ungapped_seq: str, gap_positions: list[int], aln_length: int) -> str:
    """
    Insert '-' characters into *ungapped_seq* at *gap_positions* (0-based
    indices in the final aligned sequence) so that the result has length
    *aln_length*.

    This mirrors exactly how the reference looks after alignment: the
    non-gap characters of the reference occupy all positions that are NOT in
    gap_positions, in order.

    Parameters
    ----------
    ungapped_seq : str
        The raw (unaligned) sequence — must have length == aln_length - len(gap_positions).
    gap_positions : list[int]
        Sorted list of 0-based positions where '-' should be placed in the
        output string.
    aln_length : int
        Expected total length of the returned string.

    Returns
    -------
    str
        Aligned sequence of length *aln_length*.
    """
    result = list("-" * aln_length)
    gap_set = set(gap_positions)
    seq_iter = iter(ungapped_seq)
    for i in range(aln_length):
        if i not in gap_set:
            result[i] = next(seq_iter)
    return "".join(result)


def get_gap_positions(aligned_seq: str) -> list[int]:
    """Return sorted list of 0-based indices where *aligned_seq* contains '-'."""
    return [i for i, ch in enumerate(aligned_seq) if ch == "-"]


def find_altref_gene_file(busco_dir: Path, altref_species: str, gene_id: str) -> Path | None:
    """
    Locate the unaligned gene FASTA for *altref_species* / *gene_id*.

    BUSCO stores single-copy sequences at:
        <busco_dir>/<species>/busco_sequences/single_copy_busco_sequences/<gene_id>.fna

    The *gene_id* wildcard in Snakemake has the form "<busco_id>.merged", so
    we strip ".merged" when looking for the file (BUSCO files are named by
    bare busco_id).
    """
    bare_id = gene_id.removesuffix(".merged")
    candidate = (
        busco_dir
        / altref_species
        / "busco_sequences"
        / "single_copy_busco_sequences"
        / f"{bare_id}.fna"
    )
    return candidate if candidate.is_file() else None


def add_altref_sequences(
    raw_aln: dict[str, str],
    busco_dir: Path,
    ref_to_altrefs: dict[str, list[str]],
    gene_id: str,
) -> dict[str, str]:
    """
    Return a copy of *raw_aln* extended with all available AltRef sequences.

    For each (ref_prefix, altref_species_list) pair the function:
      - Locates the aligned reference in *raw_aln* by scanning for a key
        that matches *ref_prefix* (exact or as the first token of a longer
        header, e.g. "ref1" matches ">ref1 some annotation").
      - Derives gap positions from that aligned reference sequence.
      - Loads each AltRef's unaligned gene, validates its length, inserts gaps,
        and appends it to the output alignment.

    Missing AltRef gene files are logged as warnings and skipped — the
    alignment is still written so the pipeline does not fail hard.
    """
    aln_length = len(next(iter(raw_aln.values())))  # all seqs must be the same length
    extended = dict(raw_aln)  # will be mutated below

    for ref_prefix, altref_species_list in ref_to_altrefs.items():
        # --- Find the reference sequence in the alignment ---
        ref_key = next(
            (k for k in raw_aln if k == ref_prefix or k.startswith(ref_prefix + " ")),
            None,
        )
        if ref_key is None:
            log.warning(
                "Reference '%s' not found in raw alignment for gene '%s'. "
                "Skipping AltRef species: %s.",
                ref_prefix,
                gene_id,
                altref_species_list,
            )
            continue

        aligned_ref = raw_aln[ref_key]
        gap_positions = get_gap_positions(aligned_ref)
        ungapped_ref_len = aln_length - len(gap_positions)

        log.info(
            "Reference '%s': aligned length %d, gaps %d, ungapped %d.",
            ref_prefix,
            aln_length,
            len(gap_positions),
            ungapped_ref_len,
        )

        # --- Process each AltRef species ---
        for altref_species in altref_species_list:
            gene_file = find_altref_gene_file(busco_dir, altref_species, gene_id)

            if gene_file is None:
                log.warning(
                    "Gene file not found for AltRef species '%s', gene '%s'. "
                    "This gene may be missing or fragmentary in the VCF-reconstructed set. "
                    "Skipping.",
                    altref_species,
                    gene_id,
                )
                continue

            altref_seqs = read_fasta(gene_file)
            if not altref_seqs:
                log.warning("Empty FASTA for AltRef species '%s', gene '%s'. Skipping.", altref_species, gene_id)
                continue

            # take the first (and normally only) sequence
            altref_seq = next(iter(altref_seqs.values())).replace("-", "")

            if len(altref_seq) != ungapped_ref_len:
                log.warning(
                    "Length mismatch for AltRef species '%s', gene '%s': "
                    "expected ungapped length %d (from ref '%s'), got %d. "
                    "Skipping — this species will be absent from the alignment.",
                    altref_species,
                    gene_id,
                    ungapped_ref_len,
                    ref_prefix,
                    len(altref_seq),
                )
                continue

            aligned_altref = insert_gaps(altref_seq, gap_positions, aln_length)

            # Match the case convention of the alignment (determined from the
            # reference sequence, which comes from the same merged_sequences_raw).
            non_gap_ref = [ch for ch in aligned_ref if ch != "-"]
            if non_gap_ref:
                upper_count = sum(1 for ch in non_gap_ref if ch.isupper())
                if upper_count >= len(non_gap_ref) / 2:
                    aligned_altref = aligned_altref.upper()
                else:
                    aligned_altref = aligned_altref.lower()

            extended[altref_species] = aligned_altref
            log.info("Inserted AltRef species '%s' into alignment.", altref_species)

    return extended


def remove_ref_sequences(
    alignment: dict[str, str],
    ref_prefixes: list[str],
) -> dict[str, str]:
    """
    Remove reference sequences from *alignment*.

    A sequence is considered a reference if its header equals one of
    *ref_prefixes* or starts with "<ref_prefix> " (space-separated annotation).
    """
    filtered = {}
    for key, seq in alignment.items():
        is_ref = any(
            key == ref or key.startswith(ref + " ")
            for ref in ref_prefixes
        )
        if not is_ref:
            filtered[key] = seq
        else:
            log.info("Removed reference sequence '%s' from final alignment.", key)
    return filtered


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--raw_aln",
        required=True,
        type=Path,
        metavar="FASTA",
        help="Raw alignment FASTA (refs + non-AltRef species, without AltRef).",
    )
    p.add_argument(
        "--busco_dir",
        required=True,
        type=Path,
        metavar="DIR",
        help="Root BUSCO output directory (<busco_dir>/<species>/busco_sequences/…).",
    )
    p.add_argument(
        "--ref_to_altrefs",
        required=True,
        type=str,
        metavar="JSON",
        help=(
            'JSON mapping ref_prefix → list of AltRef species names, e.g. '
            '\'{"ref1": ["A1.ref1.AltRef", "A2.ref1.AltRef"]}\''
        ),
    )
    p.add_argument(
        "--gene_id",
        required=True,
        type=str,
        metavar="ID",
        help="BUSCO gene wildcard value (e.g. '1234at7742.merged').",
    )
    p.add_argument(
        "--keep_refs",
        action="store_true",
        default=False,
        help=(
            "Keep reference sequences in the output alignment. "
            "If not set, ref sequences are removed after AltRef insertion "
            "(use when vcf_reconstruct_ref_as_species: False)."
        ),
    )
    p.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="FASTA",
        help="Output aligned FASTA path.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    # --- Parse ref → altrefs mapping ---
    try:
        ref_to_altrefs: dict[str, list[str]] = json.loads(args.ref_to_altrefs)
    except json.JSONDecodeError as exc:
        log.error("Failed to parse --ref_to_altrefs JSON: %s", exc)
        sys.exit(1)

    # --- Load raw alignment ---
    if not args.raw_aln.is_file():
        log.error("Raw alignment file not found: %s", args.raw_aln)
        sys.exit(1)

    raw_aln = read_fasta(args.raw_aln)
    if not raw_aln:
        log.error("Raw alignment is empty: %s", args.raw_aln)
        sys.exit(1)

    log.info(
        "Loaded raw alignment: %d sequences, length %d.",
        len(raw_aln),
        len(next(iter(raw_aln.values()))),
    )

    # --- Insert AltRef sequences ---
    extended_aln = add_altref_sequences(
        raw_aln=raw_aln,
        busco_dir=args.busco_dir,
        ref_to_altrefs=ref_to_altrefs,
        gene_id=args.gene_id,
    )

    # --- Optionally remove reference sequences ---
    if not args.keep_refs:
        log.info("--keep_refs not set: removing reference sequences from output.")
        extended_aln = remove_ref_sequences(
            alignment=extended_aln,
            ref_prefixes=list(ref_to_altrefs.keys()),
        )

    # --- Validate non-empty output ---
    if not extended_aln:
        log.error(
            "Output alignment is empty after processing gene '%s'. "
            "Check that ref sequences are present in the raw alignment and "
            "that AltRef gene files exist.",
            args.gene_id,
        )
        sys.exit(1)

    # --- Write output ---
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(extended_aln, args.output)
    log.info(
        "Written %d sequences to %s.",
        len(extended_aln),
        args.output,
    )


if __name__ == "__main__":
    main()
