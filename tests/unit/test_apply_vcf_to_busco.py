"""Unit tests for workflow/scripts/apply_vcf_to_busco.py.

Covers the VCF-facing logic: sample resolution and the SNP-loading rules
(SNP-only filtering, genotype-aware allele selection, and IUPAC het encoding).
"""
import pytest

pysam = pytest.importorskip("pysam")

import apply_vcf_to_busco as avb


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
)

VCF_RECORDS = (
    "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t1/1\n"   # hom-alt SNP
    "chr1\t15\t.\tC\tT,A\t.\t.\t.\tGT\t0/1\n"  # het, multi-allelic
    "chr1\t20\t.\tAT\tA\t.\t.\t.\tGT\t1/1\n"   # indel (ref len 2) -> ignored
    "chr1\t25\t.\tG\tC\t.\t.\t.\tGT\t0/0\n"    # hom-ref
)


def _make_vcf(tmp_path, header, records, name="s.vcf"):
    """Write, bgzip and tabix-index a VCF; return an open VariantFile."""
    path = tmp_path / name
    path.write_text(header + records)
    gz = pysam.tabix_index(str(path), preset="vcf", force=True)  # -> <path>.gz + .tbi
    return pysam.VariantFile(gz)


# ---------------------------------------------------------------------------
# Sample resolution
# ---------------------------------------------------------------------------

def test_resolve_vcf_sample_single_sample_autodetected(tmp_path):
    vcf = _make_vcf(tmp_path, VCF_HEADER, VCF_RECORDS)
    assert avb.resolve_vcf_sample(vcf, None) == "SAMPLE1"
    assert avb.resolve_vcf_sample(vcf, "SAMPLE1") == "SAMPLE1"


def test_resolve_vcf_sample_unknown_sample_exits(tmp_path):
    vcf = _make_vcf(tmp_path, VCF_HEADER, VCF_RECORDS)
    with pytest.raises(SystemExit):
        avb.resolve_vcf_sample(vcf, "NOPE")


def test_resolve_vcf_sample_sites_only_returns_none(tmp_path):
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    vcf = _make_vcf(tmp_path, header, "chr1\t10\t.\tA\tG\t.\t.\t.\n", name="sites.vcf")
    assert avb.resolve_vcf_sample(vcf, None) is None


# ---------------------------------------------------------------------------
# SNP loading
# ---------------------------------------------------------------------------

def test_load_vcf_snps_non_iupac_uses_alt_allele(tmp_path):
    vcf = _make_vcf(tmp_path, VCF_HEADER, VCF_RECORDS)
    snps = avb.load_vcf_snps(vcf, "chr1", 1, 30, "SAMPLE1", iupac=False)
    # 10: hom-alt -> G; 15: het picks first ALT -> T; 20 indel skipped; 25 hom-ref dropped.
    assert snps == {10: "G", 15: "T"}


def test_load_vcf_snps_iupac_encodes_hets(tmp_path):
    vcf = _make_vcf(tmp_path, VCF_HEADER, VCF_RECORDS)
    snps = avb.load_vcf_snps(vcf, "chr1", 1, 30, "SAMPLE1", iupac=True)
    # 10: {G} -> G; 15: {C,T} -> Y; 20 indel skipped; 25: {G} -> G.
    assert snps == {10: "G", 15: "Y", 25: "G"}


def test_load_vcf_snps_sites_only(tmp_path):
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    records = "chr1\t10\t.\tA\tG\t.\t.\t.\n"
    vcf = _make_vcf(tmp_path, header, records, name="sites.vcf")
    assert avb.load_vcf_snps(vcf, "chr1", 1, 30, None, iupac=False) == {10: "G"}
    vcf2 = _make_vcf(tmp_path, header, records, name="sites2.vcf")
    assert avb.load_vcf_snps(vcf2, "chr1", 1, 30, None, iupac=True) == {10: "R"}
