"""Unit tests for the annotation engine."""
import pytest
from src.annotator import (
    annotate_gene_context,
    annotate_af_tier,
    annotate_consequence,
    annotate_record,
    annotate_batch,
)


class TestAnnotateGenContext:
    def test_known_gene(self):
        result = annotate_gene_context("22", 30100000)
        assert result["gene"] == "NF2"
        assert result["is_clinical_gene"] is True

    def test_intergenic(self):
        result = annotate_gene_context("22", 17000)
        assert result["gene"] == "intergenic"

    def test_none_position(self):
        result = annotate_gene_context("22", None)
        assert result["gene"] == "intergenic"

    def test_clinical_gene_flagged(self):
        result = annotate_gene_context("22", 41950000)
        assert result["is_clinical_gene"] is True
        assert result["gene"] == "SMAD4"


class TestAnnotateAfTier:
    def test_rare(self):
        assert annotate_af_tier(0.001) == "rare"

    def test_low_frequency(self):
        assert annotate_af_tier(0.02) == "low_frequency"

    def test_common(self):
        assert annotate_af_tier(0.1) == "common"

    def test_boundary_rare(self):
        assert annotate_af_tier(0.009) == "rare"

    def test_boundary_common(self):
        assert annotate_af_tier(0.05) == "common"

    def test_none_af(self):
        assert annotate_af_tier(None) == "unknown"

    def test_zero_af(self):
        assert annotate_af_tier(0.0) == "rare"


class TestAnnotateConsequence:
    def test_snp_returns_valid_consequence(self):
        result = annotate_consequence("A", "G", "SNP")
        valid = {
            "synonymous", "missense", "nonsense", "splice_region"
        }
        assert result in valid

    def test_frameshift_insertion(self):
        result = annotate_consequence("A", "AT", "insertion")
        assert result in {"frameshift", "inframe_insertion"}

    def test_frameshift_deletion(self):
        result = annotate_consequence("ATG", "A", "deletion")
        assert result in {"frameshift", "inframe_deletion"}

    def test_inframe_insertion_divisible_by_3(self):
        result = annotate_consequence("A", "ATGC", "insertion")
        assert result in {"frameshift", "inframe_insertion"}

    def test_mnp(self):
        result = annotate_consequence("AT", "GC", "MNP")
        assert result == "MNP"


class TestAnnotateRecord:
    @pytest.fixture
    def valid_record(self):
        return {
            "chrom": "22",
            "pos": 16050075,
            "variant_id": "rs123",
            "ref": "A",
            "alt": "G",
            "qual": 100.0,
            "filter_status": "PASS",
            "af": 0.001,
            "dp": 8012,
            "mq": 60.0,
            "has_rsid": True,
            "is_multiallelic": False,
            "completeness_score": 1.0,
        }

    def test_snp_annotation(self, valid_record):
        result = annotate_record(valid_record)
        assert result["variant_class"] == "SNP"
        assert result["ts_tv"] in {"transition", "transversion"}
        assert result["af_tier"] == "rare"
        assert result["dbsnp_status"] == "known"

    def test_filter_outcome_pass(self, valid_record):
        result = annotate_record(valid_record)
        assert result["filter_outcome"] == "pass"

    def test_filter_outcome_filtered(self, valid_record):
        valid_record["filter_status"] = "LowQual"
        result = annotate_record(valid_record)
        assert result["filter_outcome"] == "filtered"

    def test_novel_variant(self, valid_record):
        valid_record["has_rsid"] = False
        result = annotate_record(valid_record)
        assert result["dbsnp_status"] == "novel"

    def test_clinical_gene_detection(self, valid_record):
        valid_record["pos"] = 30100000
        result = annotate_record(valid_record)
        assert result["is_clinical_gene"] is True
        assert result["gene"] == "NF2"


class TestAnnotateBatch:
    def test_all_records_annotated(self):
        records = [
            {
                "chrom": "22", "pos": 16050075,
                "variant_id": ".", "ref": "A", "alt": "G",
                "qual": 100.0, "filter_status": "PASS",
                "af": 0.001, "dp": 100, "mq": 60.0,
                "has_rsid": False, "is_multiallelic": False,
                "completeness_score": 1.0,
            }
        ] * 5
        annotated = annotate_batch(records)
        assert len(annotated) == 5
        assert all("variant_class" in r for r in annotated)
        assert all("consequence" in r for r in annotated)
        assert all("af_tier" in r for r in annotated)
