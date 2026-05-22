"""Unit tests for the VCF parser."""
import pytest
from src.parser import (
    parse_vcf_line,
    parse_info_field,
    classify_variant_type,
    is_transition,
    parse_vcf_data,
)


class TestParseInfoField:
    def test_standard_info(self):
        result = parse_info_field("AF=0.001;DP=100;MQ=60")
        assert result["AF"] == "0.001"
        assert result["DP"] == "100"
        assert result["MQ"] == "60"

    def test_flag_field(self):
        result = parse_info_field("INDEL;AF=0.01")
        assert result["INDEL"] is True
        assert result["AF"] == "0.01"

    def test_empty_info(self):
        result = parse_info_field(".")
        assert result == {}

    def test_empty_string(self):
        result = parse_info_field("")
        assert result == {}

    def test_multiallelic_af(self):
        result = parse_info_field("AF=0.001,0.002;DP=100")
        assert result["AF"] == "0.001,0.002"


class TestParseVcfLine:
    def test_valid_snp(self):
        line = "22\t16050075\trs123\tA\tG\t100.0\tPASS\tAF=0.001;DP=8012"
        result = parse_vcf_line(line)
        assert result["parse_error"] is False
        assert result["chrom"] == "22"
        assert result["pos"] == 16050075
        assert result["ref"] == "A"
        assert result["alt"] == "G"
        assert result["qual"] == 100.0
        assert result["filter_status"] == "PASS"
        assert result["af"] == pytest.approx(0.001)
        assert result["dp"] == 8012
        assert result["has_rsid"] is True

    def test_chr_prefix_normalised(self):
        line = "chr22\t16050075\t.\tA\tG\t100.0\tPASS\tAF=0.001"
        result = parse_vcf_line(line)
        assert result["chrom"] == "22"

    def test_missing_rsid(self):
        line = "22\t16050075\t.\tA\tG\t100.0\tPASS\tAF=0.001"
        result = parse_vcf_line(line)
        assert result["has_rsid"] is False

    def test_too_few_fields(self):
        line = "22\t16050075\trs123"
        result = parse_vcf_line(line)
        assert result["parse_error"] is True

    def test_insertion(self):
        line = "22\t100\t.\tA\tAT\t200.0\tPASS\tAF=0.01"
        result = parse_vcf_line(line)
        assert result["ref"] == "A"
        assert result["alt"] == "AT"

    def test_multiallelic(self):
        line = "22\t100\t.\tA\tG,T\t200.0\tPASS\tAF=0.01"
        result = parse_vcf_line(line)
        assert result["is_multiallelic"] is True
        assert len(result["alt_alleles"]) == 2

    def test_missing_qual(self):
        line = "22\t100\t.\tA\tG\t.\tPASS\tAF=0.01"
        result = parse_vcf_line(line)
        assert result["qual"] is None
        assert result["parse_error"] is False


class TestClassifyVariantType:
    def test_snp(self):
        assert classify_variant_type("A", "G") == "SNP"

    def test_insertion(self):
        assert classify_variant_type("A", "AT") == "insertion"

    def test_deletion(self):
        assert classify_variant_type("AT", "A") == "deletion"

    def test_mnp(self):
        assert classify_variant_type("AT", "GC") == "MNP"

    def test_empty_ref(self):
        assert classify_variant_type("", "G") == "unknown"

    def test_complex(self):
        assert classify_variant_type("A", "A") == "SNP"


class TestIsTransition:
    def test_ag_transition(self):
        assert is_transition("A", "G") is True

    def test_ct_transition(self):
        assert is_transition("C", "T") is True

    def test_ga_transition(self):
        assert is_transition("G", "A") is True

    def test_tc_transition(self):
        assert is_transition("T", "C") is True

    def test_ac_transversion(self):
        assert is_transition("A", "C") is False

    def test_at_transversion(self):
        assert is_transition("A", "T") is False

    def test_gc_transversion(self):
        assert is_transition("G", "C") is False

    def test_case_insensitive(self):
        assert is_transition("a", "g") is True


class TestParseVcfData:
    def test_parses_all_lines(self):
        data = {
            "variant_lines": [
                "22\t100\t.\tA\tG\t100\tPASS\tAF=0.001",
                "22\t200\t.\tC\tT\t200\tPASS\tAF=0.002",
            ],
            "chromosome": "22",
            "reference_genome": "GRCh37",
            "dataset": "test",
            "source": "test",
        }
        records = parse_vcf_data(data)
        assert len(records) == 2

    def test_skips_empty_lines(self):
        data = {
            "variant_lines": [
                "22\t100\t.\tA\tG\t100\tPASS\tAF=0.001",
                "",
                "22\t200\t.\tC\tT\t200\tPASS\tAF=0.002",
            ],
            "chromosome": "22",
            "reference_genome": "GRCh37",
            "dataset": "test",
            "source": "test",
        }
        records = parse_vcf_data(data)
        assert len(records) == 2
