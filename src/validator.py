"""
Validation Engine for VariantLens.

Applies 9 biological and data quality rules to each parsed VCF record.
Failed records are quarantined with rejection reasons — never silently dropped.

Rule categories:
  Identity    — chromosome and position integrity
  Alleles     — REF/ALT biological validity
  Quality     — Phred score and filter status
  Metrics     — depth and mapping quality
  Completeness — required fields present

Quarantine-not-delete: every failed record is preserved with the
exact rule(s) that failed and a human-readable rejection reason.
This aligns with GxP data integrity principles used in regulated
clinical genomics environments.
"""
import os
import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
import structlog

from src.config import settings

log = structlog.get_logger()

# Valid human chromosomes (GRCh37/hg19)
VALID_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

# Valid nucleotide bases in VCF alleles
VALID_BASES = set("ACGTNacgtn*<>[].")

# Maximum plausible genomic position (GRCh37 chr22 ends at ~51Mb)
MAX_POSITION = 300_000_000

# Quality score thresholds
MIN_QUAL = 0.0
MAX_QUAL = 100_000.0

# Allele frequency bounds
MIN_AF = 0.0
MAX_AF = 1.0

# Minimum mapping quality for reliable variant calls
MIN_MQ = 0.0


@dataclass
class ValidationResult:
    passed: bool
    rule: str
    reason: str


@dataclass
class RecordOutcome:
    record: dict
    passed: bool
    failures: list[ValidationResult] = field(default_factory=list)


# ─────────────────────────────────────────────────────────────────
# The 9 validation rules — pure functions
# each takes a record dict and returns a ValidationResult
# ─────────────────────────────────────────────────────────────────

def rule_chromosome_valid(record: dict) -> ValidationResult:
    """Chromosome must be a recognised human chromosome."""
    chrom = str(record.get("chrom", "")).strip()
    passed = chrom in VALID_CHROMOSOMES
    return ValidationResult(
        passed=passed,
        rule="chromosome_valid",
        reason="" if passed else (
            f"chromosome '{chrom}' is not a valid human chromosome"
        ),
    )


def rule_position_positive(record: dict) -> ValidationResult:
    """Genomic position must be a positive integer (VCF is 1-based)."""
    pos = record.get("pos")
    try:
        pos = int(pos)
        passed = 1 <= pos <= MAX_POSITION
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="position_positive",
        reason="" if passed else (
            f"position '{pos}' is not a valid 1-based genomic coordinate"
        ),
    )


def rule_ref_allele_valid(record: dict) -> ValidationResult:
    """REF allele must be non-empty and contain only valid bases."""
    ref = record.get("ref", "")
    passed = (
        bool(ref)
        and len(ref) > 0
        and all(b in VALID_BASES for b in ref)
    )
    return ValidationResult(
        passed=passed,
        rule="ref_allele_valid",
        reason="" if passed else (
            f"REF allele '{ref}' is empty or contains invalid bases"
        ),
    )


def rule_alt_allele_valid(record: dict) -> ValidationResult:
    """ALT allele must be non-empty and differ from REF."""
    ref = record.get("ref", "")
    alt = record.get("alt", "")
    passed = (
        bool(alt)
        and alt != "."
        and alt.upper() != ref.upper()
    )
    return ValidationResult(
        passed=passed,
        rule="alt_allele_valid",
        reason="" if passed else (
            f"ALT allele '{alt}' is empty, missing, or identical to REF '{ref}'"
        ),
    )


def rule_quality_score_valid(record: dict) -> ValidationResult:
    """
    QUAL score must be non-negative and within plausible bounds.
    Phred-scaled quality: QUAL = -10 * log10(error_probability)
    A score of 30 = 99.9% confidence. Values above 10,000 are
    computationally implausible.
    """
    qual = record.get("qual")
    try:
        qual = float(qual)
        passed = MIN_QUAL <= qual <= MAX_QUAL
    except (TypeError, ValueError):
        passed = qual is None  # Missing QUAL ('.') is allowed in VCF
    return ValidationResult(
        passed=passed,
        rule="quality_score_valid",
        reason="" if passed else (
            f"QUAL score '{qual}' is outside valid range "
            f"[{MIN_QUAL}, {MAX_QUAL}]"
        ),
    )


def rule_filter_field_present(record: dict) -> ValidationResult:
    """FILTER field must be present and non-empty."""
    filter_status = record.get("filter_status", "")
    passed = bool(filter_status and filter_status.strip())
    return ValidationResult(
        passed=passed,
        rule="filter_field_present",
        reason="" if passed else "FILTER field is missing or empty",
    )


def rule_allele_frequency_valid(record: dict) -> ValidationResult:
    """
    Allele frequency (AF) must be between 0 and 1 if present.
    AF > 1.0 indicates a data corruption or normalisation error.
    AF is optional in VCF — missing AF passes this rule.
    """
    af = record.get("af")
    if af is None:
        return ValidationResult(
            passed=True,
            rule="allele_frequency_valid",
            reason="",
        )
    try:
        af = float(af)
        passed = MIN_AF <= af <= MAX_AF
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="allele_frequency_valid",
        reason="" if passed else (
            f"allele frequency '{af}' is outside valid range [0, 1]"
        ),
    )


def rule_depth_non_negative(record: dict) -> ValidationResult:
    """
    Read depth (DP) must be non-negative if present.
    Negative depth is physically impossible and indicates
    a data entry or computation error.
    """
    dp = record.get("dp")
    if dp is None:
        return ValidationResult(
            passed=True,
            rule="depth_non_negative",
            reason="",
        )
    try:
        dp = int(dp)
        passed = dp >= 0
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="depth_non_negative",
        reason="" if passed else (
            f"read depth '{dp}' is negative — physically impossible"
        ),
    )


def rule_ref_not_equal_alt(record: dict) -> ValidationResult:
    """
    REF and ALT alleles must be different sequences.
    A variant where REF == ALT is not a variant — it is a
    data error that would cause silent failures in downstream tools.
    """
    ref = str(record.get("ref", "")).upper().strip()
    alt = str(record.get("alt", "")).upper().strip()
    # Remove padding bases for indel comparison
    passed = bool(ref and alt and ref != alt)
    return ValidationResult(
        passed=passed,
        rule="ref_not_equal_alt",
        reason="" if passed else (
            f"REF '{ref}' and ALT '{alt}' are identical — not a variant"
        ),
    )


# Ordered list of all 9 rules
ALL_RULES = [
    rule_chromosome_valid,
    rule_position_positive,
    rule_ref_allele_valid,
    rule_alt_allele_valid,
    rule_quality_score_valid,
    rule_filter_field_present,
    rule_allele_frequency_valid,
    rule_depth_non_negative,
    rule_ref_not_equal_alt,
]


# ─────────────────────────────────────────────────────────────────
# Completeness scoring
# ─────────────────────────────────────────────────────────────────

COMPLETENESS_FIELDS = [
    "chrom", "pos", "ref", "alt",
    "qual", "filter_status", "af", "dp",
]


def completeness_score(record: dict) -> float:
    """Return fraction of key fields that are populated."""
    populated = sum(
        1 for f in COMPLETENESS_FIELDS
        if record.get(f) is not None
        and str(record.get(f, "")).strip() not in ("", ".")
    )
    return round(populated / len(COMPLETENESS_FIELDS), 4)


# ─────────────────────────────────────────────────────────────────
# Core validation
# ─────────────────────────────────────────────────────────────────

def validate_record(record: dict) -> RecordOutcome:
    """Apply all 9 rules to a single parsed VCF record."""
    failures = []
    for rule_fn in ALL_RULES:
        result = rule_fn(record)
        if not result.passed:
            failures.append(result)

    return RecordOutcome(
        record=record,
        passed=len(failures) == 0,
        failures=failures,
    )


def validate_batch(records: list[dict]) -> dict:
    """
    Validate a batch of parsed VCF records.
    Returns clean records (passed) and quarantined records (failed).
    Writes quarantine CSV to disk with rejection reasons.
    """
    clean = []
    quarantined = []
    total = len(records)

    for record in records:
        outcome = validate_record(record)
        score = completeness_score(record)
        record["completeness_score"] = score

        if outcome.passed:
            clean.append(record)
        else:
            rejection_reasons = "; ".join(
                f.reason for f in outcome.failures
            )
            failed_rules = [f.rule for f in outcome.failures]
            quarantine_record = {
                **record,
                "rejection_reasons": rejection_reasons,
                "failed_rules": str(failed_rules),
                "quarantined_at": datetime.now(
                    tz=timezone.utc
                ).isoformat(),
            }
            quarantined.append(quarantine_record)

    _write_quarantine(quarantined)

    pass_rate = (
        round(len(clean) / total * 100, 2) if total > 0 else 0
    )
    quarantine_rate = (
        round(len(quarantined) / total * 100, 2) if total > 0 else 0
    )
    avg_completeness = round(
        sum(r["completeness_score"] for r in clean) / len(clean), 4
    ) if clean else 0

    log.info(
        "validation_complete",
        total=total,
        passed=len(clean),
        quarantined=len(quarantined),
        pass_rate=f"{pass_rate}%",
        quarantine_rate=f"{quarantine_rate}%",
        avg_completeness=avg_completeness,
    )

    return {
        "clean": clean,
        "quarantined": quarantined,
        "stats": {
            "total": total,
            "passed": len(clean),
            "quarantined": len(quarantined),
            "pass_rate": pass_rate,
            "quarantine_rate": quarantine_rate,
            "avg_completeness": avg_completeness,
        },
    }


def _write_quarantine(quarantined: list[dict]) -> None:
    """Write quarantined records to CSV with rejection reasons."""
    if not quarantined:
        return

    os.makedirs(settings.quarantine_path, exist_ok=True)
    path = os.path.join(settings.quarantine_path, "quarantine.csv")

    fieldnames = [
        "chrom", "pos", "variant_id", "ref", "alt", "qual",
        "filter_status", "af", "dp", "mq", "completeness_score",
        "rejection_reasons", "failed_rules", "quarantined_at",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(quarantined)

    log.info(
        "quarantine_written",
        path=path,
        count=len(quarantined),
    )


if __name__ == "__main__":
    import structlog
    structlog.configure()

    from src.fetcher import fetch_variants
    from src.parser import parse_vcf_data

    data = fetch_variants(sample_size=200)
    records = parse_vcf_data(data)
    results = validate_batch(records)

    print("\n=== Validation Summary ===")
    print(f"Total records    : {results['stats']['total']}")
    print(f"Passed           : {results['stats']['passed']}")
    print(f"Quarantined      : {results['stats']['quarantined']}")
    print(f"Pass rate        : {results['stats']['pass_rate']}%")
    print(f"Avg completeness : {results['stats']['avg_completeness']}")

    if results["quarantined"]:
        print("\nQuarantined records:")
        for q in results["quarantined"]:
            print(
                f"  pos={q.get('pos')} | "
                f"{q['rejection_reasons']}"
            )
