"""
VCF Parser for VariantLens.

Parses raw VCF format lines into structured Python dictionaries.

VCF Format (v4.1) — tab-separated fields:
  1. CHROM   — chromosome
  2. POS     — 1-based genomic position
  3. ID      — variant identifier (rsID or '.')
  4. REF     — reference allele
  5. ALT     — alternate allele(s) comma-separated
  6. QUAL    — Phred-scaled quality score
  7. FILTER  — filter status (PASS or reason for failure)
  8. INFO    — semicolon-separated key=value metadata
  9. FORMAT  — genotype format (optional)
  10+         — sample genotypes (optional)

The 1000 Genomes VCF contains 2,504 sample columns per variant.
We parse only the first 8 mandatory fields plus key INFO subfields.

Key INFO fields we extract:
  AF   — allele frequency in the population
  DP   — total read depth at this position
  MQ   — RMS mapping quality
  AN   — total number of alleles in called genotypes
  AC   — allele count in genotypes
  VT   — variant type annotation (SNP, INDEL, etc.)
"""
import re
import structlog

log = structlog.get_logger()

# Valid chromosomes for human genome (GRCh37)
VALID_CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

# Nucleotide bases — valid allele characters
VALID_BASES = set("ACGTNacgtn*")

# rsID pattern — dbSNP identifier
RSID_PATTERN = re.compile(r"^rs\d+$")


def parse_info_field(info_str: str) -> dict:
    """
    Parse the VCF INFO field into a dictionary.

    INFO field format: KEY=VALUE;KEY=VALUE;FLAG
    Flags (no value) are stored as True.

    Examples:
      AF=0.001;DP=100;MQ=60    → {'AF': '0.001', 'DP': '100', 'MQ': '60'}
      INDEL;AF=0.01;DP=50      → {'INDEL': True, 'AF': '0.01', 'DP': '50'}
    """
    result = {}
    if not info_str or info_str == ".":
        return result

    for token in info_str.split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, _, value = token.partition("=")
            result[key.strip()] = value.strip()
        else:
            result[token] = True

    return result


def parse_vcf_line(line: str, line_index: int = 0) -> dict:
    """
    Parse a single VCF data line into a structured record.

    Extracts the 8 mandatory VCF fields plus key INFO subfields.
    The 1000 Genomes VCF has 2,512 tab-separated fields per line
    (8 fixed + FORMAT + 2,504 sample columns). We split on tab
    and take only what we need.

    Returns a dict with all parsed fields plus metadata.
    """
    fields = line.strip().split("\t")

    if len(fields) < 8:
        return {
            "parse_error": True,
            "error_reason": (
                f"Expected at least 8 fields, got {len(fields)}"
            ),
            "raw_line": line[:200],
            "line_index": line_index,
        }

    chrom = fields[0].lstrip("chr")  # normalise chr22 → 22
    pos_str = fields[1]
    variant_id = fields[2]
    ref = fields[3]
    alt = fields[4]
    qual_str = fields[5]
    filter_status = fields[6]
    info_str = fields[7]

    # Parse position
    try:
        pos = int(pos_str)
    except ValueError:
        pos = None

    # Parse quality score
    try:
        qual = float(qual_str) if qual_str != "." else None
    except ValueError:
        qual = None

    # Parse INFO field
    info = parse_info_field(info_str)

    # Extract key INFO subfields
    try:
        af = float(info.get("AF", "").split(",")[0]) \
            if info.get("AF") else None
    except (ValueError, IndexError):
        af = None

    try:
        dp = int(info.get("DP", "")) if info.get("DP") else None
    except ValueError:
        dp = None

    try:
        mq = float(info.get("MQ", "")) if info.get("MQ") else None
    except ValueError:
        mq = None

    try:
        an = int(info.get("AN", "")) if info.get("AN") else None
    except ValueError:
        an = None

    try:
        ac_str = info.get("AC", "")
        ac = int(ac_str.split(",")[0]) if ac_str else None
    except (ValueError, IndexError):
        ac = None

    # Variant type from INFO VT field (1000 Genomes specific)
    vt = info.get("VT", None)

    # rsID status
    has_rsid = bool(RSID_PATTERN.match(variant_id)) \
        if variant_id and variant_id != "." else False

    # Multi-allelic check
    alt_alleles = [a.strip() for a in alt.split(",") if a.strip()]
    is_multiallelic = len(alt_alleles) > 1

    return {
        "chrom": chrom,
        "pos": pos,
        "variant_id": variant_id,
        "ref": ref,
        "alt": alt,
        "alt_alleles": alt_alleles,
        "qual": qual,
        "filter_status": filter_status,
        "af": af,
        "dp": dp,
        "mq": mq,
        "an": an,
        "ac": ac,
        "vt": vt,
        "has_rsid": has_rsid,
        "is_multiallelic": is_multiallelic,
        "line_index": line_index,
        "parse_error": False,
        "error_reason": None,
    }


def parse_vcf_data(data: dict) -> list[dict]:
    """
    Parse all variant lines from fetched VCF data.
    Returns a list of parsed variant records.
    Logs progress for large datasets.
    """
    variant_lines = data.get("variant_lines", [])
    total = len(variant_lines)

    log.info("parsing_vcf_data", total_lines=total)

    records = []
    parse_errors = 0

    for i, line in enumerate(variant_lines):
        if not line.strip():
            continue

        record = parse_vcf_line(line, line_index=i)

        if record.get("parse_error"):
            parse_errors += 1
            log.warning(
                "parse_error",
                line_index=i,
                reason=record.get("error_reason"),
            )
            continue

        # Attach source metadata
        record["chromosome"] = data.get("chromosome", "22")
        record["reference_genome"] = data.get(
            "reference_genome", "GRCh37/hg19"
        )
        record["dataset"] = data.get("dataset", "1000 Genomes Project")
        record["source"] = data.get("source", "unknown")

        records.append(record)

    log.info(
        "parsing_complete",
        parsed=len(records),
        parse_errors=parse_errors,
        total=total,
    )

    return records


def classify_variant_type(ref: str, alt: str) -> str:
    """
    Classify a variant as SNP, insertion, deletion, or MNP
    based on REF and ALT allele lengths.

    Classification rules (VCF standard):
      SNP      — REF and ALT both length 1, different bases
      Insertion — ALT longer than REF
      Deletion  — REF longer than ALT
      MNP       — REF and ALT same length > 1 (multi-nucleotide)
      Complex   — does not fit above categories
    """
    if not ref or not alt:
        return "unknown"

    # Take first ALT allele for classification
    alt_first = alt.split(",")[0].strip()

    ref_len = len(ref)
    alt_len = len(alt_first)

    if ref_len == 1 and alt_len == 1:
        return "SNP"
    elif alt_len > ref_len:
        return "insertion"
    elif ref_len > alt_len:
        return "deletion"
    elif ref_len == alt_len and ref_len > 1:
        return "MNP"
    else:
        return "complex"


def is_transition(ref: str, alt: str) -> bool:
    """
    Determine if a SNP is a transition (Ts) or transversion (Tv).

    Transitions — purine/purine or pyrimidine/pyrimidine substitutions:
      A <-> G  (purine to purine)
      C <-> T  (pyrimidine to pyrimidine)

    Transversions — purine/pyrimidine substitutions:
      A <-> C, A <-> T, G <-> C, G <-> T

    Ts/Tv ratio is a key quality metric in variant calling.
    A ratio of ~2.0-2.1 indicates high-quality SNP calls.
    """
    transitions = {
        frozenset({"A", "G"}),
        frozenset({"C", "T"}),
    }
    pair = frozenset({ref.upper(), alt.upper()})
    return pair in transitions


if __name__ == "__main__":
    import structlog
    structlog.configure()

    from src.fetcher import fetch_variants

    data = fetch_variants(sample_size=200)
    records = parse_vcf_data(data)

    print("\n=== Parser Summary ===")
    print(f"Total parsed     : {len(records)}")

    # Show variant type distribution
    types = {}
    for r in records:
        vt = classify_variant_type(r["ref"], r["alt"])
        types[vt] = types.get(vt, 0) + 1
    print(f"Variant types    : {types}")

    # Ts/Tv ratio
    snps = [
        r for r in records
        if classify_variant_type(r["ref"], r["alt"]) == "SNP"
    ]
    ts = sum(1 for r in snps if is_transition(r["ref"], r["alt"]))
    tv = len(snps) - ts
    ratio = round(ts / tv, 2) if tv > 0 else 0
    print(f"Ts/Tv ratio      : {ratio}  "
          f"(Ts={ts}, Tv={tv}) — expected ~2.0 for quality calls")

    # rsID rate
    rsid_count = sum(1 for r in records if r["has_rsid"])
    print(f"Known rsIDs      : {rsid_count}/{len(records)} "
          f"({round(rsid_count/len(records)*100, 1)}%)")

    # Show first record
    print("\nFirst record:")
    first = records[0]
    for k in ["chrom", "pos", "variant_id", "ref", "alt",
              "qual", "filter_status", "af", "dp", "mq"]:
        print(f"  {k:15} {first.get(k)}")
