"""
Data Processor for VariantLens.

Takes annotated VCF records and:
  1. Structures them into an analysis-ready pandas DataFrame
  2. Computes dataset-level summary statistics
  3. Saves to Parquet format for fast analytical queries
  4. Saves a JSON summary for the API and dashboard

Parquet is the de facto standard for analytical data in genomics
pipelines — columnar storage means consequence queries only read
the consequence column, not all 30 fields.
"""
import os
import json
import pandas as pd
import structlog

from src.config import settings

log = structlog.get_logger()


def process_annotated_records(
    annotated_records: list[dict],
    validation_stats: dict,
) -> pd.DataFrame:
    """
    Convert annotated records into a clean DataFrame.
    Saves Parquet and summary JSON to disk.
    """
    log.info("processing_records", count=len(annotated_records))

    df = pd.DataFrame(annotated_records)

    # Ensure correct column types
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["qual"] = pd.to_numeric(df["qual"], errors="coerce")
    df["af"] = pd.to_numeric(df["af"], errors="coerce")
    df["dp"] = pd.to_numeric(df["dp"], errors="coerce")
    df["mq"] = pd.to_numeric(df["mq"], errors="coerce")
    df["completeness_score"] = pd.to_numeric(
        df["completeness_score"], errors="coerce"
    )

    # Sort by position for genomic order
    df = df.sort_values("pos").reset_index(drop=True)

    # Save to Parquet
    os.makedirs(settings.processed_data_path, exist_ok=True)
    parquet_path = os.path.join(
        settings.processed_data_path, "variants.parquet"
    )
    df.to_parquet(parquet_path, index=False)
    log.info("parquet_saved", path=parquet_path, rows=len(df))

    # Compute and save summary
    summary = _compute_summary(df, validation_stats)
    summary_path = os.path.join(
        settings.processed_data_path, "summary.json"
    )
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    log.info("summary_saved", path=summary_path)

    return df


def _compute_summary(df: pd.DataFrame, validation_stats: dict) -> dict:
    """Compute dataset-level summary statistics."""
    total = len(df)

    # Variant class counts
    class_counts = df["variant_class"].value_counts().to_dict()

    # Ts/Tv ratio
    snps = df[df["variant_class"] == "SNP"]
    ts = (snps["ts_tv"] == "transition").sum()
    tv = (snps["ts_tv"] == "transversion").sum()
    tstv_ratio = round(float(ts / tv), 3) if tv > 0 else 0

    # Consequence distribution
    consequence_counts = df["consequence"].value_counts().to_dict()

    # AF tier distribution
    af_tier_counts = df["af_tier"].value_counts().to_dict()

    # Filter outcome
    filter_counts = df["filter_outcome"].value_counts().to_dict()

    # Clinical variants
    clinical_df = df[df["is_clinical_gene"]]
    clinical_variants = clinical_df[
        ["variant_id", "gene", "pos", "ref", "alt",
         "consequence", "af", "disease_association"]
    ].to_dict(orient="records")

    # dbSNP status
    dbsnp_counts = df["dbsnp_status"].value_counts().to_dict()

    # Quality metrics
    qual_stats = {}
    if "qual" in df.columns and df["qual"].notna().any():
        qual_stats = {
            "mean": round(float(df["qual"].mean()), 2),
            "median": round(float(df["qual"].median()), 2),
            "min": round(float(df["qual"].min()), 2),
            "max": round(float(df["qual"].max()), 2),
        }

    # Depth stats
    dp_stats = {}
    if "dp" in df.columns and df["dp"].notna().any():
        dp_stats = {
            "mean": round(float(df["dp"].mean()), 2),
            "median": round(float(df["dp"].median()), 2),
            "min": int(df["dp"].min()),
            "max": int(df["dp"].max()),
        }

    # Top genes by variant count
    gene_counts = (
        df[df["gene"] != "intergenic"]["gene"]
        .value_counts()
        .head(10)
        .to_dict()
    )

    return {
        "dataset": df["dataset"].iloc[0] if total > 0 else "",
        "reference_genome": (
            df["reference_genome"].iloc[0] if total > 0 else ""
        ),
        "chromosome": df["chrom"].iloc[0] if total > 0 else "22",
        "source": df["source"].iloc[0] if total > 0 else "",
        "total_variants": int(total),
        "validation": {
            "total_ingested": validation_stats.get("total", 0),
            "passed": validation_stats.get("passed", 0),
            "quarantined": validation_stats.get("quarantined", 0),
            "pass_rate": validation_stats.get("pass_rate", 0),
            "quarantine_rate": validation_stats.get("quarantine_rate", 0),
            "avg_completeness": validation_stats.get(
                "avg_completeness", 0
            ),
        },
        "variant_classes": {
            k: int(v) for k, v in class_counts.items()
        },
        "ts_tv_ratio": tstv_ratio,
        "ts_count": int(ts),
        "tv_count": int(tv),
        "consequences": {
            k: int(v) for k, v in consequence_counts.items()
        },
        "af_tiers": {
            k: int(v) for k, v in af_tier_counts.items()
        },
        "filter_outcomes": {
            k: int(v) for k, v in filter_counts.items()
        },
        "dbsnp_status": {
            k: int(v) for k, v in dbsnp_counts.items()
        },
        "clinical_variants": clinical_variants,
        "indel_count": int(total - class_counts.get("SNP", 0)),
        "clinical_variant_count": len(clinical_df),
        "gene_variant_counts": {
            k: int(v) for k, v in gene_counts.items()
        },
        "quality_stats": qual_stats,
        "depth_stats": dp_stats,
    }


def load_processed_data() -> pd.DataFrame:
    """Load the processed Parquet file."""
    parquet_path = os.path.join(
        settings.processed_data_path, "variants.parquet"
    )
    if not os.path.exists(parquet_path):
        raise FileNotFoundError(
            f"Processed data not found at {parquet_path}. "
            "Run: python3 -m src.processor"
        )
    return pd.read_parquet(parquet_path)


def load_summary() -> dict:
    """Load the summary JSON."""
    summary_path = os.path.join(
        settings.processed_data_path, "summary.json"
    )
    if not os.path.exists(summary_path):
        raise FileNotFoundError(
            f"Summary not found at {summary_path}. "
            "Run: python3 -m src.processor"
        )
    with open(summary_path) as f:
        return json.load(f)


if __name__ == "__main__":
    import structlog
    structlog.configure()

    from src.fetcher import fetch_variants
    from src.parser import parse_vcf_data
    from src.validator import validate_batch
    from src.annotator import annotate_batch

    data = fetch_variants()
    records = parse_vcf_data(data)
    results = validate_batch(records)
    annotated = annotate_batch(results["clean"])
    df = process_annotated_records(annotated, results["stats"])
    summary = load_summary()

    print("\n=== Processing Summary ===")
    print(f"Total variants   : {summary['total_variants']}")
    print(f"Reference genome : {summary['reference_genome']}")
    print(f"Ts/Tv ratio      : {summary['ts_tv_ratio']}")
    print(f"SNPs             : {summary['variant_classes'].get('SNP', 0)}")
    indel_count = summary['total_variants'] - summary['variant_classes'].get('SNP', 0)
    print(f"Indels           : {indel_count}")
    print(f"Rare variants    : {summary['af_tiers'].get('rare', 0)}")
    print(f"Common variants  : {summary['af_tiers'].get('common', 0)}")
    print(f"Clinical genes   : {summary['clinical_variant_count']}")
    print(f"Pass rate        : {summary['validation']['pass_rate']}%")
    print("Parquet saved    : data/processed/variants.parquet")
