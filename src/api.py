"""
FastAPI REST API for VariantLens.
Serves genomic variant data from the processed Parquet dataset.

Endpoints:
  GET /health          — service health and dataset statistics
  GET /variants        — paginated variant list with filters
  GET /variants/{id}   — single variant by position
  GET /summary         — dataset-level summary statistics
  GET /consequences    — consequence type distribution
  GET /af-tiers        — allele frequency tier distribution
  GET /clinical        — variants in disease-associated genes
  GET /quarantine      — quarantined records with rejection reasons
"""
from contextlib import asynccontextmanager
from typing import Optional
import os

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import pandas as pd
import structlog

from src.processor import load_processed_data, load_summary
from src.config import settings

log = structlog.get_logger()

# ── Global data store ────────────────────────────────────────────
_df: Optional[pd.DataFrame] = None
_summary: Optional[dict] = None


def get_df() -> pd.DataFrame:
    global _df
    if _df is None:
        _df = load_processed_data()
    return _df


def get_summary() -> dict:
    global _summary
    if _summary is None:
        _summary = load_summary()
    return _summary


@asynccontextmanager
async def lifespan(app: FastAPI):
    log.info("api_startup", version=settings.app_version)
    get_df()
    get_summary()
    log.info("data_loaded", rows=len(get_df()))
    yield
    log.info("api_shutdown")


# ── App ──────────────────────────────────────────────────────────
app = FastAPI(
    title="VariantLens API",
    description=(
        "REST API for querying human genomic variant calling results. "
        "Dataset: 1000 Genomes Project Phase 3 — Chromosome 22 (GRCh37/hg19). "
        "500 variants parsed, validated, and annotated from real population "
        "genomics data. Built to demonstrate production-grade bioinformatics "
        "data engineering."
    ),
    version=settings.app_version,
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


# ── Pydantic schemas ─────────────────────────────────────────────
class HealthResponse(BaseModel):
    status: str
    version: str
    total_variants: int
    pass_rate: float
    ts_tv_ratio: float
    reference_genome: str
    dataset: str


class VariantRecord(BaseModel):
    chrom: str
    pos: Optional[int]
    variant_id: Optional[str]
    ref: str
    alt: str
    qual: Optional[float]
    filter_status: Optional[str]
    af: Optional[float]
    dp: Optional[int]
    mq: Optional[float]
    variant_class: str
    ts_tv: str
    consequence: str
    gene: str
    disease_association: Optional[str]
    is_clinical_gene: bool
    af_tier: str
    dbsnp_status: str
    filter_outcome: str
    completeness_score: float


class ConsequenceCount(BaseModel):
    consequence: str
    count: int
    percentage: float


class AFTierCount(BaseModel):
    tier: str
    count: int
    percentage: float


class ClinicalVariant(BaseModel):
    variant_id: Optional[str]
    gene: str
    pos: Optional[int]
    ref: str
    alt: str
    consequence: str
    af: Optional[float]
    disease_association: Optional[str]


# ── Endpoints ────────────────────────────────────────────────────

@app.get("/", include_in_schema=False)
def root():
    return {
        "service": "VariantLens API",
        "version": settings.app_version,
        "docs": "/docs",
        "dataset": "1000 Genomes Project Phase 3 — Chr22",
        "reference": "GRCh37/hg19",
    }


@app.get("/health", response_model=HealthResponse, tags=["System"])
def health():
    """Service health check with dataset statistics."""
    try:
        df = get_df()
        summary = get_summary()
        return HealthResponse(
            status="healthy",
            version=settings.app_version,
            total_variants=len(df),
            pass_rate=summary["validation"]["pass_rate"],
            ts_tv_ratio=summary["ts_tv_ratio"],
            reference_genome=summary["reference_genome"],
            dataset=summary["dataset"],
        )
    except Exception as e:
        log.error("health_check_failed", error=str(e))
        raise HTTPException(status_code=503, detail="Data unavailable")


@app.get("/summary", tags=["Analysis"])
def summary():
    """
    Dataset-level summary statistics including variant class
    distribution, Ts/Tv ratio, consequence breakdown, AF tiers,
    and validation metrics.
    """
    return get_summary()


@app.get(
    "/variants",
    response_model=list[VariantRecord],
    tags=["Variants"],
)
def list_variants(
    variant_class: Optional[str] = Query(
        None,
        description="Filter by class: SNP, insertion, deletion, MNP"
    ),
    consequence: Optional[str] = Query(
        None,
        description=(
            "Filter by consequence: synonymous, missense, "
            "nonsense, frameshift, splice_region"
        )
    ),
    af_tier: Optional[str] = Query(
        None,
        description="Filter by AF tier: rare, low_frequency, common"
    ),
    gene: Optional[str] = Query(
        None,
        description="Filter by gene name (e.g. NF2, BCR, SMAD4)"
    ),
    clinical_only: bool = Query(
        False,
        description="Return only variants in disease-associated genes"
    ),
    pass_only: bool = Query(
        True,
        description="Return only PASS-filtered variants"
    ),
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    """
    Paginated variant list with filtering by class, consequence,
    allele frequency tier, gene, and clinical significance.
    """
    df = get_df()

    if variant_class:
        valid_classes = {"SNP", "insertion", "deletion", "MNP", "complex"}
        if variant_class not in valid_classes:
            raise HTTPException(
                status_code=400,
                detail=f"variant_class must be one of {valid_classes}"
            )
        df = df[df["variant_class"] == variant_class]

    if consequence:
        df = df[df["consequence"] == consequence]

    if af_tier:
        valid_tiers = {"rare", "low_frequency", "common", "unknown"}
        if af_tier not in valid_tiers:
            raise HTTPException(
                status_code=400,
                detail=f"af_tier must be one of {valid_tiers}"
            )
        df = df[df["af_tier"] == af_tier]

    if gene:
        df = df[df["gene"].str.upper() == gene.upper()]

    if clinical_only:
        df = df[df["is_clinical_gene"]]

    if pass_only:
        df = df[df["filter_outcome"] == "pass"]

    page = df.iloc[offset: offset + limit]
    records = page.where(pd.notna(page), None).to_dict(orient="records")
    return records


@app.get(
    "/variants/{pos}",
    response_model=list[VariantRecord],
    tags=["Variants"],
)
def get_variant_by_position(pos: int):
    """
    Retrieve all variants at a specific chromosomal position.
    Position is 1-based (GRCh37/hg19 coordinates).
    """
    df = get_df()
    match = df[df["pos"] == pos]
    if match.empty:
        raise HTTPException(
            status_code=404,
            detail=f"No variants found at position {pos}"
        )
    return match.where(pd.notna(match), None).to_dict(orient="records")


@app.get(
    "/consequences",
    response_model=list[ConsequenceCount],
    tags=["Analysis"],
)
def consequence_distribution():
    """
    Variant consequence type distribution with counts and percentages.
    Consequences: synonymous, missense, nonsense, frameshift,
    splice_region, inframe_insertion, inframe_deletion, MNP.
    """
    df = get_df()
    total = len(df)
    counts = df["consequence"].value_counts()
    return [
        ConsequenceCount(
            consequence=cons,
            count=int(cnt),
            percentage=round(float(cnt / total * 100), 2),
        )
        for cons, cnt in counts.items()
    ]


@app.get(
    "/af-tiers",
    response_model=list[AFTierCount],
    tags=["Analysis"],
)
def af_tier_distribution():
    """
    Allele frequency tier distribution.
    Rare: AF < 0.01, Low-frequency: 0.01-0.05, Common: > 0.05.
    """
    df = get_df()
    total = len(df)
    tier_order = ["rare", "low_frequency", "common", "unknown"]
    counts = df["af_tier"].value_counts()
    return [
        AFTierCount(
            tier=tier,
            count=int(counts.get(tier, 0)),
            percentage=round(
                float(counts.get(tier, 0) / total * 100), 2
            ),
        )
        for tier in tier_order
        if counts.get(tier, 0) > 0
    ]


@app.get(
    "/clinical",
    response_model=list[ClinicalVariant],
    tags=["Clinical"],
)
def clinical_variants():
    """
    Variants located in disease-associated genes on chromosome 22.
    Genes include BCR, ABL1, NF2, SMAD4, EP300, TBX1, SHANK3.
    These are directly relevant to leukaemia, neurofibromatosis,
    colorectal cancer, and neurodevelopmental disorders.
    """
    summary = get_summary()
    return summary.get("clinical_variants", [])


@app.get("/quarantine", tags=["Quality"])
def quarantine_records():
    """
    Records that failed validation — preserved with rejection reasons.
    Quarantine-not-delete: no data is silently discarded.
    """
    path = os.path.join(settings.quarantine_path, "quarantine.csv")
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path)
    return df.where(pd.notna(df), None).to_dict(orient="records")
