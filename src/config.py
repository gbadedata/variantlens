"""Configuration management using pydantic-settings."""
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    ncbi_email: str = "gbadedata@gmail.com"
    ncbi_api_key: str = ""
    app_version: str = "1.0.0"
    log_level: str = "INFO"

    # VCF dataset settings
    vcf_sample_size: int = 500

    # Data paths
    raw_data_path: str = "data/raw"
    processed_data_path: str = "data/processed"
    quarantine_path: str = "data/quarantine"
    reference_path: str = "data/reference"

    # Validation thresholds
    min_quality_score: float = 0.0
    max_quality_score: float = 10000.0
    min_depth: int = 0

    class Config:
        env_file = ".env"


settings = Settings()
