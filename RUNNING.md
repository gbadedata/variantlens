# Running VariantLens

## API (port 8001)
```bash
cd ~/variantlens && source .venv/bin/activate && uvicorn src.api:app --reload --port 8001
```
Open: http://localhost:8001/docs

## Dashboard (port 8050)
```bash
cd ~/variantlens && source .venv/bin/activate && python3 -m src.dashboard
```
Open: http://localhost:8050
