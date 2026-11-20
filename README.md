# 🧬 Data & Genomics Suite

One Streamlit website with four apps (choose from **Application** in the sidebar):

| App | What it does |
|---|---|
| 📊 **Auto EDA + ML Studio** (`eda_*.py`) | Any table → 15-step EDA, 4 dimensions, cleaning, AutoML, tuning, explainability, deep learning with epochs, clustering, HTML report |
| 🧬 **Germline Clinical Genomics** (`gx_*.py`) | FASTQ → QC → cutadapt → BWA-MEM → coverage → bcftools → Ensembl VEP → ACMG → clinical PDF |
| 🩻 **Chest X-ray AI** (`mx_*.py` + `cxr_model_part1/2.bin`) | Chest radiograph (PNG/JPG/DICOM) → QC → DenseNet-121 (TorchXRayVision, 18 findings, ONNX) → CAM / occlusion heat-maps → triage worklist → CTR → radiology PDF |
| 🎗️ **Cancer Somatic Profiling** (`cx_*.py`) | Tumour (+ normal) FASTQ → somatic calls → CNV / purity / TMB → AMP/ASCO/CAP tiers → targeted therapies → oncology PDF |

`app.py` is the entry point. All files sit in one folder (no sub-folders), so they can be uploaded to GitHub with "choose your files".

## Deploy on Streamlit Cloud
1. Upload every file (including `requirements.txt` and `packages.txt`) to a GitHub repository.
2. Go to share.streamlit.io → **Create app** → main file `app.py` → **Advanced settings** → Python **3.12** → **Deploy**.

## Run locally (Linux / macOS / WSL)
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```
For research and education. The genomics apps are not certified medical devices.
