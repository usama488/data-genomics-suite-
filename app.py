"""Qasim's Data & Genomics Suite — one Streamlit site hosting three apps.

Run:  streamlit run app.py
"""
import streamlit as st

st.set_page_config(page_title="Data & Genomics Suite", page_icon="🧬", layout="wide")

import cx_app  # noqa: E402
import eda_app  # noqa: E402
import gx_app  # noqa: E402
import mx_app  # noqa: E402

APPS = {
    "🏠 Home": None,
    "📊 Auto EDA + ML Studio": eda_app.main,
    "🧬 Germline Clinical Genomics": gx_app.main,
    "🎗️ Cancer Somatic Profiling": cx_app.main,
    "🩻 Chest X-ray AI": mx_app.main,
}

CARDS = [
    ("📊 Auto EDA + ML Studio", "#2E86AB",
     "Upload any table (CSV / Excel) → 15-step EDA, four dimensions of analysis, AutoML leaderboard, tuning, explainability, "
     "neural networks trained epoch by epoch, clustering and HTML reports.",
     ["15-step EDA + data cleaning", "10+ ML models, CV, tuning, SHAP / PDP", "Deep learning with live loss curves", "Learn mode for every concept"]),
    ("🧬 Germline Clinical Genomics", "#1B4F72",
     "Hospital-style inherited-disease testing: FASTQ → QC → trimming → BWA-MEM → coverage → bcftools variant calling → "
     "Ensembl VEP → ACMG classification → signed clinical PDF.",
     ["Real GRCh38 gene panels (BRCA1/2, CFTR, HBB…)", "ClinVar / gnomAD / REVEL annotation", "ACMG points + manual override", "IGV-style read viewer"]),
    ("🎗️ Cancer Somatic Profiling", "#6C3483",
     "Molecular pathology: tumour + matched normal FASTQ → somatic mutations, copy number, purity, TMB → AMP/ASCO/CAP tiers "
     "→ targeted therapies → tumour board report.",
     ["Tumour/normal or tumour-only calling", "ERBB2-style amplification detection", "Therapy matching (KRAS G12C, EGFR, BRAF…)", "Oncology PDF report"]),
    ("🩻 Chest X-ray AI", "#0B5345",
     "Radiology AI: chest X-ray (PNG/JPG/DICOM) → image QC → DenseNet-121 (18 findings) → heat-maps → triage worklist → "
     "radiologist review → report.",
     ["Pretrained on 8 public CXR datasets", "CAM + occlusion heat-maps", "STAT / urgent worklist triage", "CTR tool + radiology PDF"]),
]


def home():
    st.title("🧬 Data & Genomics Suite")
    st.caption("Four analysis apps in one place — pick one from the sidebar or below.")
    cols = st.columns(2) + st.columns(2)
    for col, (name, color, desc, bullets) in zip(cols, CARDS):
        with col:
            st.markdown(f"<div style='border-top:6px solid {color};padding:6px 4px 0 4px'><h3 style='margin-bottom:4px'>{name}</h3></div>",
                        unsafe_allow_html=True)
            st.write(desc)
            for b in bullets:
                st.markdown(f"- {b}")
            st.button(f"Open {name.split(' ', 1)[1]}", key=f"open_{name}", type="primary",
                      on_click=lambda n=name: st.session_state.update(suite_app=n))
    st.divider()
    st.markdown("**How to start:** EDA → *🏠 Data → Load demo*. Germline → *🏠 Home & demo → ▶️ Run full demo*. "
                "Cancer → *🏠 Home & demo → choose a scenario → ▶️ Run demo*. X-ray → *🩻 Study viewer → Open demo*.")
    st.caption("Built with Streamlit, scikit-learn, PyTorch (optional), BWA, samtools/bcftools (pysam), cutadapt, Ensembl VEP. "
               "Research & education use - not a certified medical device.")


def main():
    st.session_state.setdefault("suite_app", "🏠 Home")
    st.sidebar.selectbox("Application", list(APPS), key="suite_app")
    st.sidebar.divider()
    fn = APPS[st.session_state.suite_app]
    if fn is None:
        home()
    else:
        fn()


main()
