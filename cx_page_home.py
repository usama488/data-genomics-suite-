import streamlit as st

import cx_reference as R
import cx_simulate as simulate
import cx_tools as tools
import cx_workflow as W
from cx_learn import PIPELINE_DIAGRAM
from cx_ui_common import case

DEMO_GENES = ["EGFR", "KRAS", "NRAS", "BRAF", "PIK3CA", "TP53", "IDH1", "KIT", "ERBB2"]


def run_demo(scenario: str, mode: str, purity: float):
    bar = st.progress(0.0, "Starting...")
    online = R.ensembl_online()
    name = "cancer_panel" if online else "cancer_panel_offline"
    try:
        if not (R.REF_ROOT / name / "reference.json").exists():
            R.build_panel(DEMO_GENES, name=name, offline=not online, progress=lambda f, m: bar.progress(.25 * f, m))
        c = W.new_case("TUMOUR01")
        c.ref_name, c.mode = name, mode
        c.settings["offline"] = not online
        c.patient = {"name": "Demo Patient", "dob": "1962-03-08", "sex": "Male", "mrn": "H-004512", "clinician": "Dr. Oncologist",
                     "specimen": "FFPE tumour block + EDTA blood (normal)", "block": "A1"}
        W.make_demo(c, scenario, purity, progress=lambda f, m: bar.progress(.25 + .35 * f, m))
        W.run_all(c, progress=lambda f, m: bar.progress(.6 + .4 * f, m))
        st.session_state.cx_case = c
    except Exception as e:
        st.error(f"Demo failed: {e}")
        st.stop()
    bar.empty()
    if not online:
        st.warning("Ensembl could not be reached: a **synthetic offline** reference was used, real hotspots could not be placed and "
                   "annotation/therapy matching was skipped.")
    st.success("Demo complete - open **🎯 Tiering & therapies** and **📄 Oncology report**.")


def render():
    st.title("🎗️ Cancer Somatic Profiling Pipeline")
    st.caption("Tumour (+ matched normal) FASTQ → QC → alignment → somatic mutations → CNV / purity / TMB → "
               "AMP/ASCO/CAP tiers → targeted therapies → molecular tumour board report")
    a, b = st.columns([3, 2])
    with a:
        st.code(PIPELINE_DIAGRAM, language=None)
    with b:
        st.subheader("⚡ One-click demo")
        sc = st.selectbox("Clinical scenario", list(simulate.SCENARIOS))
        mode = st.radio("Analysis", ["paired", "tumour-only"], horizontal=True,
                        format_func=lambda m: "Tumour + matched normal" if m == "paired" else "Tumour only")
        purity = st.slider("Tumour purity", .2, 1.0, .6, .05)
        st.caption("Downloads real GRCh38 sequence for 9 cancer genes, places the real hotspot mutations (mapped from HGVS via "
                   "Ensembl), simulates tumour 300x / normal 100x reads and runs every step.")
        if st.button("▶️ Run demo", type="primary"):
            run_demo(sc, mode, purity)
    c = case()
    if c:
        st.info(f"Current case **{c.case_id}** · {c.tumour_type} · {c.mode} · steps {len(c.done)}/{len(W.STEPS)}")
    with st.expander("🔧 Tools in this environment"):
        st.json(tools.versions())
    st.caption("⚠️ Research / education use. Therapy associations come from a curated snapshot - always confirm with OncoKB, CIViC, "
               "NCCN/ESMO guidelines and a molecular tumour board.")
