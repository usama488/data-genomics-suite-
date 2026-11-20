import streamlit as st

import gx_reference as R
import gx_tools as tools
import gx_workflow as W
from gx_learn import PIPELINE_DIAGRAM
from gx_ui_common import case

DEMO_GENES = R.PANELS["Mixed demo (BRCA1, TP53, CFTR, HBB)"]


def quick_demo():
    bar = st.progress(0.0, "Starting...")
    online = R.ensembl_online()
    name = "panel_demo" if online else "panel_demo_offline"
    try:
        if not (R.REF_ROOT / name / "reference.json").exists():
            R.build_panel(DEMO_GENES, name=name, offline=not online,
                          progress=lambda f, m: bar.progress(.25 * f, m))
        c = W.new_case("DEMO01")
        c.ref_name = name
        c.settings["offline"] = not online
        c.patient = {"name": "Demo Patient", "sample_id": "DEMO01", "sex": "Female", "dob": "1985-04-12",
                     "indication": "Personal and family history of early-onset breast cancer; carrier screening.",
                     "specimen": "Peripheral blood (EDTA)", "clinician": "Dr. Example", "mrn": "H-000123"}
        W.make_demo(c, coverage=100, progress=lambda f, m: bar.progress(.25 + .25 * f, m))
        W.run_all(c, progress=lambda f, m: bar.progress(.5 + .5 * f, m))
        st.session_state.gx_case = c
    except Exception as e:
        st.error(f"Demo failed: {e}")
        st.stop()
    bar.empty()
    if not online:
        st.warning("Ensembl could not be reached, so a **synthetic offline** reference was used and VEP annotation was skipped.")
    st.success("Demo complete! Open the pages in the sidebar to explore each step, then **📄 Clinical report**.")


def render():
    st.title("🧬 Clinical Germline Genomics Pipeline")
    st.caption("FASTQ → QC → trimming → alignment → coverage → variant calling → annotation → ACMG → clinical report")
    c1, c2 = st.columns([3, 2])
    with c1:
        st.code(PIPELINE_DIAGRAM, language=None)
    with c2:
        st.subheader("⚡ One-click demo")
        st.write("Downloads real GRCh38 sequence for **BRCA1, TP53, CFTR, HBB** from Ensembl, simulates a 100x paired-end "
                 "panel run with real ClinVar variants spiked in, and runs every step.")
        if st.button("▶️ Run full demo", type="primary"):
            quick_demo()
        st.subheader("Your own data")
        st.write("1. **🗂️ Reference** — build a gene panel or register a full GRCh38 genome\n"
                 "2. **🧾 Sample & FASTQ** — patient details + upload R1/R2\n"
                 "3. **▶️ Run pipeline** — run all steps or one at a time\n"
                 "4. Review QC → coverage → variants → ACMG → **📄 report**")
    c = case()
    if c:
        st.info(f"Current case **{c.sample}** · reference `{c.ref_name}` · steps done: {len(c.done)}/{len(W.STEPS)}")
    with st.expander("🔧 Tools in this environment"):
        st.json(tools.versions())
        st.caption(f"Aligner: {tools.aligner_name()}. On Streamlit Cloud the `bwa` command is installed from packages.txt; "
                   "samtools/bcftools run from the pysam library (same code as the command-line tools).")
    st.caption("⚠️ Research / education use. Clinical use requires laboratory validation, accreditation (e.g. ISO 15189 / CAP) "
               "and sign-off by a qualified clinical scientist.")
