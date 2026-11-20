"""Clinical Germline Genomics Pipeline — Streamlit app.

Run:  streamlit run app.py
"""
import streamlit as st
import gx_page_alignment as page_alignment  # noqa: E402
import gx_page_home as page_home  # noqa: E402
import gx_page_interpret as page_interpret  # noqa: E402
import gx_page_learn as page_learn  # noqa: E402
import gx_page_qc as page_qc  # noqa: E402
import gx_page_reference as page_reference  # noqa: E402
import gx_page_report as page_report  # noqa: E402
import gx_page_run as page_run  # noqa: E402
import gx_page_sample as page_sample  # noqa: E402
import gx_page_variants as page_variants  # noqa: E402
import gx_workflow as W  # noqa: E402
from gx_ui_common import case, init  # noqa: E402

PAGES = {
    "🏠 Home & demo": page_home.render,
    "🗂️ Reference": page_reference.render,
    "🧾 Sample & FASTQ": page_sample.render,
    "▶️ Run pipeline": page_run.render,
    "🔬 FASTQ QC & trimming": page_qc.render,
    "🧭 Alignment & coverage": page_alignment.render,
    "🧪 Variants": page_variants.render,
    "⚖️ Interpretation (ACMG)": page_interpret.render,
    "📄 Clinical report": page_report.render,
    "📚 Learn": page_learn.render,
}


def main():
    init()
    st.sidebar.title("🧬 Clinical Genomics")
    page = st.sidebar.radio("Navigate", list(PAGES), label_visibility="collapsed")
    st.sidebar.toggle("🎓 Learn mode", key="learn")
    c = case()
    if c:
        st.sidebar.markdown(f"**Case:** {c.sample}  \n**Reference:** {c.ref_name}")
        for s in W.STEPS:
            st.sidebar.caption(("✅ " if s in c.done else "⬜ ") + s)
    PAGES[page]()


