"""Cancer Somatic Profiling Pipeline — Streamlit app.  Run: streamlit run app.py"""
import streamlit as st
import cx_page_case as page_case  # noqa: E402
import cx_page_cnv as page_cnv  # noqa: E402
import cx_page_home as page_home  # noqa: E402
import cx_page_learn as page_learn  # noqa: E402
import cx_page_qc as page_qc  # noqa: E402
import cx_page_reference as page_reference  # noqa: E402
import cx_page_report as page_report  # noqa: E402
import cx_page_run as page_run  # noqa: E402
import cx_page_somatic as page_somatic  # noqa: E402
import cx_page_tiering as page_tiering  # noqa: E402
import cx_workflow as W  # noqa: E402
from cx_ui_common import case, init  # noqa: E402

PAGES = {
    "🏠 Home & demo": page_home.render,
    "🗂️ Reference": page_reference.render,
    "🧾 Case & FASTQ": page_case.render,
    "▶️ Run pipeline": page_run.render,
    "🔬 QC (tumour / normal)": page_qc.render,
    "🧪 Somatic calls": page_somatic.render,
    "📈 CNV, purity & TMB": page_cnv.render,
    "🎯 Tiering & therapies": page_tiering.render,
    "📄 Oncology report": page_report.render,
    "📚 Learn": page_learn.render,
}


def main():
    init()
    st.sidebar.title("🎗️ Cancer Genomics")
    page = st.sidebar.radio("Navigate", list(PAGES), label_visibility="collapsed")
    st.sidebar.toggle("🎓 Learn mode", key="learn")
    c = case()
    if c:
        st.sidebar.markdown(f"**Case:** {c.case_id}  \n**Tumour:** {c.tumour_type}  \n**Mode:** {c.mode}")
        for s in W.STEPS:
            st.sidebar.caption(("✅ " if s in c.done else "⬜ ") + s)
    PAGES[page]()


