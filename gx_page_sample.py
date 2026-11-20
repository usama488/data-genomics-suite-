import datetime as dt
import re

import streamlit as st

import gx_reference as R
import gx_workflow as W
from gx_ui_common import case, tip


def render():
    st.title("🧾 Sample, patient & FASTQ")
    refs = [p.name for p in R.list_references()]
    if not refs:
        st.info("Build a reference first on **🗂️ Reference** (or use **▶️ Run full demo** on the home page).")
        st.stop()
    c = case()
    st.subheader("1 · Case")
    a, b = st.columns(2)
    sample = a.text_input("Sample ID", c.sample if c else "PT0001")
    sample = re.sub(r"[^A-Za-z0-9_-]", "_", sample)
    ref_name = b.selectbox("Reference", refs, index=refs.index(c.ref_name) if c and c.ref_name in refs else 0)
    def fresh():
        n = W.new_case(sample)
        n.ref_name = ref_name
        st.session_state.gx_case = n
        return n

    if c is None:
        c = fresh()
    if st.button("🆕 New case (clears current results)"):
        fresh()
        st.rerun()
    if c.ref_name != ref_name:  # switching reference invalidates results
        c.ref_name, c.done, c.results = ref_name, [], {}
    c.settings["offline"] = "synthetic" in R.Reference(R.REF_ROOT / ref_name).mode

    st.subheader("2 · Patient & referral")
    p = c.patient
    x = st.columns(3)
    p["name"] = x[0].text_input("Patient name", p.get("name", ""))
    p["mrn"] = x[1].text_input("MRN / hospital number", p.get("mrn", ""))
    p["dob"] = x[2].text_input("Date of birth", p.get("dob", ""))
    x = st.columns(3)
    p["sex"] = x[0].selectbox("Sex", ["Female", "Male", "Unknown"], index=["Female", "Male", "Unknown"].index(p.get("sex", "Unknown")))
    p["specimen"] = x[1].text_input("Specimen", p.get("specimen", "Peripheral blood (EDTA)"))
    p["clinician"] = x[2].text_input("Referring clinician", p.get("clinician", ""))
    x = st.columns(2)
    p["collected"] = x[0].text_input("Collected", p.get("collected", str(dt.date.today())))
    p["received"] = x[1].text_input("Received", p.get("received", str(dt.date.today())))
    p["indication"] = st.text_area("Clinical indication / phenotype", p.get("indication", ""))
    p["sample_id"] = c.sample

    st.subheader("3 · FASTQ input")
    tip("FASTQ", "Phred quality")
    mode = st.radio("Input", ["Upload my FASTQ files", "Generate demo FASTQ (with known variants)"], horizontal=True)
    if mode.startswith("Upload"):
        f1 = st.file_uploader("R1 (FASTQ or FASTQ.gz)", type=["fastq", "fq", "gz"], key="r1")
        f2 = st.file_uploader("R2 (optional, for paired-end)", type=["fastq", "fq", "gz"], key="r2")
        if f1 and st.button("Use these files"):
            def save(f, tag):
                ext = ".fastq.gz" if f.name.endswith(".gz") else ".fastq"
                path = c.dir / f"{c.sample}_{tag}{ext}"
                path.write_bytes(f.getvalue())
                return str(path)
            c.r1 = save(f1, "R1")
            c.r2 = save(f2, "R2") if f2 else None
            c.truth, c.done, c.results = None, [], {}
            st.success("FASTQ saved. Go to **▶️ Run pipeline**.")
        st.caption("Reads from outside the reference's target genes simply won't map — that's expected for a panel reference.")
    else:
        y = st.columns(3)
        cov = y[0].slider("Target coverage", 20, 300, 100, 10)
        nrand = y[1].slider("Random novel variants", 0, 15, 6)
        clin = y[2].checkbox("Spike real ClinVar variants", True, disabled=c.settings.get("offline", False))
        if st.button("Generate demo FASTQ", type="primary"):
            bar = st.progress(0.0)
            W.make_demo(c, cov, nrand, clin, progress=lambda f, m: bar.progress(f, m))
            st.success(f"Simulated {len(c.truth)} variants. Truth set (hidden from the pipeline, used for validation):")
            st.dataframe(c.truth[["chrom", "gpos", "ref", "alt", "zygosity", "source", "id"]], hide_index=True)
    if c.r1:
        st.info(f"Input: `{c.r1.split('/')[-1]}`" + (f" + `{c.r2.split('/')[-1]}`" if c.r2 else " (single-end)"))
