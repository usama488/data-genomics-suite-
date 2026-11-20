import re

import streamlit as st

import cx_oncokb as oncokb
import cx_reference as R
import cx_simulate as simulate
import cx_workflow as W
from cx_ui_common import case, tip


def _save(c, f, role, tag):
    ext = ".fastq.gz" if f.name.endswith(".gz") else ".fastq"
    p = c.dir / "fastq" / f"{c.case_id}_{role}_{tag}{ext}"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(f.getvalue())
    return str(p)


def render():
    st.title("🧾 Case, tumour & FASTQ")
    refs = [p.name for p in R.list_references()]
    if not refs:
        st.info("Build a reference first on **🗂️ Reference**, or run the demo on 🏠 Home.")
        st.stop()
    c = case()
    a, b = st.columns(2)
    cid = re.sub(r"[^A-Za-z0-9_-]", "_", a.text_input("Case ID", c.case_id if c else "TUMOUR02"))
    ref_name = b.selectbox("Reference", refs, index=refs.index(c.ref_name) if c and c.ref_name in refs else 0)
    if c is None or st.button("🆕 New case (clears results)"):
        c = W.new_case(cid)
        c.ref_name = ref_name
        st.session_state.cx_case = c
    if c.ref_name != ref_name:
        c.ref_name, c.done, c.results = ref_name, [], {"samples": {}}
    c.settings["offline"] = "synthetic" in R.Reference(R.REF_ROOT / ref_name).mode

    st.subheader("Clinical details")
    tip("Tumour purity", "Matched normal", "Tumour-only")
    p = c.patient
    x = st.columns(4)
    p["name"] = x[0].text_input("Patient", p.get("name", ""))
    p["mrn"] = x[1].text_input("MRN", p.get("mrn", ""))
    p["dob"] = x[2].text_input("DOB", p.get("dob", ""))
    p["sex"] = x[3].selectbox("Sex", ["Female", "Male", "Unknown"], index=["Female", "Male", "Unknown"].index(p.get("sex", "Unknown")))
    x = st.columns(3)
    tt = x[0].selectbox("Tumour type", oncokb.TUMOR_TYPES, index=oncokb.TUMOR_TYPES.index(c.tumour_type)
                        if c.tumour_type in oncokb.TUMOR_TYPES else 0)
    if tt != c.tumour_type:
        c.tumour_type = tt
        if "Tiering (AMP/ASCO/CAP)" in c.done:  # therapies depend on tumour type -> re-tier
            W.run_step(c, "Tiering (AMP/ASCO/CAP)")
    pur = x[1].number_input("Pathology tumour cellularity (%)", 0, 100, int(round((c.pathology_purity or 0) * 100)))
    if (pur / 100 or None) != c.pathology_purity:
        c.pathology_purity = pur / 100 or None
        if "Copy number, purity & TMB" in c.done:
            c.reset_from("Copy number, purity & TMB")
    p["specimen"] = x[2].text_input("Specimen", p.get("specimen", "FFPE tumour + blood normal"))
    x = st.columns(2)
    p["clinician"] = x[0].text_input("Referring oncologist", p.get("clinician", ""))
    p["block"] = x[1].text_input("Block / sample ID", p.get("block", ""))

    st.subheader("FASTQ input")
    mode = st.radio("Analysis mode", ["paired", "tumour-only"], index=0 if c.mode == "paired" else 1, horizontal=True,
                    format_func=lambda m: "Tumour + matched normal" if m == "paired" else "Tumour only")
    if mode != c.mode:
        c.mode, c.done, c.results = mode, [], {"samples": {}}
    src = st.radio("Source", ["Upload FASTQ", "Simulate demo"], horizontal=True)
    if src == "Upload FASTQ":
        cols = st.columns(2 if c.mode == "paired" else 1)
        files = {}
        for i, role in enumerate(c.roles):
            with cols[i]:
                st.markdown(f"**{role.capitalize()}**")
                files[role] = (st.file_uploader(f"{role} R1", type=["fastq", "fq", "gz"], key=f"{role}1"),
                               st.file_uploader(f"{role} R2 (optional)", type=["fastq", "fq", "gz"], key=f"{role}2"))
        if all(f[0] for f in files.values()) and st.button("Use these files", type="primary"):
            c.fastq = {r: (_save(c, f1, r, "R1"), _save(c, f2, r, "R2") if f2 else None) for r, (f1, f2) in files.items()}
            c.truth, c.done, c.results = None, [], {"samples": {}}
            st.success("Saved. Go to **▶️ Run pipeline**.")
    else:
        y = st.columns(4)
        sc = y[0].selectbox("Scenario", list(simulate.SCENARIOS))
        purity = y[1].slider("Purity", .2, 1.0, .6, .05)
        tc = y[2].slider("Tumour depth", 100, 800, 300, 50)
        nc = y[3].slider("Normal depth", 30, 200, 100, 10)
        if st.button("Simulate", type="primary"):
            bar = st.progress(0.0)
            W.make_demo(c, sc, purity, tc, nc, progress=lambda f, m: bar.progress(f, m))
            st.success(f"Simulated. Truth set ({len(c.truth)} variants, hidden from the pipeline):")
            st.dataframe(c.truth[["chrom", "gpos", "ref", "alt", "kind", "label", "tumor_af", "normal_af"]], hide_index=True)
    if c.fastq:
        st.info(" · ".join(f"{r}: `{str(v[0]).split('/')[-1]}`" for r, v in c.fastq.items()))
