import pandas as pd
import streamlit as st

import cx_oncokb as oncokb
from cx_learn import PIPELINE_DIAGRAM, STEP_INFO


def render():
    st.title("📚 Learn: cancer genomics")
    st.code(PIPELINE_DIAGRAM, language=None)
    q = st.text_input("Search", placeholder="e.g. purity, TMB, tier")
    for k, (e, ex) in STEP_INFO.items():
        if not q or q.lower() in (k + e).lower():
            st.markdown(f"**{k}** — {e}" + (f"  \n*{ex}*" if ex else ""))
    st.subheader("Actionability knowledge base (snapshot)")
    rows = [{"gene": g, "alteration": a, "tumour types": ", ".join(sorted(t)).replace("ANY_SOLID", "any solid tumour (tissue-agnostic)"),
             "therapies": d, "effect": e, "evidence": n} for (g, a), lst in oncokb.KB.items() for t, d, e, n in lst]
    st.dataframe(pd.DataFrame(rows), hide_index=True)
    st.caption("Curated from FDA labels and NCCN/ESMO guidelines (mid-2026). Always confirm with OncoKB / CIViC / current guidelines.")
    st.subheader("Expected VAF calculator")
    a, b, c = st.columns(3)
    pur = a.slider("Purity", .05, 1.0, .6, .05)
    ccf = b.slider("Cancer cell fraction (clonality)", .05, 1.0, 1.0, .05)
    cn = c.slider("Tumour copy number at locus", 1, 8, 2)
    mut = st.slider("Mutated copies", 1, cn, 1)
    vaf = pur * ccf * mut / (pur * cn + (1 - pur) * 2)
    st.write(f"Expected VAF = purity·CCF·mutated copies / (purity·CN + (1−purity)·2) = **{vaf:.1%}**")
    dp = st.slider("Sequencing depth", 50, 2000, 300, 50)
    st.write(f"≈ **{vaf * dp:.0f}** alt reads expected at {dp}x " + ("✅ detectable (≥ 4 reads)" if vaf * dp >= 4 and vaf >= .05 else "⚠️ near/below the limit of detection"))
