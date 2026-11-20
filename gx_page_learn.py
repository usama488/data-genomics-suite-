import pandas as pd
import streamlit as st

import gx_acmg as acmg
from gx_learn import PIPELINE_DIAGRAM, STEP_INFO


def render():
    st.title("📚 Learn: how clinical sequencing analysis works")
    st.code(PIPELINE_DIAGRAM, language=None)
    q = st.text_input("Search a term", placeholder="e.g. VAF, MAPQ, ClinVar")
    for term, (exp, ex) in STEP_INFO.items():
        if q and q.lower() not in (term + exp).lower():
            continue
        st.markdown(f"**{term}** — {exp}" + (f"  \n*{ex}*" if ex else ""))
    st.subheader("ACMG/AMP evidence codes")
    st.dataframe(pd.DataFrame([{"code": k, "direction": "Pathogenic" if d == "P" else "Benign", "default strength": s,
                                "meaning": m} for k, (d, s, m) in acmg.CRITERIA.items()]), hide_index=True)
    st.markdown("**Points:** Very strong 8 · Strong 4 · Moderate 2 · Supporting 1 (benign codes subtract).  \n"
                "**≥ 10** Pathogenic · **6–9** Likely pathogenic · **0–5** VUS · **−1 to −6** Likely benign · **≤ −7** Benign · "
                "BA1 alone = Benign.")
    st.subheader("Phred quality calculator")
    qv = st.slider("Phred Q", 0, 50, 30)
    st.write(f"Q{qv} → error probability **{10 ** (-qv / 10):.5f}** (1 in {int(round(10 ** (qv / 10))):,} bases), accuracy {100 - 100 * 10 ** (-qv / 10):.3f}%")
    st.subheader("How many reads do I need?")
    a, b, c = st.columns(3)
    tgt = a.number_input("Target size (bp)", 1_000, 3_200_000_000, 50_000_000)
    cov = b.number_input("Desired mean depth", 1, 2000, 100)
    rl = c.number_input("Read length", 50, 300, 150)
    reads = tgt * cov / rl
    st.write(f"≈ **{reads / 1e6:,.1f} M reads** ({reads / 2e6:,.1f} M pairs) of usable on-target data. With ~70% on-target and "
             f"~10% duplicates, sequence ≈ **{reads / .7 / .9 / 1e6:,.1f} M reads**.")
