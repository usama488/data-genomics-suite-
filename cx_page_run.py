import pandas as pd
import streamlit as st

import cx_workflow as W
from cx_ui_common import case, tip


def render():
    st.title("▶️ Run pipeline")
    c = case()
    if c is None or not c.fastq:
        st.info("Load FASTQ on **🧾 Case & FASTQ** first.")
        st.stop()
    st.write(f"Case **{c.case_id}** · {c.tumour_type} · {'tumour + normal' if c.mode == 'paired' else 'tumour-only'} · reference `{c.ref_name}`")
    with st.expander("⚙️ Somatic calling thresholds"):
        s = c.settings
        a = st.columns(4)
        s["min_tumor_dp"] = a[0].number_input("Min tumour depth", 5, 1000, int(s["min_tumor_dp"]))
        s["min_tumor_alt"] = a[1].number_input("Min alt reads", 2, 100, int(s["min_tumor_alt"]))
        s["min_tumor_vaf"] = a[2].number_input("Min tumour VAF", .01, .5, float(s["min_tumor_vaf"]), .01)
        s["max_normal_vaf"] = a[3].number_input("Max normal VAF", 0.0, .2, float(s["max_normal_vaf"]), .01)
        b = st.columns(4)
        s["fisher_p"] = b[0].number_input("Fisher p <", 1e-6, .1, float(s["fisher_p"]), format="%.4f")
        s["min_normal_dp"] = b[1].number_input("Min normal depth", 1, 500, int(s["min_normal_dp"]))
        s["min_bq"] = b[2].number_input("Min base quality", 0, 40, int(s["min_bq"]))
        s["offline"] = b[3].checkbox("Offline (skip VEP)", s.get("offline", False))
    tip("Somatic vs germline", "Fisher exact test", "Strand bias", "Limit of detection")
    st.dataframe(pd.DataFrame([{"#": i + 1, "step": s_, "status": "✅" if s_ in c.done else "⏳",
                                "seconds": c.timings.get(s_ + " (total)", "")} for i, s_ in enumerate(W.STEPS)]), hide_index=True)
    a, b, d = st.columns([1.3, 2, 1])
    go_all = a.button("▶️ Run all remaining", type="primary")
    start = b.selectbox("…or re-run from", W.STEPS, index=None, placeholder="choose a step")
    go_from = d.button("Re-run", disabled=start is None)
    if go_all or go_from:
        bar = st.progress(0.0)
        try:
            if go_from:
                c.reset_from(start)
            W.run_all(c, progress=lambda f, m: bar.progress(f, m))
        except Exception as e:
            st.error(f"Step failed: {e}")
            st.code("\n".join(c.log[-50:]))
            st.stop()
        st.rerun()
    if c.done:
        cols = st.columns(4)
        for i, (k, v) in enumerate(W.qc_summary(c).items()):
            cols[i % 4].metric(k, v)
        v = c.results.get("validation")
        if v:
            st.success(f"Validation vs truth: somatic sensitivity **{v['sensitivity']:.0f}%**, precision **{v['precision']:.0f}%** "
                       f"({v['TP']} TP / {v['FN']} FN / {v['FP']} FP); germline variants mis-called as somatic: {v['germline_miscalled_somatic']}")
    with st.expander("📜 Command log"):
        st.code("\n".join(c.log) or "-", language="bash")
