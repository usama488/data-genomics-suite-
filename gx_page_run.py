import pandas as pd
import streamlit as st

import gx_pipeline as P
import gx_workflow as W
from gx_ui_common import case, tip


def render():
    st.title("▶️ Run pipeline")
    c = case()
    if c is None or not c.r1:
        st.info("Load FASTQ on **🧾 Sample & FASTQ** first.")
        st.stop()
    st.write(f"Case **{c.sample}** · reference `{c.ref_name}` · input `{c.r1.split('/')[-1]}`"
             + (f" + `{c.r2.split('/')[-1]}`" if c.r2 else ""))
    with st.expander("⚙️ Settings", expanded=False):
        s = c.settings
        a = st.columns(4)
        s["quality"] = a[0].number_input("Trim bases below Q", 5, 40, int(s["quality"]))
        s["min_len"] = a[1].number_input("Discard reads shorter than (bp)", 15, 150, int(s["min_len"]))
        s["adapter1"] = a[2].text_input("Adapter R1", s["adapter1"])
        s["adapter2"] = a[3].text_input("Adapter R2", s["adapter2"])
        b = st.columns(4)
        s["min_qual"] = b[0].number_input("Min QUAL", 0.0, 1000.0, float(s["min_qual"]))
        s["min_dp"] = b[1].number_input("Min depth", 1, 500, int(s["min_dp"]))
        s["min_vaf"] = b[2].number_input("Min VAF", 0.0, 1.0, float(s["min_vaf"]), .05)
        s["min_gq"] = b[3].number_input("Min GQ", 0, 99, int(s["min_gq"]))
        s["offline"] = st.checkbox("Offline annotation (skip Ensembl VEP)", s.get("offline", False))
        st.caption("Changing settings? Re-run from the affected step below.")
    tip("Trimming", "Alignment", "Variant calling", "Annotation", "ACMG classification")

    rows = [{"#": i + 1, "step": s, "status": "✅ done" if s in c.done else "⏳ pending",
             "seconds": c.timings.get(s + " (total)", "")} for i, s in enumerate(W.STEPS)]
    st.dataframe(pd.DataFrame(rows), hide_index=True)
    a, b, d = st.columns([1.3, 2, 1])
    go_all = a.button("▶️ Run all remaining steps", type="primary")
    start = b.selectbox("…or re-run from step", W.STEPS, index=None, placeholder="choose a step")
    go_from = d.button("Re-run from here", disabled=start is None)
    if go_all or go_from:
        bar = st.progress(0.0)
        try:
            if go_from:
                c.reset_from(start)
            W.run_all(c, progress=lambda f, m: bar.progress(f, m))
            W.save_case(c)
        except Exception as e:
            st.error(f"Step failed: {e}")
            with st.expander("Log"):
                st.code("\n".join(c.log[-60:]))
            st.stop()
        st.success("Done.")
        st.rerun()
    if c.done:
        q = W.qc_summary(c)
        st.subheader("Run summary")
        cols = st.columns(4)
        for i, (k, v) in enumerate(q.items()):
            cols[i % 4].metric(k, v)
        if "validation" in c.results:
            v = c.results["validation"]
            st.success(f"Validation vs truth set: sensitivity **{v['sensitivity']:.1f}%**, precision **{v['precision']:.1f}%** "
                       f"({v['TP']} TP, {v['FN']} FN, {v['FP']} FP)")
            tip("Validation (truth set)")
    with st.expander("📜 Command log (every command that was run)"):
        st.code("\n".join(c.log) or "nothing yet", language="bash")
