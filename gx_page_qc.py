import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import gx_fastq_qc as Q
from gx_ui_common import STATUS_ICON, findings, need_case, show, tip


def render():
    c = need_case("Raw FASTQ QC")
    st.title("🔬 FASTQ quality control & trimming")
    tip("FASTQ", "Phred quality", "Adapter", "Duplication")
    raw, trim = c.results["qc_raw"], c.results.get("qc_trim")
    # status grid (FastQC style)
    tbl = {}
    items = [(f"Raw {k}", v) for k, v in raw.items()] + ([(f"Trimmed {k}", v) for k, v in trim.items()] if trim else [])
    for label, res in items:
        tbl[label] = {m: f"{STATUS_ICON[s]} {s}" for m, s in res["status"].items()}
    st.dataframe(pd.DataFrame(tbl), width="stretch")
    met = pd.DataFrame({k: v["metrics"] for k, v in raw.items()}).T
    if trim:
        met = pd.concat([met.rename(index=lambda x: f"raw {x}"), pd.DataFrame({k: v["metrics"] for k, v in trim.items()}).T.rename(index=lambda x: f"trimmed {x}")])
    st.dataframe(met.round(2))
    if "trim" in c.results:
        t = c.results["trim"]
        st.subheader("✂️ Trimming (cutadapt)")
        tip("Trimming")
        m = st.columns(5)
        m[0].metric("Read pairs in", f"{t['reads_in']:,}")
        m[1].metric("Read pairs kept", f"{t['reads_out']:,}", f"{t['kept_pct']:.1f}%")
        m[2].metric("Too short (removed)", f"{t['too_short']:,}")
        m[3].metric("R1 with adapter", f"{t['with_adapter_r1']:,}")
        m[4].metric("Bases removed", f"{(t['bp_in'] - t['bp_out']) / max(t['bp_in'], 1) * 100:.1f}%")
    view = st.radio("Show", [k for k in (["raw"] + (["trimmed"] if trim else []))], horizontal=True)
    src = raw if view == "raw" else trim
    mate = st.radio("Read", list(src), horizontal=True)
    res = src[mate]
    findings(Q.interpret(res))
    f = Q.figs(res, f"{view} {mate}")
    show(f["quality"], key="q1")
    if trim and st.checkbox("Overlay raw vs trimmed quality", True):
        o = go.Figure()
        for lbl, rr in (("raw", raw[mate]), ("trimmed", trim[mate])):
            o.add_scatter(x=rr["per_base"].position, y=rr["per_base"]["mean"], name=f"{lbl} mean")
            o.add_scatter(x=rr["per_base"].position, y=rr["per_base"]["q10"], name=f"{lbl} 10th pct", line=dict(dash="dot"))
        o.update_layout(title="Mean & 10th-percentile quality: raw vs trimmed", xaxis_title="position", yaxis_title="Phred")
        show(o, key="q_over")
    a, b = st.columns(2)
    with a:
        show(f["read_q"], key="q2")
        show(f["gc"], key="q4")
        show(f["adapters"], key="q6")
    with b:
        show(f["composition"], key="q3")
        show(f["length"], key="q5")
        if res["overrep"]:
            st.markdown("**Overrepresented sequences**")
            st.dataframe(pd.DataFrame(res["overrep"], columns=["sequence", "count"]), hide_index=True)
