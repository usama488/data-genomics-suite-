import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

import cx_fastq_qc as Q
from cx_ui_common import STATUS_ICON, findings, need_case, show


def render():
    c = need_case("FASTQ QC")
    st.title("🔬 Sequencing & alignment QC")
    S = c.results["samples"]
    status = {}
    for role, d in S.items():
        for mate, res in d.get("qc_raw", {}).items():
            status[f"{role} {mate}"] = {k: f"{STATUS_ICON[v]} {v}" for k, v in res["status"].items()}
    st.dataframe(pd.DataFrame(status))
    role = st.radio("Sample", list(S), horizontal=True)
    d = S[role]
    res = d["qc_raw"]["R1"]
    findings(Q.interpret(res))
    f = Q.figs(res, f"{role} R1")
    if "qc_trim" in d:
        o = go.Figure()
        for lbl, rr in (("raw", res), ("trimmed", d["qc_trim"]["R1"])):
            o.add_scatter(x=rr["per_base"].position, y=rr["per_base"]["mean"], name=f"{lbl} mean")
            o.add_scatter(x=rr["per_base"].position, y=rr["per_base"]["q10"], name=f"{lbl} 10th pct", line=dict(dash="dot"))
        o.update_layout(title=f"{role}: quality before vs after trimming", xaxis_title="position", yaxis_title="Phred")
        show(o, key="qo")
    a, b = st.columns(2)
    with a:
        show(f["quality"], key="q1")
        show(f["gc"], key="q2")
    with b:
        show(f["adapters"], key="q3")
        show(f["composition"], key="q4")
    if "flagstat" in d:
        st.subheader("Alignment")
        fl = d["flagstat"]
        m = st.columns(5)
        m[0].metric("Reads", f"{fl['total']:,}")
        m[1].metric("Mapped", f"{fl['mapped_pct']:.2f}%")
        m[2].metric("Duplicates", f"{fl['dup_pct']:.1f}%")
        m[3].metric("On-target", f"{d['on_target']:.1f}%")
        m[4].metric("Mean depth", f"{d['coverage']['overall']['mean_depth']:.0f}x")
        if len(d["bamstats"]["insert_size"]):
            show(px.histogram(x=d["bamstats"]["insert_size"], nbins=50, title="Insert size", labels={"x": "bp"}), key="is")
