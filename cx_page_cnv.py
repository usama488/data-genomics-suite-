import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from cx_ui_common import need_case, show, tip


def render():
    c = need_case("Copy number, purity & TMB")
    st.title("📈 Copy number, purity & TMB")
    R = c.results
    tip("Copy number (CNV)", "Tumour purity", "TMB")
    m = st.columns(4)
    est = R.get("purity_estimate")
    m[0].metric("Purity (VAF estimate)", f"{est:.0%}" if est else "n/a")
    m[1].metric("Purity (pathology)", f"{c.pathology_purity:.0%}" if c.pathology_purity else "not given")
    tm = R["tmb"]
    m[2].metric("TMB", f"{tm['tmb']:.1f} mut/Mb", help=f"{tm['mutations']} mutations / {tm['megabases']:.3f} Mb ({tm['basis']})")
    m[3].metric("TMB class", tm["class"] if tm["reliable"] else "not reportable")
    if not tm["reliable"]:
        st.warning(f"This panel covers only {tm['megabases'] * 1000:.0f} kb. TMB needs ≥ 1 Mb of sequence to be meaningful, so the "
                   "value is shown for teaching only.")
    if est and c.pathology_purity and abs(est - c.pathology_purity) > .2:
        st.warning("VAF-based purity differs from the pathology estimate by > 20% - check the sample or the tumour content.")
    cnv = R.get("cnv_tiered", R["cnv"])
    if c.mode != "paired":
        st.info("Tumour-only: copy number is relative to the sample's own median depth (no normal) - a panel of normals is needed "
                "for reliable calls.")
    col = {"AMPLIFICATION": "#C0392B", "gain": "#E59866", "neutral": "#95A5A6", "loss": "#5DADE2", "DEEP DELETION": "#1B4F72"}
    fig = px.bar(cnv, x="gene", y="log2_ratio", color="call", color_discrete_map=col, error_y="sd",
                 hover_data=["copy_number", "exons"], title="Gene-level log2(tumour / normal) depth ratio")
    for y, t in ((1, "amplification"), (.4, "gain"), (-.4, "loss"), (-1, "deep loss")):
        fig.add_hline(y=y, line_dash="dot", line_color="grey", annotation_text=t)
    show(fig, key="cnv")
    show_cols = [x for x in ["gene", "exons", "log2_ratio", "copy_number", "call", "tier", "therapies_here", "rationale"] if x in cnv]
    st.dataframe(cnv[show_cols].round(2), hide_index=True)
    ex = R["cnv_exons"].copy()
    ex["log2"] = np.log2(ex.ratio.clip(lower=1e-3))
    show(px.scatter(ex, x="target", y="log2", color="gene", title="Exon-level log2 ratio"), key="ex")
    som = R["annotated"]
    som = som[som.status.str.startswith("somatic") & (som["filter"] == "PASS")] if len(som) else som
    if len(som):
        f = go.Figure()
        f.add_histogram(x=som.t_vaf, nbinsx=30, name="somatic VAF")
        if est:
            f.add_vline(x=est / 2, line_dash="dash", line_color="red", annotation_text=f"expected clonal VAF (purity {est:.0%} / 2)")
        f.update_layout(title="Somatic VAF spectrum - clonal vs subclonal", xaxis_title="tumour VAF")
        show(f, key="vaf")
        tip("Clonal / subclonal")
