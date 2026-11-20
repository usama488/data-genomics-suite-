import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from gx_ui_common import need_case, show, tip


def render():
    c = need_case("Alignment + duplicates", "BAM QC + coverage")
    st.title("🧭 Alignment & coverage")
    R = c.results
    f, cov = R["flagstat"], R["coverage"]
    st.caption(f"Aligner: **{R['align']['aligner']}**" + ("" if R["align"]["paired"] else
               " — library mode aligns mates independently, so duplicates are marked single-end (rate looks higher)."))
    tip("Alignment", "BAM", "MAPQ", "Duplication", "On-target rate")
    m = st.columns(6)
    m[0].metric("Total reads", f"{f['total']:,}")
    m[1].metric("Mapped", f"{f['mapped_pct']:.2f}%")
    m[2].metric("Properly paired", f"{f['properly_paired'] / max(f['total'], 1) * 100:.1f}%")
    m[3].metric("Duplicates", f"{f['dup_pct']:.1f}%")
    m[4].metric("On-target", f"{R['on_target']:.1f}%")
    m[5].metric("Mean target depth", f"{cov['overall']['mean_depth']:.0f}x")
    with st.expander("samtools flagstat"):
        st.code(f["raw"])
    a, b = st.columns(2)
    bs = R["bamstats"]
    with a:
        show(px.histogram(x=bs["mapq"], nbins=30, title="Mapping quality (MAPQ)", labels={"x": "MAPQ"}), key="mapq")
    with b:
        if len(bs["insert_size"]):
            show(px.histogram(x=bs["insert_size"], nbins=50, title=f"Insert size (median {int(np.median(bs['insert_size']))} bp)",
                              labels={"x": "fragment length (bp)"}), key="isize")
        else:
            st.info("Insert size needs properly paired reads.")
        tip("Insert size")

    st.header("Coverage")
    tip("Coverage / depth")
    o = cov["overall"]
    m = st.columns(5)
    for i, th in enumerate((1, 10, 20, 30, 100)):
        m[i].metric(f"Target bases ≥ {th}x", f"{o[f'pct_{th}x']:.2f}%")
    pg = cov["per_gene"]
    st.dataframe(pg.round(2), hide_index=True)
    show(px.bar(pg, x="gene", y="pct_20x", color="mean_depth", range_y=[0, 100], title="% target bases ≥ 20x per gene",
                color_continuous_scale="Blues"), key="g20")
    pt = cov["per_target"]
    pt["status"] = np.where(pt.pct_20x >= 100, "complete", np.where(pt.pct_20x >= 95, "minor gap", "GAP"))
    show(px.bar(pt, x="target", y="mean_depth", color="status", color_discrete_map={"complete": "#27AE60", "minor gap": "#F39C12", "GAP": "#C0392B"},
                title="Mean depth per exon", hover_data=["pct_20x", "min_depth"]), key="exon")
    gaps = cov["low"]
    if len(gaps):
        st.warning(f"{len(gaps)} region(s) below 20x — variants here could be missed. Labs usually fill these by Sanger sequencing.")
        st.dataframe(gaps, hide_index=True)
    else:
        st.success("Every target base is covered at ≥ 20x.")
    gene = st.selectbox("Depth track for gene", sorted(pt.gene.unique()))
    sub = pt[pt.gene == gene]
    dp = cov["depth"]
    d = dp[dp.contig == sub.contig.iloc[0]]
    fig = go.Figure()
    fig.add_scatter(x=d.pos, y=d.depth, mode="lines", fill="tozeroy", name="depth")
    for t in sub.itertuples():
        fig.add_vrect(x0=t.start, x1=t.end, fillcolor="green", opacity=.08, line_width=0)
    fig.add_hline(y=20, line_dash="dash", line_color="red", annotation_text="20x")
    fig.update_layout(title=f"Per-base depth across {gene} targets (green = exons)", xaxis_title=f"{sub.contig.iloc[0]} position",
                      yaxis_title="depth", height=330)
    show(fig, key="track")
    bam = R["align"]["bam"]
    st.download_button("⬇️ BAM", open(bam, "rb").read(), bam.name)
    st.download_button("⬇️ BAM index (.bai)", open(str(bam) + ".bai", "rb").read(), bam.name + ".bai")
