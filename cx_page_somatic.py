import plotly.express as px
import streamlit as st

from cx_ui_common import need_case, show, tip

COLS = ["gene", "protein", "hgvsc", "consequence", "chrom", "gpos", "ref", "alt", "type", "t_vaf", "t_alt", "t_dp", "n_vaf", "n_dp",
        "fisher_p", "status", "filter", "cosmic", "gnomad_af"]


def render():
    c = need_case("Somatic calling")
    st.title("🧪 Somatic & germline variant calls")
    tip("Somatic vs germline", "VAF", "Clonal / subclonal", "Strand bias")
    R = c.results
    df = R.get("tiered") if "tiered" in R and len(R["tiered"]) else None
    allv = R.get("annotated", R["calls"])
    if allv is None or allv.empty:
        st.warning("No candidate variants found.")
        return
    if "protein" not in allv:
        import cx_oncokb as oncokb
        allv = allv.assign(protein=[oncokb.short_protein(h) for h in allv.get("hgvsp", [None] * len(allv))])
    m = st.columns(5)
    som = allv[allv.status.str.startswith("somatic") & (allv["filter"] == "PASS")]
    m[0].metric("Somatic PASS", len(som))
    m[1].metric("Germline", int(allv.status.str.contains("germline").sum()))
    m[2].metric("Filtered / uncertain", int((allv["filter"] != "PASS").sum()))
    m[3].metric("Median somatic VAF", f"{som.t_vaf.median():.1%}" if len(som) else "-")
    m[4].metric("Analysis", "T + N" if c.mode == "paired" else "Tumour-only")
    v = R.get("validation")
    if v:
        st.success(f"🎯 Truth set: sensitivity {v['sensitivity']:.0f}%, precision {v['precision']:.0f}% — missed {v['missed'] or 'none'}; "
                   f"false positives {v['false_pos'] or 'none'}")
    view = st.radio("Show", ["Somatic PASS", "Germline", "All candidates"], horizontal=True)
    sub = som if view == "Somatic PASS" else allv[allv.status.str.contains("germline")] if view == "Germline" else allv
    st.dataframe(sub[[x for x in COLS if x in sub.columns]], hide_index=True,
                 column_config={"t_vaf": st.column_config.ProgressColumn("tumour VAF", min_value=0, max_value=1, format="%.3f"),
                                "n_vaf": st.column_config.NumberColumn("normal VAF", format="%.3f"),
                                "fisher_p": st.column_config.NumberColumn("Fisher p", format="%.1e")})
    if view == "Germline" and len(sub):
        st.info("Germline variants are inherited. Pathogenic germline findings in cancer-predisposition genes need genetic counselling "
                "and are not part of the somatic interpretation.")
    a, b = st.columns(2)
    with a:
        show(px.histogram(allv, x="t_vaf", color="status", nbins=50, title="Tumour VAF distribution (clonal peak ≈ purity/2)"), key="v1")
    with b:
        if c.mode == "paired":
            show(px.scatter(allv, x="n_vaf", y="t_vaf", color="status", hover_data=["gene", "gpos", "ref", "alt"],
                            title="Tumour VAF vs normal VAF", range_x=[-.02, 1.02], range_y=[-.02, 1.02]), key="v2")
        else:
            show(px.scatter(allv, x="t_dp", y="t_vaf", color="status", hover_data=["gene", "gpos"], title="Tumour VAF vs depth"), key="v2")
    st.caption("Somatic mutations sit on the left of the scatter (absent in normal); germline variants lie on the diagonal at 0.5 / 1.0.")
    if R.get("vcf"):
        st.download_button("⬇️ Somatic VCF", open(R["vcf"], "rb").read(), f"{c.case_id}.somatic.vcf.gz")
    st.download_button("⬇️ All calls (TSV)", allv.to_csv(sep="\t", index=False).encode(), f"{c.case_id}_calls.tsv")
