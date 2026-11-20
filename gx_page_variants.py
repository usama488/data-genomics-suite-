import pandas as pd
import plotly.express as px
import streamlit as st

import gx_pipeline as P
from gx_ui_common import need_case, show, tip

COLS = ["gene", "chrom", "gpos", "ref", "alt", "type", "gt", "zygosity", "dp", "alt_reads", "vaf", "qual", "gq", "filter",
        "consequence", "hgvsc", "hgvsp", "rsid", "gnomad_af", "clinvar", "classification"]


def render():
    c = need_case("Variant calling", "Filtering")
    st.title("🧪 Variant calls")
    tip("Variant calling", "VCF", "Genotype (GT)", "VAF", "Normalisation", "QUAL / GQ")
    R = c.results
    df = R.get("classified", R.get("annotated", R["calls"]))
    m = st.columns(5)
    m[0].metric("Total calls", len(df))
    m[1].metric("PASS", int((df["filter"] == "PASS").sum()))
    m[2].metric("SNVs / indels", f"{(df.type == 'SNV').sum()} / {(df.type != 'SNV').sum()}")
    m[3].metric("Het / Hom", f"{(df.zygosity == 'het').sum()} / {(df.zygosity == 'hom').sum()}")
    tv = P.titv(df[df["filter"] == "PASS"])
    m[4].metric("Ti/Tv", f"{tv:.2f}" if tv else "n/a", help="Few variants in a small panel → Ti/Tv is noisy")
    tip("Ti/Tv")
    if "validation" in R:
        v = R["validation"]
        st.success(f"🎯 Validation against the spiked truth set: **sensitivity {v['sensitivity']:.1f}%**, **precision {v['precision']:.1f}%** "
                   f"— {v['TP']} found, {v['FN']} missed, {v['FP']} false positives.")
        if v["missed"]:
            st.warning(f"Missed: {v['missed']}")
        if v["false_pos"]:
            st.warning(f"False positives: {v['false_pos']}")
    a, b, d = st.columns(3)
    only_pass = a.checkbox("PASS only", True)
    genes = b.multiselect("Genes", sorted(df.gene.dropna().unique()) if "gene" in df else [])
    q = d.text_input("Search (HGVS, rsID, position)")
    view = df[df["filter"] == "PASS"] if only_pass else df
    if genes:
        view = view[view.gene.isin(genes)]
    if q:
        view = view[view.astype(str).apply(lambda r: r.str.contains(q, case=False)).any(axis=1)]
    st.dataframe(view[[c_ for c_ in COLS if c_ in view.columns]], hide_index=True,
                 column_config={"vaf": st.column_config.ProgressColumn("VAF", min_value=0, max_value=1, format="%.2f"),
                                "gnomad_af": st.column_config.NumberColumn("gnomAD AF", format="%.2e")})
    x, y = st.columns(2)
    with x:
        show(px.histogram(df, x="vaf", color="filter", nbins=40, title="Variant allele fraction (het ≈ 0.5, hom ≈ 1.0)"), key="vaf")
        show(px.histogram(df, x="dp", color="filter", nbins=40, title="Depth at variant sites"), key="dp")
    with y:
        show(px.scatter(df, x="dp", y="qual", color="filter", hover_data=["gene", "gpos", "ref", "alt"], title="QUAL vs depth"), key="qd")
        show(px.pie(df, names="type", title="Variant types"), key="types")
    if "annotated" in R and "consequence" in df and df.consequence.notna().any():
        cons = df[df["filter"] == "PASS"].consequence.str.split(",").explode().value_counts().reset_index()
        cons.columns = ["consequence", "count"]
        show(px.bar(cons, x="count", y="consequence", orientation="h", title="Consequences (PASS)"), key="cons")
    st.download_button("⬇️ Filtered VCF (.vcf.gz)", open(R["vcf"], "rb").read(), R["vcf"].name)
    st.download_button("⬇️ Variant table (TSV)", df.to_csv(sep="\t", index=False).encode(), f"{c.sample}_variants.tsv")
