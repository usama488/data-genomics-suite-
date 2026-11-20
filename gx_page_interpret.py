import pandas as pd
import plotly.express as px
import streamlit as st

import gx_acmg as acmg
import gx_viewer as viewer
import gx_workflow as W
from gx_ui_common import need_case, show, tip


def _badge(cls):
    return f"<span style='background:{acmg.COLORS[cls]};color:white;padding:3px 10px;border-radius:12px;font-weight:600'>{cls}</span>"


def render():
    c = need_case("ACMG classification")
    st.title("⚖️ Variant interpretation (ACMG/AMP)")
    tip("Annotation", "HGVS", "MANE transcript", "gnomAD", "ClinVar", "In-silico predictors", "ACMG classification", "VUS")
    df = c.results["classified"]
    df_pass = df[df["filter"] == "PASS"]
    counts = df_pass.classification.value_counts().reindex(acmg.ORDER, fill_value=0).reset_index()
    counts.columns = ["classification", "n"]
    show(px.bar(counts, x="classification", y="n", color="classification", color_discrete_map=acmg.COLORS,
                title="PASS variants by ACMG class"), key="cls")
    if "annotated" not in df_pass or not df_pass["annotated"].fillna(False).astype(bool).any():
        st.warning("Variants are not annotated (offline / synthetic reference) - population and ClinVar evidence are missing, "
                   "so most variants will be VUS.")
    st.dataframe(df_pass[["gene", "hgvsc", "hgvsp", "consequence", "zygosity", "vaf", "gnomad_af", "clinvar", "revel",
                          "criteria", "points", "classification", "auto_or_manual"]], hide_index=True)

    st.subheader("Review one variant")
    labels = [f"{acmg.SHORT[r.classification]} · {r.gene} {r.hgvsc or ''} {r.hgvsp or ''} ({r.chrom}:{r.gpos} {r.ref}>{r.alt})"
              for r in df_pass.itertuples()]
    if not labels:
        st.info("No PASS variants.")
        return
    i = st.selectbox("Variant", range(len(labels)), format_func=lambda k: labels[k])
    v = df_pass.iloc[i]
    st.markdown(_badge(v.classification) + f" &nbsp; **{v.points} points** ({v.auto_or_manual})", unsafe_allow_html=True)
    a, b = st.columns(2)
    with a:
        info = {"Gene": v.gene, "Transcript": f"{v.transcript or ''} {('(MANE ' + v.mane + ')') if isinstance(v.mane, str) else ''}",
                "HGVS c.": v.hgvsc, "HGVS p.": v.hgvsp, "Exon / intron": v.exon if isinstance(v.exon, str) else v.intron,
                "Consequence": v.consequence, "Impact": v.impact, "Genomic (GRCh38)": f"chr{v.chrom}:{v.gpos} {v.ref}>{v.alt}",
                "Genotype": f"{v.gt} ({v.zygosity}), DP {v.dp}, VAF {v.vaf}", "dbSNP": v.rsid}
        st.table(pd.DataFrame(info.items(), columns=["field", "value"]).astype(str))
    with b:
        evid = {"gnomAD AF": v.gnomad_af, "ClinVar": v.clinvar, "REVEL": v.revel, "AlphaMissense": v.alphamissense,
                "SpliceAI Δmax": v.spliceai, "CADD": v.cadd, "SIFT": v.sift, "PolyPhen": v.polyphen}
        st.table(pd.DataFrame(evid.items(), columns=["evidence", "value"]).astype(str))
        if isinstance(v.rsid, str) and v.rsid.startswith("rs"):
            st.markdown(f"[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/?term={v.rsid}) · "
                        f"[gnomAD](https://gnomad.broadinstitute.org/variant/{v.chrom}-{v.gpos}-{v.ref}-{v.alt}?dataset=gnomad_r4) · "
                        f"[Ensembl](https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={v.rsid}) · "
                        f"[Franklin](https://franklin.genoox.com/clinical-db/variant/snp/chr{v.chrom}-{v.gpos}-{v.ref}-{v.alt})")
    st.markdown("**Evidence applied**")
    for part in str(v.evidence).split(" | "):
        if part:
            st.write("•", part)

    with st.expander("✏️ Clinical scientist override — edit ACMG criteria"):
        auto = acmg.auto_criteria(v.to_dict())
        current = c.overrides.get(v.variant_id, auto)
        chosen = st.multiselect("Criteria", list(acmg.CRITERIA), default=list(current),
                                format_func=lambda k: f"{k} — {acmg.CRITERIA[k][2][:70]}", key=f"crit_{v.variant_id}")
        new = {}
        cols = st.columns(4)
        for j, code in enumerate(chosen):
            default = current.get(code, (acmg.CRITERIA[code][1], "manual"))[0]
            opts = ["StandAlone"] if code == "BA1" else ["VeryStrong", "Strong", "Moderate", "Supporting"]
            st_ = cols[j % 4].selectbox(code, opts, index=opts.index(default) if default in opts else 0, key=f"str_{v.variant_id}_{code}")
            new[code] = (st_, current.get(code, (None, "added manually"))[1])
        preview, pts = acmg.classify(new)
        st.markdown(f"Preview: {_badge(preview)} ({pts} points)", unsafe_allow_html=True)
        x, y = st.columns(2)
        if x.button("💾 Save override", key=f"save_{v.variant_id}"):
            c.overrides[v.variant_id] = new
            c.results["classified"] = W.classify_all(c)
            st.rerun()
        if y.button("↩️ Revert to automatic", key=f"rev_{v.variant_id}"):
            c.overrides.pop(v.variant_id, None)
            c.results["classified"] = W.classify_all(c)
            st.rerun()

    st.subheader("🔍 Read evidence (IGV-style)")
    w = st.slider("Window (bp each side)", 15, 100, 40)
    fig, info = viewer.pileup_fig(str(c.results["align"]["bam"]), str(c.ref().fasta), v.contig, int(v.pos), window=w,
                                  title=f"{v.gene} {v.hgvsc or ''} — reads at chr{v.chrom}:{v.gpos}")
    show(fig, key="pile")
    st.caption("Blue = forward reads, pink = reverse, grey = duplicates. Coloured squares = mismatches, black lines = deletions, "
               "purple bars = insertions. A real heterozygous variant appears in ~half the reads on both strands.")
