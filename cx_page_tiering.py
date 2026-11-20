import pandas as pd
import plotly.express as px
import streamlit as st

import cx_oncokb as oncokb
import cx_viewer as viewer
import cx_workflow as W
from cx_ui_common import need_case, show, tip


def _badge(t):
    return f"<span style='background:{oncokb.TIER_COLORS.get(t, '#888')};color:white;padding:3px 10px;border-radius:12px;font-weight:600'>{t}</span>"


def render():
    c = need_case("Tiering (AMP/ASCO/CAP)")
    st.title("🎯 Tiering & therapy matching")
    st.caption(f"Tumour type: **{c.tumour_type}** (change it on 🧾 Case & FASTQ — tiers update automatically)")
    tip("AMP/ASCO/CAP tiers", "Actionable", "Hotspot", "Resistance mutation", "COSMIC")
    t = c.results["tiered"]
    cnv = c.results["cnv_tiered"]
    acts = []
    if len(t):
        acts += [{"biomarker": f"{r.gene} {r.protein or r.hgvsc or ''}", "type": "mutation", "tier": r.tier, "VAF": f"{r.t_vaf:.1%}",
                  "therapies (this tumour)": r.therapies_here, "other tumour types": r.therapies_other, "rationale": r.rationale}
                 for r in t.itertuples() if r.tier in ("Tier I", "Tier II")]
    acts += [{"biomarker": f"{r.gene} {r.call.lower()}", "type": "copy number", "tier": r.tier, "VAF": f"CN≈{r.copy_number}",
              "therapies (this tumour)": r.therapies_here, "other tumour types": r.therapies_other, "rationale": r.rationale}
             for r in cnv.itertuples() if r.tier in ("Tier I", "Tier II")]
    if acts:
        st.subheader("Clinically relevant biomarkers")
        st.dataframe(pd.DataFrame(acts), hide_index=True)
        for a in acts:
            if a["tier"] == "Tier I" and a["therapies (this tumour)"]:
                st.success(f"💊 **{a['biomarker']}** → {a['therapies (this tumour)']}")
    else:
        st.info("No Tier I / II biomarkers detected.")
    if len(t) == 0:
        st.warning("No somatic PASS variants to tier.")
        return
    cnt = t.tier.value_counts().reindex(["Tier I", "Tier II", "Tier III", "Tier IV"], fill_value=0).reset_index()
    cnt.columns = ["tier", "n"]
    show(px.bar(cnt, x="tier", y="n", color="tier", color_discrete_map=oncokb.TIER_COLORS, title="Somatic variants by tier"), key="tiers")
    if "annotated" in t and not t.annotated.fillna(False).astype(bool).any():
        st.warning("Variants are not annotated (offline) - protein changes are unknown, so therapies cannot be matched.")
    st.dataframe(t[["gene", "protein", "hgvsc", "consequence", "t_vaf", "tier", "level", "therapies_here", "rationale", "tier_source"]],
                 hide_index=True, column_config={"t_vaf": st.column_config.NumberColumn("VAF", format="%.3f")})

    st.subheader("Review a variant")
    labels = [f"{r.tier} · {r.gene} {r.protein or r.hgvsc or ''} ({r.chrom}:{r.gpos} {r.ref}>{r.alt}) VAF {r.t_vaf:.1%}" for r in t.itertuples()]
    i = st.selectbox("Variant", range(len(labels)), format_func=lambda k: labels[k])
    v = t.iloc[i]
    st.markdown(_badge(v.tier) + f" &nbsp; Level {v.level} · {v.tier_source}", unsafe_allow_html=True)
    a, b = st.columns(2)
    with a:
        st.table(pd.DataFrame({"field": ["Gene", "Protein", "HGVS c.", "Consequence", "Transcript", "Genomic", "COSMIC", "gnomAD AF"],
                               "value": [v.gene, v.protein, v.hgvsc, v.consequence, v.transcript, f"chr{v.chrom}:{v.gpos} {v.ref}>{v.alt}",
                                         v.get("cosmic"), v.gnomad_af]}).astype(str))
    with b:
        st.table(pd.DataFrame({"field": ["Tumour VAF", "Tumour alt / depth", "Normal VAF / depth", "Fisher p", "Status"],
                               "value": [f"{v.t_vaf:.1%}", f"{v.t_alt} / {v.t_dp}",
                                         f"{v.n_vaf:.1%} / {v.n_dp}" if "n_vaf" in v and pd.notna(v.n_vaf) else "no normal",
                                         f"{v.fisher_p:.1e}" if "fisher_p" in v and pd.notna(v.fisher_p) else "-", v.status]}).astype(str))
        if isinstance(v.protein, str):
            st.markdown(f"[OncoKB](https://www.oncokb.org/gene/{v.gene}/{v.protein}) · [CIViC](https://civicdb.org/search?q={v.gene}%20{v.protein}) · "
                        f"[COSMIC](https://cancer.sanger.ac.uk/cosmic/search?q={v.gene}+{v.protein}) · "
                        f"[ClinicalTrials.gov](https://clinicaltrials.gov/search?term={v.gene}%20{v.protein})")
    st.write("**Rationale:**", v.rationale)
    if v.therapies_here:
        st.write("**Therapies in this tumour type:**", v.therapies_here)
    if v.therapies_other:
        st.write("**Evidence in other tumour types:**", v.therapies_other)
    with st.expander("✏️ Override tier (molecular tumour board decision)"):
        opts = ["Tier I", "Tier II", "Tier III", "Tier IV"]
        new = st.selectbox("Tier", opts, index=opts.index(v.tier), key=f"t_{v.variant_id}")
        note = st.text_input("Reason", key=f"n_{v.variant_id}")
        x, y = st.columns(2)
        if x.button("💾 Save", key=f"s_{v.variant_id}"):
            c.overrides[v.variant_id] = {"tier": new, "note": f"Manual: {note}" if note else "Manual override"}
            W.run_step(c, "Tiering (AMP/ASCO/CAP)")
            st.rerun()
        if y.button("↩️ Revert", key=f"r_{v.variant_id}"):
            c.overrides.pop(v.variant_id, None)
            W.run_step(c, "Tiering (AMP/ASCO/CAP)")
            st.rerun()
    st.subheader("🔍 Read evidence — tumour vs normal")
    S = c.results["samples"]
    for role in c.roles:
        fig, info = viewer.pileup_fig(str(S[role]["align"]["bam"]), str(c.ref().fasta), v.contig, int(v.pos), window=35, max_rows=40,
                                      title=f"{role.upper()}: {v.gene} {v.protein or ''} at chr{v.chrom}:{v.gpos}")
        show(fig, key=f"pile_{role}")
    st.caption("A true somatic mutation shows coloured mismatches in the tumour but none in the normal.")
