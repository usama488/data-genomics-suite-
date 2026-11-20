import io
import zipfile

import streamlit as st

import cx_report_pdf as report_pdf
import cx_workflow as W
from cx_ui_common import need_case


def render():
    c = need_case("Tiering (AMP/ASCO/CAP)")
    st.title("📄 Molecular oncology report")
    st.session_state.setdefault("cx_lab", {"name": "Molecular Pathology Laboratory", "address": "", "analyst": "", "reviewer": "",
                                        "test_name": "Somatic tumour profiling — targeted NGS panel"})
    lab = st.session_state.cx_lab
    with st.expander("🏥 Laboratory details", expanded=True):
        a, b = st.columns(2)
        lab["name"] = a.text_input("Laboratory / hospital", lab["name"])
        lab["address"] = b.text_input("Address / accreditation", lab["address"])
        lab["test_name"] = a.text_input("Test name", lab["test_name"])
        lab["analyst"] = b.text_input("Analysed by", lab["analyst"])
        lab["reviewer"] = a.text_input("Authorised by", lab["reviewer"])
    t, cnv = c.results["tiered"], c.results["cnv_tiered"]
    t1 = t[t.tier == "Tier I"] if len(t) else t
    lines = []
    for r in t1.itertuples():
        lines.append(f"{r.gene} {r.protein or r.hgvsc} (VAF {r.t_vaf:.0%}) is a Tier I biomarker in {c.tumour_type}; "
                     f"associated therapies: {r.therapies_here}.")
    for r in cnv[cnv.tier == "Tier I"].itertuples():
        lines.append(f"{r.gene} amplification (~{r.copy_number} copies) - {r.therapies_here}. Confirm by IHC/ISH.")
    if not lines:
        lines.append("No Tier I biomarkers were identified. Consider clinical-trial eligibility for Tier II findings.")
    interp = st.text_area("Interpretation (edit before authorising)", "\n".join(lines), height=140)
    include3 = st.checkbox("Include Tier III variants", True)
    an = c.results["annotated"]
    germ = an[an.status.str.contains("germline")] if len(an) else an
    pdf = report_pdf.build(c, lab, t, cnv, germ if c.mode == "paired" else None, W.qc_summary(c), W.methods(c), interp, include3)
    st.download_button("⬇️ Download oncology PDF report", pdf, f"{c.case_id}_oncology_report.pdf", "application/pdf", type="primary")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{c.case_id}_oncology_report.pdf", pdf)
        if len(t):
            z.writestr(f"{c.case_id}_somatic_tiered.tsv", t.to_csv(sep="\t", index=False))
        z.writestr(f"{c.case_id}_cnv.tsv", cnv.to_csv(sep="\t", index=False))
        z.writestr(f"{c.case_id}_all_calls.tsv", an.to_csv(sep="\t", index=False))
        z.writestr("pipeline_commands.log", "\n".join(c.log))
        if c.results.get("vcf"):
            z.write(c.results["vcf"], f"{c.case_id}.somatic.vcf.gz")
        for role, d in c.results["samples"].items():
            bam = d["align"]["bam"]
            z.write(bam, f"{role}.bam")
            z.write(str(bam) + ".bai", f"{role}.bam.bai")
    st.download_button("⬇️ Everything (ZIP: PDF, VCF, BAMs, TSVs, log)", buf.getvalue(), f"{c.case_id}_results.zip")
