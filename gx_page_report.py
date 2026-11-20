import io
import zipfile

import streamlit as st

import gx_report_pdf as report_pdf
import gx_workflow as W
from gx_ui_common import need_case


def render():
    c = need_case("ACMG classification")
    st.title("📄 Clinical report")
    st.session_state.setdefault("gx_lab", {"name": "Department of Clinical Genomics", "address": "", "analyst": "", "reviewer": "",
                                        "test_name": "Germline DNA sequencing — targeted gene panel"})
    lab = st.session_state.gx_lab
    with st.expander("🏥 Laboratory details", expanded=True):
        a, b = st.columns(2)
        lab["name"] = a.text_input("Laboratory / hospital name", lab["name"])
        lab["address"] = b.text_input("Address / accreditation", lab["address"])
        lab["test_name"] = a.text_input("Test name", lab["test_name"])
        lab["analyst"] = b.text_input("Analysed by", lab["analyst"])
        lab["reviewer"] = a.text_input("Reviewed / authorised by", lab["reviewer"])
    df = c.results["classified"]
    passing = df[df["filter"] == "PASS"]
    include_vus = st.checkbox("Report variants of uncertain significance", True)
    reportable = passing[passing.classification.isin(["Pathogenic", "Likely pathogenic"])]
    default = ""
    if len(reportable):
        r = reportable.iloc[0]
        default = (f"A {('heterozygous' if r.zygosity == 'het' else 'homozygous')} {r.classification.lower()} variant was detected in "
                   f"{r.gene} ({r.hgvsc or ''} {r.hgvsp or ''}). This result is consistent with the clinical indication. "
                   f"Pathogenic variants in {r.gene} are associated with the relevant inherited condition; clinical correlation and "
                   f"genetic counselling are recommended.")
    interp = st.text_area("Interpretation (edit before signing)", default, height=120)
    lab["recommendation"] = st.text_area("Recommendations", lab.get("recommendation") or
                                         "Genetic counselling is recommended. Cascade testing of at-risk relatives may be offered.", height=80)
    patient = {**c.patient, "sample_id": c.sample}
    pdf = report_pdf.build(patient, lab, passing, W.qc_summary(c), c.results["coverage"]["per_gene"], W.methods_text(c),
                           interp, include_vus)
    st.download_button("⬇️ Download clinical PDF report", pdf, f"{c.sample}_clinical_report.pdf", "application/pdf", type="primary")

    st.subheader("Summary")
    st.write(f"**{len(reportable)}** pathogenic / likely pathogenic · **{(passing.classification == 'Uncertain significance').sum()}** VUS · "
             f"**{passing.classification.isin(['Likely benign', 'Benign']).sum()}** benign / likely benign")
    st.dataframe(passing[["gene", "hgvsc", "hgvsp", "zygosity", "classification", "criteria"]], hide_index=True)

    st.subheader("All result files")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{c.sample}_clinical_report.pdf", pdf)
        z.writestr(f"{c.sample}_variants_classified.tsv", df.to_csv(sep="\t", index=False))
        z.writestr("pipeline_commands.log", "\n".join(c.log))
        for f in [c.results["vcf"], c.results["align"]["bam"]]:
            z.write(f, f.name)
        z.write(str(c.results["align"]["bam"]) + ".bai", c.results["align"]["bam"].name + ".bai")
        z.writestr(f"{c.sample}_coverage_per_exon.tsv", c.results["coverage"]["per_target"].to_csv(sep="\t", index=False))
    st.download_button("⬇️ Download everything (ZIP: PDF, VCF, BAM, TSV, log)", buf.getvalue(), f"{c.sample}_results.zip")
