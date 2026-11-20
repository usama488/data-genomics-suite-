"""Molecular oncology (somatic) PDF report."""
from __future__ import annotations

import datetime as dt
import io

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

SS = getSampleStyleSheet()
H1 = ParagraphStyle("h1", parent=SS["Heading1"], fontSize=16, textColor=colors.HexColor("#6C3483"), spaceAfter=4)
H2 = ParagraphStyle("h2", parent=SS["Heading2"], fontSize=12, textColor=colors.HexColor("#6C3483"), spaceBefore=8, spaceAfter=4)
BODY = ParagraphStyle("b", parent=SS["BodyText"], fontSize=9, leading=12)
SMALL = ParagraphStyle("s", parent=BODY, fontSize=7.5, leading=9.5, textColor=colors.HexColor("#555555"))
CELL = ParagraphStyle("c", parent=BODY, fontSize=7.5, leading=9)
CELLH = ParagraphStyle("ch", parent=CELL, textColor=colors.white, fontName="Helvetica-Bold")
BANNER = ParagraphStyle("ban", parent=BODY, fontSize=11.5, leading=15, alignment=TA_CENTER, textColor=colors.white)
TIER_BG = {"Tier I": "#C0392B", "Tier II": "#E67E22", "Tier III": "#B7950B", "Tier IV": "#27AE60"}


def _p(x, st=CELL):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        x = ""
    return Paragraph(str(x), st)


def _table(rows, widths, head="#6C3483"):
    rows = [[Paragraph(str(h), CELLH) for h in rows[0]]] + rows[1:]
    t = Table(rows, colWidths=widths, repeatRows=1)
    st = [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(head)), ("FONTSIZE", (0, 0), (-1, -1), 7.5),
          ("VALIGN", (0, 0), (-1, -1), "TOP"), ("GRID", (0, 0), (-1, -1), .25, colors.HexColor("#BBBBBB")),
          ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    for i in range(2, len(rows), 2):
        st.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#F5EEF8")))
    t.setStyle(TableStyle(st))
    return t


def build(case, lab: dict, tiered: pd.DataFrame, cnv: pd.DataFrame, germline: pd.DataFrame, qc: dict, methods: dict,
          interpretation: str = "", include_tier3=True) -> bytes:
    p = case.patient
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=14 * mm, rightMargin=14 * mm, topMargin=14 * mm, bottomMargin=16 * mm,
                            title=f"Molecular oncology report {case.case_id}")

    def footer(c, d):
        c.saveState()
        c.setFont("Helvetica", 7)
        c.setFillColor(colors.HexColor("#777777"))
        c.drawString(14 * mm, 9 * mm, f"{lab.get('name', '')} · Case {case.case_id} · {dt.date.today()}")
        c.drawRightString(A4[0] - 14 * mm, 9 * mm, f"Page {d.page}")
        c.restoreState()

    el = [Paragraph(lab.get("name") or "Molecular Pathology Laboratory", H1), Paragraph(lab.get("address", ""), SMALL),
          Paragraph(f"<b>{lab.get('test_name') or 'Somatic tumour profiling - targeted NGS panel'}</b>", BODY), Spacer(1, 4)]
    info = [["Patient", p.get("name", ""), "Case ID", case.case_id],
            ["DOB / Sex", f"{p.get('dob', '')} / {p.get('sex', '')}", "MRN", p.get("mrn", "")],
            ["Tumour type", case.tumour_type, "Specimen", p.get("specimen", "")],
            ["Tumour cellularity", f"{case.pathology_purity * 100:.0f}% (pathology)" if case.pathology_purity else "not provided",
             "Analysis", "Tumour vs matched normal" if case.mode == "paired" else "Tumour-only"],
            ["Referring clinician", p.get("clinician", ""), "Block / sample", p.get("block", "")]]
    t = Table(info, colWidths=[32 * mm, 60 * mm, 30 * mm, 60 * mm])
    t.setStyle(TableStyle([("FONTSIZE", (0, 0), (-1, -1), 8), ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#F5EEF8")),
                           ("BACKGROUND", (2, 0), (2, -1), colors.HexColor("#F5EEF8")), ("GRID", (0, 0), (-1, -1), .25, colors.grey)]))
    el.append(t)

    t1 = tiered[tiered.tier == "Tier I"] if len(tiered) else tiered
    t2 = tiered[tiered.tier == "Tier II"] if len(tiered) else tiered
    t3 = tiered[tiered.tier == "Tier III"] if len(tiered) else tiered
    c1 = cnv[cnv.tier == "Tier I"] if len(cnv) else cnv
    c2 = cnv[cnv.tier.isin(["Tier II", "Tier III"])] if len(cnv) else cnv
    n_t1 = len(t1) + len(c1)
    if n_t1:
        txt, col = f"{n_t1} biomarker(s) with approved / guideline therapy in {case.tumour_type} (Tier I)", "#C0392B"
    elif len(t2) + len(c2):
        txt, col = f"No Tier I biomarkers; {len(t2) + len(c2)} alteration(s) of potential significance (Tier II)", "#CA6F1E"
    else:
        txt, col = "No clinically actionable alterations detected", "#1E8449"
    el.append(Spacer(1, 8))
    b = Table([[Paragraph(f"<b>SUMMARY - {txt}</b>", BANNER)]], colWidths=[182 * mm])
    b.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(col)), ("TOPPADDING", (0, 0), (-1, -1), 6),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
    el.append(b)

    def vtab(df):
        rows = [["Gene", "Alteration", "VAF / depth", "Tier", "Therapies (this tumour type)", "Other evidence"]]
        for r in df.itertuples():
            alt = f"<b>{r.protein or ''}</b><br/>{r.hgvsc or ''}<br/>chr{r.chrom}:{r.gpos} {r.ref}>{r.alt}"
            rows.append([_p(f"<b>{r.gene}</b>"), _p(alt), _p(f"{r.t_vaf:.1%}<br/>{int(r.t_dp)}x"), _p(f"<b>{r.tier}</b> ({r.level})"),
                         _p(r.therapies_here or "-"), _p((r.therapies_other or "") + ("<br/>" if r.therapies_other else "") + (r.rationale or ""))])
        return _table(rows, [16 * mm, 38 * mm, 18 * mm, 18 * mm, 46 * mm, 46 * mm])

    el.append(Paragraph("Tier I - strong clinical significance", H2))
    if len(t1):
        el.append(vtab(t1))
    if len(c1):
        el.append(Spacer(1, 3))
        el.append(_table([["Gene", "Copy-number change", "Est. copies", "Therapies (this tumour type)"]] +
                         [[_p(f"<b>{r.gene}</b>"), _p(f"{r.call} (log2 {r.log2_ratio:.2f})"), _p(r.copy_number), _p(r.therapies_here)]
                          for r in c1.itertuples()], [20 * mm, 50 * mm, 22 * mm, 90 * mm]))
    if not n_t1:
        el.append(Paragraph("None detected.", BODY))
    el.append(Paragraph("Tier II - potential clinical significance", H2))
    if len(t2):
        el.append(vtab(t2))
    if len(c2):
        el.append(_table([["Gene", "Copy-number change", "Est. copies", "Evidence"]] +
                         [[_p(r.gene), _p(f"{r.call} (log2 {r.log2_ratio:.2f})"), _p(r.copy_number), _p((r.therapies_other or "") + " " + (r.rationale or ""))]
                          for r in c2.itertuples()], [20 * mm, 50 * mm, 22 * mm, 90 * mm]))
    if not len(t2) and not len(c2):
        el.append(Paragraph("None detected.", BODY))
    if include_tier3:
        el.append(Paragraph("Tier III - variants of unknown significance", H2))
        el.append(vtab(t3) if len(t3) else Paragraph("None detected.", BODY))

    el.append(Paragraph("Genomic signatures", H2))
    tm = case.results.get("tmb", {})
    el.append(Paragraph(f"Tumour mutational burden: <b>{tm.get('tmb', 0):.1f} mutations/Mb</b> ({tm.get('mutations', 0)} {tm.get('basis', '')} "
                        f"mutations over {tm.get('megabases', 0):.3f} Mb). " + ("" if tm.get("reliable") else
                        "<i>Not reportable: the panel covers less than 1 Mb (validated TMB needs ≥ 1 Mb, ideally ≥ 1.5 Mb).</i>"), BODY))
    if germline is not None and len(germline):
        el.append(Paragraph("Germline findings (matched normal)", H2))
        el.append(Paragraph(f"{len(germline)} inherited variant(s) were detected in the normal sample and excluded from somatic analysis. "
                            "Pathogenic germline variants in cancer-predisposition genes (e.g. TP53, BRCA1/2) require referral to clinical "
                            "genetics with patient consent; they are not interpreted in this report.", BODY))
    if interpretation:
        el.append(Paragraph("Interpretation", H2))
        el.append(Paragraph(interpretation.replace("\n", "<br/>"), BODY))
    el.append(Paragraph("Quality metrics", H2))
    items = list(qc.items())
    rows = [["Metric", "Value", "Metric", "Value"]]
    for i in range(0, len(items), 2):
        a = items[i]
        bb = items[i + 1] if i + 1 < len(items) else ("", "")
        rows.append([_p(a[0]), _p(a[1]), _p(bb[0]), _p(bb[1])])
    el.append(_table(rows, [52 * mm, 39 * mm, 52 * mm, 39 * mm]))
    el.append(Paragraph("Methodology", H2))
    el.append(Paragraph(methods.get("text", "") + f" Genes analysed: {methods.get('genes', '')}.", BODY))
    el.append(Paragraph("Limitations", H2))
    el.append(Paragraph("Only the coding exons of the listed genes were analysed. Gene fusions, structural variants, MSI and variants below "
                        "the limit of detection (~5% VAF, depending on depth and tumour content) are not assessed. Copy-number estimates "
                        "depend on tumour purity and should be confirmed (e.g. ERBB2 by IHC/ISH). Therapy associations reflect a curated "
                        "knowledge-base snapshot and regulatory approvals that may differ between countries; consult current guidelines "
                        "and a molecular tumour board.", BODY))
    el.append(Paragraph("<b>Disclaimer:</b> produced by an automated research pipeline; must be reviewed and authorised by a qualified "
                        "molecular pathologist / clinical scientist before clinical use.", SMALL))
    el.append(Spacer(1, 12))
    sig = Table([["Analysed by:", lab.get("analyst", ""), "Authorised by:", lab.get("reviewer", "")], ["Signature:", "", "Signature:", ""],
                 ["Date:", str(dt.date.today()), "Date:", ""]], colWidths=[28 * mm, 60 * mm, 30 * mm, 64 * mm])
    sig.setStyle(TableStyle([("FONTSIZE", (0, 0), (-1, -1), 8), ("LINEBELOW", (1, 1), (1, 1), .5, colors.black),
                             ("LINEBELOW", (3, 1), (3, 1), .5, colors.black), ("TOPPADDING", (0, 0), (-1, -1), 5)]))
    el.append(KeepTogether(sig))
    if len(tiered):
        el.append(PageBreak())
        el.append(Paragraph("Appendix - all somatic variants (PASS)", H2))
        rows = [["Gene", "Protein / HGVS c.", "GRCh38", "VAF", "Tumour DP", "Normal VAF", "COSMIC", "Tier"]]
        for r in tiered.itertuples():
            nv = getattr(r, "n_vaf", None)
            rows.append([_p(r.gene), _p(f"{r.protein or ''} {r.hgvsc or ''}"), _p(f"{r.chrom}:{r.gpos} {r.ref}>{r.alt}"), f"{r.t_vaf:.1%}",
                         int(r.t_dp), "-" if nv is None or pd.isna(nv) else f"{nv:.1%}", _p(getattr(r, "cosmic", "") or ""), r.tier])
        el.append(_table(rows, [16 * mm, 40 * mm, 38 * mm, 14 * mm, 16 * mm, 16 * mm, 24 * mm, 18 * mm]))
        el.append(Paragraph(methods.get("versions", ""), SMALL))
    doc.build(el, onFirstPage=footer, onLaterPages=footer)
    return buf.getvalue()
