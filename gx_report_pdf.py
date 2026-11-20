"""Clinical-style PDF report (reportlab)."""
from __future__ import annotations

import datetime as dt
import io

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle)

import gx_acmg as acmg

SS = getSampleStyleSheet()
H1 = ParagraphStyle("h1", parent=SS["Heading1"], fontSize=16, textColor=colors.HexColor("#1B4F72"), spaceAfter=4)
H2 = ParagraphStyle("h2", parent=SS["Heading2"], fontSize=12, textColor=colors.HexColor("#1B4F72"), spaceBefore=8, spaceAfter=4)
BODY = ParagraphStyle("b", parent=SS["BodyText"], fontSize=9, leading=12)
SMALL = ParagraphStyle("s", parent=BODY, fontSize=7.5, leading=9.5, textColor=colors.HexColor("#555555"))
CELL = ParagraphStyle("c", parent=BODY, fontSize=7.5, leading=9)
CELLH = ParagraphStyle("ch", parent=BODY, fontSize=7.5, leading=9, textColor=colors.white, fontName="Helvetica-Bold")
BANNER = ParagraphStyle("ban", parent=BODY, fontSize=12, leading=15, alignment=TA_CENTER, textColor=colors.white)


def _p(x, st=CELL):
    return Paragraph("" if x is None or (isinstance(x, float) and pd.isna(x)) else str(x), st)


def _table(data, widths, header_bg="#1B4F72", zebra=True):
    data = [[Paragraph(str(h), CELLH) for h in data[0]]] + data[1:]
    t = Table(data, colWidths=widths, repeatRows=1)
    style = [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(header_bg)), ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
             ("FONTSIZE", (0, 0), (-1, -1), 7.5), ("VALIGN", (0, 0), (-1, -1), "TOP"),
             ("GRID", (0, 0), (-1, -1), .25, colors.HexColor("#BBBBBB")), ("TOPPADDING", (0, 0), (-1, -1), 2),
             ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
    if zebra:
        for i in range(1, len(data)):
            if i % 2 == 0:
                style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#F4F6F7")))
    t.setStyle(TableStyle(style))
    return t


def _fmt_af(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "n/a"
    return "absent" if x == 0 else f"{x:.2e}" if x < .001 else f"{x:.4f}"


def build(patient: dict, lab: dict, variants: pd.DataFrame, qc: dict, coverage_gene: pd.DataFrame, methods: dict,
          interpretation: str = "", include_vus=True) -> bytes:
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=14 * mm, rightMargin=14 * mm, topMargin=14 * mm, bottomMargin=16 * mm,
                            title=f"Genomic report {patient.get('sample_id', '')}")

    def footer(c, d):
        c.saveState()
        c.setFont("Helvetica", 7)
        c.setFillColor(colors.HexColor("#777777"))
        c.drawString(14 * mm, 9 * mm, f"{lab.get('name', '')} · Sample {patient.get('sample_id', '')} · Report generated {dt.date.today()}")
        c.drawRightString(A4[0] - 14 * mm, 9 * mm, f"Page {d.page}")
        c.restoreState()

    el = []
    el.append(Paragraph(lab.get("name") or "Genomics Laboratory", H1))
    el.append(Paragraph(lab.get("address", ""), SMALL))
    el.append(Paragraph(f"<b>{lab.get('test_name') or 'Germline DNA sequencing - targeted gene panel'}</b>", BODY))
    el.append(Spacer(1, 4))
    pinfo = [["Patient name", patient.get("name", ""), "Sample ID", patient.get("sample_id", "")],
             ["Date of birth", patient.get("dob", ""), "Sex", patient.get("sex", "")],
             ["MRN / Hospital no.", patient.get("mrn", ""), "Specimen", patient.get("specimen", "")],
             ["Referring clinician", patient.get("clinician", ""), "Collected / received", f"{patient.get('collected', '')} / {patient.get('received', '')}"],
             ["Clinical indication", Paragraph(patient.get("indication", ""), CELL), "Genes analysed", Paragraph(methods.get("genes", ""), CELL)]]
    t = Table(pinfo, colWidths=[32 * mm, 58 * mm, 34 * mm, 58 * mm])
    t.setStyle(TableStyle([("FONTSIZE", (0, 0), (-1, -1), 8), ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#EAF2F8")),
                           ("BACKGROUND", (2, 0), (2, -1), colors.HexColor("#EAF2F8")), ("GRID", (0, 0), (-1, -1), .25, colors.grey),
                           ("VALIGN", (0, 0), (-1, -1), "TOP")]))
    el.append(t)

    # ---- result banner
    reportable = variants[variants["classification"].isin(["Pathogenic", "Likely pathogenic"])] if len(variants) else variants
    vus = variants[variants["classification"] == "Uncertain significance"] if len(variants) else variants
    if len(reportable):
        res_txt, col = f"POSITIVE: {len(reportable)} pathogenic / likely pathogenic variant(s) identified", "#C0392B"
    elif len(vus) and include_vus:
        res_txt, col = f"INCONCLUSIVE: {len(vus)} variant(s) of uncertain significance identified", "#B7950B"
    else:
        res_txt, col = "NEGATIVE: no pathogenic or likely pathogenic variants identified", "#1E8449"
    el.append(Spacer(1, 8))
    ban = Table([[Paragraph(f"<b>RESULT - {res_txt}</b>", BANNER)]], colWidths=[182 * mm])
    ban.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(col)), ("TOPPADDING", (0, 0), (-1, -1), 6),
                             ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
    el.append(ban)

    def vtable(df):
        head = ["Gene", "Variant (HGVS)", "Genomic (GRCh38)", "Zygosity", "Consequence", "gnomAD AF", "ClinVar", "Classification"]
        rows = [head]
        for r in df.itertuples():
            hg = f"{r.transcript or ''}<br/>{r.hgvsc or ''}<br/>{r.hgvsp or ''}"
            rows.append([_p(f"<b>{r.gene}</b>"), _p(hg), _p(f"chr{r.chrom}:{r.gpos}<br/>{r.ref}>{r.alt}"),
                         _p(f"{'Heterozygous' if r.zygosity == 'het' else 'Homozygous'}<br/>VAF {r.vaf}"),
                         _p((r.consequence or "").replace("_", " ").replace(",", ", ")), _p(_fmt_af(r.gnomad_af)),
                         _p(r.clinvar or "not reported"), _p(f"<b>{r.classification}</b><br/>{r.criteria}")])
        return _table(rows, [16 * mm, 38 * mm, 25 * mm, 20 * mm, 22 * mm, 15 * mm, 20 * mm, 26 * mm])

    el.append(Paragraph("Primary findings", H2))
    if len(reportable):
        el.append(vtable(reportable))
    else:
        el.append(Paragraph("No pathogenic or likely pathogenic variants were detected in the genes analysed.", BODY))
    if include_vus:
        el.append(Paragraph("Variants of uncertain significance", H2))
        el.append(vtable(vus) if len(vus) else Paragraph("None detected.", BODY))
        if len(vus):
            el.append(Paragraph("VUS should not be used for clinical decision-making. They may be reclassified as evidence accumulates.", SMALL))

    if interpretation:
        el.append(Paragraph("Interpretation", H2))
        el.append(Paragraph(interpretation.replace("\n", "<br/>"), BODY))
    el.append(Paragraph("Recommendations", H2))
    el.append(Paragraph(lab.get("recommendation") or "Genetic counselling is recommended to discuss these results, their implications "
                        "for the patient and at-risk relatives, and options for cascade testing.", BODY))

    # ---- QC
    el.append(Paragraph("Quality metrics", H2))
    qrows = [["Metric", "Value", "Metric", "Value"]]
    items = list(qc.items())
    for i in range(0, len(items), 2):
        a = items[i]
        b = items[i + 1] if i + 1 < len(items) else ("", "")
        qrows.append([a[0], a[1], b[0], b[1]])
    el.append(_table(qrows, [50 * mm, 41 * mm, 50 * mm, 41 * mm]))
    if coverage_gene is not None and len(coverage_gene):
        el.append(Spacer(1, 4))
        crow = [["Gene", "Target bp", "Mean depth", "% ≥ 20x", "% ≥ 30x"]]
        for r in coverage_gene.itertuples():
            crow.append([r.gene, f"{int(r.bp):,}", f"{r.mean_depth:.0f}x", f"{r.pct_20x:.1f}%", f"{r.pct_30x:.1f}%"])
        el.append(_table(crow, [36 * mm, 36 * mm, 36 * mm, 37 * mm, 37 * mm]))

    # ---- methods & limitations
    el.append(Paragraph("Methodology", H2))
    el.append(Paragraph(methods.get("text", ""), BODY))
    el.append(Paragraph("Limitations", H2))
    el.append(Paragraph("This test analyses only the coding exons (±20 bp) of the listed genes. It does not reliably detect large "
                        "deletions/duplications (CNVs), structural variants, repeat expansions, deep intronic variants, low-level mosaicism, "
                        "or variants in regions with insufficient coverage (&lt;20x). Variant classification follows ACMG/AMP 2015 with the "
                        "ClinGen Bayesian point system and reflects knowledge on the report date.", BODY))
    el.append(Paragraph("<b>Disclaimer:</b> generated by an automated analysis pipeline. Results must be reviewed and signed by a "
                        "qualified clinical scientist before use in patient care, and confirmed by an orthogonal method where required.", SMALL))
    el.append(Spacer(1, 14))
    sig = Table([["Analysed by:", lab.get("analyst", ""), "Reviewed / authorised by:", lab.get("reviewer", "")],
                 ["Signature:", "", "Signature:", ""], ["Date:", str(dt.date.today()), "Date:", ""]],
                colWidths=[28 * mm, 60 * mm, 38 * mm, 56 * mm])
    sig.setStyle(TableStyle([("FONTSIZE", (0, 0), (-1, -1), 8), ("LINEBELOW", (1, 1), (1, 1), .5, colors.black),
                             ("LINEBELOW", (3, 1), (3, 1), .5, colors.black), ("TOPPADDING", (0, 0), (-1, -1), 5)]))
    el.append(KeepTogether(sig))

    # ---- appendix
    if len(variants):
        el.append(PageBreak())
        el.append(Paragraph("Appendix - all variants passing quality filters", H2))
        rows = [["Gene", "HGVS c. / p.", "GRCh38", "GT", "DP", "VAF", "gnomAD", "Class"]]
        for r in variants.itertuples():
            rows.append([_p(r.gene), _p(f"{r.hgvsc or ''} {r.hgvsp or ''}"), _p(f"{r.chrom}:{r.gpos} {r.ref}>{r.alt}"),
                         r.gt, r.dp, r.vaf, _fmt_af(r.gnomad_af), acmg.SHORT.get(r.classification, r.classification)])
        el.append(_table(rows, [18 * mm, 48 * mm, 40 * mm, 12 * mm, 12 * mm, 14 * mm, 20 * mm, 18 * mm]))
        el.append(Paragraph(methods.get("versions", ""), SMALL))
    doc.build(el, onFirstPage=footer, onLaterPages=footer)
    return buf.getvalue()
