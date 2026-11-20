"""Automated ACMG/AMP 2015 evidence + ClinGen Bayesian point system (Tavtigian 2020).

Only criteria that can be judged from the data in this app are applied automatically.
Clinical scientists must review and can add/remove criteria in the UI (manual override).
Points: Very strong 8, Strong 4, Moderate 2, Supporting 1 (benign criteria negative).
Total >= 10 Pathogenic, 6-9 Likely pathogenic, 0-5 VUS, -1..-6 Likely benign, <= -7 Benign.
"""
from __future__ import annotations

import math

import gx_reference as R

POINTS = {"VeryStrong": 8, "Strong": 4, "Moderate": 2, "Supporting": 1}
LOF = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost", "transcript_ablation"}

CRITERIA = {  # code: (direction, default strength, description)
    "PVS1": ("P", "VeryStrong", "Null variant (nonsense, frameshift, canonical ±1/2 splice, start loss) in a gene where loss of function causes disease"),
    "PS1": ("P", "Strong", "Same amino-acid change as an established pathogenic variant"),
    "PS2": ("P", "Strong", "De novo (parentage confirmed)"),
    "PS3": ("P", "Strong", "Well-established functional studies show a damaging effect"),
    "PS4": ("P", "Strong", "Prevalence in affected individuals significantly increased vs controls"),
    "PM1": ("P", "Moderate", "Located in a mutational hot spot / critical functional domain"),
    "PM2": ("P", "Supporting", "Absent or extremely rare in population databases (gnomAD) - ClinGen: supporting"),
    "PM3": ("P", "Moderate", "Recessive disorder: detected in trans with a pathogenic variant"),
    "PM4": ("P", "Moderate", "Protein length change (in-frame indel in non-repeat region, stop-loss)"),
    "PM5": ("P", "Moderate", "Novel missense at a residue where a different pathogenic missense is known"),
    "PM6": ("P", "Moderate", "Assumed de novo"),
    "PP1": ("P", "Supporting", "Co-segregation with disease in affected family members"),
    "PP2": ("P", "Supporting", "Missense in a gene with low benign missense variation"),
    "PP3": ("P", "Supporting", "Computational evidence supports a deleterious effect (REVEL / SpliceAI / SIFT+PolyPhen)"),
    "PP4": ("P", "Supporting", "Patient phenotype highly specific for the gene"),
    "PP5": ("P", "Supporting", "Reputable source (ClinVar) reports pathogenic"),
    "BA1": ("B", "StandAlone", "Allele frequency > 5% in gnomAD"),
    "BS1": ("B", "Strong", "Allele frequency greater than expected for the disorder (> 1%)"),
    "BS2": ("B", "Strong", "Observed in healthy adult with full penetrance expected"),
    "BS3": ("B", "Strong", "Functional studies show no damaging effect"),
    "BS4": ("B", "Strong", "Lack of segregation in affected family members"),
    "BP1": ("B", "Supporting", "Missense in a gene where only truncating variants cause disease"),
    "BP2": ("B", "Supporting", "Observed in trans/cis with a pathogenic variant"),
    "BP3": ("B", "Supporting", "In-frame indel in a repetitive region without known function"),
    "BP4": ("B", "Supporting", "Computational evidence suggests no impact"),
    "BP5": ("B", "Supporting", "Found in a case with an alternate molecular basis for disease"),
    "BP6": ("B", "Supporting", "Reputable source (ClinVar) reports benign"),
    "BP7": ("B", "Supporting", "Synonymous variant with no predicted splice impact"),
}


def _num(x):
    try:
        f = float(x)
        return None if math.isnan(f) else f
    except (TypeError, ValueError):
        return None


def auto_criteria(v: dict) -> dict[str, tuple[str, str]]:
    """Return {code: (strength, reason)} for one annotated variant (dict-like row)."""
    c: dict[str, tuple[str, str]] = {}
    cons = set(str(v.get("consequence") or "").split(","))
    gene = v.get("gene")
    af = _num(v.get("gnomad_af"))
    revel, sai = _num(v.get("revel")), _num(v.get("spliceai"))
    am = _num(v.get("alphamissense"))
    clin = str(v.get("clinvar") or "").lower()

    if cons & LOF and gene in R.LOF_GENES:
        c["PVS1"] = ("VeryStrong", f"{', '.join(sorted(cons & LOF))} in {gene}; loss of function is a known disease mechanism")
    annotated = bool(v.get("annotated", True)) and v.get("annotated") == v.get("annotated")  # False/NaN when offline
    if annotated and (af is None or af == 0):
        c["PM2"] = ("Supporting", "Not found in gnomAD")
    elif af is not None and 0 < af < 1e-4:
        c["PM2"] = ("Supporting", f"gnomAD AF {af:.2e} < 0.0001")
    if af is not None and af > .05:
        c["BA1"] = ("StandAlone", f"gnomAD AF {af:.3f} > 5%")
    elif af is not None and af > .01:
        c["BS1"] = ("Strong", f"gnomAD AF {af:.3f} > 1%")
    if cons & {"inframe_insertion", "inframe_deletion", "stop_lost"}:
        c["PM4"] = ("Moderate", ", ".join(sorted(cons & {"inframe_insertion", "inframe_deletion", "stop_lost"})))
    if "missense_variant" in cons:
        if revel is not None:
            if revel >= .932:
                c["PP3"] = ("Strong", f"REVEL {revel:.3f} ≥ 0.932")
            elif revel >= .773:
                c["PP3"] = ("Moderate", f"REVEL {revel:.3f} ≥ 0.773")
            elif revel >= .644:
                c["PP3"] = ("Supporting", f"REVEL {revel:.3f} ≥ 0.644")
            elif revel <= .016:
                c["BP4"] = ("Strong", f"REVEL {revel:.3f} ≤ 0.016")
            elif revel <= .183:
                c["BP4"] = ("Moderate", f"REVEL {revel:.3f} ≤ 0.183")
            elif revel <= .290:
                c["BP4"] = ("Supporting", f"REVEL {revel:.3f} ≤ 0.290")
        elif am is not None:
            if am >= .564:
                c["PP3"] = ("Supporting", f"AlphaMissense {am:.3f} (likely pathogenic)")
            elif am <= .34:
                c["BP4"] = ("Supporting", f"AlphaMissense {am:.3f} (likely benign)")
        else:
            sift, pp = str(v.get("sift") or ""), str(v.get("polyphen") or "")
            if sift.startswith("deleterious") and "damaging" in pp:
                c["PP3"] = ("Supporting", f"SIFT {sift}; PolyPhen {pp}")
            elif sift.startswith("tolerated") and pp.startswith("benign"):
                c["BP4"] = ("Supporting", f"SIFT {sift}; PolyPhen {pp}")
    if sai is not None and "PVS1" not in c:
        if sai >= .2 and "PP3" not in c:
            c["PP3"] = ("Supporting", f"SpliceAI Δ {sai:.2f} ≥ 0.2")
        elif sai <= .1 and "synonymous_variant" in cons:
            c["BP7"] = ("Supporting", f"Synonymous, SpliceAI Δ {sai:.2f} ≤ 0.1")
    if "synonymous_variant" in cons and "BP7" not in c and not cons & {"splice_region_variant", "splice_donor_5th_base_variant"}:
        c["BP7"] = ("Supporting", "Synonymous, outside the splice region")
    if clin:
        if "pathogenic" in clin and "conflict" not in clin and "benign" not in clin:
            c["PP5"] = ("Supporting", f"ClinVar: {v.get('clinvar')}")
        elif "benign" in clin and "pathogenic" not in clin:
            c["BP6"] = ("Supporting", f"ClinVar: {v.get('clinvar')}")
    return c


def classify(crit: dict[str, tuple[str, str]]) -> tuple[str, int]:
    if "BA1" in crit:
        return "Benign", -8
    pts = 0
    for code, (strength, _) in crit.items():
        p = POINTS.get(strength, 0)
        pts += p if CRITERIA.get(code, ("P",))[0] == "P" else -p
    if pts >= 10:
        cls = "Pathogenic"
    elif pts >= 6:
        cls = "Likely pathogenic"
    elif pts >= 0:
        cls = "Uncertain significance"
    elif pts >= -6:
        cls = "Likely benign"
    else:
        cls = "Benign"
    return cls, pts


ORDER = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]
COLORS = {"Pathogenic": "#C0392B", "Likely pathogenic": "#E67E22", "Uncertain significance": "#F1C40F",
          "Likely benign": "#5DADE2", "Benign": "#27AE60"}
SHORT = {"Pathogenic": "P", "Likely pathogenic": "LP", "Uncertain significance": "VUS", "Likely benign": "LB", "Benign": "B"}
