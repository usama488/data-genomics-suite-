"""Clinical actionability knowledge base + AMP/ASCO/CAP 2017 somatic tiering.

The knowledge base is a small, curated snapshot of well-established, regulator-approved
biomarker–therapy associations (FDA/EMA labels and NCCN/ESMO guidelines, as of mid-2026).
It is NOT exhaustive and must be checked against OncoKB / CIViC / current guidelines before use.

Tiers (Li et al., J Mol Diagn 2017):
  Tier I   - strong clinical significance (Level A: approved therapy / guideline in THIS tumour type)
  Tier II  - potential clinical significance (Level C: approved in a different tumour type,
             Level D: clinical-trial / preclinical evidence, known oncogenic hotspot)
  Tier III - unknown clinical significance (VUS)
  Tier IV  - benign / likely benign
"""
from __future__ import annotations

import re

import pandas as pd

TUMOR_TYPES = ["Non-small cell lung cancer", "Colorectal cancer", "Melanoma", "Breast cancer", "Anaplastic thyroid cancer",
               "Gastric / GEJ cancer", "Cholangiocarcinoma", "Glioma", "Acute myeloid leukaemia", "Myelodysplastic syndrome",
               "Systemic mastocytosis", "GIST", "Pancreatic cancer", "Other solid tumour"]
SOLID = set(TUMOR_TYPES) - {"Acute myeloid leukaemia", "Myelodysplastic syndrome", "Systemic mastocytosis"}

# (gene, alteration pattern) -> list of (tumour types | {"ANY_SOLID"}, therapies, effect, evidence note)
KB: dict[tuple[str, str], list[tuple[set, str, str, str]]] = {
    ("BRAF", "V600E"): [
        ({"Melanoma"}, "dabrafenib + trametinib; vemurafenib + cobimetinib; encorafenib + binimetinib", "sensitive", "FDA approved"),
        ({"Non-small cell lung cancer"}, "dabrafenib + trametinib; encorafenib + binimetinib", "sensitive", "FDA approved"),
        ({"Colorectal cancer"}, "encorafenib + cetuximab (± mFOLFOX6)", "sensitive", "FDA approved"),
        ({"Anaplastic thyroid cancer"}, "dabrafenib + trametinib", "sensitive", "FDA approved"),
        ({"ANY_SOLID"}, "dabrafenib + trametinib (tissue-agnostic, unresectable/metastatic solid tumours)", "sensitive", "FDA approved (tumour-agnostic)"),
    ],
    ("KRAS", "G12C"): [
        ({"Non-small cell lung cancer"}, "sotorasib; adagrasib", "sensitive", "FDA approved"),
        ({"Colorectal cancer"}, "adagrasib + cetuximab; sotorasib + panitumumab", "sensitive", "FDA approved"),
    ],
    ("KRAS", "RAS_HOTSPOT"): [
        ({"Colorectal cancer"}, "cetuximab; panitumumab (anti-EGFR antibodies)", "resistant", "Guideline: RAS-mutant CRC should not receive anti-EGFR therapy"),
    ],
    ("NRAS", "RAS_HOTSPOT"): [
        ({"Colorectal cancer"}, "cetuximab; panitumumab (anti-EGFR antibodies)", "resistant", "Guideline: RAS-mutant CRC should not receive anti-EGFR therapy"),
    ],
    ("EGFR", "SENSITISING"): [
        ({"Non-small cell lung cancer"}, "osimertinib (± chemotherapy); amivantamab + lazertinib; erlotinib; gefitinib; afatinib; dacomitinib",
         "sensitive", "FDA approved"),
    ],
    ("EGFR", "T790M"): [
        ({"Non-small cell lung cancer"}, "osimertinib", "sensitive", "FDA approved"),
        ({"Non-small cell lung cancer"}, "erlotinib; gefitinib; afatinib (1st/2nd-generation EGFR TKIs)", "resistant", "Acquired resistance mutation"),
    ],
    ("PIK3CA", "PIK3CA_HOTSPOT"): [
        ({"Breast cancer"}, "alpelisib + fulvestrant; inavolisib + palbociclib + fulvestrant; capivasertib + fulvestrant (HR+/HER2-)",
         "sensitive", "FDA approved"),
    ],
    ("IDH1", "R132"): [
        ({"Acute myeloid leukaemia"}, "ivosidenib; olutasidenib", "sensitive", "FDA approved"),
        ({"Myelodysplastic syndrome"}, "ivosidenib (relapsed/refractory)", "sensitive", "FDA approved"),
        ({"Cholangiocarcinoma"}, "ivosidenib", "sensitive", "FDA approved"),
        ({"Glioma"}, "vorasidenib (grade 2 astrocytoma / oligodendroglioma)", "sensitive", "FDA approved"),
    ],
    ("KIT", "D816V"): [
        ({"Systemic mastocytosis"}, "avapritinib; midostaurin", "sensitive", "FDA approved"),
        ({"GIST"}, "imatinib", "resistant", "D816V confers imatinib resistance"),
    ],
    ("ERBB2", "AMPLIFICATION"): [
        ({"Breast cancer"}, "trastuzumab; pertuzumab; trastuzumab emtansine; trastuzumab deruxtecan; tucatinib; neratinib; lapatinib",
         "sensitive", "FDA approved (HER2-positive; confirm by IHC/ISH)"),
        ({"Gastric / GEJ cancer"}, "trastuzumab + chemotherapy (± pembrolizumab); trastuzumab deruxtecan", "sensitive", "FDA approved (HER2-positive)"),
        ({"Colorectal cancer"}, "tucatinib + trastuzumab", "sensitive", "FDA approved (HER2-positive, RAS wild-type)"),
    ],
}

KNOWN_ONCOGENIC = {  # well-characterised hotspots without a matching approved therapy here
    ("TP53", r"R(175|248|249|273|282)\w|G245\w|Y220C"): "TP53 DNA-binding domain hotspot; prognostic in several cancers; Y220C has trials of reactivators",
    ("KRAS", r"G1[23]\w|Q61\w|A146\w|K117\w"): "KRAS hotspot (MAPK activation); KRAS G12D / pan-RAS inhibitors in clinical trials",
    ("NRAS", r"Q61\w|G1[23]\w"): "NRAS hotspot (MAPK activation); MEK inhibitor trials in melanoma",
    ("BRAF", r"V600[KDR]"): "BRAF V600 non-E mutation; BRAF/MEK inhibitors used per guidelines in melanoma",
    ("EGFR", r"G719\w|L861Q|S768I"): "Uncommon EGFR sensitising mutation (afatinib/osimertinib data)",
    ("KIT", r"exon11|V560\w|W557\w"): "KIT exon 11 mutation (imatinib-sensitive in GIST)",
}
TSG = {"TP53", "PTEN", "RB1", "NF1", "APC", "CDKN2A", "BRCA1", "BRCA2", "STK11", "SMAD4", "ARID1A", "ATM"}
TRUNC = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
AA = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
      "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
      "Ter": "*", "Sec": "U"}


def short_protein(hgvsp: str | None) -> str | None:
    """p.Gly12Cys -> G12C ; p.Glu746_Ala750del -> E746_A750del."""
    if not isinstance(hgvsp, str) or not hgvsp.startswith("p."):
        return None
    s = hgvsp[2:].replace("(", "").replace(")", "")
    for k, v in AA.items():
        s = s.replace(k, v)
    return s.replace("%3D", "=")


def alteration_keys(gene: str, prot: str | None, consequence: str, exon: str | None) -> list[str]:
    keys = []
    if prot:
        keys.append(prot)
    cons = set(str(consequence or "").split(","))
    exon_no = str(exon or "").split("/")[0]
    if gene == "EGFR" and (prot in ("L858R",) or (exon_no == "19" and "inframe_deletion" in cons)):
        keys.append("SENSITISING")
    if gene in ("KRAS", "NRAS") and prot and re.match(r"(G12|G13|Q61|K117|A146|A59)", prot):
        keys.append("RAS_HOTSPOT")
    if gene == "PIK3CA" and prot and re.match(r"(E542K|E545[KAGDQ]|H1047[RLY]|C420R|Q546[EKPR]|E726K|N345K|G1049R|M1043[IV])", prot):
        keys.append("PIK3CA_HOTSPOT")
    if gene == "IDH1" and prot and prot.startswith("R132"):
        keys.append("R132")
    return keys


def evidence_for(gene: str, keys: list[str]) -> list[dict]:
    out = []
    for k in keys:
        for types, drugs, effect, note in KB.get((gene, k), []):
            out.append({"alteration": k, "tumour_types": types, "therapies": drugs, "effect": effect, "note": note})
    return out


def _matches(types: set, tumour: str) -> bool:
    return tumour in types or ("ANY_SOLID" in types and tumour in SOLID)


def tier_variant(v: dict, tumour: str) -> dict:
    gene = v.get("gene")
    cons = str(v.get("consequence") or "")
    prot = short_protein(v.get("hgvsp"))
    af = v.get("gnomad_af")
    try:
        af = float(af)
    except (TypeError, ValueError):
        af = None
    keys = alteration_keys(gene, prot, cons, v.get("exon"))
    ev = evidence_for(gene, keys)
    here = [e for e in ev if _matches(e["tumour_types"], tumour)]
    other = [e for e in ev if not _matches(e["tumour_types"], tumour)]
    other_sens = [e for e in other if e["effect"] == "sensitive"]
    res = {"protein": prot, "tier": None, "level": None, "therapies_here": "; ".join(f"{e['therapies']} [{e['effect']}]" for e in here),
           "therapies_other": "; ".join(f"{', '.join(sorted(e['tumour_types'] - {'ANY_SOLID'}) or ['tumour-agnostic'])}: {e['therapies']} [{e['effect']}]"
                                        for e in other), "rationale": ""}
    coding = any(c in cons for c in ("missense", "stop_", "frameshift", "inframe", "splice_acceptor", "splice_donor", "start_lost", "protein_altering"))
    if af is not None and af >= .01:
        res.update(tier="Tier IV", level="Benign", rationale=f"Common in the population (gnomAD AF {af:.3f}) - likely germline polymorphism")
    elif here:
        res.update(tier="Tier I", level="A", rationale=f"Approved/guideline therapy association in {tumour}")
    elif other_sens:
        res.update(tier="Tier II", level="C", rationale="Therapy approved for this biomarker in a different tumour type")
    elif any(k[0] == gene and prot and re.fullmatch(k[1], prot) for k in KNOWN_ONCOGENIC):
        note = next(n for k, n in KNOWN_ONCOGENIC.items() if k[0] == gene and prot and re.fullmatch(k[1], prot))
        res.update(tier="Tier II", level="D", rationale=note)
    elif gene in TSG and any(t in cons for t in TRUNC):
        res.update(tier="Tier II", level="D", rationale=f"Truncating variant in tumour suppressor {gene} - likely oncogenic (loss of function)")
    elif coding:
        res.update(tier="Tier III", level="VUS", rationale="Coding change of unknown significance in a cancer gene")
    elif cons:
        res.update(tier="Tier IV", level="Benign", rationale="Non-coding / synonymous - no expected protein impact")
    else:
        res.update(tier="Tier III", level="VUS", rationale="Not annotated (offline) - cannot be assessed")
    if res["tier"] in ("Tier I", "Tier II", "Tier III") and isinstance(v.get("cosmic"), str):
        res["rationale"] += f"; reported in COSMIC ({v['cosmic']})"
    return res


def tier_cnv(gene: str, call: str, tumour: str) -> dict:
    if call != "AMPLIFICATION":
        return {"tier": "Tier III" if call != "neutral" else None, "therapies_here": "", "therapies_other": "", "level": None,
                "rationale": "Copy-number change of unknown significance" if call != "neutral" else ""}
    ev = evidence_for(gene, ["AMPLIFICATION"])
    here = [e for e in ev if _matches(e["tumour_types"], tumour)]
    other = [e for e in ev if not _matches(e["tumour_types"], tumour)]
    if here:
        return {"tier": "Tier I", "level": "A", "therapies_here": "; ".join(e["therapies"] for e in here),
                "therapies_other": "", "rationale": f"{gene} amplification is an approved biomarker in {tumour} (confirm by IHC/ISH)"}
    if other:
        return {"tier": "Tier II", "level": "C", "therapies_here": "",
                "therapies_other": "; ".join(f"{', '.join(sorted(e['tumour_types']))}: {e['therapies']}" for e in other),
                "rationale": f"{gene} amplification is actionable in other tumour types"}
    return {"tier": "Tier II", "level": "D", "therapies_here": "", "therapies_other": "", "rationale": f"{gene} amplification (oncogene)"}


def tier_all(df: pd.DataFrame, tumour: str, overrides: dict | None = None) -> pd.DataFrame:
    out = df.copy()
    res = [tier_variant(r.to_dict(), tumour) for _, r in out.iterrows()]
    for col in ("protein", "tier", "level", "therapies_here", "therapies_other", "rationale"):
        out[col] = [r[col] for r in res]
    out["tier_source"] = "auto"
    for vid, ov in (overrides or {}).items():
        m = out.variant_id == vid
        if m.any():
            out.loc[m, "tier"] = ov["tier"]
            out.loc[m, "rationale"] = ov.get("note") or out.loc[m, "rationale"]
            out.loc[m, "tier_source"] = "manual"
    order = {"Tier I": 0, "Tier II": 1, "Tier III": 2, "Tier IV": 3}
    out["_o"] = out.tier.map(order).fillna(9)
    return out.sort_values(["_o", "t_vaf"], ascending=[True, False]).drop(columns="_o").reset_index(drop=True)


TIER_COLORS = {"Tier I": "#C0392B", "Tier II": "#E67E22", "Tier III": "#F1C40F", "Tier IV": "#27AE60"}
