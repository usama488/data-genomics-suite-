"""Variant annotation with the Ensembl VEP REST API (GRCh38).

Returns gene, transcript (MANE / canonical), consequence, HGVS c./p., exon, gnomAD frequency,
ClinVar significance, SIFT, PolyPhen and - when the server provides them - REVEL, AlphaMissense,
CADD and SpliceAI. Results are cached on disk so re-runs are instant and work offline.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd
import requests

import gx_reference as R

URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
CACHE = R.APP_DIR / "references" / "vep_cache.json"
BASE_PARAMS = {"hgvs": 1, "canonical": 1, "mane": 1, "numbers": 1, "af_gnomade": 1, "af_gnomadg": 1, "variant_class": 1}
MINIMAL_PARAMS = {"hgvs": 1, "canonical": 1}
PLUGINS = {"REVEL": 1, "AlphaMissense": 1, "SpliceAI": 2, "CADD": 1}
SEVERITY = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
            "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
            "missense_variant", "protein_altering_variant", "splice_donor_5th_base_variant", "splice_region_variant",
            "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant",
            "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
            "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant",
            "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
            "downstream_gene_variant", "intergenic_variant"]

_session = requests.Session()


def _load_cache() -> dict:
    try:
        return json.loads(CACHE.read_text())
    except Exception:
        return {}


def _save_cache(c: dict):
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(c))


def _vep_input(chrom, pos, ref, alt) -> str:
    return f"{chrom} {pos} . {ref} {alt} . . ."


def _post(variants: list[str], params: dict, tries=4) -> list[dict]:
    body = {"variants": variants, **params}
    for i in range(tries):
        r = _session.post(URL, json=body, headers={"Content-Type": "application/json", "Accept": "application/json"}, timeout=180)
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 2)) + 1)
            continue
        if r.ok:
            return r.json()
        if r.status_code == 400:
            raise ValueError(r.text[:300])
        time.sleep(2 + 2 * i)
    r.raise_for_status()
    return []


def query_vep(rows: list[tuple], log=None) -> dict:
    """rows: (variant_id, chrom, gpos, ref, alt). Returns {variant_id: raw VEP json}."""
    cache = _load_cache()
    todo = [r for r in rows if r[0] not in cache]
    if log is not None:
        log.append(f"VEP: {len(rows) - len(todo)} cached, {len(todo)} to query")
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        inputs = [_vep_input(c, p, rf, a) for _, c, p, rf, a in chunk]
        res = None
        for params in ({**BASE_PARAMS, **PLUGINS}, BASE_PARAMS, MINIMAL_PARAMS):
            try:  # a plugin/option not offered by the server gives HTTP 400 -> retry with fewer options
                res = _post(inputs, params)
                break
            except ValueError as e:
                if log is not None:
                    log.append(f"VEP rejected options ({str(e)[:120]}); retrying with fewer")
        if res is None:
            raise RuntimeError("VEP rejected every request")
        by_input = {x.get("input"): x for x in res}
        for (vid, *_), inp in zip(chunk, inputs):
            cache[vid] = by_input.get(inp, {"input": inp, "not_annotated": True})
    if todo:
        _save_cache(cache)
    return {r[0]: cache.get(r[0], {}) for r in rows}


def _pick_transcript(tcs: list[dict]) -> dict:
    if not tcs:
        return {}
    def rank(t):
        sev = min((SEVERITY.index(c) for c in t.get("consequence_terms", []) if c in SEVERITY), default=99)
        return (0 if t.get("mane_select") else 1, 0 if t.get("canonical") else 1, 0 if t.get("biotype") == "protein_coding" else 1, sev)
    return sorted(tcs, key=rank)[0]


def _find(d: dict, *keys):
    for k, v in d.items():
        lk = k.lower()
        if any(key in lk for key in keys):
            return v
    return None


def parse(vid: str, alt: str, js: dict) -> dict:
    out = {"variant_id": vid, "gene": None, "transcript": None, "mane": None, "consequence": None, "impact": None,
           "hgvsc": None, "hgvsp": None, "exon": None, "intron": None, "gnomad_af": None, "clinvar": None, "rsid": None,
           "sift": None, "polyphen": None, "revel": None, "alphamissense": None, "cadd": None, "spliceai": None,
           "all_consequences": None, "annotated": bool(js) and not js.get("not_annotated")}
    if not out["annotated"]:
        return out
    t = _pick_transcript(js.get("transcript_consequences", []))
    terms = t.get("consequence_terms", [js.get("most_severe_consequence")]) if t else [js.get("most_severe_consequence")]
    out.update({"gene": t.get("gene_symbol"), "transcript": t.get("transcript_id"), "mane": t.get("mane_select"),
                "consequence": ",".join(x for x in terms if x), "impact": t.get("impact"),
                "hgvsc": (t.get("hgvsc") or "").split(":")[-1] or None, "hgvsp": (t.get("hgvsp") or "").split(":")[-1] or None,
                "exon": t.get("exon"), "intron": t.get("intron"),
                "sift": f"{t['sift_prediction']} ({t.get('sift_score')})" if t.get("sift_prediction") else None,
                "polyphen": f"{t['polyphen_prediction']} ({t.get('polyphen_score')})" if t.get("polyphen_prediction") else None})
    rv = _find(t, "revel")
    try:
        out["revel"] = float(rv) if rv is not None else None
    except (TypeError, ValueError):
        pass
    am = _find(t, "am_pathogenicity", "alphamissense")
    if isinstance(am, dict):
        am = am.get("am_pathogenicity")
    try:
        out["alphamissense"] = float(am) if am is not None else None
    except (TypeError, ValueError):
        pass
    cadd = _find(t, "cadd_phred")
    out["cadd"] = float(cadd) if isinstance(cadd, (int, float)) else None
    sp = _find(t, "spliceai")
    if isinstance(sp, dict):
        ds = [v for k, v in sp.items() if k.upper().startswith("DS_") and isinstance(v, (int, float))]
        out["spliceai"] = max(ds) if ds else None
    # colocated variants: rsID, ClinVar, gnomAD
    afs, clin, ids = [], [], []
    for cv in js.get("colocated_variants", []):
        if cv.get("id"):
            ids.append(cv["id"])
        csa = cv.get("clin_sig_allele")
        if csa:
            for part in csa.split(";"):
                if ":" in part:
                    a, s = part.split(":", 1)
                    if a == alt or alt.endswith(a) or a == "-" and len(alt) == 1:
                        clin.append(s)
        elif cv.get("clin_sig"):
            clin += cv["clin_sig"]
        fr = (cv.get("frequencies") or {})
        for a, f in fr.items():
            if a == alt or alt.endswith(a) or len(fr) == 1:
                afs += [v for k, v in f.items() if k.startswith("gnomad") and isinstance(v, (int, float)) and k in ("gnomade", "gnomadg", "gnomad")]
    rs = [i for i in ids if i.startswith("rs")]
    out["rsid"] = rs[0] if rs else (ids[0] if ids else None)
    out["clinvar"] = ", ".join(sorted(set(clin))) if clin else None
    out["gnomad_af"] = max(afs) if afs else (0.0 if js.get("colocated_variants") is not None else None)
    out["all_consequences"] = len(js.get("transcript_consequences", []))
    return out


def annotate(df: pd.DataFrame, offline=False, log=None, ref: R.Reference | None = None) -> pd.DataFrame:
    if df.empty:
        return df
    rows = list(zip(df.variant_id, df.chrom, df.gpos, df.ref, df.alt))
    raw = {}
    if not offline:
        try:
            raw = query_vep(rows, log)
        except Exception as e:
            if log is not None:
                log.append(f"VEP unavailable ({e}); falling back to offline gene annotation")
            raw = {}
    ann = pd.DataFrame([parse(vid, alt, raw.get(vid, {})) for vid, alt in zip(df.variant_id, df.alt)])
    out = df.merge(ann, on="variant_id", how="left")
    if ref is not None:  # offline gene name from the panel BED
        g = pd.read_csv(ref.genes_bed, sep="\t", header=None, names=["contig", "s", "e", "gene"])
        def gene_of(r):
            hit = g[(g.contig == r.contig) & (g.s < r.pos) & (g.e >= r.pos)]
            return hit.gene.iloc[0] if len(hit) else None
        out["gene"] = [r.gene if isinstance(r.gene, str) else gene_of(r) for r in out.itertuples()]
    return out
