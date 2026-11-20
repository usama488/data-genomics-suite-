"""Somatic variant detection (tumour vs matched normal, or tumour-only), copy number, purity and TMB.

Calling approach (similar in spirit to VarScan2 somatic): bcftools mpileup is run jointly on the tumour
and normal BAMs to count reference / alternate reads per strand at every target position. Each candidate
allele is then classified with explicit, auditable rules plus a one-sided Fisher exact test.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact

import cx_pipeline as P
import cx_reference as R
import cx_tools as tools

DEFAULTS = {"min_tumor_dp": 20, "min_tumor_alt": 4, "min_tumor_vaf": 0.05, "min_normal_dp": 8,
            "max_normal_vaf": 0.02, "germline_vaf": 0.25, "fisher_p": 0.01, "min_bq": 20, "min_mq": 20}
NONSYN = {"missense_variant", "stop_gained", "frameshift_variant", "inframe_insertion", "inframe_deletion", "stop_lost",
          "start_lost", "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant"}


def joint_pileup(run: P.Run, ref: R.Reference, tumor_bam, normal_bam=None, s: dict = DEFAULTS, pad=20):
    padded = run.p("targets_padded.bed")
    t = ref.targets()
    t["start"] = (t.start - pad).clip(lower=0)
    t["end"] = t.end + pad
    t.sort_values(["contig", "start"]).to_csv(padded, sep="\t", header=False, index=False)
    out = run.p("joint_pileup.bcf")
    bams = [tumor_bam] + ([normal_bam] if normal_bam else [])
    with run.step("Joint pileup (bcftools mpileup, tumour + normal)"):
        tools.bcftools("mpileup", "-f", ref.fasta, "-a", "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP", "-q", s["min_mq"],
                       "-Q", s["min_bq"], "-d", 20000, "-m", 3, "-F", 0.01, "-R", padded, "-Ob", "-o", out, *bams, log=run.log)
    return out


def _counts(rec, i, si):
    smp = rec.samples[si]
    ad, adf, adr = smp.get("AD"), smp.get("ADF"), smp.get("ADR")
    if not ad or ad[0] is None:
        return 0, 0, 0, 0
    dp = sum(a for a in ad if a)
    alt = ad[i] or 0
    f = (adf[i] or 0) if adf else 0
    r = (adr[i] or 0) if adr else 0
    return dp, alt, f, r


def classify_sites(run: P.Run, ref: R.Reference, pile, paired: bool, s: dict = DEFAULTS) -> pd.DataFrame:
    rows = []
    with run.step("Somatic classification (rules + Fisher exact test)"):
        with pysam.VariantFile(str(pile)) as v:
            for rec in v:
                alts = rec.alts or ()
                for i, alt in enumerate(alts, start=1):
                    if alt == "<*>":
                        continue
                    tdp, talt, tf, tr = _counts(rec, i, 0)
                    if talt < 2 or tdp == 0:
                        continue
                    tvaf = talt / tdp
                    row = {"contig": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": alt, "t_dp": tdp, "t_alt": talt,
                           "t_vaf": round(tvaf, 4), "t_fwd": tf, "t_rev": tr}
                    if paired:
                        ndp, nalt, _, _ = _counts(rec, i, 1)
                        nvaf = nalt / ndp if ndp else 0.0
                        p = fisher_exact([[talt, tdp - talt], [nalt, ndp - nalt]], alternative="greater")[1] if ndp else 1.0
                        row.update({"n_dp": ndp, "n_alt": nalt, "n_vaf": round(nvaf, 4), "fisher_p": p})
                    rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    status, filt = [], []
    for r in df.itertuples():
        fl = []
        if r.t_dp < s["min_tumor_dp"]:
            fl.append("LowTumorDP")
        if r.t_alt < s["min_tumor_alt"]:
            fl.append("LowAltReads")
        if r.t_vaf < s["min_tumor_vaf"]:
            fl.append("LowVAF")
        if r.t_alt >= 6 and (r.t_fwd == 0 or r.t_rev == 0):
            fl.append("StrandBias")
        if paired:
            if r.n_dp < s["min_normal_dp"]:
                fl.append("LowNormalDP")
            if r.n_vaf >= s["germline_vaf"]:
                st = "germline"
            elif r.n_vaf <= s["max_normal_vaf"] or r.n_alt <= 1:
                st = "somatic" if r.fisher_p < s["fisher_p"] else "uncertain"
                if r.fisher_p >= s["fisher_p"]:
                    fl.append("FisherNS")
            else:
                st = "ambiguous (alt in normal)"
                fl.append("AltInNormal")
        else:  # no normal: VAF near 50% / 100% suggests an inherited variant
            st = "possible germline (tumour-only)" if (.4 <= r.t_vaf <= .6 or r.t_vaf >= .9) else "somatic (tumour-only)"
        if st == "germline":
            fl = [f for f in fl if f in ("LowTumorDP", "LowNormalDP", "StrandBias")]
        status.append(st)
        filt.append(";".join(fl) or "PASS")
    df["status"], df["filter"] = status, filt
    # keep candidates that are at least plausible (drop pure noise)
    keep = (df["filter"] == "PASS") | (df.t_vaf >= 0.03) | (df.status == "germline")
    return df[keep].reset_index(drop=True)


def normalise(run: P.Run, ref: R.Reference, df: pd.DataFrame) -> pd.DataFrame:
    """Write candidates to VCF, left-align/split with bcftools norm, read back."""
    if df.empty:
        return df
    raw, norm = run.p("somatic_candidates.vcf"), run.p(f"{run.sample}.somatic.vcf.gz")
    fai = pd.read_csv(str(ref.fasta) + ".fai", sep="\t", header=None)
    info_cols = [c for c in df.columns if c not in ("contig", "pos", "ref", "alt", "filter")]
    with open(raw, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, l in fai[[0, 1]].values:
            fh.write(f"##contig=<ID={c},length={l}>\n")
        for f in sorted({x for fs in df["filter"] for x in fs.split(";")} - {"PASS"}):
            fh.write(f'##FILTER=<ID={f},Description="{f}">\n')
        for c in info_cols:
            typ = "String" if c == "status" else ("Float" if df[c].dtype.kind == "f" else "Integer")
            fh.write(f'##INFO=<ID={c.upper()},Number=1,Type={typ},Description="{c}">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in df.sort_values(["contig", "pos"]).itertuples(index=False):
            d = r._asdict()
            info = ";".join(f"{c.upper()}={str(d[c]).replace(' ', '_').replace(';', ',')}" for c in info_cols)
            fh.write(f"{d['contig']}\t{d['pos']}\t.\t{d['ref']}\t{d['alt']}\t.\t{d['filter']}\t{info}\n")
    tools.bcftools("norm", "-f", ref.fasta, "-m", "-any", "-Oz", "-o", norm, raw, log=run.log)
    tools.bcftools("index", "-f", norm, log=run.log)
    out = []
    with pysam.VariantFile(str(norm)) as v:
        for rec in v:
            d = {"contig": rec.chrom, "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0], "filter": ";".join(rec.filter.keys()) or "PASS"}
            for c in info_cols:
                val = rec.info.get(c.upper())
                d[c] = val.replace("_", " ") if isinstance(val, str) else val
            out.append(d)
    res = pd.DataFrame(out).drop_duplicates(["contig", "pos", "ref", "alt"])
    res["chrom"], res["gpos"] = zip(*[ref.to_genomic(c, p) for c, p in zip(res.contig, res.pos)])
    res["type"] = ["SNV" if len(a) == len(b) == 1 else "deletion" if len(a) > len(b) else "insertion" if len(b) > len(a) else "MNV"
                   for a, b in zip(res.ref, res.alt)]
    res["variant_id"] = [f"{c}-{p}-{a}-{b}" for c, p, a, b in zip(res.chrom, res.gpos, res.ref, res.alt)]
    res.attrs["vcf"] = norm
    return res


# ------------------------------------------------------------------ copy number
def copy_number(cov_t: pd.DataFrame, cov_n: pd.DataFrame | None, purity: float | None) -> pd.DataFrame:
    """Gene-level log2 ratio tumour/normal (median-normalised) and estimated absolute copy number."""
    t = cov_t.set_index("target")["mean_depth"]
    if cov_n is not None:
        n = cov_n.set_index("target")["mean_depth"].reindex(t.index)
        ratio = (t / t.median()) / (n / n.median()).replace(0, np.nan)
    else:
        ratio = t / t.median()
    ex = pd.DataFrame({"target": t.index, "ratio": ratio.to_numpy()})
    ex["gene"] = ex.target.str.split("_").str[0]
    g = ex.groupby("gene").agg(exons=("ratio", "size"), ratio=("ratio", "median"), sd=("ratio", "std")).reset_index()
    g["log2_ratio"] = np.log2(g.ratio.clip(lower=1e-3))
    p = purity if purity and purity > .1 else 1.0
    g["copy_number"] = (2 * (g.ratio - (1 - p)) / p).clip(lower=0).round(1)
    def call(r):
        if r.log2_ratio >= 1 and r.copy_number >= 5:
            return "AMPLIFICATION"
        if r.log2_ratio >= .4:
            return "gain"
        if r.log2_ratio <= -1:
            return "DEEP DELETION" if r.copy_number < .5 else "loss"
        if r.log2_ratio <= -.4:
            return "loss"
        return "neutral"
    g["call"] = g.apply(call, axis=1)
    return g, ex


def estimate_purity(som: pd.DataFrame, cnv: pd.DataFrame | None = None) -> float | None:
    """Purity from clonal heterozygous somatic SNVs: purity ≈ 2 × VAF of the top VAF cluster."""
    if som.empty:
        return None
    s = som[(som["filter"] == "PASS") & som.status.str.startswith("somatic") & (som.type == "SNV")]
    if cnv is not None and "gene" in s:
        amp = set(cnv[cnv.call != "neutral"].gene)
        s = s[~s.gene.isin(amp)]
    v = s.t_vaf.dropna().to_numpy()
    if len(v) < 1:
        return None
    top = np.sort(v)[len(v) // 2:] if len(v) >= 3 else v
    return float(min(1.0, 2 * np.median(top)))


def tmb(som: pd.DataFrame, target_bp: int) -> dict:
    s = som[(som["filter"] == "PASS") & som.status.str.startswith("somatic")]
    if "consequence" in s and s.consequence.notna().any():
        cons = s.consequence.fillna("").str.split(",")
        s = s[[bool(set(c) & NONSYN) for c in cons]]
        basis = "non-synonymous coding"
    else:
        basis = "all PASS somatic (not annotated)"
    mb = target_bp / 1e6
    val = len(s) / mb if mb else float("nan")
    return {"mutations": len(s), "megabases": mb, "tmb": val, "basis": basis,
            "class": "TMB-High (≥10 mut/Mb)" if val >= 10 else "TMB-Low",
            "reliable": mb >= 1.0}


def compare_truth(som: pd.DataFrame, truth: pd.DataFrame) -> dict:
    tru_s = truth[truth.kind.str.startswith("somatic")]
    called = som[(som["filter"] == "PASS") & som.status.str.startswith("somatic")]
    tkey = {(r.chrom, r.gpos, r.ref, r.alt) for r in tru_s.itertuples()}
    # indels in the truth are pre-normalisation: match on position window for indels
    ckey = {(r.chrom, r.gpos, r.ref, r.alt) for r in called.itertuples()}
    tp = {k for k in tkey if k in ckey or any(c[0] == k[0] and abs(c[1] - k[1]) <= 10 and len(c[2]) - len(c[3]) == len(k[2]) - len(k[3])
                                              and len(k[2]) != len(k[3]) for c in ckey)}
    fn = tkey - tp
    matched_c = {c for c in ckey if c in tkey or any(k[0] == c[0] and abs(c[1] - k[1]) <= 10 and len(c[2]) - len(c[3]) == len(k[2]) - len(k[3])
                                                    and len(k[2]) != len(k[3]) for k in tkey)}
    fp = ckey - matched_c
    germ = truth[truth.kind == "germline"]
    gkey = {(r.chrom, r.gpos) for r in germ.itertuples()}
    germ_as_som = {c for c in ckey if (c[0], c[1]) in gkey}
    return {"somatic_truth": len(tkey), "called_somatic": len(ckey), "TP": len(tp), "FN": len(fn), "FP": len(fp),
            "sensitivity": len(tp) / max(len(tkey), 1) * 100, "precision": len(tp) / max(len(ckey), 1) * 100,
            "germline_miscalled_somatic": len(germ_as_som), "missed": sorted(fn), "false_pos": sorted(fp)}
