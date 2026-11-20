"""Simulate a matched tumour / normal panel sequencing run with real cancer hotspot mutations.

Hotspots are defined by their HGVS coding change (e.g. KRAS c.34G>T = p.G12C) and mapped to GRCh38
through the Ensembl REST API, then checked against the downloaded reference sequence, so the demo
mutations are at their true genomic positions. A truth set is kept for validation.
"""
from __future__ import annotations

import gzip
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import cx_reference as R

ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")

# gene, c.start, c.end, coding ref, coding alt, protein label  (alt "" = deletion)
HOTSPOTS = {
    "KRAS G12C": ("KRAS", 34, 34, "G", "T", "p.G12C"),
    "KRAS G12D": ("KRAS", 35, 35, "G", "A", "p.G12D"),
    "KRAS G13D": ("KRAS", 38, 38, "G", "A", "p.G13D"),
    "NRAS Q61R": ("NRAS", 182, 182, "A", "G", "p.Q61R"),
    "BRAF V600E": ("BRAF", 1799, 1799, "T", "A", "p.V600E"),
    "EGFR L858R": ("EGFR", 2573, 2573, "T", "G", "p.L858R"),
    "EGFR T790M": ("EGFR", 2369, 2369, "C", "T", "p.T790M"),
    "EGFR exon 19 deletion": ("EGFR", 2235, 2249, None, "", "p.E746_A750del"),
    "PIK3CA H1047R": ("PIK3CA", 3140, 3140, "A", "G", "p.H1047R"),
    "PIK3CA E545K": ("PIK3CA", 1633, 1633, "G", "A", "p.E545K"),
    "TP53 R273H": ("TP53", 818, 818, "G", "A", "p.R273H"),
    "TP53 R175H": ("TP53", 524, 524, "G", "A", "p.R175H"),
    "IDH1 R132H": ("IDH1", 395, 395, "G", "A", "p.R132H"),
    "KIT D816V": ("KIT", 2447, 2447, "A", "T", "p.D816V"),
}

SCENARIOS = {
    "Lung adenocarcinoma — KRAS G12C + TP53": {"tumor": "Non-small cell lung cancer", "drivers": ["KRAS G12C", "TP53 R273H"],
                                                "subclonal": ["PIK3CA E545K"], "amp": {}},
    "Lung adenocarcinoma — EGFR L858R + T790M (resistance)": {"tumor": "Non-small cell lung cancer",
                                                               "drivers": ["EGFR L858R", "TP53 R175H"], "subclonal": ["EGFR T790M"], "amp": {}},
    "Melanoma — BRAF V600E": {"tumor": "Melanoma", "drivers": ["BRAF V600E"], "subclonal": ["TP53 R273H"], "amp": {}},
    "Colorectal — KRAS G12D + PIK3CA + TP53": {"tumor": "Colorectal cancer", "drivers": ["KRAS G12D", "TP53 R175H"],
                                                "subclonal": ["PIK3CA E545K"], "amp": {}},
    "Breast — PIK3CA H1047R + ERBB2 amplification": {"tumor": "Breast cancer", "drivers": ["PIK3CA H1047R", "TP53 R273H"],
                                                      "subclonal": [], "amp": {"ERBB2": 8}},
}


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def hotspot_to_vcf(name: str, ref: R.Reference, fa: pysam.FastaFile, log=None):
    gene, cs, ce, cref, calt, prot = HOTSPOTS[name]
    g = next((x for x in ref.genes if x["symbol"] == gene), None)
    if g is None:
        return None
    m = R.map_cds(g["transcript"], cs, ce)
    if not m:
        if log is not None:
            log.append(f"{name}: could not map c.{cs} on {g['transcript']}")
        return None
    loc = ref.to_local(m["chrom"], m["start"])
    if not loc:
        return None
    ctg, lpos = loc
    if calt == "":  # deletion of c.cs..ce -> VCF with left anchor base
        anchor = fa.fetch(ctg, lpos - 2, lpos - 1).upper()
        deleted = fa.fetch(ctg, lpos - 1, lpos - 1 + (ce - cs + 1)).upper()
        return ctg, lpos - 1, anchor + deleted, anchor
    gref = cref if m["strand"] == 1 else cref.translate(COMP)
    galt = calt if m["strand"] == 1 else calt.translate(COMP)
    actual = fa.fetch(ctg, lpos - 1, lpos).upper()
    if actual != gref:
        if log is not None:
            log.append(f"{name}: reference base {actual} != expected {gref} - skipped")
        return None
    return ctg, lpos, gref, galt


def design(ref: R.Reference, scenario: str, purity=0.6, n_germline=6, n_passenger=2, seed=5, log=None) -> tuple[pd.DataFrame, dict]:
    """Return the variant truth table with expected VAF in tumour and normal, and copy-number changes."""
    rng = random.Random(seed)
    sc = SCENARIOS[scenario]
    fa = pysam.FastaFile(str(ref.fasta))
    rows = []
    online = "synthetic" not in ref.mode
    for kind, names, ccf in (("somatic driver", sc["drivers"], 1.0), ("somatic subclonal", sc["subclonal"], .4)):
        for n in names:
            v = hotspot_to_vcf(n, ref, fa, log) if online else None
            if v:
                rows.append((*v, kind, n, purity * ccf / 2, 0.0))
    tgt = ref.targets()
    n_som = n_passenger + (4 if not rows else 0)  # offline/no hotspots: add more random somatic mutations
    for k in range(n_som + n_germline):
        t = tgt.sample(1, random_state=seed + 31 * k).iloc[0]
        if t.end - t.start < 40:
            continue
        pos = rng.randint(t.start + 15, t.end - 15)
        base = fa.fetch(t.contig, pos - 1, pos).upper()
        alt = rng.choice([b for b in "ACGT" if b != base])
        if k < n_som:
            rows.append((t.contig, pos, base, alt, "somatic passenger", "random", purity / 2 * rng.choice([1, 1, .5]), 0.0))
        else:
            af = rng.choice([.5, .5, .5, 1.0])
            rows.append((t.contig, pos, base, alt, "germline", "random SNP", af, af))
    df = pd.DataFrame(rows, columns=["contig", "pos", "ref", "alt", "kind", "label", "tumor_af", "normal_af"])
    df = df.drop_duplicates(["contig", "pos"]).sort_values(["contig", "pos"]).reset_index(drop=True)
    if len(df):
        df["chrom"], df["gpos"] = zip(*[ref.to_genomic(c, p) for c, p in zip(df.contig, df.pos)])
    cn = {g: c for g, c in sc["amp"].items() if any(x["symbol"] == g for x in ref.genes)}
    return df, cn


def _apply(seq, seg_start, variants):
    for pos, vref, valt in sorted(variants, key=lambda v: -v[0]):
        i = pos - 1 - seg_start
        if 0 <= i and i + len(vref) <= len(seq):
            seq = seq[:i] + valt + seq[i + len(vref):]
    return seq


def _mutate(seq, qual, rng):
    perr = 10 ** (-qual / 10)
    flip = rng.random(len(seq)) < perr
    if not flip.any():
        return seq
    s = list(seq)
    for i in np.where(flip)[0]:
        s[i] = "ACGT"[(("ACGT".find(s[i]) if s[i] in "ACGT" else 0) + rng.integers(1, 4)) % 4]
    return "".join(s)


def simulate_sample(ref: R.Reference, variants: pd.DataFrame, af_col: str, out_dir: Path, sample: str, coverage=200,
                    gene_factor: dict | None = None, read_len=150, dup_rate=.08, seed=1, progress=None):
    rng = np.random.default_rng(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fa = pysam.FastaFile(str(ref.fasta))
    tgt = ref.targets()
    tgt["gene"] = tgt["name"].str.split("_").str[0]
    w = np.array((tgt.end - tgt.start), dtype=float)
    w *= np.array([(gene_factor or {}).get(g, 1.0) for g in tgt.gene])
    n_frag = int(w.sum() * coverage / (2 * read_len) * 1.35)
    tidx = rng.choice(len(tgt), n_frag, p=w / w.sum())
    vlist = {}
    for _, v in variants.iterrows():
        if v[af_col] > 0:
            vlist.setdefault(v.contig, []).append((int(v.pos), v.ref, v.alt, float(v[af_col])))
    ctg_len = {c: fa.get_reference_length(c) for c in fa.references}
    r1p, r2p = out_dir / f"{sample}_R1.fastq.gz", out_dir / f"{sample}_R2.fastq.gz"
    n = 0
    with gzip.open(r1p, "wt", compresslevel=3) as f1, gzip.open(r2p, "wt", compresslevel=3) as f2:
        for k, ti in enumerate(tidx):
            t = tgt.iloc[ti]
            flen = int(np.clip(rng.normal(250, 55), 90, 600))
            centre = int(rng.integers(t.start - 60, t.end + 60))
            s = max(0, centre - flen // 2)
            e = min(ctg_len[t.contig], s + flen + 20)
            seg = fa.fetch(t.contig, s, e).upper()
            vs = [(p, a, b) for p, a, b, af in vlist.get(t.contig, []) if s < p <= e and rng.random() < af]
            frag = _apply(seg, s, vs)[:flen]
            if len(frag) < 60:
                continue
            if rng.random() < .5:
                frag = revcomp(frag)
            q = np.clip(np.linspace(38, 31, read_len) + rng.normal(0, 3, (2, read_len)), 2, 41).astype(int)
            for _ in range(2 if rng.random() < dup_rate else 1):
                n += 1
                r1 = _mutate((frag + ADAPTER1 + "A" * read_len)[:read_len], q[0], rng)
                r2 = _mutate((revcomp(frag) + ADAPTER2 + "A" * read_len)[:read_len], q[1], rng)
                qs1, qs2 = "".join(chr(x + 33) for x in q[0]), "".join(chr(x + 33) for x in q[1])
                f1.write(f"@{sample}:{n} 1:N:0:ACGT\n{r1}\n+\n{qs1}\n")
                f2.write(f"@{sample}:{n} 2:N:0:ACGT\n{r2}\n+\n{qs2}\n")
            if progress and k % 3000 == 0:
                progress(k / n_frag, f"{sample}: simulating reads {k:,}/{n_frag:,}")
    return r1p, r2p


def tumor_gene_factor(cn: dict, purity: float) -> dict:
    """Read-depth multiplier for amplified genes in the tumour: (purity*CN + (1-purity)*2) / 2."""
    return {g: (purity * c + (1 - purity) * 2) / 2 for g, c in cn.items()}
