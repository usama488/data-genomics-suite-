"""Simulate a realistic paired-end Illumina panel sequencing run (FASTQ.gz) with known spiked variants.

Gives the app a ready-made test sample and a truth set, so the pipeline's sensitivity and
precision can be measured at the end (like a laboratory validation run).
"""
from __future__ import annotations

import gzip
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import gx_reference as R

ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def _ensembl_to_vcf(v: dict, ref: R.Reference, fa: pysam.FastaFile):
    a = v["alleles"]
    if len(a) != 2:
        return None
    r, alt = a
    if v.get("strand", 1) == -1:
        r, alt = (revcomp(r) if r != "-" else r), (revcomp(alt) if alt != "-" else alt)
    loc = ref.to_local(str(v.get("chrom")), v["start"]) if "chrom" in v else None
    if loc is None:
        return None
    ctg, lpos = loc
    if r != "-" and alt != "-" and len(r) == len(alt) == 1:
        pos, vref, valt = lpos, r, alt
    elif alt == "-" and r != "-":  # deletion
        pos = lpos - 1
        anchor = fa.fetch(ctg, pos - 1, pos)
        vref, valt = anchor + r, anchor
    elif r == "-" and alt != "-":  # insertion between end and start
        pos = lpos - 1
        anchor = fa.fetch(ctg, pos - 1, pos)
        vref, valt = anchor, anchor + alt
    else:
        return None
    if fa.fetch(ctg, pos - 1, pos - 1 + len(vref)).upper() != vref:
        return None
    return ctg, pos, vref, valt


def pick_variants(ref: R.Reference, n_random: int = 6, seed: int = 7, use_clinvar: bool = True, log=None) -> pd.DataFrame:
    """Choose variants to spike: real ClinVar variants (online) + random novel ones."""
    rng = random.Random(seed)
    fa = pysam.FastaFile(str(ref.fasta))
    tgt = ref.targets()
    rows = []
    if use_clinvar and "synthetic" not in ref.mode:
        for g in ref.genes:
            # query only a few exons per gene to keep the request small
            exons = sorted(g["exons"], key=lambda e: e[1] - e[0], reverse=True)[:4]
            cands = []
            for s, e in exons:
                for v in R.clinvar_variants(g["chrom"], s, e):
                    v["chrom"] = g["chrom"]
                    cands.append(v)
            path = [v for v in cands if any("pathogenic" in c and "conflict" not in c for c in v["clin_sig"])]
            benign = [v for v in cands if any(c in ("benign", "likely benign") for c in v["clin_sig"])]
            rng.shuffle(path)
            rng.shuffle(benign)
            chosen = path[:2] + benign[:1]
            for v in chosen:
                conv = _ensembl_to_vcf(v, ref, fa)
                if conv:
                    rows.append((*conv, "ClinVar " + "/".join(v["clin_sig"]), v["id"]))
            if log is not None:
                log.append(f"{g['symbol']}: {len(cands)} ClinVar variants in largest exons, spiked {sum(1 for r in rows if r[0] == g['contig'])}")
    # random novel variants inside targets (SNVs + small indels)
    for k in range(n_random):
        t = tgt.sample(1, random_state=seed + k).iloc[0]
        if t.end - t.start < 40:
            continue
        pos = rng.randint(t.start + 15, t.end - 15)
        base = fa.fetch(t.contig, pos - 1, pos).upper()
        kind = ["snv", "snv", "snv", "del", "ins"][k % 5]
        if kind == "snv":
            rows.append((t.contig, pos, base, rng.choice([b for b in "ACGT" if b != base]), "novel (random)", "."))
        elif kind == "del":
            ln = rng.randint(1, 4)
            rows.append((t.contig, pos, fa.fetch(t.contig, pos - 1, pos + ln).upper(), base, "novel (random)", "."))
        else:
            rows.append((t.contig, pos, base, base + "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3))), "novel (random)", "."))
    df = pd.DataFrame(rows, columns=["contig", "pos", "ref", "alt", "source", "id"]).drop_duplicates(["contig", "pos"])
    # avoid overlapping variants (keep it simple for the simulator)
    df = df.sort_values(["contig", "pos"])
    keep, last = [], {}
    for i, r in df.iterrows():
        if r.contig in last and r.pos <= last[r.contig] + 5:
            continue
        keep.append(i)
        last[r.contig] = r.pos + len(r.ref)
    df = df.loc[keep].reset_index(drop=True)
    df["zygosity"] = ["het" if rng.random() < .75 else "hom" for _ in range(len(df))]
    df["chrom"], df["gpos"] = zip(*[ref.to_genomic(c, p) for c, p in zip(df.contig, df.pos)]) if len(df) else ([], [])
    return df


def _apply(seq: str, seg_start: int, variants) -> str:
    """Apply variants (1-based local pos) to a 0-based segment starting at seg_start."""
    for pos, vref, valt in sorted(variants, key=lambda v: -v[0]):
        i = pos - 1 - seg_start
        if 0 <= i and i + len(vref) <= len(seq):
            seq = seq[:i] + valt + seq[i + len(vref):]
    return seq


def _quals(n: int, rng: np.random.Generator, read_len: int) -> np.ndarray:
    base = np.linspace(38, 31, read_len)
    q = base + rng.normal(0, 3, (n, read_len))
    drop = rng.random(n) < .03  # a few poor reads
    q[drop] -= 15
    return np.clip(q, 2, 41).astype(int)


def _mutate(seq: str, qual: np.ndarray, rng: np.random.Generator) -> str:
    perr = 10 ** (-qual / 10)
    flip = rng.random(len(seq)) < perr
    if not flip.any():
        return seq
    s = list(seq)
    for i in np.where(flip)[0]:
        s[i] = "ACGT"[(("ACGT".find(s[i]) if s[i] in "ACGT" else 0) + rng.integers(1, 4)) % 4]
    return "".join(s)


def simulate(ref: R.Reference, variants: pd.DataFrame, out_dir: Path, sample="DEMO01", coverage=100, read_len=150,
             frag_mean=260, frag_sd=60, dup_rate=.06, seed=11, progress=None) -> tuple[Path, Path]:
    rng = np.random.default_rng(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fa = pysam.FastaFile(str(ref.fasta))
    tgt = ref.targets()
    lens = (tgt.end - tgt.start).to_numpy()
    n_frag = int(lens.sum() * coverage / (2 * read_len) * 1.35)  # extra for off-target / flanks
    hap_vars = {0: {}, 1: {}}
    for _, v in variants.iterrows():
        for h in ([0, 1] if v.zygosity == "hom" else [int(rng.integers(0, 2))]):
            hap_vars[h].setdefault(v.contig, []).append((int(v.pos), v.ref, v.alt))
    r1p, r2p = out_dir / f"{sample}_R1.fastq.gz", out_dir / f"{sample}_R2.fastq.gz"
    probs = lens / lens.sum()
    tidx = rng.choice(len(tgt), n_frag, p=probs)
    ctg_len = {c: fa.get_reference_length(c) for c in fa.references}
    with gzip.open(r1p, "wt", compresslevel=3) as f1, gzip.open(r2p, "wt", compresslevel=3) as f2:
        n = 0
        for k, ti in enumerate(tidx):
            t = tgt.iloc[ti]
            flen = int(np.clip(rng.normal(frag_mean, frag_sd), 90, 600))
            centre = int(rng.integers(t.start - 60, t.end + 60))
            s = max(0, centre - flen // 2)
            e = min(ctg_len[t.contig], s + flen + 10)
            hap = int(rng.integers(0, 2))
            segment = fa.fetch(t.contig, s, e).upper()
            vs = [v for v in hap_vars[hap].get(t.contig, []) if s < v[0] <= e]
            frag = _apply(segment, s, vs)[:flen]
            if len(frag) < 60:
                continue
            if rng.random() < .5:
                frag = revcomp(frag)
            copies = 2 if rng.random() < dup_rate else 1
            q = _quals(2, rng, read_len)
            for c in range(copies):
                n += 1
                r1 = (frag + ADAPTER1 + "A" * read_len)[:read_len]
                r2 = (revcomp(frag) + ADAPTER2 + "A" * read_len)[:read_len]
                r1, r2 = _mutate(r1, q[0], rng), _mutate(r2, q[1], rng)
                name = f"@{sample}:{n}"
                qs1 = "".join(chr(x + 33) for x in q[0])
                qs2 = "".join(chr(x + 33) for x in q[1])
                f1.write(f"{name} 1:N:0:ACGTACGT\n{r1}\n+\n{qs1}\n")
                f2.write(f"{name} 2:N:0:ACGTACGT\n{r2}\n+\n{qs2}\n")
            if progress and k % 2000 == 0:
                progress(k / n_frag, f"Simulating reads {k:,}/{n_frag:,}")
    if progress:
        progress(1.0, f"Wrote {n:,} read pairs")
    return r1p, r2p
