"""FastQC-style quality control of FASTQ files (computed in Python on a read sample)."""
from __future__ import annotations

import gzip
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

ADAPTERS = {"Illumina universal": "AGATCGGAAGAGC", "Nextera": "CTGTCTCTTATACACATCT", "PolyA": "AAAAAAAAAAAAAAA",
            "PolyG (NovaSeq)": "GGGGGGGGGGGGGGG"}


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fastq(path, max_reads=200_000):
    seqs, quals = [], []
    with _open(path) as fh:
        while len(seqs) < max_reads:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            q = fh.readline().strip()
            if not h.startswith("@") or len(s) != len(q):
                raise ValueError(f"{Path(path).name}: not a valid FASTQ (record {len(seqs) + 1})")
            seqs.append(s)
            quals.append(q)
    return seqs, quals


def count_reads(path) -> int:
    n = 0
    with _open(path) as fh:
        for n, _ in enumerate(fh, 1):
            pass
    return n // 4


def analyse(path, max_reads=200_000) -> dict:
    seqs, quals = read_fastq(path, max_reads)
    if not seqs:
        raise ValueError("Empty FASTQ")
    n = len(seqs)
    lens = np.array([len(s) for s in seqs])
    L = int(lens.max())
    Q = np.full((n, L), np.nan)
    for i, q in enumerate(quals):
        Q[i, :len(q)] = np.frombuffer(q.encode(), dtype=np.uint8) - 33
    per_base = pd.DataFrame({
        "position": np.arange(1, L + 1), "mean": np.nanmean(Q, 0), "median": np.nanmedian(Q, 0),
        "q10": np.nanpercentile(Q, 10, 0), "q25": np.nanpercentile(Q, 25, 0), "q75": np.nanpercentile(Q, 75, 0),
        "q90": np.nanpercentile(Q, 90, 0)})
    mean_q_read = np.nanmean(Q, 1)
    # base composition
    comp = np.zeros((L, 5))
    for s in seqs:
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        for j, b in enumerate(b"ACGTN"):
            comp[:len(arr), j] += arr == b
    comp = comp / comp.sum(1, keepdims=True) * 100
    comp_df = pd.DataFrame(comp, columns=list("ACGTN"))
    comp_df["position"] = np.arange(1, L + 1)
    gc = np.array([(s.count("G") + s.count("C")) / max(len(s), 1) * 100 for s in seqs])
    # duplication (exact sequence) on first 50k reads
    cnt = Counter(seqs[:50_000])
    dup_pct = (1 - len(cnt) / min(n, 50_000)) * 100
    over = [(s, c) for s, c in cnt.most_common(10) if c / min(n, 50_000) > .001]
    # adapter content: % reads containing adapter by position
    step = max(1, L // 50)
    ad = {}
    for name, a in ADAPTERS.items():
        first = np.array([s.find(a) for s in seqs[:50_000]])
        ad[name] = [float(((first >= 0) & (first <= p)).mean() * 100) for p in range(0, L, step)]
    ad_df = pd.DataFrame(ad)
    ad_df["position"] = list(range(1, L + 1, step))
    q30 = float((Q >= 30).sum() / np.isfinite(Q).sum() * 100)
    n_pct = float(comp_df["N"].mean())
    m = {"reads_sampled": n, "length_min": int(lens.min()), "length_max": L, "length_mean": float(lens.mean()),
         "gc_mean": float(gc.mean()), "q30_pct": q30, "mean_quality": float(np.nanmean(Q)), "dup_pct": dup_pct,
         "n_pct": n_pct, "adapter_max_pct": float(ad_df.drop(columns="position").max().max())}
    return {"metrics": m, "per_base": per_base, "read_q": mean_q_read, "comp": comp_df, "gc": gc, "lens": lens,
            "adapters": ad_df, "overrep": over, "status": status(m, per_base)}


def status(m: dict, per_base: pd.DataFrame) -> dict:
    """PASS / WARN / FAIL per module, using FastQC-like thresholds."""
    s = {}
    lowq = per_base["q10"].min()
    s["Per-base quality"] = "FAIL" if per_base["median"].min() < 20 or lowq < 5 else "WARN" if lowq < 10 or per_base["median"].min() < 25 else "PASS"
    s["Q30 bases"] = "PASS" if m["q30_pct"] >= 80 else "WARN" if m["q30_pct"] >= 70 else "FAIL"
    s["GC content"] = "PASS" if 35 <= m["gc_mean"] <= 65 else "WARN"
    s["N content"] = "PASS" if m["n_pct"] < 1 else "WARN" if m["n_pct"] < 5 else "FAIL"
    s["Duplication"] = "PASS" if m["dup_pct"] < 20 else "WARN" if m["dup_pct"] < 50 else "FAIL"
    s["Adapter content"] = "PASS" if m["adapter_max_pct"] < 5 else "WARN" if m["adapter_max_pct"] < 10 else "FAIL"
    return s


def interpret(res: dict) -> list[dict]:
    m, s = res["metrics"], res["status"]
    out = [{"level": "info", "text": f"{m['reads_sampled']:,} reads sampled; length {m['length_min']}–{m['length_max']} bp; "
                                     f"mean quality Q{m['mean_quality']:.1f}; {m['q30_pct']:.1f}% of bases ≥ Q30."}]
    if s["Per-base quality"] != "PASS":
        out.append({"level": "warn", "text": "Quality drops towards the 3' end - quality trimming is recommended."})
    if s["Adapter content"] != "PASS":
        out.append({"level": "warn", "text": f"Adapter found in up to {m['adapter_max_pct']:.1f}% of reads (short inserts) - adapter trimming needed."})
    if s["Duplication"] != "PASS":
        out.append({"level": "warn", "text": f"{m['dup_pct']:.1f}% duplicate sequences - expected for targeted panels; duplicates are marked after alignment."})
    if s["GC content"] != "PASS":
        out.append({"level": "warn", "text": f"Mean GC {m['gc_mean']:.1f}% is unusual for human DNA (~40-50%) - possible contamination or capture bias."})
    if all(v == "PASS" for v in s.values()):
        out.append({"level": "good", "text": "All QC modules pass."})
    return out


def figs(res: dict, title: str) -> dict:
    pb = res["per_base"]
    f = {}
    q = go.Figure()
    q.add_hrect(y0=28, y1=42, fillcolor="green", opacity=.08, line_width=0)
    q.add_hrect(y0=20, y1=28, fillcolor="orange", opacity=.08, line_width=0)
    q.add_hrect(y0=0, y1=20, fillcolor="red", opacity=.08, line_width=0)
    q.add_scatter(x=pb.position, y=pb.q90, line=dict(width=0), showlegend=False, hoverinfo="skip")
    q.add_scatter(x=pb.position, y=pb.q10, fill="tonexty", line=dict(width=0), name="10–90%", fillcolor="rgba(46,134,171,.2)")
    q.add_scatter(x=pb.position, y=pb.q75, line=dict(width=0), showlegend=False, hoverinfo="skip")
    q.add_scatter(x=pb.position, y=pb.q25, fill="tonexty", line=dict(width=0), name="25–75%", fillcolor="rgba(46,134,171,.45)")
    q.add_scatter(x=pb.position, y=pb["median"], name="median", line_color="#E4572E")
    q.add_scatter(x=pb.position, y=pb["mean"], name="mean", line=dict(color="#2E86AB", dash="dot"))
    q.update_layout(title=f"Per-base sequence quality — {title}", xaxis_title="position in read (bp)", yaxis_title="Phred quality", yaxis_range=[0, 42])
    f["quality"] = q
    f["read_q"] = px.histogram(x=res["read_q"], nbins=40, title="Per-sequence mean quality", labels={"x": "mean Phred quality"})
    comp = res["comp"].melt("position", var_name="base", value_name="percent")
    f["composition"] = px.line(comp, x="position", y="percent", color="base", title="Per-base sequence content")
    f["gc"] = px.histogram(x=res["gc"], nbins=50, title="Per-sequence GC content", labels={"x": "GC %"})
    f["length"] = px.histogram(x=res["lens"], nbins=30, title="Read length distribution", labels={"x": "length (bp)"})
    ad = res["adapters"].melt("position", var_name="adapter", value_name="% reads")
    f["adapters"] = px.line(ad, x="position", y="% reads", color="adapter", title="Adapter content (cumulative)")
    return f
