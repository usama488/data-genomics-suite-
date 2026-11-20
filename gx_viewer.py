"""IGV-style read pileup viewer drawn with Plotly."""
from __future__ import annotations

import plotly.graph_objects as go
import pysam
from plotly.subplots import make_subplots

BASE_COL = {"A": "#2ECC71", "C": "#3498DB", "G": "#F39C12", "T": "#E74C3C", "N": "#7F8C8D"}


def pileup_fig(bam: str, fasta: str, contig: str, pos: int, window: int = 45, max_rows: int = 60, title: str = ""):
    start, end = max(0, pos - window - 1), pos + window
    fa = pysam.FastaFile(fasta)
    refseq = fa.fetch(contig, start, end).upper()
    rows_end: list[int] = []
    reads = []
    with pysam.AlignmentFile(bam) as f:
        cov = [0] * (end - start)
        alt_cov = [0] * (end - start)
        for r in f.fetch(contig, start, end):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            rs, re_ = r.reference_start, r.reference_end
            # greedy row packing
            row = next((i for i, e in enumerate(rows_end) if e < rs - 1), None)
            if row is None:
                if len(rows_end) >= max_rows:
                    row = -1
                else:
                    rows_end.append(re_)
                    row = len(rows_end) - 1
            else:
                rows_end[row] = re_
            seq = r.query_sequence or ""
            mism, dels, ins = [], [], []
            for (qp, rp) in r.get_aligned_pairs():
                if rp is None or rp < start or rp >= end:
                    continue
                if qp is None:
                    dels.append(rp)
                    continue
                cov[rp - start] += 1
                b = seq[qp]
                if b != refseq[rp - start]:
                    mism.append((rp, b))
                    alt_cov[rp - start] += 1
            # insertions: anchor at the reference base before the inserted bases
            last_rp = None
            for qp, rp in r.get_aligned_pairs():
                if rp is not None:
                    last_rp = rp
                elif qp is not None and last_rp is not None and start <= last_rp < end:
                    ins.append(last_rp + .5)
            if row >= 0:
                reads.append({"row": row, "s": max(rs, start), "e": min(re_, end), "rev": r.is_reverse, "mism": mism,
                              "dels": dels, "ins": sorted(set(i for i in ins if i)), "dup": r.is_duplicate,
                              "name": r.query_name, "mapq": r.mapping_quality})
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[.18, .82], vertical_spacing=.02)
    xs = list(range(start + 1, end + 1))
    fig.add_bar(x=xs, y=cov, marker_color="#95A5A6", name="depth", row=1, col=1, hovertemplate="pos %{x}<br>depth %{y}<extra></extra>")
    fig.add_bar(x=xs, y=alt_cov, marker_color="#E74C3C", name="mismatch", row=1, col=1)
    shapes = []
    for rd in reads:
        y = -rd["row"]
        col = "#AFC6E9" if not rd["rev"] else "#F5B7B1"
        if rd["dup"]:
            col = "#DDDDDD"
        shapes.append(dict(type="rect", x0=rd["s"] + .5, x1=rd["e"] + .5, y0=y - .35, y1=y + .35, fillcolor=col,
                           line=dict(width=0), xref="x2", yref="y2", layer="below"))
        for d in rd["dels"]:
            shapes.append(dict(type="line", x0=d + .5, x1=d + 1.5, y0=y, y1=y, line=dict(color="black", width=2), xref="x2", yref="y2"))
    fig.update_layout(shapes=shapes)
    # mismatches as coloured squares with letters
    for b, c in BASE_COL.items():
        pts = [(rp + 1, -rd["row"]) for rd in reads for rp, bb in rd["mism"] if bb == b]
        if pts:
            fig.add_scatter(x=[p[0] for p in pts], y=[p[1] for p in pts], mode="markers+text", text=[b] * len(pts),
                            textfont=dict(size=8, color="white"), marker=dict(symbol="square", size=10, color=c),
                            name=f"mismatch {b}", row=2, col=1, hoverinfo="skip")
    ins_pts = [(i + 1, -rd["row"]) for rd in reads for i in rd["ins"]]
    if ins_pts:
        fig.add_scatter(x=[p[0] for p in ins_pts], y=[p[1] for p in ins_pts], mode="markers", name="insertion",
                        marker=dict(symbol="line-ns", size=12, line=dict(width=3, color="#8E44AD")), row=2, col=1)
    # reference sequence along the top of the read panel
    fig.add_scatter(x=xs, y=[1.2] * len(xs), mode="text", text=list(refseq), textfont=dict(size=9,
                    color=[BASE_COL.get(b, "#555") for b in refseq]), showlegend=False, row=2, col=1, hoverinfo="skip")
    fig.add_vline(x=pos, line=dict(color="black", dash="dot", width=1))
    n_rows = max((rd["row"] for rd in reads), default=0) + 1
    fig.update_yaxes(visible=False, row=2, col=1, range=[-n_rows, 2])
    fig.update_yaxes(title_text="depth", row=1, col=1)
    fig.update_layout(height=260 + 9 * n_rows, barmode="overlay", title=title or f"{contig}:{pos}", showlegend=True,
                      legend=dict(orientation="h", y=-0.05), margin=dict(t=50, l=40, r=10))
    fig.update_xaxes(range=[start + .5, end + .5])
    return fig, {"reads_shown": len(reads), "depth_at_pos": cov[pos - start - 1] if 0 <= pos - start - 1 < len(cov) else 0}
