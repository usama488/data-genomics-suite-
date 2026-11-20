"""Pipeline steps: trimming -> alignment -> duplicate marking -> BAM QC/coverage -> variant calling -> filtering."""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import gx_reference as R
import gx_simulate as simulate
import gx_tools as tools

THREADS = max(1, min(4, os.cpu_count() or 1))


class Run:
    """A sample's working directory + a log of every command executed."""

    def __init__(self, root: Path, sample: str):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.sample = sample
        self.log: list[str] = []
        self.timings: dict[str, float] = {}

    def p(self, name: str) -> Path:
        return self.root / name

    def step(self, name):
        run = self

        class _T:
            def __enter__(self):
                self.t = time.time()
                run.log.append(f"### {name}")

            def __exit__(self, *a):
                run.timings[name] = round(time.time() - self.t, 2)
        return _T()


# ------------------------------------------------------------------ 1. trimming
def trim(run: Run, r1: Path, r2: Path | None, adapter1=simulate.ADAPTER1, adapter2=simulate.ADAPTER2,
         qual=20, min_len=36, nextseq=False) -> dict:
    with run.step("Trimming (cutadapt)"):
        o1 = run.p(f"{run.sample}_R1.trimmed.fastq.gz")
        rep = run.p("cutadapt.json")
        args = ["-j", THREADS, "-a", adapter1, ("--nextseq-trim" if nextseq else "-q"), (str(qual) if nextseq else f"{qual},{qual}"),
                "-m", min_len, "--json", rep, "-o", o1]
        o2 = None
        if r2:
            o2 = run.p(f"{run.sample}_R2.trimmed.fastq.gz")
            args += ["-A", adapter2, "-p", o2, r1, r2]
        else:
            args += [r1]
        if nextseq:
            args[args.index("--nextseq-trim")] = f"--nextseq-trim={qual}"
            args.remove(str(qual))
        tools.cutadapt(args, run.log)
        js = json.loads(rep.read_text())
    rc, bc = js["read_counts"], js["basepair_counts"]
    pairs_in = rc["input"]
    out = {"r1": o1, "r2": o2, "reads_in": pairs_in, "reads_out": rc["output"],
           "too_short": (rc.get("filtered") or {}).get("too_short") or 0,
           "with_adapter_r1": (rc.get("read1_with_adapter") or 0), "with_adapter_r2": (rc.get("read2_with_adapter") or 0),
           "bp_in": bc["input"], "bp_out": bc["output"], "quality_trimmed_bp": (bc.get("quality_trimmed") or 0)}
    out["kept_pct"] = out["reads_out"] / max(pairs_in, 1) * 100
    return out


# ------------------------------------------------------------------ 2. alignment
def align(run: Run, ref: R.Reference, r1: Path, r2: Path | None, platform="ILLUMINA", progress=None) -> dict:
    rg = f"@RG\\tID:{run.sample}\\tSM:{run.sample}\\tPL:{platform}\\tLB:{run.sample}_lib"
    sam = run.p("aligned.sam")
    with run.step("Alignment (BWA-MEM)"):
        if tools.have("bwa"):
            cmd = ["bwa", "mem", "-t", str(THREADS), "-R", rg, str(ref.fasta), str(r1)] + ([str(r2)] if r2 else [])
            run.log.append("$ " + " ".join(cmd) + f" > {sam.name}")
            with open(sam, "w") as fh:
                p = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
            if p.returncode != 0:
                raise tools.ToolError(p.stderr[-2000:])
            paired = r2 is not None
        else:
            paired = False
            _align_pybwa(run, ref, r1, r2, sam, progress)
    with run.step("Sort + mark duplicates (samtools)"):
        fx, srt, md = run.p("fixmate.bam"), run.p("sorted.bam"), run.p(f"{run.sample}.bam")
        tools.samtools("fixmate", "-m", "-O", "bam", "-@", THREADS, sam, fx, log=run.log)
        tools.samtools("sort", "-@", THREADS, "-o", srt, fx, log=run.log)
        tools.samtools("markdup", "-s", "-f", run.p("markdup.txt"), "-@", THREADS, srt, md, log=run.log)
        tools.samtools("index", md, log=run.log)
        for f in (sam, fx, srt):
            f.unlink(missing_ok=True)
    return {"bam": md, "paired": paired, "aligner": tools.aligner_name()}


def _align_pybwa(run, ref, r1, r2, sam, progress=None, chunk=20_000):
    from pybwa import BwaIndex, BwaMem, BwaMemOptions
    run.log.append(f"$ pybwa BwaMem.align({ref.fasta.name}, reads) # library mode, pairs aligned as single-end reads")
    idx = BwaIndex(prefix=str(ref.fasta))
    aligner = BwaMem(index=idx)
    opt = BwaMemOptions(threads=THREADS)
    hdr = idx.header.to_dict()
    hdr["RG"] = [{"ID": run.sample, "SM": run.sample, "PL": "ILLUMINA", "LB": f"{run.sample}_lib"}]
    header = pysam.AlignmentHeader.from_dict(hdr)
    with pysam.AlignmentFile(str(sam), "w", header=header) as out:
        for mate, path in ((1, r1), (2, r2)):
            if path is None:
                continue
            batch = []
            with pysam.FastxFile(str(path)) as fq:
                for rec in fq:
                    batch.append(rec)
                    if len(batch) >= chunk:
                        _write_batch(out, header, aligner, opt, batch, mate, run.sample)
                        batch = []
                if batch:
                    _write_batch(out, header, aligner, opt, batch, mate, run.sample)


def _write_batch(out, header, aligner, opt, batch, mate, sample):
    res = aligner.align(batch, opt)
    for rec, alns in zip(batch, res):
        for a in alns:
            seg = pysam.AlignedSegment.fromstring(a.to_string(), header) if a.header is None or a.header != header else a
            seg.query_name = f"{rec.name}/{mate}"
            seg.set_tag("RG", sample)
            out.write(seg)


# ------------------------------------------------------------------ 3. BAM QC + coverage
def flagstat(bam: Path) -> dict:
    txt = tools.samtools("flagstat", bam)
    d = {}
    for line in txt.splitlines():
        m = re.match(r"(\d+) \+ \d+ (.+)", line)
        if m:
            d[m.group(2).split("(")[0].strip()] = int(m.group(1))
    total = d.get("in total", 0)
    mapped = d.get("mapped", d.get("primary mapped", 0))
    return {"total": total, "mapped": mapped, "mapped_pct": mapped / max(total, 1) * 100,
            "duplicates": d.get("duplicates", d.get("primary duplicates", 0)),
            "dup_pct": d.get("duplicates", 0) / max(total, 1) * 100,
            "properly_paired": d.get("properly paired", 0), "secondary": d.get("secondary", 0),
            "supplementary": d.get("supplementary", 0), "raw": txt}


def bam_stats(bam: Path, max_reads=200_000) -> dict:
    mapq, isize, ontarget = [], [], 0
    with pysam.AlignmentFile(str(bam)) as f:
        for i, r in enumerate(f.fetch(until_eof=True)):
            if i >= max_reads:
                break
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            mapq.append(r.mapping_quality)
            if r.is_proper_pair and r.is_read1 and 0 < abs(r.template_length) < 2000:
                isize.append(abs(r.template_length))
    return {"mapq": np.array(mapq), "insert_size": np.array(isize)}


def coverage(run: Run, ref: R.Reference, bam: Path, thresholds=(1, 10, 20, 30, 100)) -> dict:
    with run.step("Coverage (samtools depth)"):
        txt = tools.samtools("depth", "-a", "-b", ref.targets_bed, "-Q", 1, bam, log=run.log)
    from io import StringIO
    dp = pd.read_csv(StringIO(txt), sep="\t", header=None, names=["contig", "pos", "depth"]) if txt.strip() else \
        pd.DataFrame(columns=["contig", "pos", "depth"])
    tgt = ref.targets()
    rows, low = [], []
    for _, t in tgt.iterrows():
        d = dp[(dp.contig == t.contig) & (dp.pos > t.start) & (dp.pos <= t.end)]
        vals = d.depth.to_numpy() if len(d) else np.zeros(int(t.end - t.start))
        row = {"target": t["name"], "contig": t.contig, "start": int(t.start), "end": int(t.end), "length": int(t.end - t.start),
               "mean_depth": float(vals.mean()) if len(vals) else 0.0, "min_depth": int(vals.min()) if len(vals) else 0}
        for th in thresholds:
            row[f"pct_{th}x"] = float((vals >= th).mean() * 100) if len(vals) else 0.0
        rows.append(row)
        # contiguous low-coverage stretches
        if len(d):
            lowmask = d.depth.to_numpy() < 20
            if lowmask.any():
                pos = d.pos.to_numpy()
                runs_ = np.split(pos[lowmask], np.where(np.diff(pos[lowmask]) != 1)[0] + 1)
                for rr in runs_:
                    if len(rr) >= 5:
                        g = ref.to_genomic(t.contig, int(rr[0]))
                        low.append({"target": t["name"], "contig": t.contig, "start": int(rr[0]), "end": int(rr[-1]),
                                    "bp": len(rr), "genomic": f"{g[0]}:{g[1]}-{g[1] + len(rr) - 1}"})
    per_target = pd.DataFrame(rows)
    per_target["gene"] = per_target.target.str.split("_").str[0]
    w = per_target.length
    per_gene = per_target.groupby("gene").apply(
        lambda g: pd.Series({"targets": len(g), "bp": int(g.length.sum()),
                             "mean_depth": float(np.average(g.mean_depth, weights=g.length)),
                             **{f"pct_{th}x": float(np.average(g[f'pct_{th}x'], weights=g.length)) for th in thresholds}}),
        include_groups=False).reset_index()
    overall = {"mean_depth": float(np.average(per_target.mean_depth, weights=w)) if len(per_target) else 0,
               **{f"pct_{th}x": float(np.average(per_target[f'pct_{th}x'], weights=w)) if len(per_target) else 0 for th in thresholds}}
    return {"per_target": per_target, "per_gene": per_gene, "overall": overall, "low": pd.DataFrame(low), "depth": dp}


def on_target_rate(bam: Path, ref: R.Reference) -> float:
    tgt = ref.targets()
    total = on = 0
    with pysam.AlignmentFile(str(bam)) as f:
        for _, t in tgt.iterrows():
            on += f.count(t.contig, max(0, int(t.start)), int(t.end), read_callback="all")
        total = f.mapped
    return on / max(total, 1) * 100


# ------------------------------------------------------------------ 4. variant calling
def call_variants(run: Run, ref: R.Reference, bam: Path, min_bq=20, min_mq=20, max_depth=5000, pad=50) -> Path:
    tgt = ref.targets()
    padded = run.p("targets_padded.bed")
    t = tgt.copy()
    t["start"] = (t.start - pad).clip(lower=0)
    t["end"] = t.end + pad
    t.sort_values(["contig", "start"]).to_csv(padded, sep="\t", header=False, index=False)
    with run.step("Variant calling (bcftools mpileup + call)"):
        pile, raw = run.p("pileup.bcf"), run.p("raw.vcf.gz")
        tools.bcftools("mpileup", "-f", ref.fasta, "-a", "FORMAT/AD,FORMAT/DP,INFO/AD", "-q", min_mq, "-Q", min_bq,
                       "-d", max_depth, "-R", padded, "-Ob", "-o", pile, bam, log=run.log)
        tools.bcftools("call", "-mv", "--ploidy", "GRCh38", "-Oz", "-o", raw, pile, log=run.log)
    with run.step("Normalise (bcftools norm)"):
        norm = run.p("normalized.vcf.gz")
        tools.bcftools("norm", "-f", ref.fasta, "-m", "-any", "-Oz", "-o", norm, raw, log=run.log)
        tools.bcftools("index", "-f", norm, log=run.log)
        run.p("pileup.bcf").unlink(missing_ok=True)
    return norm


DEFAULT_FILTERS = {"min_qual": 30.0, "min_dp": 10, "min_vaf": 0.20, "min_gq": 20}


def filter_variants(run: Run, vcf_in: Path, f: dict = DEFAULT_FILTERS) -> Path:
    out = run.p(f"{run.sample}.filtered.vcf.gz")
    with run.step("Filter variants"):
        run.log.append(f"# soft filters: QUAL<{f['min_qual']} -> LowQual; DP<{f['min_dp']} -> LowDP; "
                       f"VAF<{f['min_vaf']} -> LowVAF; GQ<{f['min_gq']} -> LowGQ")
        vin = pysam.VariantFile(str(vcf_in))
        hdr = vin.header
        for fid, desc in (("LowQual", f"QUAL < {f['min_qual']}"), ("LowDP", f"Depth < {f['min_dp']}"),
                          ("LowVAF", f"Variant allele fraction < {f['min_vaf']}"), ("LowGQ", f"Genotype quality < {f['min_gq']}")):
            if fid not in hdr.filters:
                hdr.filters.add(fid, None, None, desc)
        with pysam.VariantFile(str(out), "wz", header=hdr) as vout:
            for rec in vin:
                s = rec.samples[0]
                dp = s.get("DP") or rec.info.get("DP") or 0
                ad = s.get("AD")
                vaf = (ad[1] / sum(ad)) if ad and len(ad) > 1 and sum(a or 0 for a in ad) else None
                gq = s.get("GQ")
                fl = []
                if rec.qual is not None and rec.qual < f["min_qual"]:
                    fl.append("LowQual")
                if dp < f["min_dp"]:
                    fl.append("LowDP")
                if vaf is not None and vaf < f["min_vaf"]:
                    fl.append("LowVAF")
                if gq is not None and gq < f["min_gq"]:
                    fl.append("LowGQ")
                if fl:
                    for x in fl:
                        rec.filter.add(x)
                else:
                    rec.filter.add("PASS")
                vout.write(rec)
        tools.bcftools("index", "-f", out, log=run.log)
    return out


def vcf_to_df(vcf: Path, ref: R.Reference) -> pd.DataFrame:
    rows = []
    with pysam.VariantFile(str(vcf)) as v:
        for rec in v:
            s = rec.samples[0]
            ad = s.get("AD") or (None, None)
            dp = s.get("DP") or rec.info.get("DP") or 0
            tot = sum(a for a in ad if a) if ad and ad[0] is not None else 0
            alt_ad = ad[1] if ad and len(ad) > 1 and ad[1] is not None else 0
            gt = s.get("GT")
            gts = "/".join("." if x is None else str(x) for x in gt) if gt else "./."
            chrom, gpos = ref.to_genomic(rec.chrom, rec.pos)
            alt = rec.alts[0] if rec.alts else "."
            kind = "SNV" if len(rec.ref) == len(alt) == 1 else ("deletion" if len(rec.ref) > len(alt) else "insertion" if len(alt) > len(rec.ref) else "MNV")
            rows.append({"contig": rec.chrom, "pos": rec.pos, "chrom": chrom, "gpos": gpos, "ref": rec.ref, "alt": alt,
                         "type": kind, "qual": round(rec.qual or 0, 1), "filter": ";".join(rec.filter.keys()) or "PASS",
                         "gt": gts, "zygosity": "hom" if gts in ("1/1",) else "het" if "1" in gts else "ref",
                         "dp": int(dp), "alt_reads": int(alt_ad), "vaf": round(alt_ad / tot, 3) if tot else None,
                         "gq": s.get("GQ"), "variant_id": f"{chrom}-{gpos}-{rec.ref}-{alt}"})
    return pd.DataFrame(rows)


def titv(df: pd.DataFrame) -> float | None:
    snv = df[df.type == "SNV"]
    if snv.empty:
        return None
    ts = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    n_ts = sum((r, a) in ts for r, a in zip(snv.ref, snv.alt))
    n_tv = len(snv) - n_ts
    return n_ts / n_tv if n_tv else None


# ------------------------------------------------------------------ 5. validation vs truth set
def compare_truth(run: Run, ref: R.Reference, truth: pd.DataFrame, calls: pd.DataFrame) -> dict:
    """Normalise the truth set the same way as the calls, then compute sensitivity / precision."""
    tv = run.p("truth.vcf")
    fai = pd.read_csv(str(ref.fasta) + ".fai", sep="\t", header=None)
    with open(tv, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c, l in fai[[0, 1]].values:
            fh.write(f"##contig=<ID={c},length={l}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, r in truth.sort_values(["contig", "pos"]).iterrows():
            fh.write(f"{r.contig}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\tPASS\t.\n")
    tn = run.p("truth.norm.vcf")
    tools.bcftools("norm", "-f", ref.fasta, "-m", "-any", "-o", tn, tv, log=run.log)
    tset = set()
    with pysam.VariantFile(str(tn)) as v:
        for rec in v:
            tset.add((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
    passing = calls[calls["filter"] == "PASS"]
    cset = {(r.contig, r.pos, r.ref, r.alt) for r in passing.itertuples()}
    tp, fn, fp = tset & cset, tset - cset, cset - tset
    return {"truth": len(tset), "called_pass": len(cset), "TP": len(tp), "FN": len(fn), "FP": len(fp),
            "sensitivity": len(tp) / max(len(tset), 1) * 100, "precision": len(tp) / max(len(cset), 1) * 100,
            "missed": sorted(fn), "false_pos": sorted(fp)}
