"""End-to-end orchestration: one `Case` holds everything about a sample from FASTQ to report."""
from __future__ import annotations

import json
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

import gx_acmg as acmg
import gx_annotate as annotate
import gx_fastq_qc as fastq_qc
import gx_pipeline as P
import gx_reference as R
import gx_simulate as simulate
import gx_tools as tools

RUNS = R.APP_DIR / "runs_germline"

STEPS = ["Raw FASTQ QC", "Trimming", "Trimmed FASTQ QC", "Alignment + duplicates", "BAM QC + coverage",
         "Variant calling", "Filtering", "Annotation (VEP)", "ACMG classification"]


@dataclass
class Case:
    sample: str
    ref_name: str | None = None
    r1: str | None = None
    r2: str | None = None
    patient: dict = field(default_factory=dict)
    truth: pd.DataFrame | None = None
    results: dict = field(default_factory=dict)
    done: list = field(default_factory=list)
    log: list = field(default_factory=list)
    timings: dict = field(default_factory=dict)
    overrides: dict = field(default_factory=dict)  # variant_id -> {code: strength}
    settings: dict = field(default_factory=lambda: {"quality": 20, "min_len": 36, "adapter1": simulate.ADAPTER1,
                                                    "adapter2": simulate.ADAPTER2, **P.DEFAULT_FILTERS, "offline": False})

    @property
    def dir(self) -> Path:
        d = RUNS / self.sample
        d.mkdir(parents=True, exist_ok=True)
        return d

    def ref(self) -> R.Reference:
        return R.Reference(R.REF_ROOT / self.ref_name)

    def run(self) -> P.Run:
        r = P.Run(self.dir, self.sample)
        r.log = self.log
        r.timings = self.timings
        return r

    def reset_from(self, step: str):
        i = STEPS.index(step)
        self.done = [s for s in self.done if STEPS.index(s) < i]


def new_case(sample: str) -> Case:
    d = RUNS / sample
    if d.exists():
        shutil.rmtree(d)
    return Case(sample=sample)


def make_demo(case: Case, coverage=100, n_random=6, use_clinvar=True, progress=None):
    ref = case.ref()
    if progress:
        progress(.05, "Choosing variants to spike (ClinVar + random)...")
    truth = simulate.pick_variants(ref, n_random=n_random, use_clinvar=use_clinvar and not case.settings.get("offline"), log=case.log)
    r1, r2 = simulate.simulate(ref, truth, case.dir, sample=case.sample, coverage=coverage,
                               progress=(lambda f, m: progress(.1 + .9 * f, m)) if progress else None)
    case.r1, case.r2, case.truth = str(r1), str(r2), truth
    case.done, case.results = [], {}
    case.log.append(f"Simulated demo FASTQ: {len(truth)} spiked variants, {coverage}x target coverage")


def run_step(case: Case, step: str, progress=None):
    run, ref, s = case.run(), case.ref(), case.settings
    R_ = case.results
    if step == "Raw FASTQ QC":
        R_["qc_raw"] = {"R1": fastq_qc.analyse(case.r1)}
        if case.r2:
            R_["qc_raw"]["R2"] = fastq_qc.analyse(case.r2)
    elif step == "Trimming":
        R_["trim"] = P.trim(run, Path(case.r1), Path(case.r2) if case.r2 else None, s["adapter1"], s["adapter2"],
                            s["quality"], s["min_len"])
    elif step == "Trimmed FASTQ QC":
        t = R_["trim"]
        R_["qc_trim"] = {"R1": fastq_qc.analyse(t["r1"])}
        if t["r2"]:
            R_["qc_trim"]["R2"] = fastq_qc.analyse(t["r2"])
    elif step == "Alignment + duplicates":
        t = R_["trim"]
        R_["align"] = P.align(run, ref, t["r1"], t["r2"], progress=progress)
    elif step == "BAM QC + coverage":
        bam = R_["align"]["bam"]
        R_["flagstat"] = P.flagstat(bam)
        R_["bamstats"] = P.bam_stats(bam)
        R_["coverage"] = P.coverage(run, ref, bam)
        R_["on_target"] = P.on_target_rate(bam, ref)
    elif step == "Variant calling":
        R_["raw_vcf"] = P.call_variants(run, ref, R_["align"]["bam"])
    elif step == "Filtering":
        f = {k: s[k] for k in P.DEFAULT_FILTERS}
        R_["vcf"] = P.filter_variants(run, R_["raw_vcf"], f)
        R_["calls"] = P.vcf_to_df(R_["vcf"], ref)
        if case.truth is not None and len(case.truth):
            R_["validation"] = P.compare_truth(run, ref, case.truth, R_["calls"])
    elif step == "Annotation (VEP)":
        calls = R_["calls"]
        offline = s.get("offline") or "synthetic" in ref.mode
        R_["annotated"] = annotate.annotate(calls, offline=offline, log=case.log, ref=ref)
    elif step == "ACMG classification":
        R_["classified"] = classify_all(case)
    if step not in case.done:
        case.done.append(step)


def classify_all(case: Case) -> pd.DataFrame:
    df = case.results["annotated"].copy()
    crits, classes, pts, why = [], [], [], []
    for _, v in df.iterrows():
        auto = acmg.auto_criteria(v.to_dict())
        ov = case.overrides.get(v.variant_id)
        c = ov if ov is not None else auto
        cls, p = acmg.classify(c)
        crits.append(", ".join(f"{k}" + ("" if acmg.CRITERIA.get(k, ("", st))[1] == st else f"_{st}") for k, (st, _) in c.items()))
        why.append(" | ".join(f"{k}: {r}" for k, (_, r) in c.items()))
        classes.append(cls)
        pts.append(p)
    df["criteria"], df["classification"], df["points"], df["evidence"] = crits, classes, pts, why
    df["auto_or_manual"] = ["manual" if vid in case.overrides else "auto" for vid in df.variant_id]
    rank = {c: i for i, c in enumerate(acmg.ORDER)}
    df["_r"] = df.classification.map(rank)
    return df.sort_values(["_r", "gene", "gpos"]).drop(columns="_r").reset_index(drop=True)


def run_all(case: Case, progress=None, start_at: str | None = None):
    steps = STEPS[STEPS.index(start_at):] if start_at else [s for s in STEPS if s not in case.done]
    for i, st in enumerate(steps):
        if progress:
            progress(i / len(steps), f"{st} ...")
        t = time.time()
        run_step(case, st, progress=None)
        case.timings[st + " (total)"] = round(time.time() - t, 2)
    if progress:
        progress(1.0, "Pipeline complete")


def qc_summary(case: Case) -> dict:
    R_ = case.results
    out = {}
    if "trim" in R_:
        t = R_["trim"]
        out["Read pairs (raw)"] = f"{t['reads_in']:,}"
        out["Read pairs after trimming"] = f"{t['reads_out']:,} ({t['kept_pct']:.1f}%)"
    if "qc_raw" in R_:
        out["Bases ≥ Q30 (R1)"] = f"{R_['qc_raw']['R1']['metrics']['q30_pct']:.1f}%"
    if "flagstat" in R_:
        f = R_["flagstat"]
        out["Mapped reads"] = f"{f['mapped']:,} ({f['mapped_pct']:.2f}%)"
        out["Duplicate reads"] = f"{f['dup_pct']:.1f}%"
    if "on_target" in R_:
        out["On-target reads"] = f"{R_['on_target']:.1f}%"
    if "coverage" in R_:
        o = R_["coverage"]["overall"]
        out["Mean target depth"] = f"{o['mean_depth']:.0f}x"
        out["Target bases ≥ 20x"] = f"{o['pct_20x']:.2f}%"
        out["Target bases ≥ 30x"] = f"{o['pct_30x']:.2f}%"
    if "calls" in R_:
        c = R_["calls"]
        out["Variants (PASS / total)"] = f"{(c['filter'] == 'PASS').sum()} / {len(c)}"
        tv = P.titv(c[c["filter"] == "PASS"]) if len(c) else None
        out["Ti/Tv (PASS SNVs)"] = f"{tv:.2f}" if tv else "n/a"
    return out


def methods_text(case: Case) -> dict:
    ref = case.ref()
    v = tools.versions()
    genes = ", ".join(g["symbol"] for g in ref.genes) or "see target BED"
    text = (f"DNA sequencing reads (FASTQ) were quality-checked and trimmed for adapters and bases below Q{case.settings['quality']} "
            f"with cutadapt (reads &lt; {case.settings['min_len']} bp discarded). Reads were aligned to the GRCh38 reference "
            f"({ref.mode}) with BWA-MEM; duplicates were marked with samtools markdup. Variants were called with bcftools "
            f"mpileup/call (min base quality 20, min mapping quality 20), normalised (left-aligned, multi-allelic split) and "
            f"filtered (QUAL ≥ {case.settings['min_qual']}, depth ≥ {case.settings['min_dp']}, VAF ≥ {case.settings['min_vaf']}, "
            f"GQ ≥ {case.settings['min_gq']}). Variants were annotated with Ensembl VEP (MANE/canonical transcript, gnomAD, "
            f"ClinVar, SIFT, PolyPhen, REVEL/SpliceAI where available) and classified using ACMG/AMP 2015 criteria with the "
            f"ClinGen Bayesian point system.")
    return {"genes": genes, "text": text, "versions": "Software: " + "; ".join(f"{k} {val}" for k, val in v.items())}


def save_case(case: Case):
    meta = {"sample": case.sample, "ref": case.ref_name, "r1": case.r1, "r2": case.r2, "patient": case.patient,
            "done": case.done, "timings": case.timings}
    (case.dir / "case.json").write_text(json.dumps(meta, indent=1, default=str))
