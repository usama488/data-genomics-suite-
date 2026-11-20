"""Orchestration for a tumour (+ matched normal) case, from FASTQ to tiered oncology report."""
from __future__ import annotations

import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

import cx_annotate as annotate
import cx_fastq_qc as fastq_qc
import cx_oncokb as oncokb
import cx_pipeline as P
import cx_reference as R
import cx_simulate as simulate
import cx_somatic as somatic
import cx_tools as tools

RUNS = R.APP_DIR / "runs_cancer"
STEPS = ["FASTQ QC", "Trimming", "Alignment + duplicates", "BAM QC + coverage", "Somatic calling", "Annotation (VEP)",
         "Copy number, purity & TMB", "Tiering (AMP/ASCO/CAP)"]


@dataclass
class Case:
    case_id: str
    ref_name: str | None = None
    mode: str = "paired"  # or "tumour-only"
    tumour_type: str = "Non-small cell lung cancer"
    fastq: dict = field(default_factory=dict)  # role -> (r1, r2)
    patient: dict = field(default_factory=dict)
    pathology_purity: float | None = None
    truth: pd.DataFrame | None = None
    truth_cn: dict = field(default_factory=dict)
    results: dict = field(default_factory=lambda: {"samples": {}})
    done: list = field(default_factory=list)
    log: list = field(default_factory=list)
    timings: dict = field(default_factory=dict)
    overrides: dict = field(default_factory=dict)
    settings: dict = field(default_factory=lambda: {"quality": 20, "min_len": 36, "offline": False, **somatic.DEFAULTS})

    @property
    def dir(self) -> Path:
        d = RUNS / self.case_id
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def roles(self) -> list[str]:
        return ["tumour", "normal"] if self.mode == "paired" else ["tumour"]

    def ref(self) -> R.Reference:
        return R.Reference(R.REF_ROOT / self.ref_name)

    def run(self, role: str) -> P.Run:
        r = P.Run(self.dir / role, f"{self.case_id}_{role}")
        r.log, r.timings = self.log, self.timings
        return r

    def reset_from(self, step: str):
        i = STEPS.index(step)
        self.done = [s for s in self.done if STEPS.index(s) < i]


def new_case(case_id: str) -> Case:
    d = RUNS / case_id
    if d.exists():
        shutil.rmtree(d)
    return Case(case_id=case_id)


def make_demo(case: Case, scenario: str, purity=.6, tumour_cov=300, normal_cov=100, progress=None):
    ref = case.ref()
    sc = simulate.SCENARIOS[scenario]
    case.tumour_type = sc["tumor"]
    if progress:
        progress(.02, "Mapping hotspot mutations to GRCh38 (Ensembl)...")
    truth, cn = simulate.design(ref, scenario, purity=purity, log=case.log)
    case.truth, case.truth_cn, case.pathology_purity = truth, cn, round(purity, 2)
    gf = simulate.tumor_gene_factor(cn, purity)
    t = simulate.simulate_sample(ref, truth, "tumor_af", case.dir / "fastq", f"{case.case_id}_T", tumour_cov, gf, seed=1,
                                 progress=(lambda f, m: progress(.05 + .6 * f, m)) if progress else None)
    case.fastq = {"tumour": tuple(map(str, t))}
    if case.mode == "paired":
        n = simulate.simulate_sample(ref, truth, "normal_af", case.dir / "fastq", f"{case.case_id}_N", normal_cov, seed=2,
                                     progress=(lambda f, m: progress(.65 + .35 * f, m)) if progress else None)
        case.fastq["normal"] = tuple(map(str, n))
    case.done, case.results = [], {"samples": {}}
    case.log.append(f"Demo: {scenario}; purity {purity}; {len(truth)} truth variants; CN changes {cn}")


def run_step(case: Case, step: str):
    ref, s, Rr = case.ref(), case.settings, case.results
    S = Rr["samples"]
    if step == "FASTQ QC":
        for role in case.roles:
            r1, r2 = case.fastq[role]
            S.setdefault(role, {})["qc_raw"] = {"R1": fastq_qc.analyse(r1, 100_000)} | ({"R2": fastq_qc.analyse(r2, 100_000)} if r2 else {})
    elif step == "Trimming":
        for role in case.roles:
            r1, r2 = case.fastq[role]
            S[role]["trim"] = P.trim(case.run(role), Path(r1), Path(r2) if r2 else None, simulate.ADAPTER1, simulate.ADAPTER2,
                                     s["quality"], s["min_len"])
            S[role]["qc_trim"] = {"R1": fastq_qc.analyse(S[role]["trim"]["r1"], 100_000)}
    elif step == "Alignment + duplicates":
        for role in case.roles:
            t = S[role]["trim"]
            S[role]["align"] = P.align(case.run(role), ref, t["r1"], t["r2"])
    elif step == "BAM QC + coverage":
        for role in case.roles:
            bam = S[role]["align"]["bam"]
            S[role]["flagstat"] = P.flagstat(bam)
            S[role]["bamstats"] = P.bam_stats(bam)
            S[role]["coverage"] = P.coverage(case.run(role), ref, bam)
            S[role]["on_target"] = P.on_target_rate(bam, ref)
    elif step == "Somatic calling":
        run = case.run("somatic")
        nb = S["normal"]["align"]["bam"] if case.mode == "paired" else None
        pile = somatic.joint_pileup(run, ref, S["tumour"]["align"]["bam"], nb, s)
        cand = somatic.classify_sites(run, ref, pile, case.mode == "paired", s)
        Rr["calls"] = somatic.normalise(run, ref, cand)
        Rr["vcf"] = Rr["calls"].attrs.get("vcf") if len(Rr["calls"]) else None
        Path(pile).unlink(missing_ok=True)
        if case.truth is not None and len(Rr["calls"]):
            Rr["validation"] = somatic.compare_truth(Rr["calls"], case.truth)
    elif step == "Annotation (VEP)":
        calls = Rr["calls"]
        offline = s.get("offline") or "synthetic" in ref.mode
        Rr["annotated"] = annotate.annotate(calls, offline=offline, log=case.log, ref=ref) if len(calls) else calls
        a = Rr["annotated"]
        if case.mode != "paired" and len(a) and "gnomad_af" in a:  # population filter for tumour-only
            common = a.gnomad_af.fillna(0) > .001
            a.loc[common & a.status.str.startswith("somatic"), "status"] = "likely germline (in gnomAD)"
    elif step == "Copy number, purity & TMB":
        cn_n = S["normal"]["coverage"]["per_target"] if case.mode == "paired" else None
        a = Rr["annotated"]
        cnv0, _ = somatic.copy_number(S["tumour"]["coverage"]["per_target"], cn_n, None)
        est = somatic.estimate_purity(a, cnv0) if len(a) else None
        Rr["purity_estimate"] = est
        use = case.pathology_purity or est or 1.0
        Rr["cnv"], Rr["cnv_exons"] = somatic.copy_number(S["tumour"]["coverage"]["per_target"], cn_n, use)
        tbp = int(ref.targets().pipe(lambda t: (t.end - t.start).sum()))
        Rr["tmb"] = somatic.tmb(a, tbp) if len(a) else {"mutations": 0, "megabases": tbp / 1e6, "tmb": 0, "basis": "-",
                                                        "class": "TMB-Low", "reliable": tbp >= 1e6}
    elif step == "Tiering (AMP/ASCO/CAP)":
        Rr["tiered"] = tier(case)
        cnv = Rr["cnv"].copy()
        t = [oncokb.tier_cnv(g, c, case.tumour_type) for g, c in zip(cnv.gene, cnv.call)]
        for k in ("tier", "level", "therapies_here", "therapies_other", "rationale"):
            cnv[k] = [x[k] for x in t]
        Rr["cnv_tiered"] = cnv
    if step not in case.done:
        case.done.append(step)


def tier(case: Case) -> pd.DataFrame:
    a = case.results["annotated"]
    if a is None or a.empty:
        return pd.DataFrame()
    som = a[a.status.str.startswith("somatic") & (a["filter"] == "PASS")]
    return oncokb.tier_all(som, case.tumour_type, case.overrides)


def run_all(case: Case, progress=None):
    steps = [s for s in STEPS if s not in case.done]
    for i, st in enumerate(steps):
        if progress:
            progress(i / len(steps), f"{st} ...")
        t = time.time()
        run_step(case, st)
        case.timings[st + " (total)"] = round(time.time() - t, 2)
    if progress:
        progress(1.0, "Pipeline complete")


def qc_summary(case: Case) -> dict:
    out = {}
    for role, d in case.results["samples"].items():
        R_ = role.capitalize()
        if "trim" in d:
            out[f"{R_}: read pairs"] = f"{d['trim']['reads_in']:,}"
        if "flagstat" in d:
            out[f"{R_}: mapped / duplicates"] = f"{d['flagstat']['mapped_pct']:.1f}% / {d['flagstat']['dup_pct']:.1f}%"
        if "coverage" in d:
            o = d["coverage"]["overall"]
            out[f"{R_}: mean depth"] = f"{o['mean_depth']:.0f}x"
            out[f"{R_}: bases ≥ {'100' if role == 'tumour' else '20'}x"] = f"{o['pct_100x' if role == 'tumour' else 'pct_20x']:.1f}%"
    r = case.results
    if r.get("purity_estimate") is not None:
        out["Tumour purity (VAF-based estimate)"] = f"{r['purity_estimate'] * 100:.0f}%"
    if case.pathology_purity:
        out["Tumour purity (pathology)"] = f"{case.pathology_purity * 100:.0f}%"
    if "tmb" in r:
        out["TMB"] = f"{r['tmb']['tmb']:.1f} mut/Mb" + ("" if r["tmb"]["reliable"] else " (panel < 1 Mb - not reportable)")
    return out


def methods(case: Case) -> dict:
    ref = case.ref()
    v = tools.versions()
    s = case.settings
    txt = (f"{'Tumour and matched normal' if case.mode == 'paired' else 'Tumour-only'} DNA was sequenced (targeted panel). Reads were "
           f"trimmed with cutadapt and aligned to GRCh38 ({ref.mode}) with BWA-MEM; duplicates were marked with samtools. Somatic "
           f"variants were identified from a joint bcftools mpileup with minimum tumour depth {s['min_tumor_dp']}, ≥{s['min_tumor_alt']} "
           f"alternate reads, VAF ≥ {s['min_tumor_vaf']:.0%}"
           + (f", normal VAF ≤ {s['max_normal_vaf']:.0%} and one-sided Fisher exact p < {s['fisher_p']}" if case.mode == "paired" else
              "; germline variants were excluded using VAF patterns and gnomAD frequency (tumour-only)")
           + ". Variants were annotated with Ensembl VEP and tiered according to AMP/ASCO/CAP 2017 using a curated actionability "
             "knowledge base. Copy number was estimated from tumour/normal depth ratios (median-normalised) adjusted for tumour purity.")
    return {"genes": ", ".join(g["symbol"] for g in ref.genes), "text": txt,
            "versions": "Software: " + "; ".join(f"{k} {val}" for k, val in v.items())}
