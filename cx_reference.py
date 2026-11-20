"""Reference genome management.

Two modes:
* Panel  - real GRCh38 sequence of selected genes, downloaded from the Ensembl REST API,
           one contig per gene (small enough for free cloud hosting).
* Genome - a full GRCh38 FASTA on disk (needs ~8 GB RAM and ~10 GB disk for the BWA index).

Every reference directory holds `reference.json` with the contig -> chromosome/offset map,
so variant positions can always be converted back to genomic GRCh38 coordinates.
"""
from __future__ import annotations

import json
import random
import time
from pathlib import Path

import pandas as pd
import requests

import cx_tools as tools

ENSEMBL = "https://rest.ensembl.org"
APP_DIR = Path(__file__).resolve().parent
REF_ROOT = APP_DIR / "references"

PANELS = {
    "Solid tumour hotspot panel (EGFR, KRAS, NRAS, BRAF, PIK3CA, TP53, IDH1, KIT, ERBB2)":
        ["EGFR", "KRAS", "NRAS", "BRAF", "PIK3CA", "TP53", "IDH1", "KIT", "ERBB2"],
    "Lung (EGFR, KRAS, BRAF, ERBB2, TP53)": ["EGFR", "KRAS", "BRAF", "ERBB2", "TP53"],
    "Colorectal (KRAS, NRAS, BRAF, PIK3CA, TP53)": ["KRAS", "NRAS", "BRAF", "PIK3CA", "TP53"],
    "Breast (PIK3CA, ERBB2, TP53, ESR1)": ["PIK3CA", "ERBB2", "TP53", "ESR1"],
    "Small demo (KRAS, BRAF, PIK3CA, TP53, ERBB2)": ["KRAS", "BRAF", "PIK3CA", "TP53", "ERBB2"],
}
LOF_GENES = {"TP53", "PTEN", "RB1", "NF1", "APC", "CDKN2A", "BRCA1", "BRCA2", "STK11", "SMAD4", "ARID1A", "KMT2D", "ATM"}

PAD = 300  # bp of flanking sequence around each gene


# ------------------------------------------------------------------ Ensembl REST
def _get(path: str, params: dict | None = None, accept="application/json", tries=4):
    url = ENSEMBL + path
    for i in range(tries):
        r = requests.get(url, params=params, headers={"Content-Type": accept, "Accept": accept}, timeout=60)
        if r.status_code == 429:
            time.sleep(float(r.headers.get("Retry-After", 1)) + .5)
            continue
        if r.ok:
            return r.json() if accept == "application/json" else r.text
        if i == tries - 1:
            r.raise_for_status()
        time.sleep(1 + i)
    raise RuntimeError(f"Ensembl request failed: {url}")


def ensembl_online() -> bool:
    try:
        return requests.get(ENSEMBL + "/info/ping", headers={"Content-Type": "application/json"}, timeout=8).ok
    except Exception:
        return False


def lookup_gene(symbol: str) -> dict:
    """Gene coordinates + canonical transcript exons (GRCh38)."""
    g = _get(f"/lookup/symbol/homo_sapiens/{symbol}", {"expand": 1})
    tx = [t for t in g.get("Transcript", []) if t.get("is_canonical")] or g.get("Transcript", [])[:1]
    t = tx[0]
    exons = sorted(((e["start"], e["end"]) for e in t.get("Exon", [])))
    cds = (t.get("Translation") or {})
    return {"symbol": symbol, "chrom": str(g["seq_region_name"]), "start": int(g["start"]), "end": int(g["end"]),
            "strand": int(g["strand"]), "gene_id": g["id"], "transcript": t["id"], "exons": exons,
            "cds_start": cds.get("start"), "cds_end": cds.get("end")}


def map_cds(transcript: str, start: int, end: int) -> dict | None:
    """Map coding-DNA (c.) positions of a transcript to GRCh38 genomic coordinates."""
    try:
        m = _get(f"/map/cds/{transcript.split('.')[0]}/{start}..{end}")
    except Exception:
        return None
    maps = [x for x in m.get("mappings", []) if x.get("gap", 0) == 0 and x.get("seq_region_name")]
    if not maps:
        return None
    return {"chrom": str(maps[0]["seq_region_name"]), "start": min(x["start"] for x in maps),
            "end": max(x["end"] for x in maps), "strand": int(maps[0]["strand"])}


def fetch_sequence(chrom: str, start: int, end: int) -> str:
    return _get(f"/sequence/region/human/{chrom}:{start}..{end}:1", accept="text/plain").strip().upper()


def clinvar_variants(chrom: str, start: int, end: int) -> list[dict]:
    """ClinVar-annotated short variants overlapping a region (used to spike demo data)."""
    try:
        data = _get(f"/overlap/region/human/{chrom}:{start}-{end}", {"feature": "variation"})
    except Exception:
        return []
    out = []
    for v in data:
        sig = v.get("clinical_significance") or []
        if v.get("alleles") and sig:
            out.append({"id": v.get("id"), "start": v["start"], "end": v["end"], "alleles": v["alleles"],
                        "strand": v.get("strand", 1), "clin_sig": sig, "consequence": v.get("consequence_type")})
    return out


# ------------------------------------------------------------------ reference objects
class Reference:
    def __init__(self, root: Path):
        self.root = Path(root)
        meta = json.loads((self.root / "reference.json").read_text())
        self.meta = meta
        self.fasta = Path(meta["fasta"])
        self.mode = meta["mode"]
        self.name = meta["name"]
        self.contigs = meta["contigs"]  # contig -> {chrom, offset}
        self.genes = meta.get("genes", [])

    # coordinates ------------------------------------------------------------
    def to_genomic(self, contig: str, pos: int) -> tuple[str, int]:
        c = self.contigs.get(contig)
        if not c:
            return contig, pos
        return c["chrom"], pos + c["offset"]

    def to_local(self, chrom: str, gpos: int) -> tuple[str, int] | None:
        for ctg, c in self.contigs.items():
            if c["chrom"] == chrom and c["offset"] < gpos <= c["offset"] + c["length"]:
                return ctg, gpos - c["offset"]
        return None

    @property
    def targets_bed(self) -> Path:
        return self.root / "targets.bed"

    @property
    def genes_bed(self) -> Path:
        return self.root / "genes.bed"

    def targets(self) -> pd.DataFrame:
        return pd.read_csv(self.targets_bed, sep="\t", header=None, names=["contig", "start", "end", "name"])

    def ready(self) -> bool:
        return self.fasta.exists() and tools.bwa_indexed(self.fasta) and Path(str(self.fasta) + ".fai").exists()

    def summary(self) -> dict:
        t = self.targets()
        return {"name": self.name, "mode": self.mode, "contigs": len(self.contigs),
                "reference_bp": sum(c["length"] for c in self.contigs.values()),
                "target_regions": len(t), "target_bp": int((t.end - t.start).sum()), "genes": ", ".join(g["symbol"] for g in self.genes),
                "indexed": self.ready()}


def list_references() -> list[Path]:
    REF_ROOT.mkdir(exist_ok=True)
    return sorted(p for p in REF_ROOT.iterdir() if (p / "reference.json").exists())


def _write_beds(root: Path, genes: list[dict], contig_of, local, pad_exon=20):
    rows_t, rows_g = [], []
    for g in genes:
        ctg = contig_of(g)
        for i, (s, e) in enumerate(g["exons"], 1):
            ls, le = local(g, s - pad_exon), local(g, e + pad_exon)
            exon_no = i if g["strand"] == 1 else len(g["exons"]) - i + 1
            rows_t.append((ctg, max(ls - 1, 0), le, f"{g['symbol']}_exon{exon_no}"))
        rows_g.append((ctg, max(local(g, g["start"]) - 1, 0), local(g, g["end"]), g["symbol"]))
    pd.DataFrame(rows_t).sort_values([0, 1]).to_csv(root / "targets.bed", sep="\t", header=False, index=False)
    pd.DataFrame(rows_g).to_csv(root / "genes.bed", sep="\t", header=False, index=False)


def build_panel(genes: list[str], name: str | None = None, progress=None, log=None, offline=False) -> Reference:
    """Download each gene (+ flanks) from Ensembl and build an indexed panel reference."""
    name = name or "panel_" + "_".join(genes)
    root = REF_ROOT / name
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "reference.fa"
    infos, contigs = [], {}
    with open(fasta, "w") as fh:
        for i, sym in enumerate(genes):
            if progress:
                progress(i / len(genes), f"Downloading {sym} from Ensembl (GRCh38)...")
            if offline:
                g, seq = _synthetic_gene(sym, i)
            else:
                g = lookup_gene(sym)
                seq = fetch_sequence(g["chrom"], g["start"] - PAD, g["end"] + PAD)
            ctg = f"{sym}_chr{g['chrom']}_{g['start'] - PAD}"
            offset = g["start"] - PAD - 1
            contigs[ctg] = {"chrom": g["chrom"], "offset": offset, "length": len(seq), "gene": sym}
            g["contig"] = ctg
            infos.append(g)
            fh.write(f">{ctg}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")
    if progress:
        progress(.9, "Indexing (samtools faidx, bwa index)...")
    tools.faidx(fasta, log)
    tools.bwa_index(fasta, log)
    _write_beds(root, infos, lambda g: g["contig"], lambda g, p: p - contigs[g["contig"]]["offset"])
    meta = {"name": name, "mode": "panel" + (" (synthetic offline)" if offline else ""), "build": "GRCh38",
            "fasta": str(fasta), "contigs": contigs, "genes": infos, "created": time.strftime("%Y-%m-%d %H:%M")}
    (root / "reference.json").write_text(json.dumps(meta, indent=1))
    if progress:
        progress(1.0, "Reference ready")
    return Reference(root)


def register_genome(fasta: str, genes: list[str], name="GRCh38_full", bed_path: str | None = None, build_index=False,
                    log=None, progress=None) -> Reference:
    """Register an existing full-genome FASTA. Targets come from a BED you supply or from gene lookups."""
    fasta = Path(fasta).expanduser().resolve()
    if not fasta.exists():
        raise FileNotFoundError(fasta)
    root = REF_ROOT / name
    root.mkdir(parents=True, exist_ok=True)
    if not Path(str(fasta) + ".fai").exists():
        if progress:
            progress(.05, "samtools faidx ...")
        tools.faidx(fasta, log)
    if not tools.bwa_indexed(fasta):
        if not build_index:
            raise RuntimeError("BWA index not found next to the FASTA. Tick 'build index' (needs ~5 GB RAM, ~1 h) "
                               "or run `bwa index reference.fa` yourself.")
        if progress:
            progress(.1, "bwa index (this takes ~1 hour for GRCh38) ...")
        tools.bwa_index(fasta, log)
    tools.seq_dict(fasta, log)
    fai = pd.read_csv(str(fasta) + ".fai", sep="\t", header=None, usecols=[0, 1], names=["c", "l"])
    chr_prefix = fai.c.astype(str).str.startswith("chr").any()
    contigs = {c: {"chrom": c.replace("chr", "") if chr_prefix else c, "offset": 0, "length": int(l)} for c, l in fai.values}
    infos = []
    if bed_path:
        bed = pd.read_csv(bed_path, sep="\t", header=None, comment="#", usecols=[0, 1, 2])
        bed[3] = [f"region{i + 1}" for i in range(len(bed))]
        bed.to_csv(root / "targets.bed", sep="\t", header=False, index=False)
        bed.to_csv(root / "genes.bed", sep="\t", header=False, index=False)
    else:
        for i, sym in enumerate(genes):
            if progress:
                progress(.5 + .4 * i / max(len(genes), 1), f"Looking up {sym} ...")
            g = lookup_gene(sym)
            g["contig"] = ("chr" if chr_prefix else "") + g["chrom"]
            infos.append(g)
        _write_beds(root, infos, lambda g: g["contig"], lambda g, p: p)
    meta = {"name": name, "mode": "full genome", "build": "GRCh38", "fasta": str(fasta), "contigs": contigs,
            "genes": infos, "created": time.strftime("%Y-%m-%d %H:%M")}
    (root / "reference.json").write_text(json.dumps(meta, indent=1))
    return Reference(root)


def _synthetic_gene(sym: str, i: int):
    """Offline fallback: random sequence with a fake exon structure (clearly labelled synthetic)."""
    rng = random.Random(hash(sym) % 10_000)
    length = 12_000
    seq = "".join(rng.choice("ACGT") for _ in range(length + 2 * PAD))
    start = 1_000_000 * (i + 1)
    exons, pos = [], start + 500
    while pos < start + length - 600:
        ln = rng.randint(90, 350)
        exons.append((pos, pos + ln))
        pos += ln + rng.randint(600, 1500)
    return ({"symbol": sym, "chrom": "SYN", "start": start, "end": start + length - 1, "strand": 1, "gene_id": f"SYN_{sym}",
             "transcript": f"SYN_{sym}-201", "exons": exons, "cds_start": exons[0][0], "cds_end": exons[-1][1]}, seq)
