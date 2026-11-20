"""Thin layer over the bioinformatics tools.

Each tool uses the real command-line program when it is installed (bwa, samtools,
bcftools) and otherwise falls back to the same program compiled inside pip packages:
pysam bundles samtools + bcftools, pybwa bundles BWA-MEM, cutadapt is pure pip.
The arguments are identical either way, so results are the same.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path


class ToolError(RuntimeError):
    pass


def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _pysam_cmd(tool: str, sub: str):
    import pysam
    mod = pysam.samtools if tool == "samtools" else __import__("pysam.bcftools", fromlist=["x"])
    fn = getattr(mod, sub, None)
    if fn is None:
        raise ToolError(f"{tool} {sub} is not available in pysam")
    return fn


def run(tool: str, sub: str, *args: str, log: list | None = None) -> str:
    """Run `tool sub args...` (samtools / bcftools). Returns stdout as text."""
    args = [str(a) for a in args]
    line = f"{tool} {sub} {' '.join(args)}"
    if log is not None:
        log.append("$ " + line)
    if have(tool):
        p = subprocess.run([tool, sub, *args], capture_output=True, text=True)
        if p.returncode != 0:
            raise ToolError(f"{line}\n{p.stderr[-2000:]}")
        return p.stdout
    fn = _pysam_cmd(tool, sub)
    try:
        # pysam swallows "-o file" when catching stdout, so only capture when there is no -o
        out = fn(*args, catch_stdout="-o" not in args)
    except Exception as e:  # pysam raises SamtoolsError
        raise ToolError(f"{line}\n{e}") from e
    return out or ""


def samtools(sub, *args, log=None):
    return run("samtools", sub, *args, log=log)


def bcftools(sub, *args, log=None):
    return run("bcftools", sub, *args, log=log)


def faidx(fasta: str | Path, log=None):
    samtools("faidx", str(fasta), log=log)


def seq_dict(fasta: str | Path, log=None):
    fasta = Path(fasta)
    d = fasta.with_suffix(".dict")
    if not d.exists():
        samtools("dict", "-o", str(d), str(fasta), log=log)
    # pybwa looks for <fasta>.dict as well as <prefix>.dict
    alt = Path(str(fasta) + ".dict")
    if not alt.exists():
        alt.write_text(d.read_text())


def bwa_indexed(fasta: str | Path) -> bool:
    return all(Path(str(fasta) + ext).exists() for ext in (".bwt", ".pac", ".ann", ".amb", ".sa"))


def bwa_index(fasta: str | Path, log=None):
    fasta = str(fasta)
    if log is not None:
        log.append(f"$ bwa index {fasta}")
    if have("bwa"):
        p = subprocess.run(["bwa", "index", fasta], capture_output=True, text=True)
        if p.returncode != 0:
            raise ToolError(p.stderr[-2000:])
    else:
        from pybwa import BwaIndex
        BwaIndex.index(fasta=fasta)
    seq_dict(fasta, log)


def aligner_name() -> str:
    if have("bwa"):
        return "bwa mem (command line)"
    try:
        import pybwa  # noqa: F401
        return "bwa mem (pybwa library)"
    except Exception:
        return "not available"


def versions() -> dict:
    out = {}
    try:
        import pysam
        out["pysam (samtools/bcftools/htslib)"] = pysam.__version__
    except Exception:
        out["pysam"] = "missing"
    for t in ("bwa", "samtools", "bcftools"):
        if have(t):
            try:
                p = subprocess.run([t] if t == "bwa" else [t, "--version"], capture_output=True, text=True)
                txt = (p.stdout or p.stderr)
                ver = [l for l in txt.splitlines() if "Version" in l or l.startswith(t)]
                out[f"{t} (CLI)"] = (ver[0] if ver else "installed").strip()
            except Exception:
                out[f"{t} (CLI)"] = "installed"
        else:
            out[f"{t} (CLI)"] = "not installed - using bundled library"
    try:
        import pybwa  # noqa: F401
        out["pybwa"] = "installed"
    except Exception:
        out["pybwa"] = "missing"
    try:
        import cutadapt
        out["cutadapt"] = cutadapt.__version__
    except Exception:
        out["cutadapt"] = "missing"
    return out


def cutadapt(args: list[str], log=None) -> str:
    cmd = [sys.executable, "-m", "cutadapt", *map(str, args)]
    if log is not None:
        log.append("$ cutadapt " + " ".join(map(str, args)))
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise ToolError(p.stderr[-2000:])
    return p.stdout
