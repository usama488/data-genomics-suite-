import shutil

import pandas as pd
import streamlit as st

import cx_reference as R
import cx_tools as tools
from cx_ui_common import tip


def render():
    st.title("🗂️ Reference genome (GRCh38)")
    tip("Reference genome")
    refs = R.list_references()
    if refs:
        st.subheader("Available references")
        st.dataframe(pd.DataFrame([R.Reference(p).summary() for p in refs]), hide_index=True)
        c1, c2 = st.columns([3, 1])
        rm = c1.selectbox("Delete a reference", [p.name for p in refs], index=None, placeholder="choose...")
        if rm and c2.button("🗑️ Delete"):
            shutil.rmtree(R.REF_ROOT / rm)
            st.rerun()

    t1, t2 = st.tabs(["🧪 Gene panel (works on Streamlit Cloud)", "🧬 Full genome (own PC / server)"])
    with t1:
        st.write("Downloads the real GRCh38 sequence of each gene (± 300 bp) from Ensembl and the exons of its canonical "
                 "transcript as capture targets, then builds the BWA index.")
        preset = st.selectbox("Preset panel", list(R.PANELS))
        extra = st.text_input("Genes (comma separated, any HGNC symbols)", ", ".join(R.PANELS[preset]))
        genes = [g.strip().upper() for g in extra.split(",") if g.strip()]
        name = st.text_input("Reference name", "panel_" + "_".join(genes)[:60])
        online = R.ensembl_online()
        st.caption("Ensembl REST: " + ("🟢 online" if online else "🔴 unreachable — a synthetic offline reference will be built"))
        if st.button("Build panel reference", type="primary"):
            bar = st.progress(0.0)
            log = []
            try:
                ref = R.build_panel(genes, name=name, offline=not online, progress=lambda f, m: bar.progress(f, m), log=log)
            except Exception as e:
                st.error(f"Failed: {e}")
                st.stop()
            st.success(f"Built {ref.summary()}")
            with st.expander("Gene coordinates (GRCh38)"):
                st.dataframe(pd.DataFrame([{k: g[k] for k in ("symbol", "chrom", "start", "end", "strand", "transcript")}
                                           | {"exons": len(g["exons"])} for g in ref.genes]), hide_index=True)
            st.code("\n".join(log))
    with t2:
        st.warning("A full genome needs ~8 GB RAM and ~10 GB disk. It will not run on free Streamlit Cloud — use a Linux PC, "
                   "WSL on Windows, or a cloud VM.")
        st.markdown("**1. Download GRCh38 (analysis set, no alt contigs) — run in a terminal:**")
        st.code("wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
                "seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz\n"
                "gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz\n"
                "mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38.fa\n"
                "samtools faidx GRCh38.fa\n"
                "bwa index GRCh38.fa      # ~1 hour, ~5 GB RAM (or tick 'Build BWA index' below)", language="bash")
        fasta = st.text_input("Path to GRCh38 FASTA on this machine", "/data/GRCh38/GRCh38.fa")
        bed = st.file_uploader("Capture / exome target BED (optional — otherwise exons of the genes below are used)", type=["bed", "txt"])
        genes = st.text_input("Genes to analyse (used when no BED)", "EGFR, KRAS, NRAS, BRAF, PIK3CA, TP53, IDH1, KIT, ERBB2")
        build = st.checkbox("Build BWA index if missing (slow)")
        if st.button("Register full genome", type="primary"):
            bed_path = None
            if bed:
                R.REF_ROOT.mkdir(exist_ok=True)
                bed_path = R.REF_ROOT / "uploaded_targets.bed"
                bed_path.write_bytes(bed.getvalue())
            bar = st.progress(0.0)
            try:
                ref = R.register_genome(fasta, [g.strip().upper() for g in genes.split(",") if g.strip()], bed_path=bed_path,
                                        build_index=build, progress=lambda f, m: bar.progress(f, m))
                st.success(f"Registered {ref.summary()}")
            except Exception as e:
                st.error(str(e))
    st.caption(f"Aligner in this environment: **{tools.aligner_name()}**")
