"""Plain-language explanations for cancer genomics (Learn mode)."""

STEP_INFO = {
    "Somatic vs germline": ("Germline variants are inherited and present in every cell (VAF ~50% or ~100% in normal tissue). "
                            "Somatic mutations arise in the tumour only, so they are absent from the matched normal (blood/saliva).", ""),
    "Matched normal": ("A non-tumour sample from the same patient. Subtracting it removes inherited variants and many artefacts - "
                       "the gold standard for somatic calling.", ""),
    "Tumour-only": ("Without a normal sample, germline variants must be removed using population databases (gnomAD) and VAF "
                    "patterns, so rare inherited variants can be mistaken for somatic ones.", ""),
    "VAF": ("Variant allele fraction in the tumour. A clonal heterozygous mutation in a diploid region has VAF ≈ purity / 2.",
            "Purity 60% → clonal driver at ~30% VAF; a subclone in 40% of cancer cells → ~12%."),
    "Tumour purity": ("Fraction of cells in the sample that are cancer cells (pathologist estimate or computed from VAFs). Low purity "
                      "lowers VAFs and can hide mutations and copy-number changes.", ""),
    "Clonal / subclonal": ("Clonal mutations are in every cancer cell (usually early drivers). Subclonal mutations are in a subset "
                           "(often later events, e.g. resistance mutations like EGFR T790M).", ""),
    "Limit of detection": ("The lowest VAF the assay reliably detects (~5% here). Depends on depth: at 300x a 5% VAF gives ~15 alt reads.", ""),
    "Fisher exact test": ("Tests whether the alternate allele is significantly more frequent in the tumour than in the normal.", ""),
    "Strand bias": ("A real variant is seen on reads from both DNA strands; alt reads on only one strand suggest an artefact.", ""),
    "Copy number (CNV)": ("Gains or losses of DNA. Amplification of oncogenes (ERBB2, MET, EGFR) can be a drug target. Estimated from "
                          "tumour/normal depth ratios: CN = 2·(ratio − (1 − purity)) / purity.", "log2 ratio 1.0 ≈ 2-fold more DNA than normal."),
    "TMB": ("Tumour mutational burden = somatic coding mutations per megabase. TMB ≥ 10 mut/Mb is an approved biomarker for "
            "pembrolizumab in solid tumours, but needs a panel of ≥ 1 Mb to be measured reliably.", ""),
    "Driver vs passenger": ("Drivers give cancer cells a growth advantage (KRAS G12C, BRAF V600E); passengers are bystanders.", ""),
    "Hotspot": ("A codon mutated repeatedly across many tumours (KRAS G12, BRAF V600, PIK3CA H1047) - strong evidence of oncogenicity.", ""),
    "COSMIC": ("Catalogue Of Somatic Mutations In Cancer - the largest database of mutations observed in tumours.", ""),
    "AMP/ASCO/CAP tiers": ("2017 joint guideline. Tier I = strong significance (approved therapy/guideline, Level A/B); Tier II = potential "
                           "significance (approved in another cancer = Level C, trials/preclinical = Level D); Tier III = unknown; "
                           "Tier IV = benign.", ""),
    "Actionable": ("A biomarker linked to a treatment decision: sensitivity (e.g. KRAS G12C → sotorasib in NSCLC) or resistance "
                   "(KRAS mutation → no anti-EGFR antibody in colorectal cancer).", ""),
    "Resistance mutation": ("A mutation that makes a drug stop working, e.g. EGFR T790M after erlotinib/gefitinib - but it is targetable "
                            "with osimertinib.", ""),
    "Molecular tumour board": ("Oncologists, pathologists and scientists who review genomic results together to choose therapy.", ""),
    "Validation (truth set)": ("The demo spikes known mutations at known VAFs, so the app can report sensitivity and precision "
                               "like a laboratory validation run.", ""),
}

PIPELINE_DIAGRAM = """
 TUMOUR FASTQ                         NORMAL FASTQ (blood)
      │  QC + trimming (cutadapt)           │
      ▼                                     ▼
 BWA-MEM → sort → markdup → BAM       BWA-MEM → sort → markdup → BAM
      │  coverage (≥100x)                   │  coverage (≥20x)
      └──────────────┬──────────────────────┘
                     ▼
   Joint pileup (bcftools) → somatic rules + Fisher test
   (tumour alt ≥4, VAF ≥5%, normal VAF ≤2%, both strands)
                     ▼
   Annotation (Ensembl VEP: HGVS, COSMIC, gnomAD)
                     ▼
   Copy number (T/N depth ratio) · purity (VAF) · TMB
                     ▼
   AMP/ASCO/CAP tiering + therapy matching (knowledge base)
                     ▼
   Molecular tumour board report (PDF)
"""
