"""Plain-language explanations of every step (Learn mode + Learn page)."""

STEP_INFO = {
    "FASTQ": ("Raw output of the sequencer. Every read has 4 lines: name, DNA sequence, '+', and a quality string where each "
              "character encodes a Phred score for one base.", "Paired-end = two files (R1, R2) reading the same DNA fragment from both ends."),
    "Phred quality": ("Q = -10·log10(error probability). Q20 = 1 error in 100 bases, Q30 = 1 in 1,000, Q40 = 1 in 10,000.",
                      "Clinical labs usually require ≥ 80% of bases at Q30."),
    "Adapter": ("Short synthetic DNA ligated to fragments so they stick to the flow cell. When the DNA fragment is shorter than the "
                "read length the sequencer reads into the adapter, which must be cut off before alignment.", ""),
    "Trimming": ("Removing adapters and low-quality ends (cutadapt / fastp). Reads that become too short are discarded.", ""),
    "Duplication": ("Identical reads from PCR copies of the same molecule. They are marked (not deleted) after alignment so they don't "
                    "inflate the evidence for a variant.", ""),
    "Reference genome": ("The standard human DNA sequence (GRCh38 / hg38) that reads are compared to.", ""),
    "Alignment": ("Finding where each read comes from in the reference genome. BWA-MEM is the standard aligner for Illumina reads.",
                  "Output: SAM/BAM file - one line per read with chromosome, position, CIGAR (match/insert/delete), mapping quality."),
    "BAM": ("Binary, compressed, sorted and indexed alignment file. The .bai index allows jumping to any region instantly.", ""),
    "MAPQ": ("Mapping quality: confidence that the read is placed correctly. MAPQ 60 = unique; 0 = could map equally well elsewhere.", ""),
    "Insert size": ("Length of the original DNA fragment, measured from the paired reads.", "Typical library: 200-400 bp."),
    "Coverage / depth": ("How many reads cover a base. Germline clinical panels usually require ≥ 20x (often ≥ 30x) at every target base.",
                         "Mean depth 100x but one exon at 5x is still a gap - always look at % bases ≥ 20x."),
    "On-target rate": ("Percent of reads falling inside the captured target regions - shows how efficient the capture was.", ""),
    "Variant calling": ("Comparing the pile of reads at each position with the reference to find differences, and computing genotype "
                        "likelihoods (bcftools, GATK HaplotypeCaller, DeepVariant).", ""),
    "VCF": ("Variant Call Format: one line per variant with CHROM, POS, REF, ALT, QUAL, FILTER, INFO and per-sample GT/AD/DP/GQ.", ""),
    "Genotype (GT)": ("0/1 = heterozygous (one copy), 1/1 = homozygous (both copies), 0/0 = reference.", ""),
    "VAF": ("Variant allele fraction = alt reads / total reads. Germline heterozygous ≈ 50%, homozygous ≈ 100%. "
            "Very low VAF suggests an artefact or mosaicism.", ""),
    "Normalisation": ("Left-aligning indels and splitting multi-allelic sites so the same variant is always written the same way - "
                      "essential before looking it up in databases.", ""),
    "QUAL / GQ": ("QUAL = confidence that a variant exists. GQ = confidence in the genotype. Both are Phred-scaled.", ""),
    "Ti/Tv": ("Transition (A↔G, C↔T) to transversion ratio. ~2.0-2.1 genome-wide, ~3.0 in exomes; very different values suggest false positives.", ""),
    "Annotation": ("Adding biological meaning: gene, transcript, effect on protein (HGVS), population frequency and clinical databases. "
                   "Ensembl VEP is the standard tool.", ""),
    "HGVS": ("Standard variant naming. c. = coding DNA position (c.68_69del), p. = protein change (p.Glu23ValfsTer17).", ""),
    "MANE transcript": ("Matched Annotation from NCBI and EMBL-EBI: the one agreed transcript per gene used for clinical reporting.", ""),
    "gnomAD": ("Population database of >800,000 people. Variants common in healthy people are unlikely to cause rare disease.", ""),
    "ClinVar": ("NCBI database where labs submit their classifications of variants (pathogenic, benign, VUS...).", ""),
    "In-silico predictors": ("Computer predictions of damage: REVEL and AlphaMissense (missense), SpliceAI (splicing), SIFT, PolyPhen, CADD.", ""),
    "ACMG classification": ("ACMG/AMP 2015 guideline: evidence codes (PVS1, PS, PM, PP for pathogenic; BA1, BS, BP for benign) are combined "
                            "into 5 classes: Pathogenic, Likely pathogenic, VUS, Likely benign, Benign.",
                            "ClinGen point system: very strong 8, strong 4, moderate 2, supporting 1. ≥10 P, 6-9 LP, 0-5 VUS."),
    "VUS": ("Variant of uncertain significance - not enough evidence either way; must not be used for clinical decisions.", ""),
    "Validation (truth set)": ("Labs prove their pipeline by sequencing samples with known variants and measuring sensitivity "
                               "(found / known) and precision (true / all called). The demo FASTQ has a truth set built in.", ""),
}

PIPELINE_DIAGRAM = """
FASTQ (R1, R2)
   │  ① Raw QC  ─ per-base quality, GC, adapters, duplicates
   ▼
② Trimming (cutadapt) ─ adapters + low-quality ends removed
   │  ③ Trimmed QC
   ▼
④ Alignment (BWA-MEM → GRCh38) ─ sort ─ mark duplicates ─ index  →  BAM
   │  ⑤ BAM QC ─ mapping %, duplicates, insert size, on-target, coverage ≥20x
   ▼
⑥ Variant calling (bcftools mpileup + call) ─ normalise  →  VCF
   │  ⑦ Filtering ─ QUAL, depth, VAF, GQ
   ▼
⑧ Annotation (Ensembl VEP) ─ gene, HGVS, gnomAD, ClinVar, REVEL/SpliceAI
   ▼
⑨ ACMG classification (auto + clinical scientist review)
   ▼
⑩ Clinical report (PDF) + VCF / BAM / TSV downloads
"""
