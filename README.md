# Phasing + Haplotagging + Methylation Pipeline

End-to-end pipeline for:

- VCF filtering (high-confidence heterozygous SNPs)
- Read-based phasing using WhatsHap
- Haplotagging BAM alignments
- Generating phased methylation profiles using modkit

Designed for ONT + 5-base / methylation-aware workflows.

---

## 📦 Requirements

Make sure the following tools are available in `$PATH`:

- bcftools
- tabix
- whatshap
- samtools
- modkit

---

## 🚀 Usage

```bash
./phase.sh <VCF> <BAM> <REFERENCE> <OUT_DIR> [THREADS] [MODE]
```

---

## ⚙️ Modes

### full
- Whole autosomes: chr1–chr22

### fast
- Region: chr1:1-10,000,000

---

## 🔬 Pipeline Steps

1. VCF indexing  
2. Variant filtering (SNPs, PASS, het, AF 0.2–0.8, DP ≥ 30)  
3. Phasing (WhatsHap)  
4. Haplotagging  
5. Phasing stats  
6. Phased methylation (modkit)  

---

## 📂 Outputs

- *.filtered.vcf.gz  
- *.phased.vcf.gz  
- *.phased.bam  
- *.whatshap.stats.tsv  
- *.modkit.bedmethyl.gz  

---

## 🧠 Notes

- AF filtering stabilizes phasing  
- DP ≥ 30 improves reliability  
- Designed for ONT / somatic-like data  

