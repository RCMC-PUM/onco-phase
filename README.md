# Phasing + Haplotagging + Methylation Pipeline

End-to-end pipeline for:

- VCF filtering (high-confidence heterozygous SNPs)
- Read-based phasing (WhatsHap)
- Haplotagging BAM alignments
- Phased methylation profiling (modkit)

Designed for ONT + 5-base workflows.

---

## Requirements

Tools in `$PATH`:

- bcftools
- tabix
- whatshap
- samtools
- modkit

---

## Usage

```bash
./phase.sh <VCF> <BAM> <REFERENCE> <OUT_DIR> [THREADS] [MODE] [DP] [AF_MIN] [AF_MAX]
```

### Defaults

- THREADS = 10
- MODE = full
- DP = 30
- AF_MIN = 0.30
- AF_MAX = 0.70

---

## Modes

### full
- Autosomes: chr1–chr22

### fast
- Region: chr7:130386171-130606465 (MEST - imprinted gene)

---

## Steps

1. Index VCF + BAM  
2. Filter SNPs (PASS, het, AF + DP)  
3. Phase (WhatsHap)  
4. Haplotag BAM  
5. Stats  
6. modkit pileup  

---

## Output

All results in:

```
OUT_DIR/SAMPLE_NAME/
```

Files:

- *.filtered.vcf.gz  
- *.phased.vcf.gz  
- *.phased.bam  
- *.whatshap.stats.tsv  
- modkit/*  

---

## Notes

- Single-sample VCF expected  
- Uses FORMAT/AF and FORMAT/DP  
- Tunable AF/DP for somatic data
