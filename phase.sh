#!/usr/bin/env bash
set -euo pipefail

VCF="$1"
BAM="$2"
REFERENCE="$3"
OUT_DIR="$4"
THREADS="${5:-10}"
MODE="${6:-full}"

AUTOSOMES=$(printf "chr%s," {1..22} | sed 's/,$//')
FAST_REGION="chr1:1-10000000"

mkdir -p "$OUT_DIR"

base=$(basename "$VCF" .vcf.gz)
sample="${base%%.*}"

FILTERED_VCF="${OUT_DIR}/${sample}.filtered.vcf.gz"
PHASED_VCF="${OUT_DIR}/${sample}.phased.vcf.gz"
HAPLOTAGGED_BAM="${OUT_DIR}/${sample}.phased.bam"
STATS_TSV="${OUT_DIR}/${sample}.whatshap.stats.tsv"
MODKIT_BED="${OUT_DIR}/${sample}.modkit.bedmethyl.gz"

case "$MODE" in
    full)
        REGION="$AUTOSOMES"
        MODE_DESC="full autosomal analysis"
        ;;
    fast)
        REGION="$FAST_REGION"
        MODE_DESC="fast test mode (${FAST_REGION})"
        ;;
    *)
        echo "[ERROR] MODE must be 'full' or 'fast'" >&2
        exit 1
        ;;
esac

echo "========================================"
echo "[RUN INFO]"
echo "Date:              $(date)"
echo "Sample:            $sample"
echo "Mode:              $MODE_DESC"
echo "Input VCF:         $VCF"
echo "Input BAM:         $BAM"
echo "Reference:         $REFERENCE"
echo "Output directory:  $OUT_DIR"
echo "Threads:           $THREADS"
echo
echo "[FILTER PARAMS]"
echo "Region(s):         $REGION"
echo "Variant type:      SNPs only"
echo "Filters:           PASS"
echo "Genotype:          heterozygous only"
echo "AF range:          0.20-0.80"
echo "DP >=              30"
echo
echo "[WHATSHAP PARAMS]"
echo "phase:             --merge-reads --ignore-read-groups --distrust-genotypes"
echo "haplotag:          --ignore-read-groups --output-threads=$THREADS"
echo
echo "[MODKIT PARAMS]"
echo "pileup:            --phased --cpg --combine-strands --modified-bases 5mC"
echo
echo "[OUTPUT FILES]"
echo "Filtered VCF:      $FILTERED_VCF"
echo "Phased VCF:        $PHASED_VCF"
echo "Haplotagged BAM:   $HAPLOTAGGED_BAM"
echo "Stats TSV:         $STATS_TSV"
echo "Modkit BED:        $MODKIT_BED"
echo "========================================"
echo

for cmd in bcftools tabix whatshap samtools modkit; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "[ERROR] Missing $cmd"; exit 1; }
done

echo "[INFO] Indexing input VCF"
tabix -f -p vcf "$VCF"

echo "[INFO] Filtering VCF"
bcftools view \
    -r "$REGION" \
    -m2 -M2 \
    -v snps \
    -f PASS \
    -g het \
    -i 'FORMAT/AF[0:0]>=0.20 && FORMAT/AF[0:0]<=0.80 && FORMAT/DP[0]>=30' \
    -Oz \
    -o "$FILTERED_VCF" \
    "$VCF"

tabix -f -p vcf "$FILTERED_VCF"

echo "[INFO] WhatsHap phase"
whatshap phase \
    -o "$PHASED_VCF" \
    --reference "$REFERENCE" \
    --merge-reads \
    --ignore-read-groups \
    --distrust-genotypes \
    "$FILTERED_VCF" \
    "$BAM"

tabix -f -p vcf "$PHASED_VCF"

echo "[INFO] WhatsHap haplotag"
whatshap haplotag \
    --reference "$REFERENCE" \
    --ignore-read-groups \
    --output-threads "$THREADS" \
    -o "$HAPLOTAGGED_BAM" \
    "$PHASED_VCF" \
    "$BAM"

samtools index -@ "$THREADS" "$HAPLOTAGGED_BAM"

echo "[INFO] Generating WhatsHap stats"
whatshap stats \
    --tsv "$STATS_TSV" \
    "$PHASED_VCF"

echo "[INFO] Running modkit pileup"
modkit pileup \
    "$HAPLOTAGGED_BAM" \
    "$MODKIT_BED" \
    --prefix "$sample" \
    --threads "$THREADS" \
    --phased \
    --cpg \
    --bgzf \
    --combine-strands \
    --reference "$REFERENCE" \
    --modified-bases 5mC 5hmC

echo "[DONE] Pipeline completed successfully"
echo "Phased VCF:        $PHASED_VCF"
echo "Haplotagged BAM:   $HAPLOTAGGED_BAM"
echo "Stats TSV:         $STATS_TSV"
echo "Modkit BED:        $MODKIT_BED"
