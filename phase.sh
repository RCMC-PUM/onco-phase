#!/usr/bin/env bash
set -euo pipefail

VCF="$1"
BAM="$2"
REFERENCE="$3"
OUT_DIR="$4"
THREADS="${5:-10}"
MODE="${6:-full}"
DP="${7:-30}"
AF_MIN="${8:-0.30}"
AF_MAX="${9:-0.70}"

AUTOSOMES=$(printf "chr%s," {1..22} | sed 's/,$//')
FAST_REGION="chr7:130386171-130606465"

base=$(basename "$VCF" .vcf.gz)
sample="${base%%.*}"

SAMPLE_DIR="${OUT_DIR}/${sample}"
mkdir -p "$SAMPLE_DIR"

FILTERED_VCF="${SAMPLE_DIR}/${sample}.filtered.vcf.gz"
PHASED_VCF="${SAMPLE_DIR}/${sample}.phased.vcf.gz"
HAPLOTAGGED_BAM="${SAMPLE_DIR}/${sample}.phased.bam"
STATS_TSV="${SAMPLE_DIR}/${sample}.whatshap.stats.tsv"
MODKIT_OUTPUT="${SAMPLE_DIR}/${sample}-modkit"

case "$MODE" in
    full)
        REGION="$AUTOSOMES"
        MODE_DESC="full autosomal analysis"

        HAPLOTAG_REGIONS=()
        for i in {1..22}; do
            HAPLOTAG_REGIONS+=(--regions "chr${i}")
        done
        ;;
    fast)
        REGION="$FAST_REGION"
        MODE_DESC="fast test mode (${FAST_REGION})"

        HAPLOTAG_REGIONS=(--regions "$FAST_REGION")
        ;;
    *)
        echo "[ERROR] MODE must be 'full' or 'fast'" >&2
        exit 1
        ;;
esac

echo "========================================"
echo "[RUN INFO]"
echo "Sample:            $sample"
echo "Mode:              $MODE_DESC"
echo "Output directory:  $SAMPLE_DIR"
echo "Threads:           $THREADS"
echo
echo "[FILTER PARAMS]"
echo "Region(s):         $REGION"
echo "AF range:          ${AF_MIN}-${AF_MAX}"
echo "DP >=              $DP"
echo "========================================"
echo

for cmd in bcftools tabix whatshap samtools modkit; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "[ERROR] Missing required command: $cmd" >&2
        exit 1
    }
done

mkdir -p "$MODKIT_OUTPUT"

echo "[INFO] Indexing input VCF"
tabix -f -p vcf "$VCF"

echo "[INFO] Indexing input BAM"
samtools index -@ "$THREADS" "$BAM"

echo "[INFO] Filtering VCF"
bcftools view \
    -r "$REGION" \
    -m2 -M2 \
    -v snps \
    -f PASS \
    -g het \
    -i "FORMAT/AF[0:0]>=${AF_MIN} && FORMAT/AF[0:0]<=${AF_MAX} && FORMAT/DP[0]>=${DP}" \
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
    "${HAPLOTAG_REGIONS[@]}" \
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
    "$MODKIT_OUTPUT" \
    --prefix "$sample" \
    --threads "$THREADS" \
    --phased \
    --cpg \
    --bgzf \
    --combine-strands \
    --reference "$REFERENCE" \
    --modified-bases 5mC 5hmC

echo "[DONE] Pipeline completed successfully"
