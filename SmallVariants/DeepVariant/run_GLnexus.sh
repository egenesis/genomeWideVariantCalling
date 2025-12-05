#!/bin/bash
set -euo pipefail

# === INPUTS ===
REF_NAME=$1                      # e.g. Yuc104F, Duroc
THREADS=$2                      # number of threads

# === PATHS ===
INPUT_DIR="./data/${REF_NAME}/gvcfs"
OUTPUT_BCF="./data/${REF_NAME}/cohort.bcf"
OUTPUT_VCF="./data/${REF_NAME}/cohort.vcf.gz"

# === Run GLnexus inside Docker ===
docker run --rm \
    -v "$(pwd)/data:/data" \
    ghcr.io/dnanexus-rnd/glnexus:v1.4.1 \
    glnexus_cli --threads "$THREADS" --config DeepVariant "$INPUT_DIR"/*.gvcf.gz > "$OUTPUT_BCF"

# === Convert BCF to compressed VCF + index ===
bcftools view "$OUTPUT_BCF" | bgzip -@4 -c > "$OUTPUT_VCF"
tabix -p vcf "$OUTPUT_VCF"
