#!/bin/bash

### USAGE: ./biallelify.sh path/to/cohort.vcf.gz path/to/Sscrofa11.1.fasta ###
### PURPOSE: Splits multiallelic variants into biallelics and left-aligns alleles against a reference FASTA using bcftools. ###
#   Automatically manages a conda/mamba environment ("biallelify") with bcftools & samtools, indexes the reference if missing,
#   normalizes the input VCF, and produces a bgzipped, indexed output VCF (*.biallelic.norm.vcf.gz).

# Environment management
ENV_NAME="biallelify"

# Dynamically find conda base and source conda.sh
CONDA_BASE=$(conda info --base 2>/dev/null)
CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ ! -f "$CONDA_SH" ]]; then
  echo "Error: Could not find conda.sh. Is Conda or Mamba installed?"
  exit 1
fi

source "$CONDA_SH"

# Choose environment manager: mamba preferred, fallback to conda
if command -v mamba &>/dev/null; then
  ENV_MGR="mamba"
  RUN_CMD="mamba run -n $ENV_NAME"
else
  echo "Mamba not found, falling back to conda (slower)."
  ENV_MGR="conda"
  RUN_CMD="conda run -n $ENV_NAME"
fi

# Create env if missing
if ! conda info --envs | grep -q "^${ENV_NAME}[[:space:]]"; then
  echo "Creating environment '$ENV_NAME' using $ENV_MGR..."
  $ENV_MGR create -y -n "$ENV_NAME" -c bioconda -c conda-forge bcftools samtools
fi

# Check for required arguments
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input.vcf.gz> <reference.fa>"
  echo "Splits multiallelics into biallelics and left-aligns via bcftools norm"
  exit 1
fi

# Input files
INPUT_VCF="$1"
REFERENCE="$2"

# Check that files exist
if [[ ! -f "$INPUT_VCF" ]]; then
  echo "Error: Input VCF file not found: $INPUT_VCF"
  exit 1
fi
if [[ ! -f "$REFERENCE" ]]; then
  echo "Error: Reference FASTA not found: $REFERENCE"
  exit 1
fi

# Index the reference if needed
if [[ ! -f "${REFERENCE}.fai" ]]; then
  echo "Indexing reference FASTA..."
  $RUN_CMD samtools faidx "$REFERENCE"
fi

# Derive output path from input file location
BASENAME=$(basename "$INPUT_VCF" | sed 's/\.vcf\.gz$//; s/\.gz$//')
OUTDIR=$(dirname "$INPUT_VCF")
OUTFILE="${OUTDIR}/${BASENAME}.biallelic.norm.vcf.gz"
OUTINDEX="${OUTFILE}.tbi"

# Normalize with bcftools
echo "Running bcftools norm with -m -any and reference alignment..."
$RUN_CMD bcftools norm -f "$REFERENCE" -m -any -Oz -o "$OUTFILE" "$INPUT_VCF"

# Index the output
$RUN_CMD bcftools index "$OUTFILE"

echo "Done. Output written to:"
echo " - $OUTFILE"
echo " - $OUTINDEX"
