#!/bin/bash

# ------------------------
# Script: FalseAmpHunter
# Description: End-to-end pipeline for QC, assembly, alignment, and analysis
#              of false amplicons from paired-end NGS data.
# Usage:
#   FalseAmpHunter --fastq R1.fastq --fastq2 R2.fastq --out prefix \
#                  --refdb /path/to/hg38_blast_db --primer ACTGACTGACTG --k-mer 81
# ------------------------

set -e
set -u

# ------------------------
# Parse arguments
# ------------------------
R1=""
R2=""
OUT=""
REFDB=""
PRIMER=""
KMER=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --fastq)
      R1="$2"
      shift 2
      ;;
    --fastq2)
      R2="$2"
      shift 2
      ;;
    --out)
      OUT="$2"
      shift 2
      ;;
    --refdb)
      REFDB="$2"
      shift 2
      ;;
    --primer)
      PRIMER="$2"
      shift 2
      ;;
    --k-mer)
      KMER="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# ------------------------
# Validate required args
# ------------------------
if [[ -z "$R1" || -z "$R2" || -z "$OUT" || -z "$REFDB" || -z "$PRIMER" || -z "$KMER" ]]; then
  echo -e "\nError: Missing required arguments."
  echo "Usage: FalseAmpHunter --fastq R1.fastq --fastq2 R2.fastq --out prefix \\"
  echo "       --refdb /path/to/hg38_blast_db --primer ACTGACTGACTG --k-mer 81"
  exit 1
fi

# ------------------------
# Validate k-mer value
# ------------------------
if [[ "$KMER" -gt 128 ]]; then
  echo "Warning: k-mer value $KMER exceeds ABySS maximum of 128. Setting k-mer to 81."
  KMER=81
fi

# ------------------------
# Check dependencies
# ------------------------
for tool in fastqc prinseq bwa samtools blastn abyss-pe; do
  if ! command -v $tool &>/dev/null; then
    echo "Error: $tool is not installed or not in PATH."
    exit 1
  fi
done

# ------------------------
# Step 1: QC with PRINSEQ
# ------------------------
echo "[Step 1] Running PRINSEQ for quality control..."
prinseq -fastq "$R1" \
        -fastq2 "$R2" \
        -out_good "good_$OUT" \
        -out_bad "bad_$OUT" \
        -trim_qual_left 20 \
        -trim_qual_right 20 \
        -min_qual_mean 20 \
        -min_len 20

# ------------------------
# Step 2: FastQC on QC-passed reads
# ------------------------
echo "[Step 2] Running FastQC on QC-passed reads..."
fastqc "good_${OUT}_1.fastq" "good_${OUT}_2.fastq"

# ------------------------
# Step 3: De novo assembly with ABySS
# ------------------------
echo "[Step 3] Running ABySS de novo assembly with k=$KMER..."
abyss-pe name="$OUT" k="$KMER" B=6G in="good_${OUT}_1.fastq good_${OUT}_2.fastq"
mv "${OUT}-scaffolds.fa" "${OUT}_scaffold_k${KMER}.fasta"

# Check assembly produced output
if [[ ! -s "${OUT}_scaffold_k${KMER}.fasta" ]]; then
  echo "Error: ABySS produced no scaffolds. Assembly failed or input reads are insufficient."
  exit 1
fi

echo "Assembly complete: $(grep -c '^>' ${OUT}_scaffold_k${KMER}.fasta) scaffolds produced."

# ------------------------
# Step 4: Index and map reads to assembly
# ------------------------
echo "[Step 4] Mapping QC-passed reads to assembled scaffolds..."
bwa index "${OUT}_scaffold_k${KMER}.fasta"
bwa mem -t 8 "${OUT}_scaffold_k${KMER}.fasta" \
    "good_${OUT}_1.fastq" "good_${OUT}_2.fastq" > "${OUT}.sam"
samtools view -b "${OUT}.sam" -o "${OUT}.bam"
samtools sort "${OUT}.bam" -o "${OUT}_sorted.bam"
samtools index "${OUT}_sorted.bam"

# Coverage report with mapping quality filter >= 20
samtools coverage -q 20 "${OUT}_sorted.bam" > "${OUT}_coverage.txt"

# ------------------------
# Step 5: BLAST assembled scaffolds vs reference genome
# ------------------------
echo "[Step 5] Running BLAST against reference database: $REFDB..."
blastn -query "${OUT}_scaffold_k${KMER}.fasta" \
       -db "$REFDB" \
       -out "${OUT}_vs_ref.txt" \
       -outfmt 0 \
       -task blastn

# Also output tab-delimited format for programmatic use
blastn -query "${OUT}_scaffold_k${KMER}.fasta" \
       -db "$REFDB" \
       -out "${OUT}_vs_ref_tab.txt" \
       -outfmt 6 \
       -task blastn

# ------------------------
# Step 6: Search for primer matches in BLAST results
# ------------------------
echo "[Step 6] Searching for primer matches in BLAST hits..."
grep -i -B 5 -A 20 "$PRIMER" "${OUT}_vs_ref.txt" > "${OUT}_primer_matches.txt" || true

if [[ ! -s "${OUT}_primer_matches.txt" ]]; then
  echo "Warning: No primer matches found in BLAST results."
else
  echo "Primer matches found. See ${OUT}_primer_matches.txt for details."
fi

# ------------------------
# Cleanup intermediate files
# ------------------------
retain_files=(
    "good_${OUT}_1.fastq"
    "good_${OUT}_2.fastq"
    "good_${OUT}_1_fastqc.html"
    "good_${OUT}_2_fastqc.html"
    "${OUT}_scaffold_k${KMER}.fasta"
    "${OUT}_sorted.bam"
    "${OUT}_sorted.bam.bai"
    "${OUT}_coverage.txt"
    "${OUT}_vs_ref.txt"
    "${OUT}_vs_ref_tab.txt"
    "${OUT}_primer_matches.txt"
)

echo "Cleaning up intermediate files..."

for f in *; do
    skip=0
    for keep in "${retain_files[@]}"; do
        if [[ "$f" == "$keep" ]]; then
            skip=1
            break
        fi
    done
    if [[ "$skip" -eq 0 ]]; then
        rm -rf "$f"
    fi
done

echo "Cleanup complete."

# ------------------------
# Summary
# ------------------------
echo -e "\nFalseAmpHunter pipeline complete. Final outputs:"
echo "1.  QC-passed FASTQ        : good_${OUT}_1.fastq, good_${OUT}_2.fastq"
echo "2.  FastQC reports         : good_${OUT}_1_fastqc.html, good_${OUT}_2_fastqc.html"
echo "3.  Assembly               : ${OUT}_scaffold_k${KMER}.fasta"
echo "4.  Alignment              : ${OUT}_sorted.bam (+ .bai index)"
echo "5.  Coverage report        : ${OUT}_coverage.txt (MAPQ >= 20)"
echo "6.  BLAST results          : ${OUT}_vs_ref.txt (pairwise), ${OUT}_vs_ref_tab.txt (tabular)"
echo "7.  Primer hits            : ${OUT}_primer_matches.txt"
