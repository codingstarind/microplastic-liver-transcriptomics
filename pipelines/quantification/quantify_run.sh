#!/usr/bin/env bash
# Re-quantify one GSE294180 sequencing run with the authors' reported settings.
#
#   fastp v0.23.4 : adapter trimming, poly-G and poly-X removal   (GEO data processing)
#   Salmon v1.10.3: selective alignment, -l A, --seqBias, --useVBOpt, GENCODE M35 (GEO)
#
# Input: local FASTQ pairs (FASTQ_R1/FASTQ_R2), a local .sra file (SRA_LOCAL), or the
# .sra file fetched from the NCBI SRA open-data mirror on AWS. fastp writes trimmed
# reads to temporary gzip files in <work>/trim_tmp/<label>, Salmon quantifies them,
# and the temporary files are deleted. (An earlier version streamed fastp -> two
# FIFOs -> Salmon, which deadlocked near the end of the input.) Run inside WSL:
#
#   FASTQ_R1=... FASTQ_R2=... bash pipelines/quantification/quantify_run.sh \
#       SRR33050963 /mnt/d/mpliver_work CD_1_rep1
#
# Assumptions not stated by the authors (recorded in the validation report):
# no decoy sequences, k = 31, --gencode header parsing, fastp defaults otherwise.
set -euo pipefail

RUN="$1"
WORK="${2:-/mnt/d/mpliver_work}"
LABEL="${3:-$RUN}"
THREADS="${THREADS:-8}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/mpliver-tools/mamba}"
MM="${MICROMAMBA:-$HOME/mpliver-tools/bin/micromamba}"
INDEX="$WORK/index/salmon_index_gencode_vM35_k31"
OUT="$WORK/requant/$LABEL"
SRA_DIR="$WORK/sra_tmp"
SRA_FILE="$SRA_DIR/$LABEL.sra"

mkdir -p "$OUT" "$SRA_DIR"
if [[ -f "$OUT/salmon/quant.sf" ]]; then
    echo "$LABEL already quantified; skipping"
    exit 0
fi

# Input, in order of preference:
#   FASTQ_R1/FASTQ_R2: local paired FASTQ(.gz) files (MD5-verified ENA downloads)
#   SRA_LOCAL: a .sra file already downloaded; used as is and left in place
#   otherwise: the .sra file is fetched from the AWS mirror and removed afterwards
FIFO_DIR="$(mktemp -d)"
if [[ -n "${FASTQ_R1:-}" && -n "${FASTQ_R2:-}" ]]; then
    # Trimmed reads go to (lightly compressed) temporary files, not FIFOs: with two
    # FIFOs, fastp and Salmon can deadlock near the end of the input, each blocked
    # on a different mate file.
    TRIM="$WORK/trim_tmp/$LABEL"
    mkdir -p "$TRIM"
    trap 'rm -rf "$FIFO_DIR" "$TRIM"' EXIT
    "$MM" run -n mpliver-quant fastp -i "$FASTQ_R1" -I "$FASTQ_R2" \
        -o "$TRIM/r1.fq.gz" -O "$TRIM/r2.fq.gz" -z 1 \
        --detect_adapter_for_pe --trim_poly_g --trim_poly_x \
        --thread 6 --json "$OUT/fastp.json" --html "$OUT/fastp.html" 2> "$OUT/fastp.log"
    "$MM" run -n mpliver-quant salmon quant \
        -i "$INDEX" -l A \
        -1 "$TRIM/r1.fq.gz" -2 "$TRIM/r2.fq.gz" \
        --validateMappings --seqBias --useVBOpt \
        -p "$THREADS" -o "$OUT/salmon" \
        > "$OUT/salmon.log" 2>&1
    echo "$LABEL ($RUN) done: $(grep -o '"num_processed": [0-9]*' "$OUT/salmon/aux_info/meta_info.json"),"\
" $(grep -o '"percent_mapped": [0-9.]*' "$OUT/salmon/aux_info/meta_info.json")"
    exit 0
fi
rmdir "$FIFO_DIR"

if [[ -n "${SRA_LOCAL:-}" ]]; then
    SRA_FILE="$SRA_LOCAL"
    CLEANUP=""
else
    curl -sS --retry 5 --retry-delay 10 -o "$SRA_FILE" \
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN"
    CLEANUP="$SRA_FILE"
fi

TRIM="$WORK/trim_tmp/$LABEL"
mkdir -p "$TRIM"
trap 'rm -rf "$TRIM" $CLEANUP' EXIT

# fastq-dump --split-spot writes mate 1 and mate 2 as consecutive records (interleaved).
# Trimmed reads go to temporary files (two FIFOs can deadlock fastp and Salmon).
"$MM" run -n mpliver-quant bash -c "
    fastq-dump --split-spot --skip-technical --stdout '$SRA_FILE' 2> '$OUT/fastq_dump.log' |
    fastp --stdin --interleaved_in \
        -o '$TRIM/r1.fq.gz' -O '$TRIM/r2.fq.gz' -z 1 \
        --detect_adapter_for_pe --trim_poly_g --trim_poly_x \
        --thread 6 --json '$OUT/fastp.json' --html '$OUT/fastp.html' 2> '$OUT/fastp.log'
"

"$MM" run -n mpliver-quant salmon quant \
    -i "$INDEX" -l A \
    -1 "$TRIM/r1.fq.gz" -2 "$TRIM/r2.fq.gz" \
    --validateMappings --seqBias --useVBOpt \
    -p "$THREADS" -o "$OUT/salmon" \
    > "$OUT/salmon.log" 2>&1

echo "$LABEL ($RUN) done: $(grep -o '"num_processed": [0-9]*' "$OUT/salmon/aux_info/meta_info.json"),"\
" $(grep -o '"percent_mapped": [0-9.]*' "$OUT/salmon/aux_info/meta_info.json")"
