#!/usr/bin/env bash
# Re-quantify one GSE294180 sequencing run with the authors' reported settings.
#
#   fastp v0.23.4 : adapter trimming, poly-G and poly-X removal   (GEO data processing)
#   Salmon v1.10.3: selective alignment, -l A, --seqBias, --useVBOpt, GENCODE M35 (GEO)
#
# Reads come from the NCBI SRA open-data mirror on AWS (about 5x faster than ENA here).
# The .sra file is downloaded to the work directory, streamed through fastq-dump
# (interleaved pairs) -> fastp -> two FIFOs -> Salmon, then deleted. No FASTQ is
# written to disk. Run inside WSL:
#
#   bash pipelines/quantification/quantify_run.sh SRR33050963 /mnt/d/mpliver_work [label]
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

# SRA_LOCAL: a .sra file already downloaded (e.g. in a browser). It is used as is
# and left in place; otherwise the file is fetched here and removed afterwards.
if [[ -n "${SRA_LOCAL:-}" ]]; then
    SRA_FILE="$SRA_LOCAL"
    CLEANUP=""
else
    curl -sS --retry 5 --retry-delay 10 -o "$SRA_FILE" \
        "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN"
    CLEANUP="$SRA_FILE"
fi

FIFO_DIR="$(mktemp -d)"
trap 'rm -rf "$FIFO_DIR" $CLEANUP' EXIT
mkfifo "$FIFO_DIR/r1.fq" "$FIFO_DIR/r2.fq"

# fastq-dump --split-spot writes mate 1 and mate 2 as consecutive records (interleaved).
"$MM" run -n mpliver-quant bash -c "
    fastq-dump --split-spot --skip-technical --stdout '$SRA_FILE' 2> '$OUT/fastq_dump.log' |
    fastp --stdin --interleaved_in \
        -o '$FIFO_DIR/r1.fq' -O '$FIFO_DIR/r2.fq' \
        --detect_adapter_for_pe --trim_poly_g --trim_poly_x \
        --thread 4 --json '$OUT/fastp.json' --html '$OUT/fastp.html' 2> '$OUT/fastp.log'
" &
STREAM_PID=$!

"$MM" run -n mpliver-quant salmon quant \
    -i "$INDEX" -l A \
    -1 "$FIFO_DIR/r1.fq" -2 "$FIFO_DIR/r2.fq" \
    --validateMappings --seqBias --useVBOpt \
    -p "$THREADS" -o "$OUT/salmon" \
    > "$OUT/salmon.log" 2>&1

wait "$STREAM_PID"
echo "$LABEL ($RUN) done: $(grep -o '"num_processed": [0-9]*' "$OUT/salmon/aux_info/meta_info.json"),"\
" $(grep -o '"percent_mapped": [0-9.]*' "$OUT/salmon/aux_info/meta_info.json")"
