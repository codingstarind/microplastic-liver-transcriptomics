#!/usr/bin/env bash
# Direction C: refit the factorial model under the pre-specified sensitivity variants.
# Run inside WSL from the repository root:
#   bash analysis/robustness/run_sensitivity.sh [primary_label]
# Every variant re-runs DESeq2 and the (fast) pathway analysis; outputs go to
# results/deseq2/<variant>/ and results/pathways/<variant>/.
set -euo pipefail

REPO="$(pwd)"
PRIMARY="${1:-provisional_deposited}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/mpliver-tools/mamba}"
MM="${MICROMAMBA:-$HOME/mpliver-tools/bin/micromamba}"
R() { "$MM" run -n mpliver-r Rscript "$@" > /dev/null 2>&1; }

fit() {   # fit <variant label> [extra args for run_deseq2_factorial.R]
    local label="$1"; shift
    if [[ ! -f "$REPO/results/deseq2/$label/contrast_summary.csv" ]]; then
        R analysis/deseq2/run_deseq2_factorial.R "$REPO" deposited "$label" "$@"
        R analysis/pathways/pathway_analysis.R "$REPO" "$label" fast
    fi
    echo "done: $label"
}

fit "${PRIMARY}__drop_MD_1" --drop=MD_1
for cov in muscle depth quality_proxy; do
    fit "${PRIMARY}__cov_${cov}" --covariate="$cov"
done
R analysis/deseq2/run_alternative_methods.R "$REPO" "$PRIMARY" && echo "done: edgeR/limma"

# Leave-one-mouse-out: every mouse, including those in the smallest group.
for sample in $(tail -n +2 metadata/sample_metadata.csv | cut -d, -f3 | sed 's/_S[0-9]*$//'); do
    fit "loo/${PRIMARY}__without_${sample}" --drop="$sample"
done
