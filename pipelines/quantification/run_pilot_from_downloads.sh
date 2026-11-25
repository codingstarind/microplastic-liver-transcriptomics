#!/usr/bin/env bash
# Pilot driver: quantify the four pilot runs from .sra files downloaded in a browser.
# Waits for each file in the Downloads folder to reach its full size, quantifies it
# (CD_1 twice, for the run-to-run noise floor), then deletes the downloaded file.
#
#   bash pipelines/quantification/run_pilot_from_downloads.sh /mnt/c/Users/Soham/Downloads
set -euo pipefail

DL="${1:-/mnt/c/Users/Soham/Downloads}"
WORK="${2:-/mnt/d/mpliver_work}"
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
declare -A LABELS=([SRR33050963]="CD_1_rep1 CD_1_rep2" [SRR33050950]="CD_MP_5"
                   [SRR33050959]="MD_1" [SRR33050946]="MD_MP_4")
pending=(SRR33050963 SRR33050950 SRR33050959 SRR33050946)
deadline=$(( $(date +%s) + 12 * 3600 ))

expected_size() {
    curl -sSI "https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1" |
        awk 'tolower($1)=="content-length:" {print $2+0}'
}
declare -A SIZE
for run in "${pending[@]}"; do SIZE[$run]=$(expected_size "$run"); done

while (( ${#pending[@]} )) && (( $(date +%s) < deadline )); do
    remaining=()
    for run in "${pending[@]}"; do
        f="$DL/$run"
        if [[ -f "$f" ]] && (( $(stat -c %s "$f") == ${SIZE[$run]} )); then
            ok=1
            for label in ${LABELS[$run]}; do
                SRA_LOCAL="$f" bash "$REPO/pipelines/quantification/quantify_run.sh" \
                    "$run" "$WORK" "$label" || ok=0
                [[ -s "$WORK/requant/$label/salmon/quant.sf" ]] || ok=0
            done
            if (( ok )); then rm -f "$f"; echo "$(date +%H:%M) $run quantified, download deleted"
            else echo "$(date +%H:%M) $run FAILED; file kept"; fi
        else
            remaining+=("$run")
        fi
    done
    pending=("${remaining[@]+"${remaining[@]}"}")
    (( ${#pending[@]} )) && sleep 60
done
(( ${#pending[@]} )) && echo "timed out waiting for: ${pending[*]}" || echo "pilot complete"
