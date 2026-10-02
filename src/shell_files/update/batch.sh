#!/usr/bin/env bash
# Run from the directory containing ./main:
# nohup bash run_pending_updates.sh 3 > updates_queue.log 2>&1 &

set -u

if [[ $# -ne 1 || ! $1 =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: bash $0 <k: positive integer>" >&2
    exit 2
fi

tasks=(
    'uk2005 FOREST remove 0'
    'uk2005 FOREST remove 1'
    'uk2005 PPT remove 0'
    'it2004 FOREST remove 0'
    'it2004 FOREST remove 1'
    'it2004 PPT remove 0'
    'it2004 PPT remove 1'
    'friendster FOREST remove 0'
    'friendster FOREST remove 1'
    'friendster PPT remove 0'
    'friendster PPT remove 1'
    'it2004 FOREST insert 0'
    'it2004 FOREST insert 1'
    'it2004 PPT insert 0'
    'it2004 PPT insert 1'
    'friendster FOREST insert 0'
    'friendster FOREST insert 1'
)

# There are only 17 tasks; cap large k before shell integer conversion.
k=${#tasks[@]}
if [[ ${#1} -le 2 ]] && (( $1 < k )); then
    k=$1
fi
if [[ ! -x ./main ]]; then
    echo 'Error: run this script from the directory containing executable ./main.' >&2
    exit 2
fi
mkdir -p updates/{uk2005,it2004,friendster} || exit 2

pids=()
labels=()
next=0
failed=0
while (( next < ${#tasks[@]} || ${#pids[@]} > 0 )); do
    # Reap every finished task, regardless of its launch order.
    active_pids=()
    active_labels=()
    for i in "${!pids[@]}"; do
        pid=${pids[i]}
        if kill -0 "$pid" 2>/dev/null; then
            active_pids+=("$pid")
            active_labels+=("${labels[i]}")
        else
            if wait "$pid"; then
                echo "[DONE] ${labels[i]}"
            else
                status=$?
                echo "[FAIL exit=$status] ${labels[i]}" >&2
                failed=$((failed + 1))
            fi
        fi
    done
    pids=("${active_pids[@]}")
    labels=("${active_labels[@]}")

    while (( next < ${#tasks[@]} && ${#pids[@]} < k )); do
        read -r dataset index operation mode <<< "${tasks[next]}"
        logfile="updates/$dataset/${index,,}_${operation}_${mode}.txt"
        echo "[START] ${tasks[next]} -> $logfile"
        ./main "$dataset" "$index" exp-update "$operation" "$mode" \
            > "$logfile" 2>&1 < /dev/null &
        pids+=("$!")
        labels+=("${tasks[next]}")
        next=$((next + 1))
    done
    if (( ${#pids[@]} > 0 )); then
        sleep 0.2
    fi
done

echo "Finished ${#tasks[@]} tasks; failed: $failed."
(( failed == 0 ))
