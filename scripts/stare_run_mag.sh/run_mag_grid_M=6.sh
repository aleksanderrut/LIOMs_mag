#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs
LOG="liom_mag_dane/batch_logs/run_grid_M4_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1

for DELTA in 0.8 1.5; do
for J_PRIME in 0.0 0.25; do
for FID in no yes; do
    echo
    echo "Start: d=$DELTA, Jp=$J_PRIME, FId=$FID, eig=1:10"
    echo "Czas rozpoczęcia: $(date)"
    julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
        -d "$DELTA" \
        -M 6 \
        -T both \
        -P both \
        -F yes \
        -B both \
        --include-fermion-identity "$FID" \
        --J-prime "$J_PRIME" \
        --grid-omega 1 \
        --grid-g 501 \
        --eig-first 1 \
        --eig-last 10 \
        --omega-min 0.5 \
        --omega-max 0.5 \
        --g-min 0.0 \
        --g-max 2.0
    echo "Zakończono: d=$DELTA, Jp=$J_PRIME, FId=$FID, eig=1:10"
    echo "Czas zakończenia: $(date)"
done
done
done
echo
echo "Zakończono wszystkie obliczenia: $(date)"