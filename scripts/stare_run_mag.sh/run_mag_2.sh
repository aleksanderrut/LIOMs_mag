#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs

LOG="liom_mag_dane/batch_logs/run_grid_M3_d08_d15_d2_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1

for DELTA in 0.8 1.5; do
for DELTA_2 in 0.0 0.3; do

    echo
    echo "================================================================="
    echo "Start: M=3, d=$DELTA, d2=$DELTA_2, Jp=0.0, FId=yes, grid=60x60, eig=1:20"
    echo "Czas rozpoczęcia: $(date)"
    echo "================================================================="

    julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
        -d "$DELTA" \
        --delta-2 "$DELTA_2" \
        -M 3 \
        -T both \
        -P both \
        -F yes \
        -B both \
        --include-fermion-identity yes \
        --J-prime 0.0 \
        --grid-omega 60 \
        --grid-g 60 \
        --eig-first 1 \
        --eig-last 20 \
        --omega-min 0.0 \
        --omega-max 2.0 \
        --g-min 0.0 \
        --g-max 2.0

    echo
    echo "Zakończono: M=3, d=$DELTA, d2=$DELTA_2, grid=60x60, eig=1:20"
    echo "Czas zakończenia: $(date)"

done
done


echo
echo "================================================================="
echo "Zakończono wszystkie obliczenia: $(date)"
echo "================================================================="