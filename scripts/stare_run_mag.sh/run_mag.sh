#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs
LOG="liom_mag_dane/batch_logs/run_single_M3_M4_g0_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1

for M in 3 4; do
for DELTA in 0.8 1.5; do

    echo
    echo "================================================================="
    echo "Start: M=$M, d=$DELTA, Jp=0.0, FId=yes, omega=0.5, g=0.0"
    echo "Czas rozpoczęcia: $(date)"
    echo "================================================================="

    julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
        -d "$DELTA" \
        -M "$M" \
        -T both \
        -P both \
        -F yes \
        -B both \
        --include-fermion-identity yes \
        --J-prime 0.0 \
        -w 0.5 \
        -g 0.0 \
        --grid-omega 1 \
        --grid-g 1

    echo
    echo "Zakończono: M=$M, d=$DELTA, Jp=0.0, FId=yes, g=0.0"
    echo "Czas zakończenia: $(date)"

done
done
echo
echo "================================================================="
echo "Zakończono wszystkie obliczenia: $(date)"
echo "================================================================="