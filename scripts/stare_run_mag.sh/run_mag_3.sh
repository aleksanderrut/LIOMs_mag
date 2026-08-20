#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs

LOG="liom_mag_dane/batch_logs/run_line_g_d08_d15_d2_03_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1


for DELTA in 0.8 1.5; do

    DELTA_2=0.3

    echo
    echo "================================================================="
    echo "Start: M=3, d=$DELTA, d2=$DELTA_2, Jp=0.0, FId=yes"
    echo "omega=0.5, g=0:2, 100 points, eig=1:20"
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
        --grid-omega 1 \
        --grid-g 100 \
        --eig-first 1 \
        --eig-last 20 \
        --omega-min 0.5 \
        --omega-max 0.5 \
        --g-min 0.0 \
        --g-max 2.0


    echo
    echo "Zakończono: M=3, d=$DELTA, d2=$DELTA_2"
    echo "Czas zakończenia: $(date)"

done


echo
echo "================================================================="
echo "Zakończono wszystkie obliczenia: $(date)"
echo "================================================================="