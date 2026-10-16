#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

THREADS=30

export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

mkdir -p liom_mag_dane/batch_logs

LOG="liom_mag_dane/batch_logs/run_M4_omega075_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1

echo
echo "================================================================="
echo "GRID: M=4, d=0.001, d2=0.0, grid=1x40, eig=1:200"
echo "omega: 0.75"
echo "g:     0 -> 10"
echo "Liczba punktow: 40"
echo "Czas rozpoczecia: $(date)"
echo "================================================================="

julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
    -M 4 \
    -d 0.001 \
    --delta-2 0.0 \
    --J-prime 0.0 \
    -T both \
    -P both \
    -F yes \
    -B both \
    --include-fermion-identity yes \
    --grid-omega 1 \
    --grid-g 40 \
    --omega-min 0.75 \
    --omega-max 0.75 \
    --g-min 0 \
    --g-max 10 \
    --eig-first 1 \
    --eig-last 200

echo
echo "================================================================="
echo "GRID 1x40 zakonczony"
echo "Czas zakonczenia: $(date)"
echo "================================================================="
