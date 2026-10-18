#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

THREADS=30

export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"

mkdir -p liom_mag_dane/batch_logs

LOG="liom_mag_dane/batch_logs/run_M4_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1


# ============================================================
# OBLICZENIE 1: omega_0 = 1, g = 1
# ============================================================

echo "START 1x1: $(date)"

julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
    -M 3 \
    -d 0.001 \
    --delta-2 0.0 \
    --J-prime 0.0 \
    -T both \
    -P both \
    -F both \
    -B both \
    --include-fermion-identity yes \
    --grid-omega 1 \
    --grid-g 1 \
    --omega-min 1 \
    --omega-max 1 \
    --g-min 1 \
    --g-max 1 \
    --eig-first 1 \
    --eig-last 50

echo "KONIEC 1x1: $(date)"


# ============================================================
# OBLICZENIE 2: siatka 50x50, omega_0 i g od 0 do 5
# ============================================================

echo "START 50x50: $(date)"

julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
    -M 3 \
    -d 0.001 \
    --delta-2 0.0 \
    --J-prime 0.0 \
    -T both \
    -P both \
    -F both \
    -B both \
    --include-fermion-identity yes \
    --grid-omega 50 \
    --grid-g 50 \
    --omega-min 0 \
    --omega-max 5 \
    --g-min 0 \
    --g-max 5 \
    --eig-first 1 \
    --eig-last 50

echo "KONIEC 50x50: $(date)"