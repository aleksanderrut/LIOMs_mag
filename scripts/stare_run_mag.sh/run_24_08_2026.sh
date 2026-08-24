#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs

LOG="liom_mag_dane/batch_logs/run_M3_d0001_d2_00_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1


echo
echo "================================================================="
echo "M=3, Delta=0.001, Delta2=0.0, Jp=0.0, FId=yes"
echo "Start wszystkich obliczen: $(date)"
echo "================================================================="


# =============================================================================
# 1. TEST POJEDYNCZEGO PUNKTU
# omega_0 = 0.5
# g = 0.5
# BEZ GRIDU
# =============================================================================

echo
echo "================================================================="
echo "TEST: M=3, d=0.001, d2=0.0, omega=0.5, g=0.5"
echo "Czas rozpoczecia: $(date)"
echo "================================================================="

julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
    -d 0.001 \
    --delta-2 0.0 \
    -M 3 \
    -T both \
    -P both \
    -F yes \
    -B both \
    --include-fermion-identity yes \
    --J-prime 0.0 \
    -w 0.5 \
    -g 0.5 \
    --grid-omega 1 \
    --grid-g 1

echo
echo "================================================================="
echo "TEST zakonczony poprawnie"
echo "Czas zakonczenia: $(date)"
echo "================================================================="


# =============================================================================
# 2. GRID 80 x 80
# omega_0 = 0 -> 2
# g       = 0 -> 2
# eig     = 1 -> 40
# =============================================================================

echo
echo "================================================================="
echo "GRID: M=3, d=0.001, d2=0.0, grid=80x80, eig=1:40"
echo "omega: 0 -> 2"
echo "g:     0 -> 2"
echo "Liczba punktow: 4900"
echo "Czas rozpoczecia: $(date)"
echo "================================================================="

julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
    -d 0.001 \
    --delta-2 0.0 \
    -M 3 \
    -T both \
    -P both \
    -F yes \
    -B both \
    --include-fermion-identity yes \
    --J-prime 0.0 \
    -w 0.5 \
    -g 0.5 \
    --grid-omega 80 \
    --grid-g 80 \
    --eig-first 1 \
    --eig-last 40 \
    --omega-min 0.0 \
    --omega-max 2.0 \
    --g-min 0.0 \
    --g-max 2.0

echo
echo "================================================================="
echo "GRID 80x80 zakonczony"
echo "Czas zakonczenia: $(date)"
echo "================================================================="

echo
echo "================================================================="
echo "Zakonczono wszystkie obliczenia: $(date)"
echo "================================================================="