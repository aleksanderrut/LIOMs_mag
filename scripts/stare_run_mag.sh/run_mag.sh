#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs

LOG="liom_mag_dane/batch_logs/run_M4_d08_d15_d2_03_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1


echo
echo "================================================================="
echo "M=4, Delta2=0.3, Jp=0.0, FId=yes"
echo "Start wszystkich obliczen: $(date)"
echo "================================================================="


# ================================================================
# 1. PUNKTY POJEDYNCZE
# omega_0 = 0.5
# g = 0
# ================================================================

echo
echo "================================================================="
echo "Punkty pojedyncze: omega=0.5, g=0"
echo "================================================================="


for DELTA in 0.8 1.5; do

    echo
    echo "Start punkt: Delta=$DELTA"
    echo "Czas: $(date)"

    julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
        -d "$DELTA" \
        --delta-2 0.3 \
        -M 4 \
        -T both \
        -P both \
        -F yes \
        -B both \
        --include-fermion-identity yes \
        --J-prime 0.0 \
        --grid-omega 1 \
        --grid-g 1 \
        --eig-first 1 \
        --eig-last 20 \
        --omega-min 0.5 \
        --omega-max 0.5 \
        --g-min 0.0 \
        --g-max 0.0

    echo
    echo "Zakonczono punkt: Delta=$DELTA"
    echo "Czas: $(date)"

done



# ================================================================
# 2. SCIEZKA PO g
# omega_0 = 0.5
# g = 0 -> 2
# 80 punktow
# ================================================================

echo
echo "================================================================="
echo "Sciezka g: omega=0.5, g=0-2, 80 punktow"
echo "================================================================="


for DELTA in 0.8 1.5; do

    echo
    echo "Start sciezka g: Delta=$DELTA"
    echo "Czas: $(date)"

    julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
        -d "$DELTA" \
        --delta-2 0.3 \
        -M 4 \
        -T both \
        -P both \
        -F yes \
        -B both \
        --include-fermion-identity yes \
        --J-prime 0.0 \
        --grid-omega 1 \
        --grid-g 80 \
        --eig-first 1 \
        --eig-last 30 \
        --omega-min 0.5 \
        --omega-max 0.5 \
        --g-min 0.0 \
        --g-max 2.0

    echo
    echo "Zakonczono sciezka g: Delta=$DELTA"
    echo "Czas: $(date)"

done



# ================================================================
# 3. HEATMAPA
# omega = 0 -> 2
# g = 0 -> 2
# grid 30x30
# ================================================================

echo
echo "================================================================="
echo "Heatmapa: grid 30x30"
echo "================================================================="


for DELTA in 0.8 1.5; do

    echo
    echo "Start heatmapa: Delta=$DELTA"
    echo "Czas: $(date)"

    julia -t "$THREADS" --project=. lioms_phonon_siatka_w_g.jl \
        -d "$DELTA" \
        --delta-2 0.3 \
        -M 4 \
        -T both \
        -P both \
        -F yes \
        -B both \
        --include-fermion-identity yes \
        --J-prime 0.0 \
        --grid-omega 25 \
        --grid-g 25 \
        --eig-first 1 \
        --eig-last 30 \
        --omega-min 0.0 \
        --omega-max 2.0 \
        --g-min 0.0 \
        --g-max 2.0

    echo
    echo "Zakonczono heatmapa: Delta=$DELTA"
    echo "Czas: $(date)"

done



echo
echo "================================================================="
echo "Zakonczono wszystkie obliczenia"
echo "Czas zakonczenia: $(date)"
echo "================================================================="