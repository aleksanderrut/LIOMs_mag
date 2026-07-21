#!/bin/bash

set -e

THREADS=10

echo "Start wszystkich obliczeń: $(date)"

echo "Siatka omega-g: M=3, eig-index=3"
julia -t $THREADS --project=. lioms_phonon_siatka_w_g.jl \
    -d 0.3 \
    -M 3 \
    -T both \
    -P both \
    -F yes \
    -B both \
    --grid-omega 100 \
    --grid-g 100 \
    --eig-index 3 \
    --omega-min 0.0 \
    --omega-max 2.0 \
    --g-min 0.0 \
    --g-max 2.0 \
    2>&1 | tee run_mag_M3_eig3.log

echo "Zakończono eig-index=3: $(date)"

echo "Siatka omega-g: M=3, eig-index=4"
julia -t $THREADS --project=. lioms_phonon_siatka_w_g.jl \
    -d 0.3 \
    -M 3 \
    -T both \
    -P both \
    -F yes \
    -B both \
    --grid-omega 100 \
    --grid-g 100 \
    --eig-index 4 \
    --omega-min 0.0 \
    --omega-max 2.0 \
    --g-min 0.0 \
    --g-max 2.0 \
    2>&1 | tee run_mag_M3_eig4.log

echo "Zakończono eig-index=4: $(date)"

echo "Siatka omega-g: M=3, eig-index=5"
julia -t $THREADS --project=. lioms_phonon_siatka_w_g.jl \
    -d 0.3 \
    -M 3 \
    -T both \
    -P both \
    -F yes \
    -B both \
    --grid-omega 100 \
    --grid-g 100 \
    --eig-index 5 \
    --omega-min 0.0 \
    --omega-max 2.0 \
    --g-min 0.0 \
    --g-max 2.0 \
    2>&1 | tee run_mag_M3_eig5.log

echo "Zakończono eig-index=5: $(date)"
echo "Koniec wszystkich obliczeń: $(date)"
