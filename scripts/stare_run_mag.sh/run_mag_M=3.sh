#!/bin/bash

set -euo pipefail

THREADS=10
M=3
OMEGA=0.5
G=1.0

echo "Start wszystkich obliczeń: $(date)"

for DELTA in 0.8 1.5; do
    for J_PRIME in 0.0 0.25; do
        for FERMION_IDENTITY in no yes; do

            LOG_FILE="run_M${M}_d${DELTA}_w${OMEGA}_g${G}_Jp${J_PRIME}_FId${FERMION_IDENTITY}.log"

            echo "=================================================="
            echo "Start obliczenia: $(date)"
            echo "M=${M}"
            echo "delta=${DELTA}"
            echo "omega=${OMEGA}"
            echo "g=${G}"
            echo "J-prime=${J_PRIME}"
            echo "include-fermion-identity=${FERMION_IDENTITY}"
            echo "Log: ${LOG_FILE}"
            echo "=================================================="

            julia -t "$THREADS" --project=. lioms_mag_phonon.jl \
                -d "$DELTA" \
                -w "$OMEGA" \
                -g "$G" \
                -M "$M" \
                --J-prime "$J_PRIME" \
                --include-fermion-identity "$FERMION_IDENTITY" \
                -T both \
                -P both \
                -F yes \
                -B both \
                2>&1 | tee "$LOG_FILE"

            echo "Zakończono obliczenie: $(date)"
            echo

        done
    done
done

echo "Zakończono wszystkie obliczenia: $(date)"