#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

THREADS=10

mkdir -p liom_mag_dane/batch_logs
LOG="liom_mag_dane/batch_logs/run_M5_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "$LOG") 2>&1

for DELTA in 0.8 1.5; do
    for J_PRIME in 0.0 0.25; do
        for FID in no yes; do

            echo "Start: d=$DELTA, Jp=$J_PRIME, FId=$FID"

            julia -t "$THREADS" --project=. lioms_mag_phonon.jl \
                -d "$DELTA" \
                -w 0.5 \
                -g 1.0 \
                -M 5 \
                --J-prime "$J_PRIME" \
                --include-fermion-identity "$FID" \
                -T both \
                -P both \
                -F yes \
                -B both

            echo "Zakończono: d=$DELTA, Jp=$J_PRIME, FId=$FID"
            echo "Czas zakończenia: $(date)"
            
        done
    done
done

echo "Zakończono wszystkie obliczenia: $(date)"