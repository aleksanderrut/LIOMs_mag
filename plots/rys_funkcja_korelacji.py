import numpy as np
import matplotlib.pyplot as plt
import os
import re

# =========================
# ŚCIEŻKA DO PLIKU
# =========================

plik = r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state\cor_ms_XXZ_M_14_Nup_7_J_1.0_Delta_0.3_PBC_yes_operator_J_1_R_20.txt"

# =========================
# USTAWIENIA
# =========================

show_column_2 = True   # Lanczos
show_column_3 = False  # ED
show_column_4 = True   # ED Trace T -> infinity

# =========================
# PARAMETRY Z NAZWY PLIKU
# =========================

nazwa_pliku = os.path.basename(plik)

wzorzec = (
    r"M_(\d+)_"
    r"Nup_(\d+)_"
    r"J_([-+]?\d*\.?\d+)_"
    r"Delta_([-+]?\d*\.?\d+)_"
    r"PBC_(yes|no)_"
    r"operator_(A_1|J_1)_"
    r"R_(\d+)"
)

match = re.search(wzorzec, nazwa_pliku)

if match is None:
    raise ValueError("Nie udało się odczytać parametrów z nazwy pliku.")

M = int(match.group(1))
N_up = int(match.group(2))
J = float(match.group(3))
Delta = float(match.group(4))
PBC = match.group(5)
operator = match.group(6)
R = int(match.group(7))

print("Odczytane parametry:")
print("M =", M)
print("N_up =", N_up)
print("J =", J)
print("Delta =", Delta)
print("PBC =", PBC)
print("operator =", operator)
print("R =", R)

# =========================
# NAZWY OPERATORÓW
# =========================

nazwy_operatorow = {
    "A_1": r"$A_1 = \sum_j S_j^z S_{j+1}^z$",
    "J_1": r"$J_1 = \frac{iJ}{2}\sum_j \left(S_{j+1}^{+}S_j^{-} - S_j^{+}S_{j+1}^{-}\right)$"
}

nazwa_operatora = nazwy_operatorow[operator]

# =========================
# WCZYTANIE DANYCH
# =========================

dane = np.loadtxt(plik)

t = dane[:, 0]

if show_column_2:
    C_lanczos = dane[:, 1]

if show_column_3:
    C_ED = dane[:, 2]

if show_column_4:
    C_ED_trace = dane[:, 3]

# =========================
# WYKRES
# =========================

plt.figure(figsize=(10, 6))

# Kolumna 4 - ED Trace
if show_column_4:
    plt.plot(
        t,
        C_ED_trace,
        color="orange",
        linewidth=3,
        linestyle="-",
        label=r"ED Trace $T\to\infty$"
    )

# Kolumna 2 - Lanczos
if show_column_2:
    plt.plot(
        t,
        C_lanczos,
        color="blue",
        linewidth=2,
        linestyle="--",
        label=rf"Lanczos, $R={R}$"
    )

# Kolumna 3 - ED
if show_column_3:
    plt.plot(
        t,
        C_ED,
        color="green",
        linewidth=2,
        linestyle="-",
        label="ED"
    )

plt.xlabel("Time $t$")
plt.ylabel("Correlation function $C(t)$")

plt.title(
    "Comparison of correlation functions for "
    + nazwa_operatora
    + "\n"
    + rf"$M={M},\ N_{{up}}={N_up},\ J={J},\ \Delta={Delta},\ PBC={PBC},\ R={R}$"
)

plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()