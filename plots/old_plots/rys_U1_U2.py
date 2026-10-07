import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# USTAWIENIA
filename = r"C:\Users\aleks\Desktop\praca magisterska\kod_JP\InfiniteLIOMs\scripts\liom_mag_dane\siatka_U1_U2\grid_ladder_mag_M_3_J_1.0_d_0.3_T_odd_P_even_Sz_cons_yes_U1n_100_U2n_100_eig_1.txt"

use_log_axes = True        # logarytmiczna skala osi U1 i U2
use_log_values = False      # logarytmiczna skala wartości lambda
drop_zero_lambda = False   # jeśli True, usuwa punkty z lambda = 0

save_plot = True
output_filename = r"C:\Users\aleks\Desktop\praca magisterska\kod_JP\InfiniteLIOMs\plots\heatmap_lambda_M3_eig1_log.png"

eps = 1e-30

# SPRAWDZENIE PLIKU
filename_path = Path(filename)
if not filename_path.exists():
    raise FileNotFoundError(f"Nie znaleziono pliku:\n{filename_path}")

print("Wczytuję plik:")
print(filename_path)

plot_title = filename_path.name

# WCZYTANIE DANYCH
data = np.loadtxt(filename_path)

U1 = data[:, 0]
U2 = data[:, 1]
lam = data[:, 2]

# FILTROWANIE DANYCH
mask = np.ones(len(data), dtype=bool)

# jeśli logarytmiczne osie, trzeba wyrzucić U1<=0 i U2<=0
if use_log_axes:
    mask &= (U1 > 0) & (U2 > 0)

# opcjonalnie można wyrzucić lambda = 0
if use_log_values and drop_zero_lambda:
    mask &= (lam != 0)

U1 = U1[mask]
U2 = U2[mask]
lam = lam[mask]

if len(U1) == 0:
    raise ValueError("Po odfiltrowaniu nie zostały żadne dane do narysowania.")

U1_values = np.unique(U1)
U2_values = np.unique(U2)

# lepiej wypełniać NaN niż empty
Z = np.full((len(U1_values), len(U2_values)), np.nan)

for u1, u2, value in zip(U1, U2, lam):
    i = np.where(U1_values == u1)[0][0]
    j = np.where(U2_values == u2)[0][0]
    Z[i, j] = value

# SKALA KOLORÓW / WARTOŚCI
if use_log_values:
    Z_plot = np.where(np.abs(Z) > 0, np.log10(np.abs(Z)), np.nan)
    colorbar_label = r"$\log_{10}(|\lambda|)$"
else:
    Z_plot = Z
    colorbar_label = r"$\lambda$"

# WYKRES
plt.figure(figsize=(7, 6))

mesh = plt.pcolormesh(
    U2_values,
    U1_values,
    Z_plot,
    shading="auto"
)

# SKALA LOGARYTMICZNA NA OSIACH
if use_log_axes:
    plt.xscale("log")
    plt.yscale("log")

# NAZWY OSI
plt.xlabel(r"$U_2$")
plt.ylabel(r"$U_1$")

# TYTUŁ
plt.title(plot_title, fontsize=8)

cbar = plt.colorbar(mesh)
cbar.set_label(colorbar_label)

plt.tight_layout()

# ZAPIS I WYŚWIETLENIE
if save_plot:
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    print("Zapisano obrazek do:")
    print(output_path)

plt.show()