import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# USTAWIENIA
filename = r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\liom_mag_dane\siatka_omega_g\grid_ladder_mag_M_2_J_1.0_d_0.3_T_both_P_both_Sz_cons_both_omegan_100_gn_100_eig_3.txt"

use_log_axes = False        # logarytmiczna skala osi omega_0 i g
use_log_values = False     # logarytmiczna skala wartości lambda
drop_zero_lambda = False   # jeśli True, usuwa punkty z lambda = 0

save_plot = True
output_filename = r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\plots\heatmap_lambda_M2_eig3_omega_g_log.png"

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

omega = data[:, 0]
g = data[:, 1]
lam = data[:, 2]

# FILTROWANIE DANYCH
mask = np.ones(len(data), dtype=bool)

# jeśli logarytmiczne osie, trzeba wyrzucić omega<=0 i g<=0
if use_log_axes:
    mask &= (omega > 0) & (g > 0)

# opcjonalnie można wyrzucić lambda = 0
if use_log_values and drop_zero_lambda:
    mask &= (lam != 0)

omega = omega[mask]
g = g[mask]
lam = lam[mask]

if len(omega) == 0:
    raise ValueError("Po odfiltrowaniu nie zostały żadne dane do narysowania.")

omega_values = np.unique(omega)
g_values = np.unique(g)

# lepiej wypełniać NaN niż empty
Z = np.full((len(g_values), len(omega_values)), np.nan)

for omega_value, g_value, value in zip(omega, g, lam):
    i = np.where(g_values == g_value)[0][0]
    j = np.where(omega_values == omega_value)[0][0]
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
    omega_values,
    g_values,
    Z_plot,
    shading="auto"
)

# SKALA LOGARYTMICZNA NA OSIACH
if use_log_axes:
    plt.xscale("log")
    plt.yscale("log")

# NAZWY OSI
plt.xlabel(r"$\omega_0$")
plt.ylabel(r"$g$")

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