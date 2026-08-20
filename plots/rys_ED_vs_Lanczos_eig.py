import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# Ścieżka do pliku
sciezka = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\Lanczos_XXZ_wyniki"
    r"\roznica_energii_XXZ_M_8_Nup_4_J_1.0_Delta_0.3_PBC_yes.txt"
)

# Wczytanie danych
dane = np.loadtxt(sciezka)

iteracje = dane[:, 0]
roznica = dane[:, 1]

# Ręczne ustawienie zakresu osi X
x_min = 1
x_max = 100

# Wykres
plt.figure(figsize=(9, 6))

plt.plot(iteracje, roznica, marker="o")

plt.xlabel("Iteracja Lanczosa")
plt.ylabel("ED - Lanczos (najmniejsza wartość własna)")
plt.title(sciezka.name)

plt.xlim(x_min, x_max)

plt.grid(True)
plt.tight_layout()

plt.show()