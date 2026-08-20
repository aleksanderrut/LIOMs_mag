import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# PARAMETR
# ==========================================
ile_najnizszych = 2


plik = r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state\wartosci_wlasne_wszystkie_k_XXZ_M_14_Nup_7_J_1.0_Delta_0.3.txt"

dane = np.loadtxt(plik)

k_fiz = dane[:, 0]
eigenvalues = dane[:, 1]

# Unikalne sektory momentum
unikalne_k = np.unique(k_fiz)

plt.figure(figsize=(8, 6))

# Kolejne wartości własne rysujemy osobno,
# dzięki czemu każda dostaje inny kolor
for n in range(ile_najnizszych):

    k_do_wykresu = []
    energie_do_wykresu = []

    for k in unikalne_k:

        energie_k = eigenvalues[np.isclose(k_fiz, k)]

        # Sortowanie energii w danym sektorze k
        energie_k = np.sort(energie_k)

        # Sprawdzenie, czy sektor ma co najmniej n+1 wartości własnych
        if len(energie_k) > n:
            k_do_wykresu.append(k)
            energie_do_wykresu.append(energie_k[n])

    plt.scatter(k_do_wykresu, energie_do_wykresu, label=f"{n + 1}. eigenvalue")


# Oś momentum
tick_positions = [0, np.pi / 3, 2 * np.pi / 3, np.pi, 4 * np.pi / 3, 5 * np.pi / 3, 2 * np.pi]

tick_labels = [r"$0$", r"$\frac{\pi}{3}$", r"$\frac{2\pi}{3}$", r"$\pi$", r"$\frac{4\pi}{3}$", r"$\frac{5\pi}{3}$", r"$2\pi$"]

plt.xticks(tick_positions, tick_labels)

plt.xlim(0, 2 * np.pi)
plt.ylim(-5, -3)

plt.xlabel(r"$k$")
plt.ylabel(r"$E$")

plt.title(f"XXZ Hamiltonian Spectrum — {ile_najnizszych} Lowest Eigenvalues in Each Momentum Sector")

# Tylko pionowe linie siatki
plt.grid(axis="x")

# Legenda kolorów
plt.legend()

plt.tight_layout()
plt.show()