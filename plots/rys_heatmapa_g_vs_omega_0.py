from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# PARAMETRY

# Dowolny plik z tej samej serii, numer eig w tej ścieżce zostanie automatycznie podmieniony.
filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g\grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"
)

# "single" -> jeden wybrany eig
# "sum"    -> suma wartości dla zakresu eig
mode = "single"

selected_eig = 3
eig_min = 3
eig_max = 20

# Oś X: "g" lub "g2"
# "g"  -> heatmapa omega_0 vs g
# "g2" -> heatmapa omega_0 vs g^2
x_axis_mode = "g"


# =============================================================================
# LINIE STAŁE

show_constant_lines = True

show_omega_line_1 = True
omega_line_1 = 0.483

show_omega_line_2 = True
omega_line_2 = 1.52

show_omega_line_3 = True
omega_line_3 = 0.621

show_g_line_1 = True
g_line_1 = 0.483

show_g_line_2 = True
g_line_2 = 1.52

show_g_line_3 = True
g_line_3 = 0.759

constant_line_width = 2.5
constant_line_style = "--"


# =============================================================================
# WYŚWIETLANIE / ZAPIS - FOLDER GDZIE ZAPISUJE WYKRERS

show_plot = True
save_plot = True
dpi = 300
cmap = "viridis"

vmin = None
vmax = None

output_folder = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_05.10.2026"
)


# =============================================================================
# PARAMETRY Z NAZWY PLIKU

nazwa_pliku = filename.name
wzorzec = (
    r"grid_"
    r"M(?P<M>\d+)_"
    r"Jp(?P<Jp>[-+]?\d*\.?\d+)_"
    r"d(?P<Delta>[-+]?\d*\.?\d+)_"
    r"d2(?P<Delta2>[-+]?\d*\.?\d+)_"
    r"T(?P<T>[^_]+)_"
    r"P(?P<P>[^_]+)_"
    r"F(?P<F>[^_]+)_"
    r"B(?P<B>[^_]+)_"
    r"FId(?P<FId>[^_]+)_"
    r"eig(?P<eig>\d+)\.txt"
)

match = re.match(wzorzec, nazwa_pliku)

if match is None:
    raise ValueError(
        "Nie udało się odczytać parametrów z nazwy pliku:\n"
        f"{nazwa_pliku}"
    )

parametry = match.groupdict()

M = int(parametry["M"])
Jp = float(parametry["Jp"])
Delta = float(parametry["Delta"])
Delta2 = float(parametry["Delta2"])

T = parametry["T"]
P = parametry["P"]
F = parametry["F"]
B = parametry["B"]
FId = parametry["FId"]

# =============================================================================
# FUNKCJE

def plik_dla_eig(eig):
    nowa_nazwa = re.sub(
        r"_eig\d+\.txt$",
        f"_eig{eig}.txt",
        filename.name
    )
    return filename.parent / nowa_nazwa


def wczytaj_plik(sciezka):
    if not sciezka.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku:\n{sciezka}")

    dane = np.loadtxt(sciezka)

    # KOLEJNOŚĆ KOLUMN 
    omega = dane[:, 0]
    g = dane[:, 1]
    wartosci = dane[:, 2]

    return g, omega, wartosci


# =============================================================================
# WCZYTANIE DANYCH

if mode == "single":
    sciezka = plik_dla_eig(selected_eig)

    print(f"Wczytuję eig = {selected_eig}")
    print(sciezka)

    g, omega, wartosci = wczytaj_plik(sciezka)

    opis_eig = rf"$\lambda_{{{selected_eig}}}$"
    opis_pliku = f"eig{selected_eig}"

elif mode == "sum":
    suma = None
    g_ref = None
    omega_ref = None

    for eig in range(eig_min, eig_max + 1):
        sciezka = plik_dla_eig(eig)

        print(f"Wczytuję eig = {eig}")
        print(sciezka)

        g, omega, wartosci = wczytaj_plik(sciezka)

        if suma is None:
            suma = np.zeros_like(wartosci)
            g_ref = g.copy()
            omega_ref = omega.copy()

        else:
            if not np.allclose(g, g_ref):
                raise ValueError(
                    f"Siatka g w eig={eig} różni się od poprzednich."
                )

            if not np.allclose(omega, omega_ref):
                raise ValueError(
                    f"Siatka omega_0 w eig={eig} różni się od poprzednich."
                )

        suma += wartosci

    g = g_ref
    omega = omega_ref
    wartosci = suma

    opis_eig = rf"$\sum_{{k={eig_min}}}^{{{eig_max}}}\lambda_k$"
    opis_pliku = f"sum_eig{eig_min}-{eig_max}"

else:
    raise ValueError('mode musi mieć wartość "single" albo "sum".')

# =============================================================================
# SIATKA 2D

g_values = np.unique(g)
omega_values = np.unique(omega)

heatmap = np.full((len(omega_values), len(g_values)),np.nan)

for g_i, omega_i, value_i in zip(g, omega, wartosci):

    indeks_g = np.where(np.isclose(g_values, g_i))[0][0]
    indeks_omega = np.where(np.isclose(omega_values, omega_i))[0][0]
    heatmap[indeks_omega, indeks_g] = value_i

# =============================================================================
# OŚ X

if x_axis_mode == "g":
    x_values = g_values
    x_label = r"$g$"
    opis_osi_pliku = "g_x"

elif x_axis_mode == "g2":
    x_values = g_values**2
    x_label = r"$g^2$"
    opis_osi_pliku = "g2_x"

else:
    raise ValueError('x_axis_mode musi mieć wartość "g" albo "g2".')

# =============================================================================
# TYTUŁ

title = (
    rf"$M={M}$, "
    rf"$J'={Jp}$, "
    rf"$\Delta={Delta}$, "
    rf"$\Delta_2={Delta2}$"
    "\n"
    rf"$T={T}$, "
    rf"$P={P}$, "
    rf"$F={F}$, "
    rf"$B={B}$, "
    rf"$FId={FId}$, "
    f"{opis_eig}"
)

# =============================================================================
# RYSOWANIE

fig, ax = plt.subplots(figsize=(8, 6))

mesh = ax.pcolormesh(x_values, omega_values, heatmap, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)

cbar = fig.colorbar(mesh, ax=ax)

# =============================================================================
# LINIE STAŁE

if show_constant_lines:
    # POZIOME: omega_0 = const
    if show_omega_line_1:
        ax.axhline(
            y=omega_line_1,
            color="red",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$\omega_0 = {omega_line_1:.3g}$"
        )

    if show_omega_line_2:
        ax.axhline(
            y=omega_line_2,
            color="orange",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$\omega_0 = {omega_line_2:.3g}$"
        )

    if show_omega_line_3:
        ax.axhline(
            y=omega_line_3,
            color="magenta",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$\omega_0 = {omega_line_3:.3g}$"
        )

    # PIONOWE: g = const
    if show_g_line_1:
        x_g_1 = g_line_1 if x_axis_mode == "g" else g_line_1**2

        ax.axvline(
            x=x_g_1,
            color="cyan",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$g = {g_line_1:.3g}$"
        )

    if show_g_line_2:
        x_g_2 = g_line_2 if x_axis_mode == "g" else g_line_2**2

        ax.axvline(
            x=x_g_2,
            color="white",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$g = {g_line_2:.3g}$"
        )

    if show_g_line_3:
        x_g_3 = g_line_3 if x_axis_mode == "g" else g_line_3**2

        ax.axvline(
            x=x_g_3,
            color="lime",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$g = {g_line_3:.3g}$"
        )

    ax.legend(fontsize=12, loc="best")

# =============================================================================
# COLORBAR

if mode == "single":
    cbar.set_label(rf"$\lambda_{{{selected_eig}}}$", fontsize=14)

else:
    cbar.set_label(rf"$\sum_{{k={eig_min}}}^{{{eig_max}}}\lambda_k$", fontsize=14)

# =============================================================================
# OSIE

ax.set_xlabel(x_label, fontsize=14)
ax.set_ylabel(r"$\omega_0$", fontsize=14)
ax.set_title(title, fontsize=14)
ax.tick_params(axis="both", labelsize=12)
plt.tight_layout()

# =============================================================================
# ZAPIS

output_folder.mkdir(parents=True, exist_ok=True)

opis_linii = (
    "const_lines_yes"
    if show_constant_lines
    else "const_lines_no"
)

if save_plot:
    output_filename = output_folder / (
        f"heatmap_{opis_osi_pliku}"
        f"_M{M}"
        f"_Jp{Jp}"
        f"_d{Delta}"
        f"_d2{Delta2}"
        f"_T{T}"
        f"_P{P}"
        f"_F{F}"
        f"_B{B}"
        f"_FId{FId}"
        f"_{opis_pliku}"
        f"_{opis_linii}.png"
    )

    plt.savefig(output_filename, dpi=dpi, bbox_inches="tight")

    print()
    print("Zapisano wykres:")
    print(output_filename)

# =============================================================================
# WYŚWIETLENIE

if show_plot:
    plt.show()
else:
    plt.close()