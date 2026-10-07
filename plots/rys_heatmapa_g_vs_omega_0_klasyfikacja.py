from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USTAWIENIA

filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g_liomsy\lioms_grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3_kompozycja.data"
)

# "I_n" -> |#S+ - #S-|
# "P_n" -> |#pozycja(S+) - #pozycja(S-)| tylko dla I_0
# "typ" -> b, f, m

# WYBÓR KLASYFIKACJI
klasyfikacja = "P_n"

# WYBÓR KLASY
wybrana_grupa = 1

# "g" albo "g2"
x_axis_mode = "g"


# =============================================================================
# LINIE STAŁE

show_constant_lines = False

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
# WYKRES / ZAPIS

show_plot = True
save_plot = True
dpi = 300

vmin, vmax = 0.0, 1.0
cmap = "viridis"

output_folder = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_05.10.2026"
)

# =============================================================================
# PARAMETRY Z NAZWY PLIKU

wzorzec = (
    r"lioms_grid_"
    r"M(?P<M>\d+)_"
    r"Jp(?P<Jp>[-+]?\d*\.?\d+)_"
    r"d(?P<Delta>[-+]?\d*\.?\d+)_"
    r"d2(?P<Delta2>[-+]?\d*\.?\d+)_"
    r"T(?P<T>[^_]+)_"
    r"P(?P<P>[^_]+)_"
    r"F(?P<F>[^_]+)_"
    r"B(?P<B>[^_]+)_"
    r"FId(?P<FId>[^_]+)_"
    r"eig(?P<eig>\d+)_kompozycja\.data"
)

match = re.match(wzorzec, filename.name)

if match is None:
    raise ValueError(
        f"Nie udało się odczytać parametrów z nazwy pliku:\n{filename.name}"
    )

p = match.groupdict()

M = int(p["M"])
Jp, Delta, Delta2 = float(p["Jp"]), float(p["Delta"]), float(p["Delta2"])
T, P, F, B, FId = p["T"], p["P"], p["F"], p["B"], p["FId"]
eig = int(p["eig"])

# =============================================================================
# WYBÓR KOLUMNY

with open(filename, "r") as f:
    naglowek = f.readline().strip().split()

if klasyfikacja == "I_n":
    nazwa_kolumny = f"I_n={wybrana_grupa}"

elif klasyfikacja == "P_n":
    nazwa_kolumny = f"P_n={wybrana_grupa}"

    # starsza nazwa kolumny
    if nazwa_kolumny not in naglowek:
        nazwa_kolumny = f"I_0={wybrana_grupa}"

elif klasyfikacja == "typ":
    nazwa_kolumny = f"typ={wybrana_grupa}"

else:
    raise ValueError(
        'klasyfikacja musi mieć wartość "I_n", "P_n" albo "typ".'
    )

if nazwa_kolumny not in naglowek:
    raise ValueError(
        f"Nie znaleziono kolumny {nazwa_kolumny} w pliku."
    )

indeks_kolumny = naglowek.index(nazwa_kolumny)

# =============================================================================
# WCZYTANIE DANYCH

dane = np.loadtxt(filename, skiprows=1)

omega = dane[:, 0]
g = dane[:, 1]
wartosci = dane[:, indeks_kolumny]

# =============================================================================
# SIATKA 2D

g_values = np.unique(g)
omega_values = np.unique(omega)

heatmap = np.full((len(omega_values), len(g_values)), np.nan)

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
# OPIS KLASYFIKACJI

if klasyfikacja == "I_n":
    opis_klasyfikacji = rf"$I_{{{wybrana_grupa}}}$"
    opis_klasyfikacji_pliku = f"I_{wybrana_grupa}"

elif klasyfikacja == "P_n":
    opis_klasyfikacji = rf"$P_{{{wybrana_grupa}}}$"
    opis_klasyfikacji_pliku = f"P_{wybrana_grupa}"

else:
    opis_klasyfikacji = str(wybrana_grupa)
    opis_klasyfikacji_pliku = str(wybrana_grupa)

# =============================================================================
# TYTUŁ

title = (
    rf"$M={M}$, $J'={Jp}$, $\Delta={Delta}$, $\Delta_2={Delta2}$"
    "\n"
    rf"$T={T}$, $P={P}$, $F={F}$, $B={B}$, $FId={FId}$, "
    rf"$\lambda_{{{eig}}}$, class: {opis_klasyfikacji}"
)

# =============================================================================
# RYSOWANIE

fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
ax.set_box_aspect(1)

mesh = ax.pcolormesh(
    x_values, omega_values, heatmap,
    shading="auto", cmap=cmap, vmin=vmin, vmax=vmax
)

cbar = fig.colorbar(mesh, ax=ax)

ticks = np.arange(0.0, 1.01, 0.1)
cbar.set_ticks(ticks)
cbar.set_ticklabels([f"{x:.1f}" for x in ticks])

# =============================================================================
# LINIE STAŁE

if show_constant_lines:

    if show_omega_line_1:
        ax.axhline(
            omega_line_1,
            color="red",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$\omega_0 = {omega_line_1:.3g}$"
        )

    if show_omega_line_2:
        ax.axhline(
            omega_line_2,
            color="orange",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$\omega_0 = {omega_line_2:.3g}$"
        )

    if show_omega_line_3:
        ax.axhline(
            omega_line_3,
            color="magenta",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$\omega_0 = {omega_line_3:.3g}$"
        )

    if show_g_line_1:
        x_g_1 = g_line_1 if x_axis_mode == "g" else g_line_1**2

        ax.axvline(
            x_g_1,
            color="cyan",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$g = {g_line_1:.3g}$"
        )

    if show_g_line_2:
        x_g_2 = g_line_2 if x_axis_mode == "g" else g_line_2**2

        ax.axvline(
            x_g_2,
            color="white",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$g = {g_line_2:.3g}$"
        )

    if show_g_line_3:
        x_g_3 = g_line_3 if x_axis_mode == "g" else g_line_3**2

        ax.axvline(
            x_g_3,
            color="lime",
            linewidth=constant_line_width,
            linestyle=constant_line_style,
            label=rf"$g = {g_line_3:.3g}$"
        )

    ax.legend(fontsize=12, loc="best")


# =============================================================================
# FORMATOWANIE

cbar.set_label(opis_klasyfikacji, fontsize=14)

ax.set_xlabel(x_label, fontsize=14)
ax.set_ylabel(r"$\omega_0$", fontsize=14)
ax.set_title(title, fontsize=14)
ax.tick_params(axis="both", labelsize=12)


# =============================================================================
# ZAPIS

output_folder.mkdir(parents=True, exist_ok=True)

opis_linii = "const_lines_yes" if show_constant_lines else "const_lines_no"

output_filename = output_folder / (
    f"heatmap_{opis_klasyfikacji_pliku}_{opis_osi_pliku}"
    f"_M{M}_Jp{Jp}_d{Delta}_d2{Delta2}"
    f"_T{T}_P{P}_F{F}_B{B}_FId{FId}"
    f"_eig{eig}_{opis_linii}.png"
)

if save_plot:
    plt.savefig(output_filename, dpi=dpi, bbox_inches="tight")
    print(f"\nZapisano wykres:\n{output_filename}")

# =============================================================================
# WYŚWIETLENIE

if show_plot:
    plt.show()
else:
    plt.close()