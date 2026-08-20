from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA
# =============================================================================

# Plik z rozpisanymi współczynnikami LIOM-u
filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\se_2_13.08.2026\siatka_omega_g_liomsy\lioms_grid_M3_Jp0.0_d1.5_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"
)

# Numery elementów bazy, które chcemy narysować
#wybrane_operatory = [3,8, 194, 7, 12, 202, 451]
#wybrane_operatory = [3,7,8,12,50,62,118]
wybrane_operatory = [3,7,8,12,50,53,62,113,118,166]

show_plot = True
save_plot = True

# =============================================================================
# RĘCZNY OPIS POD WYKRESEM
# =============================================================================

show_manual_legend = True
manual_legend_title = ""
manual_legend_fontsize = 14

manual_legend_entries = [

    ###########################################################################
    # M=3, FId=yes
    ###########################################################################

    (3, r"$\;S^-_{i,2}+S^+_{i,2}\quad$"),

    (4, r"$\;S^z_{i,2}S^z_{i+1,2}\quad$"),

    (7, r"$\;i\left(S^+_{i,2}S^-_{i+1,2}"r"-S^-_{i,2}S^+_{i+1,2}\right)\quad$"),

    (8, r"$\;S^+_{i,2}S^-_{i+1,2}"r"+S^-_{i,2}S^+_{i+1,2}\quad$"),

    (12, r"$\;S^-_{i,2}S^-_{i+1,2}"r"+S^+_{i,2}S^+_{i+1,2}\quad$"),

    (50, r"$\;S^z_{i,1}S^z_{i,2}\quad$"),

    (53, r"$\;S^z_{i,1}S^z_{i+1,2}\quad$"),

    (62, r"$\;S^z_{i,1}S^z_{i,2}\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    (113, r"$\;S^z_{i+1,1}S^z_{i,2}\quad$"),

    (118, r"$\;S^z_{i+1,1}S^z_{i+1,2}\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    (166, r"$\;S^z_{i,1}S^z_{i+1,1}S^z_{i,2}S^z_{i+1,2}\quad$"),

    ###########################################################################
    # M=3, FId=no
    ###########################################################################
    
    # (72, r"$\;S^z_{i+1,1}\left(S^+_{i,2}S^-_{i+1,2}+S^-_{i,2}S^+_{i+1,2}\right)\quad$"),

    # (12, r"$\;S^z_{i,1}\left(S^+_{i,2}S^-_{i+1,2}+S^-_{i,2}S^+_{i+1,2}\right)\quad$"),

    # (244, r"$\;\left(S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}\right)"r"\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    # (250, r"$\;\left(S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}\right)"r"\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    #(2, r"$\;S^z_{i,1}S^z_{i,2}\quad$"),

    #(5, r"$\;S^z_{i,1}S^z_{i+1,2}\quad$"),

    # (10, r"$\;S^z_{i,1}\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    # (67, r"$\;S^z_{i+1,1}\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    #(65, r"$\;S^z_{i+1,1}S^z_{i,2}\quad$"),

    #(116, r"$\;S^z_{i,1}S^z_{i+1,1}\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    #(122, r"$\;S^z_{i,1}S^z_{i+1,1}\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    ###########################################################################
    # M=4, FId=yes
    ###########################################################################

    # (3, r"$\;S^-_{i,2}+S^+_{i,2}\quad$"),

    # (8, r"$\;S^+_{i,2}S^-_{i+1,2}"
    # r"+S^-_{i,2}S^+_{i+1,2}\quad$"),

    # (194, r"$\;S^z_{i,1}S^z_{i,2}\quad$"),

    # (7, r"$\;i\left("
    # r"S^+_{i,2}S^-_{i+1,2}"
    # r"-S^-_{i,2}S^+_{i+1,2}"
    # r"\right)\quad$"),

    # (12, r"$\;S^-_{i,2}S^-_{i+1,2}"
    #  r"+S^+_{i,2}S^+_{i+1,2}\quad$"),

    # (206, r"$\;S^z_{i,1}S^z_{i,2}\left("
    #   r"S^-_{i+1,2}+S^+_{i+1,2}"
    #   r"\right)\quad$"),

    # (454, r"$\;S^z_{i+1,1}\left("
    #   r"S^-_{i,2}+S^+_{i,2}"
    #   r"\right)S^z_{i+1,2}\quad$"),

    # (641, r"$\;S^z_{i,1}S^z_{i+1,1}\quad$"),

    # (202, r"$\;S^z_{i,1}\left("
    #   r"S^-_{i+1,2}+S^+_{i+1,2}"
    #   r"\right)\quad$"),

    # (451, r"$\;S^z_{i+1,1}\left("
    #   r"S^-_{i,2}+S^+_{i,2}"
    #   r"\right)\quad$"),

    # (1153, r"$\;S^+_{i,1}S^-_{i+1,1}"
    #    r"+S^-_{i,1}S^+_{i+1,1}\quad$"),

    # (1164, r"$\;\left("
    #    r"S^+_{i,1}S^-_{i+1,1}"
    #    r"+S^-_{i,1}S^+_{i+1,1}"
    #    r"\right)"
    #    r"\left("
    #    r"S^+_{i,2}S^-_{i+1,2}"
    #    r"+S^-_{i,2}S^+_{i+1,2}"
    #    r"\right)\quad$"),

    ##########################################################################
    # M=4, FId=no
    ###########################################################################

    # (12353, r"$\;S^-_{i,1}S^+_{i+1,1}S^+_{i+2,1}S^-_{i+3,1}"
    #     r"+S^+_{i,1}S^-_{i+1,1}S^-_{i+2,1}S^+_{i+3,1}\quad$"),

    # (8513, r"$\;S^z_{i,1}S^+_{i+1,1}S^-_{i+2,1}S^z_{i+3,1}"
    #    r"+S^z_{i,1}S^-_{i+1,1}S^+_{i+2,1}S^z_{i+3,1}\quad$"),

    # (2881, r"$\;S^+_{i,1}S^-_{i+2,1}"
    #     r"+S^-_{i,1}S^+_{i+2,1}\quad$"),

    # (449, r"$\;S^z_{i,1}S^z_{i+1,1}\quad$"),

    # (5, r"$\;S^z_{i,1}S^z_{i+1,2}\quad$"),

    # (257, r"$\;S^z_{i+1,1}S^z_{i,2}\quad$"),

    # (4033, r"$\;i\left("
    #    r"S^+_{i,1}S^z_{i+1,1}S^-_{i+2,1}"
    #    r"-S^-_{i,1}S^z_{i+1,1}S^+_{i+2,1}"
    #    r"\right)\quad$"),

    # (2113, r"$\;i\left("
    #     r"S^+_{i,1}S^-_{i+1,1}S^z_{i+2,1}"
    #     r"-S^-_{i,1}S^+_{i+1,1}S^z_{i+2,1}"
    #     r"\right)\quad$"),

    # (3521, r"$\;i\left("
    #     r"S^z_{i,1}S^+_{i+1,1}S^-_{i+2,1}"
    #     r"-S^z_{i,1}S^-_{i+1,1}S^+_{i+2,1}"
    #     r"\right)\quad$"),


    # (12, r"$\;iS^z_{i,1}\left("
    #     r"S^+_{i,2}S^-_{i+1,2}"
    #     r"-S^-_{i,2}S^+_{i+1,2}"
    #     r"\right)\quad$"),

    # (264, r"$\;S^z_{i+1,1}\left("
    #     r"S^+_{i,2}S^-_{i+1,2}"
    #     r"+S^-_{i,2}S^+_{i+1,2}"
    #     r"\right)\quad$"),

    # (10, r"$\;S^z_{i,1}\left("
    #     r"S^-_{i+1,2}+S^+_{i+1,2}"
    #     r"\right)\quad$"),

    # (259, r"$\;S^z_{i+1,1}\left("
    #     r"S^-_{i,2}+S^+_{i,2}"
    #     r"\right)\quad$"),

    # (452, r"$\;S^z_{i,1}S^z_{i+1,1}\left("
    #     r"S^-_{i,2}+S^+_{i,2}"
    #     r"\right)\quad$"),

    # (458, r"$\;S^z_{i,1}S^z_{i+1,1}\left("
    #     r"S^-_{i+1,2}+S^+_{i+1,2}"
    #     r"\right)\quad$"),

    # (15681, r"$\;S^+_{i,1}S^z_{i+1,1}S^z_{i+2,1}S^-_{i+3,1}"
    #     r"+S^-_{i,1}S^z_{i+1,1}S^z_{i+2,1}S^+_{i+3,1}\quad$"),

    # (962, r"$\;S^z_{i,2}\left("
    #     r"S^+_{i,1}S^-_{i+1,1}"
    #     r"+S^-_{i,1}S^+_{i+1,1}"
    #     r"\right)\quad$"),

    # (965, r"$\;S^z_{i+1,2}\left("
    #     r"S^+_{i,1}S^-_{i+1,1}"
    #     r"+S^-_{i,1}S^+_{i+1,1}"
    #     r"\right)\quad$"),

    # (970, r"$\;S^+_{i,1}S^-_{i+1,1}S^-_{i+1,2}"
    #      r"+S^-_{i,1}S^+_{i+1,1}S^+_{i+1,2}\quad$"),

    # (9025, r"$\;S^+_{i,1}S^z_{i+1,1}S^-_{i+2,1}S^z_{i+3,1}"
    #    r"+S^-_{i,1}S^z_{i+1,1}S^+_{i+2,1}S^z_{i+3,1}\quad$"),

]


# =============================================================================
# USTAWIENIA WYKRESU
# =============================================================================

x_min = 0.0
x_max = 2.0

y_min = 0.0
y_max = 1.0


plot_title_fontsize = 16
dpi = 300

output_directory = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\M=3"
)


# =============================================================================
# ODCZYT PARAMETRÓW DO TYTUŁU Z NAZWY PLIKU
# =============================================================================

match_title = re.search(
    r"lioms_grid_"
    r"M(?P<M>\d+)_"
    r"Jp(?P<Jp>[^_]+)_"
    r"d(?P<Delta>[^_]+)_"
    r"T(?P<T>[^_]+)_"
    r"P(?P<P>[^_]+)_"
    r"F(?P<F>[^_]+)_"
    r"B(?P<B>[^_]+)_"
    r"FId(?P<FId>[^_]+)_"
    r"eig(?P<eig>\d+)",
    filename.stem
)

if match_title is None:
    raise ValueError(
        f"Nie udało się odczytać parametrów z nazwy pliku:\n"
        f"{filename.name}"
    )

M = match_title.group("M")
Jp = match_title.group("Jp")
Delta = match_title.group("Delta")
T = match_title.group("T")
P = match_title.group("P")
F = match_title.group("F")
B = match_title.group("B")
FId = match_title.group("FId")
eig = match_title.group("eig")


# =============================================================================
# AUTOMATYCZNA NAZWA PLIKU WYJŚCIOWEGO
# =============================================================================

output_filename = output_directory / (
    f"liom_i_"
    f"M{M}_"
    f"Jp{Jp}_"
    f"d{Delta}_"
    f"T{T}_"
    f"P{P}_"
    f"F{F}_"
    f"B{B}_"
    f"FId{FId}_"
    f"eig{eig}.png"
)


# =============================================================================
# PRZYGOTOWANIE TABLIC
# =============================================================================

g_values = []

wspolczynniki = {
    indeks: []
    for indeks in wybrane_operatory
}


# =============================================================================
# WCZYTANIE DANYCH
# =============================================================================
#
# Zakładamy:
#
# pierwsza kolumna = omega
# druga kolumna    = g
#
# dalej znajdują się wpisy:
#
# 3:0.123
# 8:0.456
# 50:0.789
#
# Jeżeli danego operatora nie ma w danej linii,
# jego współczynnik przyjmujemy jako 0.
#

with open(filename, "r", encoding="utf-8") as f:

    for linia in f:

        linia = linia.strip()

        # Pomijamy puste linie
        if not linia:
            continue

        # ---------------------------------------------------------------------
        # PODZIAŁ LINII NA KOLUMNY
        # ---------------------------------------------------------------------

        kolumny = linia.split()

        if len(kolumny) < 2:
            continue

        # ---------------------------------------------------------------------
        # g = DRUGA KOLUMNA
        # ---------------------------------------------------------------------

        try:
            g = float(kolumny[1])

        except ValueError:
            # np. nagłówek tekstowy
            continue

        # ---------------------------------------------------------------------
        # SZUKANIE WSPÓŁCZYNNIKÓW
        # ---------------------------------------------------------------------

        znalezione = re.findall(
            r"(\d+)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
            linia
        )

        wartosci_w_linii = {
            int(indeks): float(wartosc)
            for indeks, wartosc in znalezione
        }

        # ---------------------------------------------------------------------
        # ZAPIS g
        # ---------------------------------------------------------------------

        g_values.append(g)

        # ---------------------------------------------------------------------
        # ZAPIS WSPÓŁCZYNNIKÓW
        # ---------------------------------------------------------------------

        for indeks in wybrane_operatory:

            wartosc = wartosci_w_linii.get(
                indeks,
                0.0
            )

            wspolczynniki[indeks].append(
                wartosc
            )


# =============================================================================
# KONWERSJA NA NUMPY
# =============================================================================

g_values = np.array(g_values)


# =============================================================================
# SORTOWANIE PO g
# =============================================================================

kolejnosc = np.argsort(g_values)

g_values = g_values[kolejnosc]

for indeks in wybrane_operatory:

    wspolczynniki[indeks] = np.array(
        wspolczynniki[indeks]
    )

    wspolczynniki[indeks] = wspolczynniki[indeks][kolejnosc]


# =============================================================================
# KONTROLA WCZYTANYCH DANYCH
# =============================================================================

print()
print("============================================")
print("WCZYTANE DANE")
print("============================================")

print(f"Liczba punktów g: {len(g_values)}")

if len(g_values) > 0:

    print(f"g_min = {np.min(g_values)}")
    print(f"g_max = {np.max(g_values)}")

print()

for indeks in wybrane_operatory:

    liczba_niezerowych = np.count_nonzero(
        wspolczynniki[indeks]
    )

    print(
        f"Operator {indeks}: "
        f"{liczba_niezerowych} niezerowych wartości"
    )


# =============================================================================
# WYKRES
# =============================================================================

fig, ax = plt.subplots(
    figsize=(10, 6)
)


# =============================================================================
# TYTUŁ WYKRESU
# =============================================================================

tytul = (
    rf"$M = {M};\ J' = {Jp};\ \Delta = {Delta};\ "
    rf"\mathrm{{T}} = \mathrm{{{T}}};\ "
    rf"\mathrm{{P}} = \mathrm{{{P}}};\ "
    rf"\mathrm{{F}} = \mathrm{{{F}}};\ "
    rf"\mathrm{{B}} = \mathrm{{{B}}};\ "
    rf"\mathrm{{FId}} = \mathrm{{{FId}}};\ "
    rf"\mathrm{{eig}} = {eig}$"
)

ax.set_title(
    tytul,
    fontsize=plot_title_fontsize,
    pad=15
)


# =============================================================================
# RYSOWANIE
# =============================================================================
#
# Duże punkty połączone cienką linią.
# Linia i punkty mają automatycznie ten sam kolor.
#

for indeks in wybrane_operatory:

    ax.plot(
        g_values,
        wspolczynniki[indeks],

        marker="o",
        markersize=2,

        linewidth=1.0,

        label=str(indeks)
    )


# =============================================================================
# OSIE
# =============================================================================

ax.set_xlim(
    x_min,
    x_max
)

ax.set_ylim(
    y_min,
    y_max
)

ax.set_xlabel(
    r"$g$",
    fontsize=14
)

ax.set_ylabel(
    r"$|c_i|^2$",
    fontsize=14
)

ax.tick_params(
    axis="both",
    labelsize=12
)


# =============================================================================
# SIATKA
# =============================================================================

ax.grid(
    True,
    alpha=0.3
)


# =============================================================================
# LEGENDA
# =============================================================================
#
# Tylko numery operatorów.
# Bez napisu "Element bazy".
#

ax.legend(
    fontsize=12
)


# =============================================================================
# TWORZENIE RĘCZNEGO OPISU POD WYKRESEM
# =============================================================================

if show_manual_legend:

    aktywne_opisy = []

    for indeks, opis in manual_legend_entries:

        # Pokazujemy tylko operatory aktualnie rysowane
        if indeks in wybrane_operatory:

            aktywne_opisy.append(
                rf"$\mathbf{{{indeks}:}}$ "
                + opis
            )

    # =========================================================================
    # PODZIAŁ OPISU NA DWIE LINIE
    # =========================================================================

    podzial = (len(aktywne_opisy) + 1) // 2

    linia_1 = "     ".join(
        aktywne_opisy[:podzial]
    )

    linia_2 = "     ".join(
        aktywne_opisy[podzial:]
    )

    tekst_operatorow = (
        linia_1
        + "\n"
        + linia_2
    )


    # =========================================================================
    # PODŁUŻNA RAMKA POD WYKRESEM
    # =========================================================================

    fig.text(
        0.5,
        0.025,
        tekst_operatorow,

        fontsize=manual_legend_fontsize,

        horizontalalignment="center",
        verticalalignment="bottom",
        multialignment="center",

        bbox=dict(
            boxstyle="round,pad=0.8",
            facecolor="white",
            edgecolor="black",
            linewidth=1.5,
            alpha=1.0
        )
    )


# =============================================================================
# MIEJSCE NA RAMKĘ POD WYKRESEM
# =============================================================================

plt.subplots_adjust(
    bottom=0.25
)


# =============================================================================
# ZAPIS WYKRESU
# =============================================================================

if save_plot:

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    plt.savefig(
        output_filename,
        dpi=dpi,
        bbox_inches="tight"
    )

    print()
    print("Zapisano wykres:")
    print(output_filename)


# =============================================================================
# WYŚWIETLENIE
# =============================================================================

if show_plot:
    plt.show()

else:
    plt.close()