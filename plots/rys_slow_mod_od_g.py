from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA
# =============================================================================

# Plik z rozpisanymi współczynnikami LIOM-u
filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_6_25.08.2026\siatka_omega_g_liomsy\lioms_grid_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig4.txt"
)


# =============================================================================
# TRYB WYKRESU
# =============================================================================
#
# "fixed_omega"
#     -> ustalone omega_0
#     -> skrypt sam wyszukuje wszystkie odpowiednie wiersze
#     -> oś X = g
#     -> wykres |c_i|^2(g)
#
# "fixed_g"
#     -> ustalone g
#     -> skrypt sam wyszukuje wszystkie odpowiednie wiersze
#     -> oś X = omega_0
#     -> wykres |c_i|^2(omega_0)
#
# =============================================================================

plot_mode = "fixed_omega"


# =============================================================================
# STAŁE WARTOŚCI
# =============================================================================
#
# Używana jest tylko wartość odpowiadająca wybranemu plot_mode.
#
# Przykład:
#
# plot_mode = "fixed_omega"
# omega_fixed = 0.5
#
# albo:
#
# plot_mode = "fixed_g"
# g_fixed = 0.5
#
# =============================================================================

omega_fixed = 0.5063291139240507
g_fixed = 0.5


# =============================================================================
# TOLERANCJA PRZY WYSZUKIWANIU WARTOŚCI
# =============================================================================
#
# Używamy np.isclose zamiast zwykłego ==,
# ponieważ liczby w pliku są zapisane jako Float64.
#
# =============================================================================

value_tolerance = 1e-10


# =============================================================================
# WYBRANE OPERATORY
# =============================================================================

#wybrane_operatory = [3,8,194,7,12,202,451]
#wybrane_operatory = [3,7,8,12,50,62,118]
#wybrane_operatory = [3,7,8,12,50,53,62,113,118,166]
#wybrane_operatory = [1,3,4,8,12,20,34,50,181,289]

wybrane_operatory = [3,4,7,8,12,50,292,298,1121]


# =============================================================================
# ZAPIS / WYŚWIETLANIE
# =============================================================================

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

    (1, r"$\;S^z{i,2}\quad$"),

    (3, r"$\;S^-_{i,2}+S^+_{i,2}\quad$"),

    (4, r"$\;S^z_{i,2}S^z_{i+1,2}\quad$"),

    (7, r"$\;i\left(S^+_{i,2}S^-_{i+1,2}"
        r"-S^-_{i,2}S^+_{i+1,2}\right)\quad$"),

    (8, r"$\;S^+_{i,2}S^-_{i+1,2}"
        r"+S^-_{i,2}S^+_{i+1,2}\quad$"),

    (10, r"$\;S^z_{i+1,1}\left(S^-_{i,1}+S^+_{i,1}\right)\quad$"),

    (12, r"$\;S^-_{i,2}S^-_{i+1,2}"
         r"+S^+_{i,2}S^+_{i+1,2}\quad$"),

    (20, r"$\;S^+_{i,2}S^z_{i+1,2}S^-_{i+2,2}"
         r"+S^-_{i,2}S^z_{i+1,2}S^+_{i+2,2}\quad$"),

    (32, r"$\;S^+_{i,2}S^+_{i+1,2}S^-_{i+2,2}"
         r"+S^-_{i,2}S^-_{i+1,2}S^+_{i+2,2}\quad$"),

    (50, r"$\;S^z_{i,1}S^z_{i,2}\quad$"),

    (53, r"$\;S^z_{i,1}S^z_{i+1,2}\quad$"),

    (62, r"$\;S^z_{i,1}S^z_{i,2}"
         r"\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    (113, r"$\;S^z_{i+1,1}S^z_{i,2}\quad$"),

    (118, r"$\;S^z_{i+1,1}S^z_{i+1,2}"
          r"\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    (166, r"$\;S^z_{i,1}S^z_{i+1,1}"
          r"S^z_{i,2}S^z_{i+1,2}\quad$"),

    (181, r"$\;S^z_{i,1}S^z_{i+1,1}"
          r"S^z_{i+1,2}S^z_{i+2,2}\quad$"),

    (289, r"$\;S^+_{i,1}S^-_{i+1,1}"
          r"+S^-_{i,1}S^+_{i+1,1}\quad$"),

    (292, r"$\;\left(S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}\right)"
      r"\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    (298, r"$\;\left(S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}\right)"
      r"\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    (1121, r"$\;S^+_{i,1}S^z_{i+1,1}S^-_{i+2,1}"
       r"+S^-_{i,1}S^z_{i+1,1}S^+_{i+2,1}\quad$"),


    ###########################################################################
    # M=3, FId=no
    ###########################################################################

    # (72, r"$\;S^z_{i+1,1}\left("
    #      r"S^+_{i,2}S^-_{i+1,2}+S^-_{i,2}S^+_{i+1,2}"
    #      r"\right)\quad$"),

    # (12, r"$\;S^z_{i,1}\left("
    #      r"S^+_{i,2}S^-_{i+1,2}+S^-_{i,2}S^+_{i+1,2}"
    #      r"\right)\quad$"),

    # (244, r"$\;\left("
    #       r"S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}"
    #       r"\right)"
    #       r"\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),

    # (250, r"$\;\left("
    #       r"S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}"
    #       r"\right)"
    #       r"\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),


    ###########################################################################
    # M=4, FId=yes
    ###########################################################################

    # (3, r"$\;S^-_{i,2}+S^+_{i,2}\quad$"),

    # (8, r"$\;S^+_{i,2}S^-_{i+1,2}"
    #     r"+S^-_{i,2}S^+_{i+1,2}\quad$"),

    # (194, r"$\;S^z_{i,1}S^z_{i,2}\quad$"),

    # (7, r"$\;i\left("
    #     r"S^+_{i,2}S^-_{i+1,2}"
    #     r"-S^-_{i,2}S^+_{i+1,2}"
    #     r"\right)\quad$"),

    # (12, r"$\;S^-_{i,2}S^-_{i+1,2}"
    #      r"+S^+_{i,2}S^+_{i+1,2}\quad$"),

    # (206, r"$\;S^z_{i,1}S^z_{i,2}\left("
    #       r"S^-_{i+1,2}+S^+_{i+1,2}"
    #       r"\right)\quad$"),

    # (454, r"$\;S^z_{i+1,1}\left("
    #       r"S^-_{i,2}+S^+_{i,2}"
    #       r"\right)S^z_{i+1,2}\quad$"),

    # (641, r"$\;S^z_{i,1}S^z_{i+1,1}\quad$"),

    # (202, r"$\;S^z_{i,1}\left("
    #       r"S^-_{i+1,2}+S^+_{i+1,2}"
    #       r"\right)\quad$"),

    # (451, r"$\;S^z_{i+1,1}\left("
    #       r"S^-_{i,2}+S^+_{i,2}"
    #       r"\right)\quad$"),

    # (1153, r"$\;S^+_{i,1}S^-_{i+1,1}"
    #        r"+S^-_{i,1}S^+_{i+1,1}\quad$"),

    # (1164, r"$\;\left("
    #        r"S^+_{i,1}S^-_{i+1,1}"
    #        r"+S^-_{i,1}S^+_{i+1,1}"
    #        r"\right)"
    #        r"\left("
    #        r"S^+_{i,2}S^-_{i+1,2}"
    #        r"+S^-_{i,2}S^+_{i+1,2}"
    #        r"\right)\quad$"),


    ##########################################################################
    # M=4, FId=no
    ###########################################################################

    # Zostają tutaj Twoje dotychczasowe opisy M=4 FId=no.
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
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer"
)


# =============================================================================
# SPRAWDZENIE USTAWIEŃ
# =============================================================================

if plot_mode not in [
    "fixed_omega",
    "fixed_g"
]:
    raise ValueError(
        'plot_mode musi być równy "fixed_omega" albo "fixed_g".'
    )


# =============================================================================
# ODCZYT PARAMETRÓW DO TYTUŁU Z NAZWY PLIKU
# =============================================================================

match_title = re.search(
    r"lioms_grid_"
    r"M(?P<M>\d+)_"
    r"Jp(?P<Jp>[^_]+)_"
    r"d(?P<Delta>[^_]+)_"
    r"d2(?P<Delta2>[^_]+)_"
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
Delta2 = match_title.group("Delta2")
T = match_title.group("T")
P = match_title.group("P")
F = match_title.group("F")
B = match_title.group("B")
FId = match_title.group("FId")
eig = match_title.group("eig")


# =============================================================================
# PRZYGOTOWANIE TABLIC
# =============================================================================

omega_values_all = []
g_values_all = []

wspolczynniki_all = {
    indeks: []
    for indeks in wybrane_operatory
}


# =============================================================================
# WCZYTANIE CAŁEGO PLIKU
# =============================================================================
#
# Format pliku:
#
# kolumna 1 -> omega_0
# kolumna 2 -> g
#
# dalej:
#
# 3:0.123;8:0.456;50:0.789
#
# Wczytujemy CAŁY plik.
# Dopiero później wybieramy:
#
# - konkretne omega_0
# albo
# - konkretne g
#
# =============================================================================

with open(
    filename,
    "r",
    encoding="utf-8"
) as f:

    for linia in f:

        linia = linia.strip()

        if not linia:
            continue


        # =====================================================================
        # PODZIAŁ LINII
        # =====================================================================

        kolumny = linia.split()

        if len(kolumny) < 2:
            continue


        # =====================================================================
        # OMEGA_0 = PIERWSZA KOLUMNA
        # g       = DRUGA KOLUMNA
        # =====================================================================

        try:

            omega_0 = float(
                kolumny[0]
            )

            g = float(
                kolumny[1]
            )

        except ValueError:

            # np. nagłówek tekstowy
            continue


        # =====================================================================
        # SZUKANIE WSPÓŁCZYNNIKÓW LIOM-U
        # =====================================================================

        znalezione = re.findall(
            r"(\d+)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
            linia
        )

        wartosci_w_linii = {
            int(indeks): float(wartosc)
            for indeks, wartosc in znalezione
        }


        # =====================================================================
        # ZAPIS OMEGA_0 I g
        # =====================================================================

        omega_values_all.append(
            omega_0
        )

        g_values_all.append(
            g
        )


        # =====================================================================
        # ZAPIS WSPÓŁCZYNNIKÓW
        # =====================================================================

        for indeks in wybrane_operatory:

            wartosc = wartosci_w_linii.get(
                indeks,
                0.0
            )

            wspolczynniki_all[indeks].append(
                wartosc
            )


# =============================================================================
# KONWERSJA NA NUMPY
# =============================================================================

omega_values_all = np.array(
    omega_values_all
)

g_values_all = np.array(
    g_values_all
)


for indeks in wybrane_operatory:

    wspolczynniki_all[indeks] = np.array(
        wspolczynniki_all[indeks]
    )


# =============================================================================
# KONTROLA CZY WCZYTANO DANE
# =============================================================================

if len(g_values_all) == 0:

    raise ValueError(
        "Nie wczytano żadnych danych z pliku."
    )


# =============================================================================
# WYBÓR DANYCH
# =============================================================================

if plot_mode == "fixed_omega":

    # =========================================================================
    # STAŁE OMEGA_0
    # =========================================================================
    #
    # Szukamy wszystkich wierszy posiadających zadane omega_fixed.
    #
    # Następnie:
    #
    # x = g
    #
    # =========================================================================

    mask = np.isclose(
        omega_values_all,
        omega_fixed,
        rtol=0.0,
        atol=value_tolerance
    )


    # =========================================================================
    # JEŻELI NIE ZNALEZIONO OMEGA_0
    # =========================================================================

    if not np.any(mask):

        unique_omega = np.unique(
            omega_values_all
        )

        nearest_index = np.argmin(
            np.abs(
                unique_omega - omega_fixed
            )
        )

        nearest_omega = unique_omega[
            nearest_index
        ]

        difference = abs(
            nearest_omega - omega_fixed
        )

        raise ValueError(
            f"\nNie znaleziono omega_0 = {omega_fixed}.\n"
            f"Najbliższa wartość w pliku to:\n"
            f"omega_0 = {nearest_omega}\n"
            f"różnica = {difference}\n\n"
            f"Jeżeli chcesz jej użyć, ustaw:\n"
            f"omega_fixed = {nearest_omega}"
        )


    # =========================================================================
    # WYBRANE OMEGA I g
    # =========================================================================

    omega_values = omega_values_all[
        mask
    ]

    g_values = g_values_all[
        mask
    ]


    wspolczynniki = {}

    for indeks in wybrane_operatory:

        wspolczynniki[indeks] = wspolczynniki_all[indeks][
            mask
        ]


    # =========================================================================
    # STAŁA WARTOŚĆ
    # =========================================================================

    fixed_value = omega_values[0]


    # =========================================================================
    # KONTROLA
    # =========================================================================

    if not np.allclose(
        omega_values,
        fixed_value,
        rtol=0.0,
        atol=value_tolerance
    ):

        raise ValueError(
            "Wybrane wartości omega_0 nie są jednakowe.\n"
            f"omega_min = {np.min(omega_values)}\n"
            f"omega_max = {np.max(omega_values)}"
        )


    # =========================================================================
    # OŚ X
    # =========================================================================

    x_values = g_values

    x_label = r"$g$"


    # =========================================================================
    # SORTOWANIE PO g
    # =========================================================================

    kolejnosc = np.argsort(
        x_values
    )


elif plot_mode == "fixed_g":

    # =========================================================================
    # STAŁE g
    # =========================================================================
    #
    # Szukamy wszystkich wierszy posiadających zadane g_fixed.
    #
    # Następnie:
    #
    # x = omega_0
    #
    # =========================================================================

    mask = np.isclose(
        g_values_all,
        g_fixed,
        rtol=0.0,
        atol=value_tolerance
    )


    # =========================================================================
    # JEŻELI NIE ZNALEZIONO g
    # =========================================================================

    if not np.any(mask):

        unique_g = np.unique(
            g_values_all
        )

        nearest_index = np.argmin(
            np.abs(
                unique_g - g_fixed
            )
        )

        nearest_g = unique_g[
            nearest_index
        ]

        difference = abs(
            nearest_g - g_fixed
        )

        raise ValueError(
            f"\nNie znaleziono g = {g_fixed}.\n"
            f"Najbliższa wartość w pliku to:\n"
            f"g = {nearest_g}\n"
            f"różnica = {difference}\n\n"
            f"Jeżeli chcesz jej użyć, ustaw:\n"
            f"g_fixed = {nearest_g}"
        )


    # =========================================================================
    # WYBRANE OMEGA I g
    # =========================================================================

    omega_values = omega_values_all[
        mask
    ]

    g_values = g_values_all[
        mask
    ]


    wspolczynniki = {}

    for indeks in wybrane_operatory:

        wspolczynniki[indeks] = wspolczynniki_all[indeks][
            mask
        ]


    # =========================================================================
    # STAŁA WARTOŚĆ
    # =========================================================================

    fixed_value = g_values[0]


    # =========================================================================
    # KONTROLA
    # =========================================================================

    if not np.allclose(
        g_values,
        fixed_value,
        rtol=0.0,
        atol=value_tolerance
    ):

        raise ValueError(
            "Wybrane wartości g nie są jednakowe.\n"
            f"g_min = {np.min(g_values)}\n"
            f"g_max = {np.max(g_values)}"
        )


    # =========================================================================
    # OŚ X
    # =========================================================================

    x_values = omega_values

    x_label = r"$\omega_0$"


    # =========================================================================
    # SORTOWANIE PO OMEGA_0
    # =========================================================================

    kolejnosc = np.argsort(
        x_values
    )


# =============================================================================
# SORTOWANIE WSZYSTKICH DANYCH
# =============================================================================

x_values = x_values[
    kolejnosc
]

omega_values = omega_values[
    kolejnosc
]

g_values = g_values[
    kolejnosc
]


for indeks in wybrane_operatory:

    wspolczynniki[indeks] = wspolczynniki[indeks][
        kolejnosc
    ]


# =============================================================================
# KONTROLA WCZYTANYCH DANYCH
# =============================================================================

print()
print("============================================")
print("WCZYTANE DANE")
print("============================================")


if plot_mode == "fixed_omega":

    print(
        "Tryb = stałe omega_0"
    )

    print(
        f"Żądane omega_0 = {omega_fixed}"
    )

    print(
        f"Znalezione omega_0 = {fixed_value:.16g}"
    )

    print(
        f"Liczba punktów g = {len(x_values)}"
    )

    print(
        f"g_min = {np.min(x_values)}"
    )

    print(
        f"g_max = {np.max(x_values)}"
    )


elif plot_mode == "fixed_g":

    print(
        "Tryb = stałe g"
    )

    print(
        f"Żądane g = {g_fixed}"
    )

    print(
        f"Znalezione g = {fixed_value:.16g}"
    )

    print(
        f"Liczba punktów omega_0 = {len(x_values)}"
    )

    print(
        f"omega_0_min = {np.min(x_values)}"
    )

    print(
        f"omega_0_max = {np.max(x_values)}"
    )


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
# AUTOMATYCZNA NAZWA PLIKU WYJŚCIOWEGO
# =============================================================================
#
# Stała wartość jest zapisywana z trzema cyframi znaczącymi.
#
# =============================================================================

if plot_mode == "fixed_omega":

    plot_name = "liom_i_vs_g"

    fixed_tag = (
        f"_omega0_{fixed_value:.3g}"
    )


elif plot_mode == "fixed_g":

    plot_name = "liom_i_vs_omega0"

    fixed_tag = (
        f"_g_{fixed_value:.3g}"
    )


output_filename = output_directory / (

    f"{plot_name}_"

    f"M{M}_"
    f"Jp{Jp}_"
    f"d{Delta}_"
    f"d2{Delta2}_"
    f"T{T}_"
    f"P{P}_"
    f"F{F}_"
    f"B{B}_"
    f"FId{FId}_"
    f"eig{eig}"

    f"{fixed_tag}.png"
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
    rf"$M = {M};\ "
    rf"J' = {Jp};\ "
    rf"\Delta = {Delta};\ "
    rf"\Delta_2 = {Delta2};\ "
    rf"\mathrm{{T}} = \mathrm{{{T}}};\ "
    rf"\mathrm{{P}} = \mathrm{{{P}}};\ "
    rf"\mathrm{{F}} = \mathrm{{{F}}};\ "
    rf"\mathrm{{B}} = \mathrm{{{B}}};\ "
    rf"\mathrm{{FId}} = \mathrm{{{FId}}};\ "
    rf"\mathrm{{eig}} = {eig}"
)


# =============================================================================
# DOPISANIE STAŁEGO PARAMETRU DO TYTUŁU
# =============================================================================
#
# 3 cyfry znaczące:
#
# .3g
#
# =============================================================================

if plot_mode == "fixed_omega":

    tytul += (
        rf";\ \omega_0 = {fixed_value:.3g}"
    )


elif plot_mode == "fixed_g":

    tytul += (
        rf";\ g = {fixed_value:.3g}"
    )


tytul += "$"


ax.set_title(
    tytul,
    fontsize=plot_title_fontsize,
    pad=15
)


# =============================================================================
# RYSOWANIE
# =============================================================================

for indeks in wybrane_operatory:

    ax.plot(
        x_values,
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
    x_label,
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

ax.legend(
    fontsize=12
)


# =============================================================================
# TWORZENIE RĘCZNEGO OPISU POD WYKRESEM
# =============================================================================

if show_manual_legend:

    aktywne_opisy = []

    for indeks, opis in manual_legend_entries:

        if indeks in wybrane_operatory:

            aktywne_opisy.append(
                rf"$\mathbf{{{indeks}:}}$ "
                + opis
            )


    # =========================================================================
    # PODZIAŁ OPISU NA DWIE LINIE
    # =========================================================================

    podzial = (
        len(aktywne_opisy) + 1
    ) // 2


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
    # RAMKA POD WYKRESEM
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