from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# =============================================================================
# USTAWIENIA
# =============================================================================

FOLDER_DANE = Path(
    #r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_9_27_10_2026"
)

# Plik z korelacjami dla losowego stanu
PLIK = FOLDER_DANE / "Korelacje_ED_vs_g_t_PBCtrue_L_x6_N_upall_J1.0_Delta0.001_omega00.5_gmin0.0_gmax10.0_ng11_tmin0.0_tmax100.0_nt501.data"

# Plik z korelacjami sladowymi
PLIK_TRACE = PLIK.with_name(
    PLIK.name.replace("Korelacje_ED_vs_g_t_", "Korelacje_trace_ED_vs_g_t_", 1)
)

# WYBOR WARTOSCI SPRZEZENIA
G_WYBRANE = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]

# WYBOR FUNKCJI KORELACJI
WYBRANE_KORELACJE = [
    "Q_60",
]

# OPCJE WYKRESU
SKALA_Y = "linear"
POKAZ_PUNKTY = False
POKAZ_SIATKE = True

# ZAKRES OSI X
X_MIN = 0.0
X_MAX = None

# ZAKRES OSI Y
Y_MIN = None
Y_MAX = None

FIGSIZE = (12, 7)


# =============================================================================
# ODCZYTYWANIE PARAMETROW Z NAZWY PLIKU
# =============================================================================

def odczytaj_parametry(plik):

    wzor = re.compile(
        r"^Korelacje_(?:trace_)?ED_vs_g_t_"
        r"PBC(?P<PBC>true|false)"
        r"_L_x(?P<M>\d+)"
        r"_N_up(?P<N_up>\d+|all)"
        r"_J(?P<J>[-+0-9.eE]+)"
        r"_Delta(?P<Delta>[-+0-9.eE]+)"
        r"_omega0(?P<omega0>[-+0-9.eE]+)"
        r"_gmin(?P<gmin>[-+0-9.eE]+)"
        r"_gmax(?P<gmax>[-+0-9.eE]+)"
        r"_ng(?P<ng>\d+)"
        r"_tmin(?P<tmin>[-+0-9.eE]+)"
        r"_tmax(?P<tmax>[-+0-9.eE]+)"
        r"_nt(?P<nt>\d+)"
        r"\.data$"
    )

    wynik = wzor.fullmatch(plik.name)

    if wynik is None:
        raise ValueError(
            f"Nie udalo sie odczytac parametrow z nazwy:\n{plik.name}"
        )

    return wynik.groupdict()


# =============================================================================
# WCZYTYWANIE DANYCH
# =============================================================================

def wczytaj_dane(plik):

    if not plik.is_file():
        raise FileNotFoundError(f"Nie znaleziono pliku: {plik}")

    with plik.open("r", encoding="utf-8-sig") as f:
        kolumny = f.readline().strip().split()

    if len(kolumny) < 3 or kolumny[:2] != ["g", "t"]:
        raise ValueError(
            "Pierwsze dwie kolumny musza miec nazwy 'g' i 't'."
        )

    dane = np.loadtxt(plik, skiprows=1, ndmin=2)

    if dane.shape[1] != len(kolumny):
        raise ValueError(
            "Liczba kolumn naglowka nie zgadza sie z danymi."
        )

    return kolumny, dane


# =============================================================================
# WYBOR WARTOSCI g
# =============================================================================

def wybierz_dane_g(kolumny, dane, g_wybrane):

    maska = np.isclose(
        dane[:, kolumny.index("g")],
        g_wybrane,
        rtol=0,
        atol=1e-8
    )

    if not np.any(maska):
        raise ValueError(
            f"Nie znaleziono danych dla g = {g_wybrane}\n"
            f"Dostepne wartosci g: {np.unique(dane[:, 0])}"
        )

    dane_g = dane[maska]

    indeksy_sortowania = np.argsort(
        dane_g[:, kolumny.index("t")]
    )

    return dane_g[indeksy_sortowania]


# =============================================================================
# ETYKIETY FUNKCJI KORELACJI
# =============================================================================

def etykieta_korelacji(nazwa):

    if nazwa.startswith("C_slow_"):
        indeksy = ",".join(
            nazwa.removeprefix("C_slow_").split("_")
        )
        return rf"$C_{{\mathrm{{slow}},\,{indeksy}}}(t)$"

    if nazwa.startswith("Q_"):
        indeks = nazwa.removeprefix("Q_")
        return rf"$C_{{Q_{{{indeks}}}}}(t)$"

    return nazwa


# =============================================================================
# RYSOWANIE WSZYSTKICH g NA JEDNYM WYKRESIE
# =============================================================================

def rysuj_wykres(plik, plik_trace, wybrane_korelacje, g_wybrane):

    parametry = odczytaj_parametry(plik)
    odczytaj_parametry(plik_trace)

    kolumny, dane = wczytaj_dane(plik)
    kolumny_trace, dane_trace = wczytaj_dane(plik_trace)

    if not wybrane_korelacje:
        raise ValueError("Nie wybrano zadnej korelacji.")

    if not g_wybrane:
        raise ValueError("Nie wybrano zadnej wartosci g.")

    for nazwa in wybrane_korelacje:
        if nazwa not in kolumny or nazwa not in kolumny_trace:
            raise ValueError(f"Brak korelacji {nazwa} w jednym z plikow.")

    print("Dostepne wartosci g:", np.unique(dane[:, 0]))
    print("Wybrane wartosci g:", g_wybrane)

    # -------------------------------------------------------------------------
    # KOLORY - INNY KOLOR DLA KAZDEJ KOMBINACJI (g, OPERATOR)
    # -------------------------------------------------------------------------

    liczba_krzywych = len(g_wybrane) * len(wybrane_korelacje)

    if liczba_krzywych <= 10:
        kolory = plt.cm.tab10(np.arange(liczba_krzywych))
    else:
        kolory = plt.cm.turbo(
            np.linspace(0.05, 0.95, liczba_krzywych)
        )

    # -------------------------------------------------------------------------
    # WYKRES
    # -------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=FIGSIZE)

    numer_koloru = 0

    for g in g_wybrane:

        dane_g = wybierz_dane_g(kolumny, dane, g)
        dane_g_trace = wybierz_dane_g(kolumny_trace, dane_trace, g)

        t = dane_g[:, kolumny.index("t")]
        t_trace = dane_g_trace[:, kolumny_trace.index("t")]

        g_rzeczywiste = dane_g[0, kolumny.index("g")]

        for nazwa in wybrane_korelacje:

            C_t = dane_g[:, kolumny.index(nazwa)]
            C_trace_t = dane_g_trace[:, kolumny_trace.index(nazwa)]

            kolor = kolory[numer_koloru]
            numer_koloru += 1

            etykieta = etykieta_korelacji(nazwa)

            # RANDOM - LINIA CIAGLA
            ax.plot(
                t,
                C_t,
                color=kolor,
                linestyle="-",
                linewidth=1.8,
                marker="o" if POKAZ_PUNKTY else None,
                markersize=3.5,
                label=rf"{etykieta}, $g={g_rzeczywiste:g}$"
            )

            # TRACE - TEN SAM KOLOR, LINIA PRZERYWANA
            ax.plot(
                t_trace,
                C_trace_t,
                color=kolor,
                linestyle="--",
                linewidth=1.8,
                label="_nolegend_"
            )

    # -------------------------------------------------------------------------
    # OPISY OSI
    # -------------------------------------------------------------------------

    ax.set_xlabel(r"$t$", fontsize=13)
    ax.set_ylabel(r"$C_Q(t)$", fontsize=13)

    ax.set_title(
        rf"$M={parametry['M']},\; "
        rf"N_\uparrow={parametry['N_up']},\; "
        rf"J={parametry['J']},\; "
        rf"\Delta={parametry['Delta']},\; "
        rf"\omega_0={parametry['omega0']},\; "
        rf"\mathrm{{PBC}}={parametry['PBC']}$",
        fontsize=12
    )

    ax.set_yscale(SKALA_Y)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)

    ax.grid(POKAZ_SIATKE, alpha=0.3)

    # -------------------------------------------------------------------------
    # LEGENDA
    # -------------------------------------------------------------------------

    handles, labels = ax.get_legend_handles_labels()

    # Dwa dodatkowe przyklady stylu linii
    handles += [
        Line2D([0], [0], color="black", linestyle="-", linewidth=2),
        Line2D([0], [0], color="black", linestyle="--", linewidth=2)
    ]

    labels += [
        "Random state",
        "Trace"
    ]

    ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=9,
        frameon=True,
        ncol=1
    )

    fig.tight_layout()
    plt.show()
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    rysuj_wykres(
        PLIK,
        PLIK_TRACE,
        WYBRANE_KORELACJE,
        G_WYBRANE
    )
