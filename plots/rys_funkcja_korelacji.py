
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# USTAWIENIA
# =============================================================================

FOLDER_DANE = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
)

# Plik z korelacjami dla losowego stanu
PLIK = FOLDER_DANE / "Korelacje_ED_vs_g_t_PBCtrue_L_x4_N_upall_J1.0_Delta0.001_omega00.5_gmin1.0_gmax1.0_ng1_tmin0.0_tmax100.0_nt501.data"
# Plik z korelacjami sladowymi - automatycznie ta sama seria
PLIK_TRACE = PLIK.with_name(
    PLIK.name.replace("Korelacje_ED_vs_g_t_", "Korelacje_trace_ED_vs_g_t_", 1)
)

# WYBOR WARTOSCI SPRZEZENIA
G_WYBRANE = 1.0

# WYBOR FUNKCJI KORELACJI
WYBRANE_KORELACJE = [
    "C_slow_3_50",
    "C_slow_8_12_62_118",
    "C_slow_8_4_231_237_62_118_1128_1170",
]

# KOLORY LINII CIAGLYCH - LOSOWY STAN
KOLORY_RANDOM = ["#1f77b4", "#ff7f0e", "#2ca02c"]

# KOLORY LINII PRZERYWANYCH - TRACE
KOLORY_TRACE = ["#87CEFA", "#E0B000", "#90EE90"]

# OPCJE WYKRESU
SKALA_Y = "linear"       # "linear" albo "log"
POKAZ_PUNKTY = False
POKAZ_SIATKE = True

# ZAKRES OSI X
X_MIN = 0.0
X_MAX = None

# ZAKRES OSI Y
Y_MIN = None
Y_MAX = None


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

    dostepne_g = np.unique(dane[:, kolumny.index("g")])

    maska = np.isclose(
        dane[:, kolumny.index("g")],
        g_wybrane,
        rtol=0,
        atol=1e-8
    )

    if not np.any(maska):
        raise ValueError(
            f"Nie znaleziono danych dla g = {g_wybrane}\n"
            f"Dostepne wartosci g: {dostepne_g}"
        )

    dane_g = dane[maska]
    indeksy_sortowania = np.argsort(dane_g[:, kolumny.index("t")])

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

    return nazwa


# =============================================================================
# RYSOWANIE WYKRESU
# =============================================================================

def rysuj_wykres(plik, plik_trace, wybrane_korelacje, g_wybrane):

    parametry = odczytaj_parametry(plik)
    odczytaj_parametry(plik_trace)

    kolumny, dane = wczytaj_dane(plik)
    kolumny_trace, dane_trace = wczytaj_dane(plik_trace)

    if not wybrane_korelacje:
        raise ValueError("Nie wybrano zadnej funkcji korelacji.")

    for nazwa in wybrane_korelacje:
        if nazwa not in kolumny or nazwa not in kolumny_trace:
            raise ValueError(f"Brak korelacji {nazwa} w jednym z plikow.")

    if len(wybrane_korelacje) > len(KOLORY_RANDOM):
        raise ValueError("Za malo kolorow w KOLORY_RANDOM.")

    if len(wybrane_korelacje) > len(KOLORY_TRACE):
        raise ValueError("Za malo kolorow w KOLORY_TRACE.")

    print("Dostepne wartosci g:", np.unique(dane[:, kolumny.index("g")]))

    dane_g = wybierz_dane_g(kolumny, dane, g_wybrane)
    dane_g_trace = wybierz_dane_g(kolumny_trace, dane_trace, g_wybrane)

    t = dane_g[:, kolumny.index("t")]
    t_trace = dane_g_trace[:, kolumny_trace.index("t")]
    g_rzeczywiste = dane_g[0, kolumny.index("g")]

    fig, ax = plt.subplots(figsize=(10, 6))

    # RYSOWANIE OBU METOD DLA KAZDEGO OPERATORA
    for i, nazwa in enumerate(wybrane_korelacje):

        C_t = dane_g[:, kolumny.index(nazwa)]
        C_trace_t = dane_g_trace[:, kolumny_trace.index(nazwa)]
        etykieta = etykieta_korelacji(nazwa)

        # Losowy stan - linia ciagla
        ax.plot(
            t,
            C_t,
            linestyle="-",
            linewidth=1.8,
            color=KOLORY_RANDOM[i],
            marker="o" if POKAZ_PUNKTY else None,
            markersize=3.5,
            label=etykieta + " — random",
        )

        # Dokladny slad - linia przerywana
        ax.plot(
            t_trace,
            C_trace_t,
            linestyle="--",
            linewidth=2.0,
            color=KOLORY_TRACE[i],
            label=etykieta + " — trace",
        )

    # OPISY OSI
    ax.set_xlabel(r"$t$", fontsize=13)
    ax.set_ylabel(r"$C_Q(t)$", fontsize=13)

    # TYTUL
    ax.set_title(
        rf"$M={parametry['M']},\; "
        rf"N_\uparrow={parametry['N_up']},\; "
        rf"J={parametry['J']},\; "
        rf"\Delta={parametry['Delta']},\; "
        rf"\omega_0={parametry['omega0']},\; "
        rf"g={g_rzeczywiste:.4g},\; "
        rf"\mathrm{{PBC}}={parametry['PBC']}$",
        fontsize=12,
    )

    ax.set_yscale(SKALA_Y)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)

    ax.grid(POKAZ_SIATKE, alpha=0.3)
    ax.legend(fontsize=9)

    fig.tight_layout()
    plt.show()


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
