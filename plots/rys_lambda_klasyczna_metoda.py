
from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USTAWIENIA
# =============================================================================

# Folder z wynikami
FOLDER_DANE = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
)

# Wybór pliku
PLIK = FOLDER_DANE / "Tr_komutatorow_vs_g_PBCtrue_L_x6_N_upall_J1.0_Delta0.001_omega00.5_gmin0.0_gmax10.0_ng41.data"
# WYBOR OPERATOROW DO WYSWIETLENIA
WYBRANE_OPERATORY = [
    "Q_slow_3_50",
    "Q_slow_8_12_62_118",
    "Q_slow_8_4_231_237_62_118_1128_1170",
]

# WYBRANE_OPERATORY = [
#     "Q_3", "Q_4", "Q_8", "Q_12", "Q_50",
#     "Q_62", "Q_118", "Q_231", "Q_237", "Q_1128"
# ]

# OPCJE WYKRESU
SKALA_Y = "linear"  # "linear" albo "log"
POKAZ_PUNKTY = False
POKAZ_SIATKE = True

# ZAKRES OSI X
X_MIN = 0.0
X_MAX = 10.0

# ZAKRES OSI Y
Y_MIN = 0.0
Y_MAX = 2.0

# None = automatyczne dopasowanie danej granicy


# =============================================================================
# ODCZYTYWANIE PARAMETROW Z NAZWY PLIKU
# =============================================================================

def odczytaj_parametry(plik):

    wzor = re.compile(
        r"^Tr_komutatorow_vs_g_"
        r"PBC(?P<PBC>true|false)"
        r"_L_x(?P<M>\d+)"
        r"_N_up(?P<N_up>\d+|all)"
        r"_J(?P<J>[-+0-9.eE]+)"
        r"_Delta(?P<Delta>[-+0-9.eE]+)"
        r"_omega0(?P<omega0>[-+0-9.eE]+)"
        r"_gmin(?P<gmin>[-+0-9.eE]+)"
        r"_gmax(?P<gmax>[-+0-9.eE]+)"
        r"_ng(?P<ng>\d+)"
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

    # Pierwszy wiersz zawiera nazwy kolumn
    with plik.open("r", encoding="utf-8-sig") as f:
        kolumny = f.readline().strip().split()

    if not kolumny or kolumny[0] != "g":
        raise ValueError(
            "Pierwszy wiersz pliku musi zaczynac sie od 'g'."
        )

    # Wczytanie danych liczbowych
    dane = np.loadtxt(plik, skiprows=1, ndmin=2)

    if dane.shape[1] != len(kolumny):
        raise ValueError(
            "Liczba kolumn naglowka nie zgadza sie z danymi."
        )

    return kolumny, dane


# =============================================================================
# ETYKIETY OPERATOROW
# =============================================================================

def etykieta_operatora(nazwa):

    # Zwykly operator, np. Q_8
    if re.fullmatch(r"Q_\d+", nazwa):
        return rf"$Q_{{{nazwa[2:]}}}$"

    # Kombinacja operatorow, np. Q_slow_3_50
    if nazwa.startswith("Q_slow_"):

        indeksy = ",".join(
            nazwa.removeprefix("Q_slow_").split("_")
        )

        return rf"$Q_{{\mathrm{{slow}},{indeksy}}}$"

    return nazwa


# =============================================================================
# RYSOWANIE WYKRESU
# =============================================================================

def rysuj_wykres(plik, wybrane_operatory):

    # Parametry automatycznie z nazwy pliku
    parametry = odczytaj_parametry(plik)

    kolumny, dane = wczytaj_dane(plik)

    # Sprawdzenie wybranych operatorow
    for nazwa in wybrane_operatory:

        if nazwa not in kolumny:
            raise ValueError(
                f"Nie znaleziono operatora: {nazwa}\n"
                f"Dostepne kolumny: {kolumny}"
            )

    if not wybrane_operatory:
        raise ValueError("Nie wybrano zadnego operatora.")

    # Pierwsza kolumna - g
    g = dane[:, kolumny.index("g")]

    fig, ax = plt.subplots(figsize=(9, 5.5))

    # Rysowanie wybranych operatorow
    for nazwa in wybrane_operatory:

        y = dane[:, kolumny.index(nazwa)]

        # Dla jednego punktu zawsze pokaz marker
        pokaz_marker = POKAZ_PUNKTY or len(g) == 1

        ax.plot(
            g,
            y,
            label=etykieta_operatora(nazwa),
            linewidth=1.8,
            marker="o" if pokaz_marker else None,
            markersize=3.5,
        )

    # Opisy osi
    ax.set_xlabel(r"$g$", fontsize=13)

    ax.set_ylabel(
        r"$\frac{1}{D}\mathrm{Tr}\!\left([H,Q]^\dagger[H,Q]\right)$",
        fontsize=13,
    )

    # Przygotowanie zapisu N_up
    if parametry["N_up"] == "all":
        n_up_opis = r"\mathrm{all}"
    else:
        n_up_opis = parametry["N_up"]

    # Tytul z parametrami odczytanymi z nazwy pliku
    ax.set_title(
        rf"$M={parametry['M']},\; "
        rf"N_\uparrow={n_up_opis},\; "
        rf"J={parametry['J']},\; "
        rf"\Delta={parametry['Delta']},\; "
        rf"\omega_0={parametry['omega0']},\; "
        rf"\mathrm{{PBC}}=\mathrm{{{parametry['PBC']}}}$",
        fontsize=12,
    )

    # Skala osi Y
    ax.set_yscale(SKALA_Y)

    # Zakres osi X i Y
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)

    # Siatka i legenda
    ax.grid(POKAZ_SIATKE, alpha=0.3)
    ax.legend(fontsize=10)

    fig.tight_layout()
    plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    rysuj_wykres(PLIK, WYBRANE_OPERATORY)
