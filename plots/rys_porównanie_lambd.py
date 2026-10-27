
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# USTAWIENIA
# =============================================================================

FOLDER = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag"
)

# Plik z przekątną macierzy komutatorów
PLIK_DIAGONAL = FOLDER / (
    "scripts/liom_mag_dane/siatka_omega_g/"
    "diagonal_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes.txt"
)

# Plik z komutatorami wybranych operatorów Q
PLIK_Q = FOLDER / (
    "scripts/fun_cor_data/"
    "Tr_komutatorow_vs_g_PBCtrue_L_x6_N_upall_J1.0_Delta0.001_omega00.5_gmin0.0_gmax10.0_ng41.data"
)


# =============================================================================
# WYBÓR OPERATORÓW
# =============================================================================

# Indeksy operatorów O_i z przekątnej
# Przykład: [1, 3, 8]
# Zakres: list(range(1, 21))
# Pusta lista [] oznacza brak krzywych z tego pliku

#LAMBDA_WYBOR = [3, 4, 8, 12, 50, 62, 118, 231, 237, 1128, 1170]
#LAMBDA_WYBOR = list(range(49, 100))
LAMBDA_WYBOR = [8]


# Operatory Q z drugiego pliku
# Można podawać nazwy albo numery (liczone od 1)
# Przykład: ["Q_3", "Q_8", "Q_slow_3_50"]
# Przykład: [1, 3, 5]
# Zakres: list(range(1, 6))

#Q_WYBOR = ["Q_3", "Q_4", "Q_8", "Q_12", "Q_50", "Q_62",
#           "Q_118", "Q_231", "Q_237", "Q_1128", "Q_1170"]

Q_WYBOR = ["Q_8"]


# Wartość omega_0 dla danych z przekątnej
OMEGA_WYBRANA = 0.5


# =============================================================================
# USTAWIENIA WYKRESU
# =============================================================================

SKALA_Y = "linear"        # "linear" albo "log"

X_LIM = (0, 10)           # None -> automatyczny zakres

# Ograniczenia osi Y - niezależne
Y_MIN = 0.0               # None -> automatyczne minimum
Y_MAX = 10.0              # None -> automatyczne maksimum

FIGSIZE = (10, 6)
DPI = 300

ZAPISZ_WYKRES = True
POKAZ_WYKRES = True

FOLDER_WYKRESY = FOLDER / "plots"
NAZWA_WYKRESU = "porownanie_lambda_O_Q.png"


# =============================================================================
# WCZYTYWANIE PRZEKĄTNEJ
# =============================================================================

def wczytaj_diagonal(filename, omega_wybrana):

    data = np.loadtxt(filename, ndmin=2)

    if data.shape[1] < 3:
        raise ValueError(
            "Plik diagonal musi zawierać omega_0, g "
            "i co najmniej jedną lambdę."
        )

    # Kolumny:
    # 0 -> omega_0
    # 1 -> g
    # 2 -> lambda_O1
    # 3 -> lambda_O2
    # ...

    omega = data[:, 0]

    mask = np.isclose(
        omega,
        omega_wybrana,
        atol=1e-10,
        rtol=0
    )

    if not np.any(mask):
        raise ValueError(
            f"Nie znaleziono omega_0 = {omega_wybrana}. "
            f"Dostępne wartości: {np.unique(omega)}"
        )

    data = data[mask]

    # Sortowanie rosnąco po g
    data = data[np.argsort(data[:, 1])]

    g = data[:, 1]
    lambdy = data[:, 2:]

    return g, lambdy


# =============================================================================
# WCZYTYWANIE PLIKU Q
# =============================================================================

def wczytaj_Q(filename):

    # Pierwszy wiersz zawiera nazwy:
    # g Q_3 Q_4 Q_8 ...

    with open(filename, "r", encoding="utf-8-sig") as f:
        naglowek = f.readline().split()

    if not naglowek or naglowek[0] != "g":
        raise ValueError(
            "Niepoprawny nagłówek pliku Q."
        )

    nazwy_Q = naglowek[1:]

    data = np.loadtxt(
        filename,
        skiprows=1,
        ndmin=2
    )

    if data.shape[1] != len(naglowek):
        raise ValueError(
            "Liczba kolumn danych nie zgadza się "
            "z nagłówkiem pliku Q."
        )

    # Sortowanie rosnąco po g
    data = data[np.argsort(data[:, 0])]

    g = data[:, 0]
    wartosci_Q = data[:, 1:]

    return g, wartosci_Q, nazwy_Q


# =============================================================================
# WYBÓR OPERATORÓW Q
# =============================================================================

def indeksy_Q_z_wyboru(wybor, nazwy_Q):

    indeksy = []

    for element in wybor:

        if isinstance(element, str):

            if element not in nazwy_Q:
                raise ValueError(
                    f"Nie znaleziono operatora {element}.\n"
                    f"Dostępne: {nazwy_Q}"
                )

            indeks = nazwy_Q.index(element)

        else:

            # Numeracja od 1
            indeks = int(element) - 1

            if indeks < 0 or indeks >= len(nazwy_Q):
                raise ValueError(
                    f"Numer Q = {element} poza zakresem "
                    f"1...{len(nazwy_Q)}"
                )

        if indeks not in indeksy:
            indeksy.append(indeks)

    return indeksy


# =============================================================================
# RYSOWANIE
# =============================================================================

def main():

    # Wczytujemy oba pliki
    g_lambda, lambdy = wczytaj_diagonal(
        PLIK_DIAGONAL,
        OMEGA_WYBRANA
    )

    g_Q, wartosci_Q, nazwy_Q = wczytaj_Q(PLIK_Q)

    n_operatorow = lambdy.shape[1]

    # Sprawdzenie indeksów O_i
    for i in LAMBDA_WYBOR:

        if (
            not isinstance(i, (int, np.integer))
            or not 1 <= i <= n_operatorow
        ):
            raise ValueError(
                f"Nieprawidłowy indeks O_{i}. "
                f"Dostępny zakres: 1...{n_operatorow}"
            )

    indeksy_Q = indeksy_Q_z_wyboru(
        Q_WYBOR,
        nazwy_Q
    )

    if not LAMBDA_WYBOR and not indeksy_Q:
        raise ValueError(
            "Nie wybrano żadnego operatora."
        )

    print(f"Liczba elementów bazy: {n_operatorow}")
    print(f"Liczba punktów g (diagonal): {len(g_lambda)}")
    print(f"Liczba punktów g (Q): {len(g_Q)}")
    print(f"Wybrane O_i: {LAMBDA_WYBOR}")
    print(f"Wybrane Q: {[nazwy_Q[i] for i in indeksy_Q]}")

    # -------------------------------------------------------------------------
    # WYKRES
    # -------------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # -------------------------------------------------------------------------
    # WSPÓLNA PALETA KOLORÓW
    # Każda linia otrzymuje inny kolor, niezależnie od pliku
    # -------------------------------------------------------------------------

    liczba_linii = len(LAMBDA_WYBOR) + len(indeksy_Q)

    if liczba_linii <= 10:

        kolory = plt.cm.tab10(
            np.arange(liczba_linii)
        )

    elif liczba_linii <= 20:

        kolory = plt.cm.tab20(
            np.arange(liczba_linii)
        )

    else:

        kolory = plt.cm.turbo(
            np.linspace(0.05, 0.95, liczba_linii)
        )

    # -------------------------------------------------------------------------
    # OPERATORY BAZOWE O_i - linie ciągłe
    # -------------------------------------------------------------------------

    for k, i in enumerate(LAMBDA_WYBOR):

        # O_1 -> kolumna 0 tablicy lambdy
        y = lambdy[:, i - 1]

        ax.plot(
            g_lambda,
            y,
            color=kolory[k],
            linewidth=2.0,
            linestyle="-",
            label=rf"$\lambda_{{O_{{{i}}}}}$"
        )

    # -------------------------------------------------------------------------
    # OPERATORY Q - linie przerywane
    # -------------------------------------------------------------------------

    for k, indeks in enumerate(indeksy_Q):

        nazwa = nazwy_Q[indeks]
        y = wartosci_Q[:, indeks]

        # Kontynuacja wspólnej palety kolorów
        numer_koloru = len(LAMBDA_WYBOR) + k

        ax.plot(
            g_Q,
            y,
            color=kolory[numer_koloru],
            linewidth=1.8,
            linestyle="--",
            label=nazwa
        )

    # -------------------------------------------------------------------------
    # FORMATOWANIE
    # -------------------------------------------------------------------------

    ax.set_xlabel(
        r"$g$",
        fontsize=14
    )

    ax.set_ylabel(
        r"$\mathrm{Tr}([H,O]^\dagger[H,O])/D$",
        fontsize=14
    )

    ax.set_title(
        rf"$\omega_0={OMEGA_WYBRANA}$",
        fontsize=14,
        pad=12
    )

    ax.set_yscale(SKALA_Y)

    if X_LIM is not None:
        ax.set_xlim(X_LIM)

    # Niezależne ograniczenia dolne i górne osi Y
    if Y_MIN is not None:
        ax.set_ylim(bottom=Y_MIN)

    if Y_MAX is not None:
        ax.set_ylim(top=Y_MAX)

    ax.grid(
        True,
        alpha=0.25,
        linestyle="--"
    )

    ax.tick_params(
        axis="both",
        labelsize=11,
        direction="in"
    )

    # Legenda z prawej strony wykresu
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=10,
        frameon=True,
        ncol=1
    )

    fig.tight_layout()

    # -------------------------------------------------------------------------
    # ZAPIS
    # -------------------------------------------------------------------------

    if ZAPISZ_WYKRES:

        FOLDER_WYKRESY.mkdir(
            parents=True,
            exist_ok=True
        )

        sciezka = FOLDER_WYKRESY / NAZWA_WYKRESU

        fig.savefig(
            sciezka,
            dpi=DPI,
            bbox_inches="tight"
        )

        print(f"Zapisano wykres: {sciezka}")

    if POKAZ_WYKRES:
        plt.show()

    plt.close(fig)


if __name__ == "__main__":
    main()
