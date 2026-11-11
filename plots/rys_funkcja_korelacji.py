from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# =============================================================================
# WSPOLNE FORMATOWANIE - TAK JAK NA POPRZEDNICH WYKRESACH
# =============================================================================

# Wspólny rozmiar wszystkich napisów na wykresie.
FONT_SIZE = 14
plt.rcParams["font.size"] = FONT_SIZE

# =============================================================================
# USTAWIENIA
# =============================================================================

# Folder z wynikami pobranymi z serwera mm8.
FOLDER_DANE = Path(
    #r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data\doróbki"
    #r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_11_07.11.2026\jeden_elektron"
    #r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_11_07.11.2026\wszystkie_magnetyzacje"
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_12_17.11.2026\half filling"
)

# TRYB WYKRESU: "omega0", "g" albo "omega0_g".
TRYB_SKANU = "g"

# Niezalezne wlaczanie metod korelacji.
POKAZ_RANDOM_STATE = False
POKAZ_TRACE = True
POKAZ_FERMION_LOKALNY = False

# Pozycja elektronu w trzeciej metodzie (zgodna z x w nazwie pliku).
POZYCJA_FERMIONU = 5

# Mozna wskazac konkretne pliki. None = automatyczne wyszukiwanie.
PLIK = None
PLIK_TRACE = None
PLIK_FERMION_LOKALNY = None

# WYBOR WARTOSCI SPRZEZENIA (dla skanowania g)
G_WYBRANE = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 4.0, 6.0, 8.0, 10.0]

# WYBOR WARTOSCI omega_0 (dla nowego skanowania omega_0)
# None: automatycznie rysuj wszystkie wartosci omega0 z pliku.
# OMEGA0_WYBRANE = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0]
# OMEGA0_WYBRANE = [0.0, 0.4, 0.8, 1.0, 2.0]
OMEGA0_WYBRANE = None

# Stale g na potrzeby kontroli pliku i ewentualnego skanu omega0_g.
G_STALE = 0.65
omega0_STALE = 0.5

# Plik w dawnej nazwie ma ng11 dla zakresu 0..10, co zwykle oznacza krok 1.
# Wartosci 0.2, 0.4, ... moga wiec nie wystepowac w starym pliku.
# True: pomin brakujace wartosci parametru, wypisz ostrzezenie i rysuj dostepne.
# False: zatrzymaj program, gdy ktorejs wartosci parametru brakuje.
POMIN_BRAKUJACE_G = True

# WYBOR FUNKCJI KORELACJI
# Tutaj wpisujesz RZECZYWISTE operatory, ktore chcesz narysowac.
WYBRANE_KORELACJE = [
    "C_slow_289_50_1_304_300_3",
    #"C_slow_289_50_1_304_300",
    #"C_slow_289_50_1_304",
    #"C_slow_289_50_1",
    #"C_slow_289_50",
    # "Q_289",
    # "Q_50",
    # "Q_1",
    # "Q_304",
    # "Q_300",
    # "Q_3",
]

# Plik z opisem elementow bazy operatorow (indeks, sektor RR/RI/IR/II, slowo ferm|bos).
PLIK_BAZY = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\seminarium\data\operators_ladder_mag_M_3_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_1.0_g_1.0_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true.txt"
)

# True: pokaz pod wykresem rozwiniecia spinowe Q_i wystepujacych w wybranych korelacjach.
POKAZ_OPIS_OPERATOROW = True
FONT_SIZE_OPERATOROW = FONT_SIZE - 1

# Wartość omega_0 wyświetlana w tytule (nie zmienia wczytywanych danych).
OMEGA_0_TYTUL = 0.65  # Uzywane jako rezerwowy podpis dla dawnego skanu g

# -----------------------------------------------------------------------------
# POPRAWKA DLA PLIKOW ng11
# -----------------------------------------------------------------------------
# W plikach wygenerowanych pokazanym kodem Julii kolejność zapisanych danych
# nie zgadza się z nazwami kolumn. Tutaj podajemy:
#   KLUCZ   = rzeczywisty operator, który chcesz wyświetlić
#   WARTOSC = nazwa kolumny w pliku, z której trzeba te dane odczytać
#
# Dotyczy to obu plików: random state i trace.
POPRAW_BLEDNE_NAGLOWKI = False  # NOWA Julia ma poprawne naglowki; True tylko do starych plikow

KOLUMNA_DLA_OPERATORA = {
    "Q_3": "Q_60",
    "Q_8": "C_slow_3_50",
    "Q_26": "C_slow_8_12_62_118",
    "Q_60": "C_slow_8_4_231_237_62_118_1128_1170",
    "C_slow_3_50": "Q_3",
    "C_slow_8_12_62_118": "Q_8",
    "C_slow_8_4_231_237_62_118_1128_1170": "Q_26",
}

# OPCJE WYKRESU
SKALA_Y = "linear"
POKAZ_PUNKTY = False
POKAZ_SIATKE = True

# ZAKRES OSI X (t) - wpisz granice czasu, ktore maja byc widoczne.
# Przyklad: X_MIN = 0.0 i X_MAX = 20.0 pokazuje tylko t od 0 do 20.
# X_MAX = None wyswietla dane az do konca dostepnego zakresu.
X_MIN = 0.0
X_MAX = 100.0

# ZAKRES OSI Y
Y_MIN = None
Y_MAX = None

# Format poziomego wykresu do prezentacji; miejsce po prawej na legende.
FIGSIZE = (15.0, 6.5)


# =============================================================================
# WYKRES / ZAPIS
# =============================================================================

show_plot = True
save_plot = False    # Zapisuje PDF z bialym, nieprzezroczystym tlem.
save_png = False    # Dodatkowo przezroczysty PNG; False = tylko PDF.
dpi = 300

LINE_WIDTH = 2.0
LINE_WIDTH_LEGEND = 2.2

# PDF i PNG zapisujemy w folderze, w ktorym znajduje sie plik .py.
output_folder = Path(__file__).resolve().parent


# =============================================================================
# ODCZYTYWANIE PARAMETROW Z NAZWY PLIKU
# =============================================================================

def odczytaj_parametry(plik):
    # Obslugujemy pliki random state, trace oraz lokalnego elektronu.
    wzor = re.compile(
        r"^Korelacje_(?:(?P<metoda>trace|fermion_lokalny)_)?ED_vs_"
        r"(?P<skan>omega0_g|omega0|g)_t_"
        r"PBC(?P<PBC>true|false)(?P<reszta>.*)\.data$"
    )
    wynik = wzor.fullmatch(plik.name)
    if wynik is None:
        raise ValueError(f"Nie udalo sie odczytac parametrow z nazwy:\n{plik.name}")

    parametry = {
        "metoda": wynik.group("metoda") or "random",
        "skan": wynik.group("skan"),
        "PBC": wynik.group("PBC"),
    }
    wzor_pola = re.compile(
        r"_(?P<klucz>L_x|N_up|J|Delta|omega0min|omega0max|omega0|"
        r"gmin|gmax|ng|nw|x|tmin|tmax|nt)(?P<wartosc>[^_]+)"
    )
    for m in wzor_pola.finditer(wynik.group("reszta")):
        parametry[m.group("klucz")] = m.group("wartosc")

    if "L_x" not in parametry:
        raise ValueError(f"Brak dlugosci L_x w nazwie: {plik.name}")
    parametry["M"] = parametry["L_x"]
    return parametry


def aktywne_metody():
    """Zwraca wlaczone metody w stalej kolejnosci rysowania."""
    metody = []
    if POKAZ_RANDOM_STATE:
        metody.append(("random", PLIK))
    if POKAZ_TRACE:
        metody.append(("trace", PLIK_TRACE))
    if POKAZ_FERMION_LOKALNY:
        metody.append(("fermion_lokalny", PLIK_FERMION_LOKALNY))
    if not metody:
        raise ValueError("Wszystkie metody korelacji sa wylaczone.")
    return metody


def znajdz_pliki_metod(folder, skan):
    """Znajduje osobne pliki dla wlaczonych metod i sprawdza, czy opisuja ten sam skan."""
    konfiguracja = {
        "random": ("Korelacje_ED_vs_", PLIK),
        "trace": ("Korelacje_trace_ED_vs_", PLIK_TRACE),
        "fermion_lokalny": ("Korelacje_fermion_lokalny_ED_vs_", PLIK_FERMION_LOKALNY),
    }

    znalezione = {}
    parametry_metod = {}

    for metoda, plik_reczny in aktywne_metody():
        prefiks, _ = konfiguracja[metoda]

        if plik_reczny is not None:
            plik = Path(plik_reczny)
        else:
            kandydaci = list(folder.glob(f"{prefiks}{skan}_t_*.data"))

            if skan == "omega0":
                # Przy pojedynczym g w pliku zapisano gmin=gmax=G_STALE i ng=1.
                pasujace = []
                for kandydat in kandydaci:
                    p = odczytaj_parametry(kandydat)
                    if (
                        np.isclose(float(p.get("gmin", "nan")), G_STALE)
                        and np.isclose(float(p.get("gmax", "nan")), G_STALE)
                        and p.get("ng") == "1"
                    ):
                        pasujace.append(kandydat)
                kandydaci = pasujace

            if metoda == "fermion_lokalny":
                znacznik_x = f"_x{POZYCJA_FERMIONU}_tmin"
                kandydaci = [p for p in kandydaci if znacznik_x in p.name]

            if not kandydaci:
                raise FileNotFoundError(
                    f"Nie znaleziono pliku dla metody '{metoda}' i skanu '{skan}' "
                    f"w folderze:\n{folder}"
                )
            plik = max(kandydaci, key=lambda p: p.stat().st_mtime)

        if not plik.is_file():
            raise FileNotFoundError(f"Nie znaleziono pliku dla metody '{metoda}':\n{plik}")

        parametry = odczytaj_parametry(plik)
        if parametry["metoda"] != metoda:
            raise ValueError(
                f"Plik przypisany do metody '{metoda}' ma typ '{parametry['metoda']}':\n{plik}"
            )
        if parametry["skan"] != skan:
            raise ValueError(
                f"Wybrany TRYB_SKANU={skan}, ale plik ma skan {parametry['skan']}:\n{plik}"
            )
        if metoda == "fermion_lokalny":
            x_z_pliku = int(parametry.get("x", -1))
            if x_z_pliku != POZYCJA_FERMIONU:
                raise ValueError(
                    f"Plik lokalnego elektronu ma x={x_z_pliku}, "
                    f"a POZYCJA_FERMIONU={POZYCJA_FERMIONU}:\n{plik}"
                )

        znalezione[metoda] = plik
        parametry_metod[metoda] = parametry
        print(f"{metoda:18s}: {plik}")

    # Wszystkie wlaczone pliki musza dotyczyc tych samych parametrow fizycznych i skanu.
    klucze_wspolne = (
        "skan", "PBC", "L_x", "N_up", "J", "Delta", "omega0",
        "omega0min", "omega0max", "gmin", "gmax", "ng", "nw",
        "tmin", "tmax", "nt",
    )
    metoda_ref = next(iter(parametry_metod))
    ref = parametry_metod[metoda_ref]

    for metoda, parametry in parametry_metod.items():
        for klucz in klucze_wspolne:
            a = ref.get(klucz)
            b = parametry.get(klucz)
            if a is None and b is None:
                continue
            if a is None or b is None:
                raise ValueError(
                    f"Pliki metod '{metoda_ref}' i '{metoda}' maja rozne parametry ({klucz})."
                )
            try:
                zgodne = np.isclose(float(a), float(b), rtol=0, atol=1e-10)
            except ValueError:
                zgodne = (a == b)
            if not zgodne:
                raise ValueError(
                    f"Pliki metod '{metoda_ref}' i '{metoda}' nie pochodza z tego samego skanu: "
                    f"{klucz}={a} oraz {klucz}."
                )

    return znalezione


# =============================================================================
# WCZYTYWANIE DANYCH
# =============================================================================

def wczytaj_dane(plik):
    if not plik.is_file():
        raise FileNotFoundError(f"Nie znaleziono pliku: {plik}")

    with plik.open("r", encoding="utf-8-sig") as f:
        kolumny = f.readline().strip().split()

    # Stara seria: g,t,...; nowa seria: omega0,t,...;
    # skan dwuwymiarowy: omega0,g,t,...
    dostepne_naglowki = (["g", "t"], ["omega0", "t"], ["omega0", "g", "t"])
    if len(kolumny) < 3 or not any(kolumny[:len(w)] == w for w in dostepne_naglowki):
        raise ValueError(f"Nieprawidlowe poczatkowe kolumny w {plik.name}: {kolumny[:4]}")

    dane = np.loadtxt(plik, skiprows=1, ndmin=2)
    if dane.shape[1] != len(kolumny):
        raise ValueError("Liczba kolumn naglowka nie zgadza sie z danymi.")
    return kolumny, dane


# =============================================================================
# WYBOR WARTOSCI g / omega_0
# =============================================================================

def wybierz_dane_parametru(kolumny, dane, nazwa_parametru, wartosc, g_stale=None):
    maska = np.isclose(dane[:, kolumny.index(nazwa_parametru)], wartosc,
                       rtol=0, atol=1e-8)
    if g_stale is not None and "g" in kolumny:
        maska &= np.isclose(dane[:, kolumny.index("g")], g_stale,
                            rtol=0, atol=1e-8)
    if not np.any(maska):
        raise ValueError(f"Brak danych dla {nazwa_parametru}={wartosc} (g={g_stale})")
    wybrane = dane[maska]
    return wybrane[np.argsort(wybrane[:, kolumny.index("t")])]


def sprawdz_dostepne_wartosci(zbiory, wybrane, nazwa_parametru, g_stale=None):
    """Zostawia wartosci dostepne we wszystkich wlaczonych plikach."""
    poprawne = []
    brakujace = []

    for x in wybrane:
        dostepna_we_wszystkich = True
        for metoda, (nazwy, tab) in zbiory.items():
            if nazwa_parametru not in nazwy:
                raise ValueError(
                    f"Brak kolumny '{nazwa_parametru}' w pliku metody '{metoda}'."
                )
            maska = np.isclose(
                tab[:, nazwy.index(nazwa_parametru)], x, rtol=0, atol=1e-8
            )
            if g_stale is not None and "g" in nazwy:
                maska &= np.isclose(
                    tab[:, nazwy.index("g")], g_stale, rtol=0, atol=1e-8
                )
            if not np.any(maska):
                dostepna_we_wszystkich = False
                break

        if dostepna_we_wszystkich:
            poprawne.append(x)
        else:
            brakujace.append(x)

    if brakujace:
        komunikat = f"Brak wskazanych wartosci {nazwa_parametru}: {brakujace}"
        if POMIN_BRAKUJACE_G:
            print("UWAGA:", komunikat, "- pominieto")
        else:
            raise ValueError(komunikat)

    if not poprawne:
        raise ValueError(
            "Zadna wybrana wartosc parametru nie wystepuje we wszystkich wlaczonych plikach."
        )
    return poprawne


# =============================================================================
# MAPOWANIE RZECZYWISTY OPERATOR -> KOLUMNA W PLIKU
# =============================================================================

def nazwa_kolumny_w_pliku(nazwa_operatora):
    """
    Zwraca nazwę kolumny w pliku, z której należy odczytać dane dla
    rzeczywistego operatora 'nazwa_operatora'.
    """
    if POPRAW_BLEDNE_NAGLOWKI:
        return KOLUMNA_DLA_OPERATORA.get(nazwa_operatora, nazwa_operatora)
    return nazwa_operatora


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


def etykieta_operatora(nazwa):
    """Podpis operatora Q przeznaczony do tytułu wykresu."""

    if nazwa.startswith("C_slow_"):
        indeksy = ",".join(nazwa.removeprefix("C_slow_").split("_"))
        return rf"\hat Q_{{\mathrm{{slow}},\,{indeksy}}}"

    if nazwa.startswith("Q_"):
        indeks = nazwa.removeprefix("Q_")
        return rf"\hat Q_{{{indeks}}}"

    return rf"\mathrm{{{nazwa}}}"


def etykieta_do_legendy(nazwa_operatora, wartosc, liczba_operatorow, liczba_wartosci,
                        nazwa_parametru):
    """
    Dobiera podpis legendy: rozne omega_0 (lub g) dla jednego operatora,
    wzglednie wiele operatorow dla jednego parametru.
    """
    parametr = r"\omega_0" if nazwa_parametru == "omega0" else "g"
    if liczba_wartosci == 1 and liczba_operatorow > 1:
        return etykieta_korelacji(nazwa_operatora)
    if liczba_wartosci > 1 and liczba_operatorow == 1:
        return rf"${parametr}={wartosc:g}$"
    return rf"{etykieta_korelacji(nazwa_operatora)}, ${parametr}={wartosc:g}$"


# =============================================================================
# ODCZYT OPERATOROW BAZY I ROZWINIECIE NA OPERATORY SPINOWE
# =============================================================================

def indeksy_operatorow_w_korelacjach(wybrane_korelacje):
    """Zwraca Q_i z Q_i lub C_slow_i_j_..., w kolejnosci pierwszego wystapienia."""

    indeksy = []
    for nazwa in wybrane_korelacje:
        if nazwa.startswith("Q_"):
            czesci = [nazwa.removeprefix("Q_")]
        elif nazwa.startswith("C_slow_"):
            czesci = nazwa.removeprefix("C_slow_").split("_")
        else:
            continue

        for czesc in czesci:
            if not czesc.isdigit():
                raise ValueError(f"Niepoprawny indeks operatora w nazwie {nazwa!r}.")
            indeks = int(czesc)
            if indeks not in indeksy:
                indeksy.append(indeks)

    return indeksy


def wczytaj_opisy_operatorow(plik_bazy, potrzebne_indeksy):
    """
    Format wiersza (przykladowo): 1 II 130|320.
    Pierwsza kolumna NIE oznacza numeru Q_i i jest ignorowana.
    Q_1 to pierwszy wpis operatora w pliku, Q_2 to drugi itd.
    Naglowki, puste wiersze i komentarze nie sa liczone jako operatory.
    0=I, 1=S+, 2=Sz, 3=S-; RR/RI/IR/II dotyczy obu nog drabiny.
    """

    if not plik_bazy.is_file():
        raise FileNotFoundError(f"Nie znaleziono pliku z baza operatorow: {plik_bazy}")

    potrzebne = set(potrzebne_indeksy)
    znalezione = {}

    # Pierwsze pole w ogole nas nie interesuje; czytamy tylko sektor i slowo.
    wzor = re.compile(
        r"^\s*\S+\s+(RR|RI|IR|II)\s+([0-3]+)\s*\|\s*([0-3]+)(?=\s|$)"
    )

    numer_operatora = 0
    with plik_bazy.open("r", encoding="utf-8-sig") as f:
        for numer_linii, linia in enumerate(f, 1):
            m = wzor.match(linia)
            if m is None:
                continue  # naglowki, komentarze, puste wiersze i inne nie-operatorowe linie

            # WAZNE: numer operatora to jego kolejnosc w bazie, nie pierwsza kolumna.
            numer_operatora += 1
            if numer_operatora in potrzebne:
                znalezione[numer_operatora] = (m.group(1), m.group(2), m.group(3))
                print(
                    f"  Q_{numer_operatora} (linia pliku {numer_linii}): "
                    f"{m.group(1)} {m.group(2)}|{m.group(3)}"
                )

    brakujace = sorted(potrzebne - set(znalezione))
    if brakujace:
        raise ValueError(
            f"Brak operatorow o indeksach {brakujace}. "
            f"Znaleziono lacznie {numer_operatora} wpisow bazy w pliku {plik_bazy}."
        )

    return znalezione


def operatory_na_nodze(slowo, noga, sprzezenie_hermitowskie=False):
    """Zapisuje slowo 0/1/2/3 jako iloczyn spinow w LaTeX (bez normalizacji)."""

    znaki = {
        "0": None,
        "1": "-" if sprzezenie_hermitowskie else "+",
        "2": "z",
        "3": "+" if sprzezenie_hermitowskie else "-",
    }

    czynniki = []
    for pozycja, cyfra in enumerate(slowo):
        znak = znaki[cyfra]
        if znak is None:
            continue
        pozycja_l = "l" if pozycja == 0 else f"l+{pozycja}"
        czynniki.append(rf"S^{{{znak}}}_{{{pozycja_l},{noga}}}")

    return "".join(czynniki) if czynniki else r"\mathbb{1}"


def czesc_sektorowa(slowo, noga, sektor):
    """
    R: A + A†, I: i(A - A†). Pomijamy wspolne wspolczynniki normalizacji.
    Dla operatora samosprzezonego w sektorze R wystarczy A.
    """

    a = operatory_na_nodze(slowo, noga)
    a_dag = operatory_na_nodze(slowo, noga, sprzezenie_hermitowskie=True)

    if sektor == "R":
        if a == a_dag:
            return a
        return rf"\left({a}+{a_dag}\right)"
    if sektor == "I":
        if a == a_dag:
            raise ValueError(f"Sektor I jest zerowy dla slowa {slowo} na nodze {noga}.")
        return rf"\mathrm{{i}}\left({a}-{a_dag}\right)"

    raise ValueError(f"Nieznany sektor {sektor!r}.")


def wzor_spinowy_operatora(indeks, opis):
    """Buduje zapis do Matplotlib MathText: Q_i ∝ Σ_l [czesc nogi 1][czesc nogi 2]."""

    sektor, ferm, bos = opis
    if len(sektor) != 2 or len(ferm) != len(bos):
        raise ValueError(f"Niepoprawny opis Q_{indeks}: {opis}")

    pierwsza = czesc_sektorowa(ferm, 1, sektor[0])
    druga = czesc_sektorowa(bos, 2, sektor[1])
    return rf"$\hat{{Q}}_{{{indeks}}}\propto\sum_l {pierwsza}\,{druga}$"


def dodaj_opisy_pod_wykresem(fig, indeksy, opisy_bazy, wysokosc_rysunku, wysokosc_bloku):
    """Jeden wiersz na element bazy. Polozenie w calach, niezalezne od liczby g."""

    if not indeksy:
        return

    # fig.text(
    #     0.075,
    #     (wysokosc_bloku - 0.08) / wysokosc_rysunku,
    #     "Elementy bazy operatorow (z dokladnoscia do normalizacji):",
    #     fontsize=FONT_SIZE_OPERATOROW,
    #     ha="left", va="top"
    # )

    for nr, indeks in enumerate(indeksy):
        etykieta = wzor_spinowy_operatora(indeks, opisy_bazy[indeks])
        fig.text(
            0.075,
            (wysokosc_bloku - 0.52 - 0.48 * nr) / wysokosc_rysunku,
            etykieta,
            fontsize=FONT_SIZE_OPERATOROW,
            ha="left", va="top"
        )


# =============================================================================
# RYSOWANIE WSZYSTKICH omega_0 / g NA JEDNYM WYKRESIE
# =============================================================================

def rysuj_wykres(pliki, wybrane_korelacje, wartosci_wybrane=None):
    # Wczytanie wszystkich wlaczonych metod.
    parametry_metod = {}
    zbiory = {}

    for metoda, plik in pliki.items():
        parametry_metod[metoda] = odczytaj_parametry(plik)
        zbiory[metoda] = wczytaj_dane(plik)

    metoda_ref = next(iter(pliki))
    parametry = parametry_metod[metoda_ref]
    kolumny = zbiory[metoda_ref][0]
    dane = zbiory[metoda_ref][1]

    if not wybrane_korelacje:
        raise ValueError("Nie wybrano zadnej korelacji.")

    # Dla skanu 2D rysujemy krzywe omega0 przy ustalonym G_STALE.
    nazwa_parametru = "g" if TRYB_SKANU == "g" else "omega0"
    g_filtr = G_STALE if TRYB_SKANU == "omega0_g" else None

    if TRYB_SKANU == "omega0":
        g_min = float(parametry.get("gmin", "nan"))
        g_max = float(parametry.get("gmax", "nan"))
        if not np.isclose(g_min, G_STALE) or not np.isclose(g_max, G_STALE):
            raise ValueError(
                f"Wybrano plik z g={g_min:g}..{g_max:g}, a G_STALE={G_STALE:g}"
            )

    # Sprawdzamy, czy wszystkie wybrane kolumny wystepuja w kazdym wlaczonym pliku.
    for metoda, (kolumny_metody, _) in zbiory.items():
        for nazwa in wybrane_korelacje:
            nazwa_w_pliku = nazwa_kolumny_w_pliku(nazwa)
            if nazwa_w_pliku not in kolumny_metody:
                raise ValueError(
                    f"Brak korelacji {nazwa} (kolumna: {nazwa_w_pliku}) "
                    f"w pliku metody '{metoda}'.\n"
                    f"Dostepne kolumny: {kolumny_metody}"
                )

    wartosci_dostepne = np.unique(dane[:, kolumny.index(nazwa_parametru)])
    if wartosci_wybrane is None:
        wartosci_wybrane = list(wartosci_dostepne)
    print(f"Dostepne wartosci {nazwa_parametru}:", wartosci_dostepne)
    print(f"Wybrane wartosci {nazwa_parametru}:", wartosci_wybrane)

    wartosci_do_rysowania = sprawdz_dostepne_wartosci(
        zbiory, wartosci_wybrane, nazwa_parametru, g_stale=g_filtr
    )
    print("Rysowane wartosci:", wartosci_do_rysowania)

    if POPRAW_BLEDNE_NAGLOWKI:
        print("\nMapowanie operator -> kolumna w pliku:")
        for nazwa in wybrane_korelacje:
            print(f"  {nazwa} -> {nazwa_kolumny_w_pliku(nazwa)}")
        print()

    # Wypisujemy wszystkie Q_i wystepujace w wybranych korelacjach,
    # rowniez wszystkie skladniki operatorow C_slow_i_j_... .
    indeksy_bazy = indeksy_operatorow_w_korelacjach(wybrane_korelacje) if POKAZ_OPIS_OPERATOROW else []
    opisy_bazy = wczytaj_opisy_operatorow(PLIK_BAZY, indeksy_bazy) if indeksy_bazy else {}

    # -------------------------------------------------------------------------
    # KOLORY - TEN SAM KOLOR DLA TEGO SAMEGO OPERATORA I WARTOSCI PARAMETRU.
    # Metody rozrozniamy stylem linii.
    # -------------------------------------------------------------------------
    liczba_krzywych = len(wartosci_do_rysowania) * len(wybrane_korelacje)
    if liczba_krzywych <= 10:
        kolory = plt.cm.tab10(np.arange(liczba_krzywych))
    else:
        kolory = plt.cm.turbo(np.linspace(0.05, 0.95, liczba_krzywych))

    style_metod = {
        "random": ("Random state", "--"),
        "trace": ("Trace", "-"),
        "fermion_lokalny": ("Lokalny elektron", ":"),
    }

    # -------------------------------------------------------------------------
    # WYKRES
    # -------------------------------------------------------------------------
    # Zwiekszamy figure w pionie tylko o tyle, ile potrzeba na opisy Q_i.
    # Sam obszar danych zachowuje dotychczasowa wysokosc.
    wysokosc_bloku = 0.65 + 0.48 * len(indeksy_bazy) if indeksy_bazy else 0.0
    wysokosc_rysunku = FIGSIZE[1] + wysokosc_bloku
    fig, ax = plt.subplots(figsize=(FIGSIZE[0], wysokosc_rysunku))

    # Biale tlo figury i osi.
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Miejsce na legende po prawej, proporcje jak w wykresach do seminarium.
    if indeksy_bazy:
        fig.subplots_adjust(
            left=0.075,
            right=0.715,
            bottom=(0.14 * FIGSIZE[1] + wysokosc_bloku) / wysokosc_rysunku,
            top=1.0 - (0.14 * FIGSIZE[1]) / wysokosc_rysunku
        )
    else:
        fig.subplots_adjust(left=0.075, right=0.715, bottom=0.14, top=0.86)

    numer_koloru = 0
    liczba_operatorow = len(wybrane_korelacje)
    liczba_wartosci = len(wartosci_do_rysowania)

    for wartosc in wartosci_do_rysowania:
        # Dane dla kazdej wlaczonej metody przy tej samej wartosci parametru.
        dane_metod_x = {}
        for metoda, (kolumny_metody, dane_metody) in zbiory.items():
            dane_metod_x[metoda] = wybierz_dane_parametru(
                kolumny_metody, dane_metody, nazwa_parametru, wartosc, g_stale=g_filtr
            )

        wartosc_rzeczywista = dane_metod_x[metoda_ref][
            0, zbiory[metoda_ref][0].index(nazwa_parametru)
        ]

        for nazwa in wybrane_korelacje:
            nazwa_w_pliku = nazwa_kolumny_w_pliku(nazwa)
            kolor = kolory[numer_koloru]
            numer_koloru += 1

            podpis = etykieta_do_legendy(
                nazwa, wartosc_rzeczywista, liczba_operatorow,
                liczba_wartosci, nazwa_parametru
            )

            # Kazda metoda ma ten sam kolor dla danego operatora i parametru,
            # ale inny styl linii.
            for nr_metody, (metoda, (kolumny_metody, _)) in enumerate(zbiory.items()):
                dane_x = dane_metod_x[metoda]
                t = dane_x[:, kolumny_metody.index("t")]
                C_t = dane_x[:, kolumny_metody.index(nazwa_w_pliku)]
                etykieta_metody, styl_linii = style_metod[metoda]

                ax.plot(
                    t,
                    C_t,
                    color=kolor,
                    linestyle=styl_linii,
                    linewidth=LINE_WIDTH,
                    marker="o" if (POKAZ_PUNKTY and metoda == "random") else None,
                    markersize=3.5,
                    label=podpis if nr_metody == 0 else "_nolegend_"
                )

    # -------------------------------------------------------------------------
    # OPISY OSI
    # -------------------------------------------------------------------------
    ax.set_xlabel(r"$t$", fontsize=FONT_SIZE, labelpad=3)
    ax.set_ylabel(r"$C_Q(t)$", fontsize=FONT_SIZE, labelpad=4)

    # Parametr L_x w nazwie oznacza dlugosc ukladu L.
    # Lista wybranych operatorow jest w tytule, a nie w legendzie.
    opis_operatorow = r",\; ".join(etykieta_operatora(n) for n in wybrane_korelacje)

    if nazwa_parametru == "omega0":
        # Przy skanowaniu omega_0, g jest stale; omega0 jest podpisane w legendzie.
        g_tytul = G_STALE if TRYB_SKANU == "omega0_g" else float(parametry["gmin"])
        opis_stalego_parametru = rf"g={g_tytul:g}"
    else:
        omega_tytul = float(parametry.get("omega0", OMEGA_0_TYTUL))
        opis_stalego_parametru = rf"\omega_0={omega_tytul:g}"

    if "fermion_lokalny" in pliki:
        x_fermionu = parametry_metod["fermion_lokalny"].get("x", str(POZYCJA_FERMIONU))
        opis_pozycji = rf",\; x_{{\mathrm{{ferm}}}}={x_fermionu}"
    else:
        opis_pozycji = ""

    tytul = (
        rf"$L={parametry['M']},\; "
        rf"N_{{\uparrow}}={parametry['N_up']},\; "
        rf"J={parametry['J']},\; "
        rf"\Delta={parametry['Delta']},\; "
        rf"{opis_stalego_parametru}"
        rf"{opis_pozycji},\; "
        rf"{opis_operatorow}$"
    )

    # Jeden wspolny tytul, ten sam rozmiar czcionki co osie i legenda.
    fig.suptitle(tytul, fontsize=FONT_SIZE, y=0.975)

    ax.set_yscale(SKALA_Y)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.tick_params(axis="both", labelsize=FONT_SIZE, pad=2)
    ax.grid(POKAZ_SIATKE, alpha=0.20, linewidth=0.7)
    ax.set_axisbelow(True)

    # -------------------------------------------------------------------------
    # LEGENDA
    # -------------------------------------------------------------------------
    handles, labels = ax.get_legend_handles_labels()

    # Style linii rozrozniaja wlaczone metody.
    for metoda in pliki:
        etykieta_metody, styl_linii = style_metod[metoda]
        handles.append(
            Line2D([0], [0], color="black", linestyle=styl_linii,
                   linewidth=LINE_WIDTH_LEGEND)
        )
        labels.append(etykieta_metody)

    ax.legend(
        handles, labels, loc="upper left", bbox_to_anchor=(1.025, 1.0),
        fontsize=FONT_SIZE, frameon=False, ncol=1, labelspacing=0.45,
        handlelength=2.0, borderaxespad=0.0
    )

    # -------------------------------------------------------------------------
    # ROZWINIECIA SPINOWE ELEMENTOW BAZY - POD WYKRESEM
    # -------------------------------------------------------------------------
    if indeksy_bazy:
        dodaj_opisy_pod_wykresem(fig, indeksy_bazy, opisy_bazy, wysokosc_rysunku, wysokosc_bloku)

    # -------------------------------------------------------------------------
    # ZAPIS - PDF Z BIALYM TLEM (PNG pozostaje opcjonalny)
    # -------------------------------------------------------------------------
    # Nazwa pliku zawiera wybrane metody, operator i zakres parametru.
    opis_korelacji = "_".join(wybrane_korelacje)
    opis_param = (
        f"{nazwa_parametru}{wartosci_do_rysowania[0]:g}-"
        f"{wartosci_do_rysowania[-1]:g}_n{len(wartosci_do_rysowania)}"
    )
    dopisek_kolumn = "_kolumny_poprawione" if POPRAW_BLEDNE_NAGLOWKI else ""
    opis_stalego_parametru_plik = opis_stalego_parametru.replace("\\omega_0", "omega0")
    opis_metod = "_".join(pliki.keys())
    dopisek_pozycji = f"_x{POZYCJA_FERMIONU}" if "fermion_lokalny" in pliki else ""
    nazwa_wyniku = (
        f"korelacje_{opis_metod}_{opis_korelacji}"
        f"_L{parametry['M']}_{opis_stalego_parametru_plik}"
        f"{dopisek_pozycji}_{opis_param}{dopisek_kolumn}"
    )

    if save_plot or save_png:
        output_folder.mkdir(parents=True, exist_ok=True)
    if save_plot:
        output_filename = output_folder / f"{nazwa_wyniku}.pdf"
        fig.savefig(
            output_filename, dpi=dpi, format="pdf", bbox_inches="tight",
            pad_inches=0.015, transparent=False, facecolor="white", edgecolor="white"
        )
        print(f"\nZapisano PDF z bialym tlem:\n{output_filename}")
    if save_png:
        output_filename_png = output_folder / f"{nazwa_wyniku}.png"
        fig.savefig(
            output_filename_png, dpi=dpi, format="png", bbox_inches="tight",
            pad_inches=0.015, transparent=True, facecolor="none", edgecolor="none"
        )
        print(f"Zapisano przezroczysty PNG:\n{output_filename_png}")

    # -------------------------------------------------------------------------
    # WYSWIETLENIE
    # -------------------------------------------------------------------------
    if show_plot:
        plt.show()
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    pliki = znajdz_pliki_metod(FOLDER_DANE, TRYB_SKANU)
    wartosci = G_WYBRANE if TRYB_SKANU == "g" else OMEGA0_WYBRANE
    rysuj_wykres(
        pliki,
        WYBRANE_KORELACJE,
        wartosci
    )
