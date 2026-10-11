from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# PARAMETRY

# Przykladowy plik z serii - numer eig na koncu nazwy bedzie zamieniany
# automatycznie na kolejne wartosci z zakresu [eig_min, eig_max].
# Format linii: omega_0 \t g \t indeks:wartosc;indeks:wartosc;...
filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g_liomsy\lioms_grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"
)

# ZAKRES LIOM-ow (oba konce WLACZNIE), np. eig3, eig4, ..., eig20.
eig_min = 3
eig_max = 6

# Numery elementow bazy, ktorych udzialy sumujemy we WSZYSTKICH tych LIOM-ach.
# selected_basis_indices = [1,1153,1164,1168]
selected_basis_indices = [4,8,200,458]

# Plik bazy DLA M=4, zgodnie z przykladowym skryptem.
basis_filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\bazy\operators_ladder_mag_M_4_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_0.5_g_0.5_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true.txt"
)

# =============================================================================
# WYSWIETLANIE / ZAPIS

show_plot = True
save_plot = True
show_operator_descriptions = True

# Ustawienia ramki pod wykresem (takie jak w skrypcie referencyjnym)
operator_description_fontsize = 14
operator_description_lines = 3

plot_title_fontsize = 13
dpi = 300
cmap = "viridis"
vmin = None
vmax = None

# Nie pokazuj pierwszego punktu siatki na kazdej osi
# (czyli np. g=0 i omega_0=0, jesli to sa pierwsze wartosci).
skip_first_g_point = True
skip_first_omega_point = True

output_folder = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_12.10.2026"
)

# =============================================================================
# PARAMETRY Z NAZWY PLIKU

nazwa_pliku = filename.name
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
    r"eig(?P<eig>\d+)\.txt"
)

match = re.match(wzorzec, nazwa_pliku)
if match is None:
    raise ValueError(f"Nie udało się odczytać parametrów z nazwy pliku:\n{nazwa_pliku}")

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
# Numer eig w zmiennej filename jest tylko wzorcem nazwy;
# faktyczny zakres przeszukiwania ustalaja eig_min i eig_max.

# =============================================================================
# FUNKCJE WCZYTUJACE DANE

def parse_coeff_string(coeff_string):
    """Zamienia '7:0.3;8:0.4' na {7: 0.3, 8: 0.4}."""
    coeffs = {}
    if not coeff_string.strip():
        return coeffs

    for entry in coeff_string.strip().split(";"):
        entry = entry.strip()
        if not entry:
            continue
        idx_str, value_str = entry.split(":")
        coeffs[int(idx_str)] = float(value_str)
    return coeffs

def filename_for_eig(numer_eig):
    """Zamienia koncowe 'eig3.txt' na 'eig<numer_eig>.txt'."""
    nowa_nazwa = re.sub(r"eig\d+\.txt$", f"eig{numer_eig}.txt", filename.name)
    return filename.with_name(nowa_nazwa)

def wczytaj_plik_liomow(sciezka, selected_indices):
    """
    Wczytuje udzialy wybranych operatorow dla wszystkich (omega_0, g).

    W plikach 'indeks:wartosc' wartosc oznacza juz |c_i|^2,
    zatem sumujemy wartosci BEZ ponownego podnoszenia do kwadratu.
    """
    wyniki = {}

    with open(sciezka, "r", encoding="utf-8") as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split(maxsplit=2)
            if len(parts) < 3:
                raise ValueError(
                    f"Zly format w {sciezka.name}, linia {line_number}: {line}"
                )

            omega_0 = float(parts[0])
            g = float(parts[1])
            coeffs = parse_coeff_string(parts[2])
            udzial = sum(coeffs.get(i, 0.0) for i in selected_indices)

            # Zaokraglenie chroni przed nieistotnymi roznicami zapisu float.
            punkt = (round(omega_0, 12), round(g, 12))
            if punkt in wyniki:
                raise ValueError(
                    f"Powtorzony punkt (omega_0, g)={punkt} "
                    f"w pliku {sciezka.name}, linia {line_number}."
                )
            wyniki[punkt] = udzial

    if not wyniki:
        raise ValueError(f"Plik nie zawiera danych: {sciezka}")

    return wyniki

def zsumuj_wszystkie_eig(eig_od, eig_do, selected_indices):
    """
    Dla kazdego (omega_0, g) liczy:
        W = sum_{k=eig_od}^{eig_do} sum_{i in selected_indices} |c_i^(k)|^2.
    Wymaga tej samej siatki w kazdym pliku, zeby nie zgubic danych.
    """
    if eig_od < 0 or eig_do < eig_od:
        raise ValueError("Wymagane: 0 <= eig_min <= eig_max.")
    if not selected_indices or len(set(selected_indices)) != len(selected_indices):
        raise ValueError("Podaj niepusta liste unikatowych indeksow operatorow.")

    pliki = [filename_for_eig(k) for k in range(eig_od, eig_do + 1)]
    nieobecne = [str(path) for path in pliki if not path.is_file()]
    if nieobecne:
        raise FileNotFoundError(
            "Nie znaleziono wymaganych plikow z zakresu eig:\n"
            + "\n".join(nieobecne)
        )

    suma = None
    punkty_referencyjne = None

    for k, path in zip(range(eig_od, eig_do + 1), pliki):
        print(f"Wczytuje eig{k}: {path.name}")
        dane = wczytaj_plik_liomow(path, selected_indices)

        if suma is None:
            suma = dane.copy()
            punkty_referencyjne = set(dane)
        else:
            punkty = set(dane)
            if punkty != punkty_referencyjne:
                brak = punkty_referencyjne - punkty
                nadmiar = punkty - punkty_referencyjne
                raise ValueError(
                    f"Niezgodna siatka w pliku {path.name}. "
                    f"Brakuje {len(brak)} punktow, "
                    f"dodatkowych punktow: {len(nadmiar)}. "
                    f"Przykladowy brak: {next(iter(brak), None)}; "
                    f"przykladowy dodatkowy: {next(iter(nadmiar), None)}."
                )

            for punkt, udzial in dane.items():
                suma[punkt] += udzial

    print(f"Wczytano {len(pliki)} plikow, {len(suma)} punktow siatki.")
    # Klucze: (omega_0, g). Na koncu funkcji zwracamy wektory
    # w kolejnosci wymaganej w obecnej czesci rysujacej.
    punkty = sorted(suma)
    omega = np.array([punkt[0] for punkt in punkty])
    g = np.array([punkt[1] for punkt in punkty])
    wartosci = np.array([suma[punkt] for punkt in punkty])
    return g, omega, wartosci

# =============================================================================
# FUNKCJE OPISUJACE BAZE - TAK JAK W PRZYKLADOWYM SKRYPCIE

def hc_code(code):
    """Kod operatora sprzezonego hermitowsko: S+ <-> S-, Sz -> Sz."""
    mapa = {"0": "0", "1": "3", "2": "2", "3": "1"}
    return "".join(mapa[x] for x in code)


def site_index(offset, leg):
    site = "i" if offset == 0 else f"i+{offset}"
    return rf"_{{{site},{leg}}}"


def code_to_latex(code, leg):
    mapa = {"1": "S^+", "2": "S^z", "3": "S^-"}
    wynik = []
    for offset, digit in enumerate(code):
        if digit != "0":
            wynik.append(mapa[digit] + site_index(offset, leg))
    return "".join(wynik) if wynik else "1"


def leg_to_latex(code, sector, leg):
    code_hc = hc_code(code)
    op = code_to_latex(code, leg)
    op_hc = code_to_latex(code_hc, leg)

    if sector == "R":
        if code == code_hc:
            return op, 1
        return rf"\left({op}+{op_hc}\right)", 1

    if sector == "I":
        if code == code_hc:
            raise ValueError(
                f"Operator {code} jest samosprzężony, ale ma sektor I."
            )
        return rf"\left({op}-{op_hc}\right)", 1j

    raise ValueError(f"Nieznany sektor: {sector}")


def operator_to_latex(sector, code):
    upper_code, lower_code = code.split("|")
    upper, pref_upper = leg_to_latex(upper_code, sector[0], 1)
    lower, pref_lower = leg_to_latex(lower_code, sector[1], 2)

    pref = pref_upper * pref_lower
    factors = []
    if upper != "1":
        factors.append(upper)
    if lower != "1":
        factors.append(lower)

    body = "".join(factors) if factors else "1"
    prefix = {1: "", 1j: "i", -1: "-", -1j: "-i"}[pref]
    return prefix + body


def read_basis(path):
    """Indeksy odpowiadaja kolejnosci poprawnych rekordow bazy (od 1)."""
    if not path.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku bazy:\n{path}")

    baza = {}
    indeks = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            sector, code = parts[1], parts[2]
            if len(sector) != 2 or any(x not in "RI" for x in sector) or "|" not in code:
                continue

            indeks += 1
            baza[indeks] = (sector, code)
    return baza

def build_operator_descriptions(path, selected):
    baza = read_basis(path)
    entries = []

    print("\n============================================")
    print("OPERATORY Z BAZY")
    print("============================================")

    for indeks in selected:
        if indeks not in baza:
            raise ValueError(
                f"Operator {indeks} nie istnieje w bazie. "
                f"Baza ma {len(baza)} operatorów."
            )
        sector, code = baza[indeks]
        opis = operator_to_latex(sector, code)
        entries.append((indeks, opis))
        print(f"{indeks}\t{sector}\t{code}")

    return entries

def rows_of_operator_descriptions(entries, requested_lines):
    """Ten sam sposob rownowazenia wierszy co w dostarczonym przykladzie."""
    opisy = [rf"$\mathbf{{{indeks}:}}\;{opis}$" for indeks, opis in entries]
    if not opisy:
        return []

    line_count = max(1, requested_lines)
    if len(opisy) <= line_count:
        linie = [[x] for x in opisy]
    elif line_count == 3:
        best = None
        best_score = float("inf")
        n = len(opisy)

        for i in range(1, n - 1):
            for j in range(i + 1, n):
                candidate = [opisy[:i], opisy[i:j], opisy[j:]]
                lengths = [sum(len(x) for x in row) for row in candidate]
                mean = sum(lengths) / 3
                score = sum((x - mean) ** 2 for x in lengths)
                if score < best_score:
                    best_score = score
                    best = candidate
        linie = best
    else:
        linie = [[] for _ in range(line_count)]
        lengths = [0] * line_count
        for opis in opisy:
            i = int(np.argmin(lengths))
            linie[i].append(opis)
            lengths[i] += len(opis)

    return ["     ".join(row) for row in linie if row]

# =============================================================================
# WCZYTANIE DANYCH I BAZY

print(f"Sumuje udzialy operatorow {selected_basis_indices} w eig{eig_min}-eig{eig_max}:")
g, omega, wartosci = zsumuj_wszystkie_eig(eig_min, eig_max, selected_basis_indices)

operator_descriptions = (
    build_operator_descriptions(basis_filename, selected_basis_indices)
    if show_operator_descriptions
    else []
)

# =============================================================================
# SIATKA 2D

g_values_all = np.unique(g)
omega_values_all = np.unique(omega)

if skip_first_g_point and len(g_values_all) < 2:
    raise ValueError("Nie mozna pominac pierwszego punktu osi g: siatka ma mniej niz 2 wartosci g.")
if skip_first_omega_point and len(omega_values_all) < 2:
    raise ValueError("Nie mozna pominac pierwszego punktu osi omega_0: siatka ma mniej niz 2 wartosci omega_0.")

g_values = g_values_all[1:] if skip_first_g_point else g_values_all
omega_values = omega_values_all[1:] if skip_first_omega_point else omega_values_all

mask = np.isin(g, g_values) & np.isin(omega, omega_values)
g_plot = g[mask]
omega_plot = omega[mask]
wartosci_plot = wartosci[mask]

heatmap = np.full((len(omega_values), len(g_values)), np.nan)

for g_i, omega_i, value_i in zip(g_plot, omega_plot, wartosci_plot):
    indeks_g = np.where(np.isclose(g_values, g_i))[0][0]
    indeks_omega = np.where(np.isclose(omega_values, omega_i))[0][0]
    heatmap[indeks_omega, indeks_g] = value_i

# =============================================================================
# TYTUL

selected_indices_str = ", ".join(str(x) for x in selected_basis_indices)

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
    rf"$\mathrm{{eig}}={eig_min}\ldots {eig_max}$"
    "\n"
    rf"selected basis indices: {selected_indices_str}"
)

# =============================================================================
# RYSOWANIE HEATMAPY BEZ LINII STALYCH

fig, ax = plt.subplots(figsize=(10, 8))
mesh = ax.pcolormesh(
    g_values,
    omega_values,
    heatmap,
    shading="auto",
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
)
cbar = fig.colorbar(mesh, ax=ax)

ax.set_xlabel(r"$g$", fontsize=14)
ax.set_ylabel(r"$\omega_0$", fontsize=14)
ax.set_title(title, fontsize=plot_title_fontsize, pad=15)
ax.tick_params(axis="both", labelsize=12)

# Kwadratowe pole heatmapy
ax.set_box_aspect(1)

# =============================================================================
# ZAOKRAGLONA RAMKA Z WZORAMI OPERATOROW POD WYKRESEM

if show_operator_descriptions and operator_descriptions:
    linie = rows_of_operator_descriptions(
        operator_descriptions, operator_description_lines
    )
    tekst = "\n".join(linie)

    fig.text(
        0.5,
        0.018,
        tekst,
        fontsize=operator_description_fontsize,
        horizontalalignment="center",
        verticalalignment="bottom",
        multialignment="center",
        linespacing=1.5,
        bbox=dict(
            boxstyle="round,pad=0.8",
            facecolor="white",
            edgecolor="black",
            linewidth=1.5,
        ),
    )
    plt.subplots_adjust(left=0.12, right=0.88, top=0.79, bottom=0.34)
else:
    plt.subplots_adjust(left=0.12, right=0.88, top=0.79, bottom=0.12)

# =============================================================================
# ZAPIS

output_folder.mkdir(parents=True, exist_ok=True)
selected_tag = "-".join(str(x) for x in selected_basis_indices)
output_filename = output_folder / (
    f"heatmap_selected_basis_sum_"
    f"M{M}_Jp{Jp}_d{Delta}_d2{Delta2}_"
    f"T{T}_P{P}_F{F}_B{B}_FId{FId}_"
    f"eig{eig_min}-{eig_max}_basis_{selected_tag}.png"
)

if save_plot:
    plt.savefig(output_filename, dpi=dpi, bbox_inches="tight")
    print("\nZapisano wykres:")
    print(output_filename)

# =============================================================================
# WYSWIETLENIE

if show_plot:
    plt.show()
else:
    plt.close()
