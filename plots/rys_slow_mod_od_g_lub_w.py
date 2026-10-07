from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USTAWIENIA

filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g_liomsy\lioms_grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"
)

basis_filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\bazy\operators_ladder_mag_M_4_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_0.5_g_0.5_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true.txt"
)

plot_mode = "fixed_g"  # "fixed_omega" lub "fixed_g"

# omega_fixed = 1.5172413793103448
# g_fixed = 1.5172413793103448
# omega_fixed = 0.4827586206896552
# g_fixed = 0.4827586206896552
omega_fixed = 0.62068965517241381
g_fixed = 0.75862068965517238
value_tolerance = 1e-10

wybrane_operatory = [3, 4, 8, 194, 197, 458, 449, 1154, 1156, 1162]

show_plot = True
save_plot = True
show_operator_descriptions = True

x_min, x_max = 0.0, 2.0
y_min, y_max = 0.0, 1.0

plot_title_fontsize = 16
operator_description_fontsize = 14
operator_description_lines = 3
dpi = 300

output_directory = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_05.10.2026"
)

# =============================================================================
# FUNKCJE OPISUJĄCE BAZĘ

def hc_code(code):
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
            raise ValueError(f"Operator {code} jest samosprzężony, ale ma sektor I.")
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

            sector = parts[1]
            code = parts[2]

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
            raise ValueError(f"Operator {indeks} nie istnieje w bazie. Baza ma {len(baza)} operatorów.")

        sector, code = baza[indeks]
        opis = operator_to_latex(sector, code)

        entries.append((indeks, opis))
        print(f"{indeks}\t{sector}\t{code}")

    return entries

# =============================================================================
# KONTROLA

if plot_mode not in ["fixed_omega", "fixed_g"]:
    raise ValueError('plot_mode musi być "fixed_omega" albo "fixed_g".')

if not filename.exists():
    raise FileNotFoundError(filename)

if not basis_filename.exists():
    raise FileNotFoundError(basis_filename)

operator_descriptions = build_operator_descriptions(basis_filename, wybrane_operatory)

# =============================================================================
# PARAMETRY Z NAZWY PLIKU

match = re.search(
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

if match is None:
    raise ValueError(f"Nie udało się odczytać parametrów z {filename.name}")

p = match.groupdict()

M, Jp = p["M"], p["Jp"]
Delta, Delta2 = p["Delta"], p["Delta2"]
T, P, F, B, FId = p["T"], p["P"], p["F"], p["B"], p["FId"]
eig = p["eig"]

# =============================================================================
# WCZYTANIE DANYCH

omega_all = []
g_all = []
wspolczynniki_all = {indeks: [] for indeks in wybrane_operatory}

pattern = re.compile(
    r"(\d+)\s*:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
)

with open(filename, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        try:
            omega = float(parts[0])
            g = float(parts[1])
        except ValueError:
            continue

        znalezione = {int(i): float(v) for i, v in pattern.findall(line)}

        omega_all.append(omega)
        g_all.append(g)

        for indeks in wybrane_operatory:
            wspolczynniki_all[indeks].append(znalezione.get(indeks, 0.0))

omega_all = np.array(omega_all)
g_all = np.array(g_all)

for indeks in wybrane_operatory:
    wspolczynniki_all[indeks] = np.array(wspolczynniki_all[indeks])

if len(g_all) == 0:
    raise ValueError("Nie wczytano danych.")

# =============================================================================
# WYBÓR PRZEKROJU

if plot_mode == "fixed_omega":
    mask = np.isclose(omega_all, omega_fixed, rtol=0.0, atol=value_tolerance)

    if not np.any(mask):
        nearest = np.unique(omega_all)
        nearest = nearest[np.argmin(abs(nearest - omega_fixed))]
        raise ValueError(f"Nie znaleziono omega_0 = {omega_fixed}\nNajbliższe: {nearest}")

    omega = omega_all[mask]
    g = g_all[mask]

    fixed_value = omega[0]
    x_values = g
    x_label = r"$g$"

else:
    mask = np.isclose(g_all, g_fixed, rtol=0.0, atol=value_tolerance)

    if not np.any(mask):
        nearest = np.unique(g_all)
        nearest = nearest[np.argmin(abs(nearest - g_fixed))]
        raise ValueError(f"Nie znaleziono g = {g_fixed}\nNajbliższe: {nearest}")

    omega = omega_all[mask]
    g = g_all[mask]

    fixed_value = g[0]
    x_values = omega
    x_label = r"$\omega_0$"

wspolczynniki = {
    indeks: wspolczynniki_all[indeks][mask]
    for indeks in wybrane_operatory
}

kolejnosc = np.argsort(x_values)

x_values = x_values[kolejnosc]
omega = omega[kolejnosc]
g = g[kolejnosc]

for indeks in wybrane_operatory:
    wspolczynniki[indeks] = wspolczynniki[indeks][kolejnosc]

# =============================================================================
# INFORMACJE

print("\n============================================")
print("WCZYTANE DANE")
print("============================================")

if plot_mode == "fixed_omega":
    print(f"omega_0 = {fixed_value:.16g}")
    print(f"g = {x_values.min()} -> {x_values.max()}")
else:
    print(f"g = {fixed_value:.16g}")
    print(f"omega_0 = {x_values.min()} -> {x_values.max()}")

for indeks in wybrane_operatory:
    print(f"Operator {indeks}: {np.count_nonzero(wspolczynniki[indeks])} niezerowych wartości")

# =============================================================================
# NAZWA PLIKU

if plot_mode == "fixed_omega":
    plot_name = "liom_i_vs_g"
    fixed_tag = f"_omega0_{fixed_value:.3g}"
else:
    plot_name = "liom_i_vs_omega0"
    fixed_tag = f"_g_{fixed_value:.3g}"

output_filename = output_directory / (
    f"{plot_name}_"
    f"M{M}_Jp{Jp}_d{Delta}_d2{Delta2}_"
    f"T{T}_P{P}_F{F}_B{B}_FId{FId}_"
    f"eig{eig}{fixed_tag}.png"
)

# =============================================================================
# WYKRES

fig, ax = plt.subplots(figsize=(10, 6))

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

if plot_mode == "fixed_omega":
    tytul += rf";\ \omega_0 = {fixed_value:.3g}"
else:
    tytul += rf";\ g = {fixed_value:.3g}"

tytul += "$"

ax.set_title(tytul, fontsize=plot_title_fontsize, pad=15)

for indeks in wybrane_operatory:
    ax.plot(
        x_values,
        wspolczynniki[indeks],
        marker="o",
        markersize=2,
        linewidth=1.0,
        label=str(indeks)
    )

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_xlabel(x_label, fontsize=14)
ax.set_ylabel(r"$|c_i|^2$", fontsize=14)
ax.tick_params(axis="both", labelsize=12)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)

# =============================================================================
# OPIS OPERATORÓW POD WYKRESEM

if show_operator_descriptions:
    opisy = [
        rf"$\mathbf{{{indeks}:}}\;{opis}$"
        for indeks, opis in operator_descriptions
    ]

    if len(opisy) <= operator_description_lines:
        linie = [[x] for x in opisy]

    elif operator_description_lines == 3:
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
        linie = [[] for _ in range(operator_description_lines)]
        lengths = [0] * operator_description_lines

        for opis in opisy:
            i = int(np.argmin(lengths))
            linie[i].append(opis)
            lengths[i] += len(opis)

    tekst = "\n".join("     ".join(row) for row in linie if row)

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
            linewidth=1.5
        )
    )

    plt.subplots_adjust(bottom=0.34)

else:
    plt.subplots_adjust(bottom=0.12)

# =============================================================================
# ZAPIS / WYŚWIETLENIE

if save_plot:
    output_directory.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_filename, dpi=dpi, bbox_inches="tight")
    print(f"\nZapisano wykres:\n{output_filename}")

if show_plot:
    plt.show()
else:
    plt.close()