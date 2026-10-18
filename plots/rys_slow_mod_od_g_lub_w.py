
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA

# Plik wzorcowy: numer eig w nazwie jest automatycznie podmieniany.
filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\liom_mag_dane\siatka_omega_g_liomsy\lioms_grid_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig19.txt"
)

basis_filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\liom_mag_dane\bazy\operators_ladder_mag_M_3_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_1.0_g_1.0_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true.txt"
)

# PRZELACZNIK:
# "single" - wspolczynniki z jednego pliku eig
# "sum"    - suma wspolczynnikow tego samego operatora z eig_min ... eig_max
eig_mode = "sum"
selected_eig = 19

eig_min = 3
eig_max = 19

# "fixed_omega": omega_0 ustalone, na osi X jest g
# "fixed_g":     g ustalone, na osi X jest omega_0
plot_mode = "fixed_omega"
omega_fixed = 0.5
g_fixed = 0.25
value_tolerance = 1e-10

# Kazdy operator jest rysowany jako osobna linia.
# wybrane_operatory = [3,4,7,8,12,13,16,19,20,25]
#wybrane_operatory = [26,32,33,34,36,37,38,44,48,58]
#wybrane_operatory = [60,64,115,120,124,138,142,144,148,156]
#wybrane_operatory = [225,289,292,300,304,1121]

#wybrane_operatory = [4,7,13,16,19,20,25,33,34,37,38,225,1057,1121]
#wybrane_operatory = [3,12,32,36,44,48,58,60,64,115,120,124,138,142,144,148,156,292,300,304]
#wybrane_operatory = [8,26]
#wybrane_operatory = [49,1,52,289]

wybrane_operatory = [8,26]


show_plot = True
save_plot = True
show_operator_descriptions = True

x_min, x_max = 0.0, 10.0
y_min, y_max = 0.0, None  # None = automatyczny zakres

plot_title_fontsize = 16
operator_description_fontsize = 14
operator_description_lines = 3
dpi = 300

output_directory = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_19.10.2026\M=3"
)


# =============================================================================
# FUNKCJE OPISUJACE BAZE

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
            raise ValueError(
                f"Operator {code} jest samosprzezony, ale ma sektor I."
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

            if (
                len(sector) != 2
                or any(x not in "RI" for x in sector)
                or "|" not in code
            ):
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
                f"Baza ma {len(baza)} operatorow."
            )

        sector, code = baza[indeks]
        entries.append((indeks, operator_to_latex(sector, code)))

        print(f"{indeks}\t{sector}\t{code}")

    return entries


# =============================================================================
# ODCZYT PLIKOW eig

pattern = re.compile(
    r"(\d+)\s*:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
)


def path_for_eig(eig):
    """Podmiana numeru eig w nazwie pliku."""
    new_name, replacements = re.subn(
        r"_eig\d+\.txt$",
        f"_eig{eig}.txt",
        filename.name
    )

    if replacements != 1:
        raise ValueError(
            f"Brak koncowki _eigNN.txt w nazwie: {filename.name}"
        )

    return filename.with_name(new_name)


def read_eig(path, selected):
    """
    Zwraca:
    {(omega_0, g): {indeks_operatora: |c_i|^2}}
    """

    if not path.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku:\n{path}"
        )

    data = {}

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()

            if len(parts) < 2:
                continue

            try:
                omega, g = float(parts[0]), float(parts[1])
            except ValueError:
                continue

            key = (omega, g)

            if key in data:
                raise ValueError(
                    f"Powtorzony punkt (omega_0, g) = {key} "
                    f"w {path.name}"
                )

            found = {
                int(i): float(v)
                for i, v in pattern.findall(line)
            }

            # Brak operatora oznacza wspolczynnik rowny 0.
            data[key] = {
                i: found.get(i, 0.0)
                for i in selected
            }

    if not data:
        raise ValueError(
            f"Nie wczytano zadnych punktow z pliku: {path}"
        )

    return data


def load_selected_eigs(selected):

    if eig_mode == "single":
        eig_list = [selected_eig]
        eig_description = f"eig {selected_eig}"
        file_tag = f"eig{selected_eig}"

    elif eig_mode == "sum":
        if eig_min > eig_max:
            raise ValueError("eig_min musi byc <= eig_max")

        eig_list = list(range(eig_min, eig_max + 1))
        eig_description = f"sum eig {eig_min}-{eig_max}"
        file_tag = f"sum_eig{eig_min}-{eig_max}"

    else:
        raise ValueError(
            'eig_mode musi byc "single" albo "sum".'
        )

    totals = None
    reference_points = None

    for eig in eig_list:
        path = path_for_eig(eig)

        print(f"Wczytuje eig = {eig}: {path}")

        current = read_eig(path, selected)
        points = set(current)

        if totals is None:
            reference_points = points
            totals = {
                point: {i: 0.0 for i in selected}
                for point in points
            }

        elif points != reference_points:
            missing = reference_points - points
            extra = points - reference_points

            raise ValueError(
                f"Siatka (omega_0, g) w eig={eig} jest inna. "
                f"Brakuje {len(missing)} punktow, "
                f"dodatkowych {len(extra)}."
            )

        # SUMOWANIE WSPOLCZYNNIKOW
        for point, values in current.items():
            for i in selected:
                totals[point][i] += values[i]

    return totals, eig_description, file_tag


# =============================================================================
# RYSOWANIE

def main():

    if plot_mode not in ("fixed_omega", "fixed_g"):
        raise ValueError(
            'plot_mode musi byc "fixed_omega" albo "fixed_g".'
        )

    if not wybrane_operatory:
        raise ValueError(
            "Podaj przynajmniej jeden operator."
        )

    if len(set(wybrane_operatory)) != len(wybrane_operatory):
        raise ValueError(
            "Lista wybrane_operatory zawiera powtarzajace sie indeksy."
        )

    if not basis_filename.exists():
        raise FileNotFoundError(basis_filename)

    descriptions = build_operator_descriptions(
        basis_filename,
        wybrane_operatory
    )

    # =========================================================================
    # PARAMETRY Z NAZWY PLIKU

    match = re.fullmatch(
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
        r"eig\d+\.txt",
        filename.name,
    )

    if match is None:
        raise ValueError(
            f"Nie mozna odczytac parametrow z nazwy: {filename.name}"
        )

    params = match.groupdict()

    # =========================================================================
    # WCZYTANIE DANYCH

    totals, eig_description, file_tag = load_selected_eigs(
        wybrane_operatory
    )

    points = sorted(totals)

    omega_all = np.array(
        [point[0] for point in points],
        dtype=float
    )

    g_all = np.array(
        [point[1] for point in points],
        dtype=float
    )

    # =========================================================================
    # WYBOR PRZEKROJU

    if plot_mode == "fixed_omega":

        mask = np.isclose(
            omega_all,
            omega_fixed,
            rtol=0.0,
            atol=value_tolerance
        )

        if not np.any(mask):
            available = np.unique(omega_all)
            nearest = available[
                np.argmin(np.abs(available - omega_fixed))
            ]

            raise ValueError(
                f"Nie ma omega_0 = {omega_fixed}; "
                f"najblizsze: {nearest}"
            )

        fixed_value = omega_all[mask][0]
        x_values = g_all[mask]

        x_label = r"$g$"
        fixed_tag = f"_omega0_{fixed_value:.3g}"
        plot_name = "liom_i_vs_g"
        fixed_title = rf"\omega_0 = {fixed_value:.3g}"

    else:

        mask = np.isclose(
            g_all,
            g_fixed,
            rtol=0.0,
            atol=value_tolerance
        )

        if not np.any(mask):
            available = np.unique(g_all)
            nearest = available[
                np.argmin(np.abs(available - g_fixed))
            ]

            raise ValueError(
                f"Nie ma g = {g_fixed}; "
                f"najblizsze: {nearest}"
            )

        fixed_value = g_all[mask][0]
        x_values = omega_all[mask]

        x_label = r"$\omega_0$"
        fixed_tag = f"_g_{fixed_value:.3g}"
        plot_name = "liom_i_vs_omega0"
        fixed_title = rf"g = {fixed_value:.3g}"

    # =========================================================================
    # SORTOWANIE PUNKTOW

    order = np.argsort(x_values)
    x_values = x_values[order]

    selected_points = [
        point
        for point, keep in zip(points, mask)
        if keep
    ]

    selected_points = [
        selected_points[j]
        for j in order
    ]

    # =========================================================================
    # INFORMACJE

    print("\n============================================")
    print(f"TRYB: {eig_description}")
    print(f"PRZEKROJ: {fixed_title}")
    print("============================================")

    for indeks in wybrane_operatory:

        y = np.array([
            totals[point][indeks]
            for point in selected_points
        ])

        print(
            f"Operator {indeks}: "
            f"{np.count_nonzero(y)} niezerowych punktow; "
            f"max = {y.max():.6g}"
        )

    # =========================================================================
    # WYKRES

    fig, ax = plt.subplots(figsize=(10, 6))

    p = params

    title = (
        rf"$M={p['M']};\ J'={p['Jp']};\ "
        rf"\Delta={p['Delta']};\ \Delta'={p['Delta2']}$"
        "\n"
        rf"$T={p['T']};\ P={p['P']};\ F={p['F']};\ "
        rf"B={p['B']};\ FId={p['FId']};\ {fixed_title}$"
        f"\n{eig_description}"
    )

    ax.set_title(
        title,
        fontsize=plot_title_fontsize,
        pad=15
    )

    # =========================================================================
    # KRZYWE OPERATOROW

    for indeks in wybrane_operatory:

        y = np.array([
            totals[point][indeks]
            for point in selected_points
        ])

        ax.plot(
            x_values,
            y,
            marker="o",
            markersize=2,
            linewidth=1.0,
            label=str(indeks)
        )

    # =========================================================================
    # OSIE

    if x_min is not None or x_max is not None:
        ax.set_xlim(
            left=x_min,
            right=x_max
        )

    if y_min is not None or y_max is not None:
        ax.set_ylim(
            bottom=y_min,
            top=y_max
        )

    ax.set_xlabel(
        x_label,
        fontsize=14
    )

    if eig_mode == "single":
        ax.set_ylabel(
            r"$|c_i|^2$",
            fontsize=14
        )

    else:
        ax.set_ylabel(
            rf"$\sum_{{k={eig_min}}}^{{{eig_max}}}|c_{{i,k}}|^2$",
            fontsize=14
        )

    ax.tick_params(
        axis="both",
        labelsize=12
    )

    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)

    # =========================================================================
    # OPIS OPERATOROW POD WYKRESEM

    if show_operator_descriptions:

        opisy = [
            rf"$\mathbf{{{i}:}}\;{latex}$"
            for i, latex in descriptions
        ]

        number_of_lines = max(
            1,
            operator_description_lines
        )

        if len(opisy) <= number_of_lines:
            linie = [[x] for x in opisy]

        elif number_of_lines == 3:

            best = None
            best_score = float("inf")
            n = len(opisy)

            for i in range(1, n - 1):
                for j in range(i + 1, n):

                    candidate = [
                        opisy[:i],
                        opisy[i:j],
                        opisy[j:]
                    ]

                    lengths = [
                        sum(len(x) for x in row)
                        for row in candidate
                    ]

                    mean = sum(lengths) / 3

                    score = sum(
                        (x - mean) ** 2
                        for x in lengths
                    )

                    if score < best_score:
                        best_score = score
                        best = candidate

            linie = best

        else:

            linie = [
                [] for _ in range(number_of_lines)
            ]

            lengths = [0] * number_of_lines

            for opis in opisy:

                j = int(np.argmin(lengths))

                linie[j].append(opis)
                lengths[j] += len(opis)

        tekst = "\n".join(
            "     ".join(row)
            for row in linie
            if row
        )

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

    # =========================================================================
    # NAZWA PLIKU

    output_filename = output_directory / (
        f"{plot_name}_"
        f"M{p['M']}_Jp{p['Jp']}_d{p['Delta']}_d2{p['Delta2']}_"
        f"T{p['T']}_P{p['P']}_F{p['F']}_B{p['B']}_FId{p['FId']}_"
        f"{file_tag}{fixed_tag}.png"
    )

    # =========================================================================
    # ZAPIS / WYSWIETLENIE

    if save_plot:

        output_directory.mkdir(
            parents=True,
            exist_ok=True
        )

        fig.savefig(
            output_filename,
            dpi=dpi,
            bbox_inches="tight"
        )

        print(
            f"\nZapisano wykres:\n{output_filename}"
        )

    if show_plot:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
