from dataclasses import dataclass
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize


# =============================================================================
# USTAWIENIA
# =============================================================================

input_filename = Path(r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\se_2_13.08.2026\siatka_omega_g_liomsy\lioms_grid_M3_Jp0.0_d0.8_Tboth_Pboth_Fyes_Bboth_FIdno_eig3.txt")

output_directory = Path(r"C:\Users\aleks\Desktop\praca magisterska\M=3")

save_plot = False
show_plot = True

# "g" -> y = g
# "omega" -> y = omega_0
y_axis_parameter = "g"
fixed_parameter_value = None

# Zakresy osi X:
# None -> cała baza
# Można podać kilka rozłącznych zakresów, np. [(1, 100), (200, 300)]
x_ranges = [(1, 1000)]

# Zakres osi Y; None -> automatycznie
y_min, y_max = None, None

# None -> rozmiar bazy wyznaczany automatycznie
basis_size = None

# Skala kolorów
use_log_values = False

# Stały zakres skali kolorów
color_min, color_max = 0.0, 1.0

# Wartości poniżej progu będą białe
white_below_threshold = True
threshold_white = 1e-2

# TYLKO te indeksy dostaną tick i podpis na osi X
labeled_basis_indices = [1, 100, 3, 7, 8, 12, 50]

# Ręczny opis pod wykresem
show_manual_legend = True
manual_legend_title = "Basis elements:"
manual_legend_fontsize = 14
plot_title_fontsize = 16
dpi = 300

manual_legend_entries = [
    #(1,  r"RR$\;1_{i}1_{i+1}\;|\;S^{z}_{i}1_{i+1}\quad$"),
    #(3,  r"RR$\;1_{i}1_{i+1}\;|\;S^{-}_{i}1_{i+1}\quad$"),
    #(4,  r"RR$\;1_{i}1_{i+1}\;|\;S^{z}_{i}S^{z}_{i+1}\quad$"),
    #(7,  r"RI$\;1_{i}1_{i+1}\;|\;S^{+}_{i}S^{-}_{i+1}\quad$"),
    #(8,  r"RR$\;1_{i}1_{i+1}\;|\;S^{+}_{i}S^{-}_{i+1}\quad$"),
    #(12, r"RR$\;S_{i}1_{i+1}\;|\;S^{-}_{i}S^{-}_{i+1}\quad$"),
    #(14, r"RR$\;S^{z}_{i}1_{i+1}\;|\;S^{z}_{i}1_{i+1}\quad$"),
    #(41, r"RR$\;S^{z}_{i}S^{z}_{i+1}\;|\;1^{z}_{i}1^{z}_{i+1}\quad$"),
    #(46, r"RR$\;S^{z}_{i}S^{z}_{i+1}\;|\;S^{z}_{i}S^{z}_{i+1}\quad$"),
    #(73, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;1_{i}1_{i+1}\quad$"),
    #(84, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;S^{+}_{i}S^{-}_{i+1}\quad$"),

    ###########################################################################

    #(2,  r"RR$\;S^{z}_{i}1_{i+1}\;|\;S^{z}_{i}1_{i+1}\quad$"),
    #(4,  r"RR$\;S^{z}_{i}1_{i+1}\;|\;S^{-}_{i}1_{i+1}\quad$"),
    #(5,  r"RR$\;S^{z}_{i}1_{i+1}\;|\;1_{i}S^{z}_{i+1}\quad$"),
    #(12, r"RR$\;S^{z}_{i}1_{i+1}\;|\;S^{-}_{i}S^{+}_{i+1}\quad$"),
    #(17, r"RR$\;1_{i}S^{z}_{i+1}\;|\;S^{z}_{i}1_{i+1}\quad$"),
    #(24, r"RR$\;1_{i}S^{z}_{i+1}\;|\;S^{-}_{i}S^{+}_{i+1}\quad$"),
    #(29, r"RR$\;S^{z}_{i}S^{z}_{i+1}\;|\;1_{i}1_{i+1}\quad$"),
    #(45, r"IR$\;S^{+}_{i}S^{-}_{i+1}\;|\;1_{i}1_{i+1}\quad$"),
    #(61, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;1_{i}1_{i+1}\quad$"),
    #(62, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;S^{z}_{i}1_{i+1}\quad$"),
    #(64, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;S^{-}_{i}1_{i+1}\quad$"),
    #(65, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;1_{i}S^{z}_{i+1}\quad$"),
    #(70, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;1_{i}S^{-}_{i+1}\quad$"),
    #(72, r"RR$\;S^{+}_{i}S^{-}_{i+1}\;|\;S^{+}_{i}S^{-}_{i+1}\quad$"),

    ###########################################################################

    #(3, r"$\;S^{-}_{i,2}1_{i+1,2}+S^{+}_{i,2}1_{i+1,2}\quad$"),
    #(8, r"$\;S^{+}_{i,2}S^{-}_{i+1,2}+S^{-}_{i,2}S^{+}_{i+1,2}\quad$"),

    ###########################################################################

    #(4,  r"$\;S^{z}_{i,1}S^{-}_{i,2}+S^{z}_{i,1}S^{+}_{i,2}\quad$"),
    #(29, r"$\;S^{z}_{i,1}S^{z}_{i+1,1}\quad$"),
    #(61, r"$\;S^{+}_{i,1}S^{-}_{i+1,1}+S^{-}_{i,1}S^{+}_{i+1,1}\quad$"),


    ###########################################################################
    # M=3, FId=no
    ###########################################################################
    #(4,   r"$\;S^{z}_{i,1}\left(S^{-}_{i,2}+S^{+}_{i,2}\right)\quad$"),
    #(113, r"$\;S^{z}_{i,1}S^{z}_{i+1,1}\quad$"),
    #(241, r"$\;S^{+}_{i,1}S^{-}_{i+1,1}+S^{-}_{i,1}S^{+}_{i+1,1}\quad$"),

    # (12, r"$\;S^z_{i,1}\left(S^+_{i,2}S^-_{i+1,2}+S^-_{i,2}S^+_{i+1,2}\right)\quad$"),
    # (72, r"$\;S^z_{i+1,1}\left(S^+_{i+1,2}S^-_{i+2,2}+S^-_{i+1,2}S^+_{i+2,2}\right)\quad$"),
    # (244, r"$\;\left(S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}\right)\left(S^-_{i,2}+S^+_{i,2}\right)\quad$"),
    # (250, r"$\;\left(S^+_{i,1}S^-_{i+1,1}+S^-_{i,1}S^+_{i+1,1}\right)\left(S^-_{i+1,2}+S^+_{i+1,2}\right)\quad$"),

    ###########################################################################
    # M=3, FId=yes
    ###########################################################################
    (3, r"$\;S^-_{i,2}+S^+_{i,2}\quad$"),
    (7, r"$\;i\left(S^+_{i,2}S^-_{i+1,2}-S^-_{i,2}S^+_{i+1,2}\right)\quad$"),
    (8, r"$\;S^+_{i,2}S^-_{i+1,2}+S^-_{i,2}S^+_{i+1,2}\quad$"),
    (12, r"$\;S^-_{i,2}S^-_{i+1,2}+S^+_{i,2}S^+_{i+1,2}\quad$"),
    (50, r"$\;S^z_{i,1}S^z_{i,2}\quad$"),



]


# =============================================================================
# STRUKTURA DANYCH
# =============================================================================

@dataclass
class GridRow:
    omega_0: float
    g: float
    coefficients: dict[int, float]


# =============================================================================
# WCZYTYWANIE
# =============================================================================

def parse_sparse_coefficients(text: str, *, line_number: int) -> dict[int, float]:
    coefficients = {}

    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue

        try:
            index_text, value_text = item.split(":", maxsplit=1)
            basis_index, squared_coefficient = int(index_text), float(value_text)
        except ValueError as error:
            raise ValueError(f"Niepoprawna para indeks:wartość w linii {line_number}: {item!r}") from error

        if basis_index < 1:
            raise ValueError(f"Indeks bazy musi być dodatni, linia {line_number}: {basis_index}")

        if basis_index in coefficients:
            raise ValueError(f"Indeks {basis_index} występuje dwa razy w linii {line_number}.")

        coefficients[basis_index] = squared_coefficient

    return coefficients


def load_grid_file(filename: Path) -> list[GridRow]:
    if not filename.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku:\n{filename}")

    rows = []

    with filename.open("r", encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = raw_line.rstrip("\n\r").split("\t", maxsplit=2)

            if len(parts) != 3:
                raise ValueError(f"Niepoprawny format linii {line_number}:\n{raw_line}")

            try:
                omega_0, g = float(parts[0]), float(parts[1])
            except ValueError as error:
                raise ValueError(f"Nie udało się wczytać omega_0 lub g w linii {line_number}:\n{raw_line}") from error

            coefficients = parse_sparse_coefficients(parts[2], line_number=line_number)
            rows.append(GridRow(omega_0, g, coefficients))

    if not rows:
        raise ValueError("Plik nie zawiera żadnych danych.")

    return rows


# =============================================================================
# BUDOWANIE MACIERZY
# =============================================================================

def unique_sorted(values: list[float]) -> np.ndarray:
    return np.array(sorted(set(values)), dtype=float)


def select_rows(rows: list[GridRow]) -> tuple[list[GridRow], str, str, float]:
    if y_axis_parameter == "g":
        fixed_name, y_name = "omega_0", "g"
        fixed_values = unique_sorted([row.omega_0 for row in rows])
        get_fixed, get_y = lambda row: row.omega_0, lambda row: row.g

    elif y_axis_parameter == "omega":
        fixed_name, y_name = "g", "omega_0"
        fixed_values = unique_sorted([row.g for row in rows])
        get_fixed, get_y = lambda row: row.g, lambda row: row.omega_0

    else:
        raise ValueError('y_axis_parameter musi mieć wartość "g" albo "omega".')

    if fixed_parameter_value is None:
        if len(fixed_values) != 1:
            available = ", ".join(f"{value:.16g}" for value in fixed_values)
            raise ValueError(f"W pliku znaleziono kilka wartości {fixed_name}:\n{available}\nUstaw fixed_parameter_value.")

        selected_fixed_value = float(fixed_values[0])

    else:
        selected_fixed_value = float(fixed_parameter_value)

        if not np.any(np.isclose(fixed_values, selected_fixed_value, rtol=1e-10, atol=1e-12)):
            available = ", ".join(f"{value:.16g}" for value in fixed_values)
            raise ValueError(f"Nie znaleziono {fixed_name} = {selected_fixed_value:.16g}.\nDostępne: {available}")

    selected_rows = [row for row in rows if np.isclose(get_fixed(row), selected_fixed_value, rtol=1e-10, atol=1e-12)]
    selected_rows.sort(key=get_y)

    y_values = [get_y(row) for row in selected_rows]

    if len(y_values) != len(set(y_values)):
        raise ValueError(f"Występują powtórzone wartości parametru {y_name}.")

    return selected_rows, y_name, fixed_name, selected_fixed_value


def build_heatmap_matrix(rows: list[GridRow], y_name: str) -> tuple[np.ndarray, np.ndarray, int]:
    if not rows:
        raise ValueError("Brak danych do utworzenia heat mapy.")

    largest_index = max((max(row.coefficients, default=0) for row in rows), default=0)
    final_basis_size = largest_index if basis_size is None else int(basis_size)

    if final_basis_size < largest_index:
        raise ValueError(f"basis_size = {final_basis_size}, ale w pliku występuje indeks {largest_index}.")

    if final_basis_size < 1:
        raise ValueError("Nie znaleziono żadnych współczynników.")

    matrix = np.zeros((len(rows), final_basis_size), dtype=float)
    y_values = np.zeros(len(rows), dtype=float)

    for row_number, row in enumerate(rows):
        y_values[row_number] = row.g if y_name == "g" else row.omega_0

        for basis_index, value in row.coefficients.items():
            matrix[row_number, basis_index - 1] = value

    return matrix, y_values, final_basis_size


# =============================================================================
# PARAMETRY Z NAZWY PLIKU
# =============================================================================

def extract_plot_parameters(filename: Path) -> dict[str, str]:
    patterns = {
        "M": r"_M([^_]+)",
        "Jp": r"_Jp([^_]+)",
        "delta": r"_d([^_]+)",
        "T": r"_T([^_]+)",
        "P": r"_P([^_]+)",
        "F": r"_F(yes|no|both)",
        "B": r"_B([^_]+)",
        "FId": r"_FId([^_]+)",
        "eig_index": r"_eig(\d+)$",
    }

    parameters = {}

    for key, pattern in patterns.items():
        match = re.search(pattern, filename.stem)
        if match:
            parameters[key] = match.group(1)

    return parameters


def build_plot_title(filename: Path) -> str:
    p = extract_plot_parameters(filename)
    parts = []

    if "M" in p: parts.append(f"M = {p['M']}")
    if "Jp" in p: parts.append(rf"$J' = {p['Jp']}$")
    if "delta" in p: parts.append(rf"$\Delta = {p['delta']}$")
    if "T" in p: parts.append(f"T = {p['T']}")
    if "P" in p: parts.append(f"P = {p['P']}")
    if "F" in p: parts.append(f"F = {p['F']}")
    if "B" in p: parts.append(f"B = {p['B']}")
    if "FId" in p: parts.append(f"FId = {p['FId']}")
    if "eig_index" in p: parts.append(f"eig = {p['eig_index']}")

    return "; ".join(parts)


# =============================================================================
# WYKRES
# =============================================================================

def centers_to_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)

    if len(values) == 1:
        half_width = max(abs(values[0]) * 0.05, 0.5)
        return np.array([values[0] - half_width, values[0] + half_width])

    differences = np.diff(values)

    if np.any(differences <= 0):
        raise ValueError("Wartości osi y muszą być rosnące i unikalne.")

    return np.concatenate((
        [values[0] - differences[0] / 2],
        values[:-1] + differences / 2,
        [values[-1] + differences[-1] / 2],
    ))


def plot_heatmap(matrix: np.ndarray, y_values: np.ndarray, final_basis_size: int, y_name: str) -> Path | None:

    # -------------------------------------------------------------------------
    # WYBÓR ZAKRESÓW OSI X
    # -------------------------------------------------------------------------

    ranges = [(1, final_basis_size)] if x_ranges is None else x_ranges
    selected_indices = []

    for start, end in ranges:
        if start < 1 or end > final_basis_size or start > end:
            raise ValueError(f"Niepoprawny zakres osi X: ({start}, {end}). Rozmiar bazy = {final_basis_size}.")

        selected_indices.extend(range(start, end + 1))

    if len(selected_indices) != len(set(selected_indices)):
        raise ValueError("Zakresy w x_ranges nachodzą na siebie.")

    selected_indices = np.array(selected_indices, dtype=int)

    # Wybieramy tylko kolumny należące do podanych zakresów
    plot_matrix = matrix[:, selected_indices - 1]

    x_edges = np.arange(0.5, len(selected_indices) + 1.5)
    y_edges = centers_to_edges(y_values)

    fig, ax = plt.subplots(figsize=(16, 9))

    colormap = plt.get_cmap().copy()
    colormap.set_bad("white")

    # -------------------------------------------------------------------------
    # PRÓG BIAŁEGO KOLORU
    # -------------------------------------------------------------------------

    if white_below_threshold:
        plot_matrix = np.ma.masked_less(plot_matrix, threshold_white)

    # -------------------------------------------------------------------------
    # STAŁA SKALA KOLORÓW
    # -------------------------------------------------------------------------

    if use_log_values:
        if color_max <= 0:
            raise ValueError("color_max musi być dodatni dla skali logarytmicznej.")

        log_min = max(threshold_white, 1e-15)
        plot_matrix = np.ma.masked_less_equal(plot_matrix, 0.0)
        norm = LogNorm(vmin=log_min, vmax=color_max)

    else:
        norm = Normalize(vmin=color_min, vmax=color_max)

    heatmap = ax.pcolormesh(x_edges, y_edges, plot_matrix, shading="flat", norm=norm, cmap=colormap)

    colorbar = fig.colorbar(heatmap, ax=ax, pad=0.02)
    colorbar.set_label(r"$|a_i|^2$")

    ax.set_ylabel(r"$g$" if y_name == "g" else r"$\omega_0$")

    # -------------------------------------------------------------------------
    # TICKI X — TYLKO WYBRANE
    # -------------------------------------------------------------------------

    index_to_position = {original_index: position for position, original_index in enumerate(selected_indices, start=1)}

    shown_labels = [index for index in labeled_basis_indices if index in index_to_position]
    tick_positions = [index_to_position[index] for index in shown_labels]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels([str(index) for index in shown_labels])

    # -------------------------------------------------------------------------
    # LINIE ODDZIELAJĄCE ROZŁĄCZNE ZAKRESY X
    # -------------------------------------------------------------------------

    current_position = 0

    for start, end in ranges[:-1]:
        current_position += end - start + 1
        ax.axvline(current_position + 0.5, color="black", linestyle="--", linewidth=1)

    # -------------------------------------------------------------------------
    # OŚ Y
    # -------------------------------------------------------------------------

    if len(y_values) <= 25:
        ax.set_yticks(y_values)
    else:
        idx = np.linspace(0, len(y_values) - 1, min(15, len(y_values)), dtype=int)
        ax.set_yticks(y_values[idx])
        ax.set_yticklabels([f"{value:.4g}" for value in y_values[idx]])

    ax.set_xlim(0.5, len(selected_indices) + 0.5)
    ax.set_ylim(y_edges[0] if y_min is None else y_min, y_edges[-1] if y_max is None else y_max)

    # -------------------------------------------------------------------------
    # TYTUŁ
    # -------------------------------------------------------------------------

    title = build_plot_title(input_filename)

    # -------------------------------------------------------------------------
    # LEGENDA
    # -------------------------------------------------------------------------

    if show_manual_legend and manual_legend_entries:
        legend_text = manual_legend_title + " " + "; ".join(
            f"{i}: {description}" for i, description in manual_legend_entries
        )

        fig.tight_layout(rect=(0, 0.10, 1, 0.92))
        fig.text(0.5, 0.965, title, ha="center", va="top", fontsize=plot_title_fontsize)
        fig.text(0.5, 0.030, legend_text, ha="center", va="bottom", fontsize=manual_legend_fontsize)

    else:
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        fig.text(0.5, 0.965, title, ha="center", va="top", fontsize=plot_title_fontsize)

    # -------------------------------------------------------------------------
    # ZAPIS
    # -------------------------------------------------------------------------

    output_path = None

    if save_plot:
        output_directory.mkdir(parents=True, exist_ok=True)
        p = extract_plot_parameters(input_filename)

        output_name = (
            f"heatmap_M_{p.get('M', 'unknown')}_"
            f"Jp_{p.get('Jp', 'unknown')}_"
            f"d_{p.get('delta', 'unknown')}_"
            f"FId_{p.get('FId', 'unknown')}_"
            f"eig_{p.get('eig_index', 'unknown')}.png"
        )

        output_path = output_directory / output_name
        fig.savefig(output_path, dpi=dpi, facecolor="white", bbox_inches="tight")

        print(f"Zapisano wykres:\n{output_path}")

    plt.show() if show_plot else plt.close(fig)

    return output_path


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print(f"Wczytuję plik:\n{input_filename}")

    rows = load_grid_file(input_filename)
    selected_rows, y_name, fixed_name, fixed_value = select_rows(rows)
    matrix, y_values, final_basis_size = build_heatmap_matrix(selected_rows, y_name)

    print(f"Liczba wczytanych punktów: {len(rows)}")
    print(f"Liczba punktów na osi y: {len(y_values)}")
    print(f"Liczba elementów bazy: {final_basis_size}")
    print(f"Wybrany przekrój: {fixed_name} = {fixed_value:.16g}")
    print(f"Zakresy osi X: {x_ranges if x_ranges is not None else 'cała baza'}")
    print(f"Zakres skali kolorów: {color_min} - {color_max}")

    if white_below_threshold:
        print(f"Wartości mniejsze niż {threshold_white:.3e} będą białe.")

    plot_heatmap(matrix, y_values, final_basis_size, y_name)

    print("Gotowe.")


if __name__ == "__main__":
    main()