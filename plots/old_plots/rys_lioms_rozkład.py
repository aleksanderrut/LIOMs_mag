import re
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA
# =============================================================================

# Katalog, w którym znajduje się ten skrypt.
SCRIPT_DIR = Path(__file__).resolve().parent

lioms_filename = Path(r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\liom_mag_dane\LIOMs_found\LIOMs_ladder_mag_M_3_J_1.0_d_0.3_w_1.0_g_1.0_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both.txt")

basis_filename = Path(r"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\liom_mag_dane\bazy\operators_ladder_mag_M_3_J_1.0_d_0.3_w_1.0_g_1.0_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both.txt")

output_directory = SCRIPT_DIR / "lioms_wykresy" / "lioms_rozklad_w_bazie"

# None oznacza: narysuj wszystkie LIOM-y i wszystkie NONZERO MODE.
# Przykłady wyboru:
# selected_modes = ["LIOM 2"]
# selected_modes = ["LIOM 1", "NONZERO MODE 1", "NONZERO MODE 3"]
selected_modes = None

save_plots = True
show_plots = True

# False: zwykła skala liniowa.
# True: logarytmiczna skala osi y; zera nie będą widoczne jako słupki.
use_log_y = False

# Wszystkie indeksy bazy pozostają na osi x.
# Ta opcja określa tylko, czy numer ma być wypisany przy każdym indeksie.
show_every_basis_index = True

# Rozdzielczość zapisywanych obrazów.
dpi = 300

# =============================================================================
# STRUKTURY DANYCH
# =============================================================================

@dataclass
class ModeData:
    name: str
    eigenvalue: float | None = None
    coefficients: dict[str, float] = field(default_factory=dict)


# =============================================================================
# WCZYTYWANIE BAZY
# =============================================================================

def load_basis(filename: Path) -> tuple[list[str], dict[str, int]]:
    """
    Wczytuje plik bazy.
    Z linii:
        1    RR    02|30
    tworzy etykietę:
        RR 02|30
    Indeksy bazy zaczynają się od 1, zgodnie z numerem linii
    niepustego elementu w pliku.
    """
    if not filename.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku bazy:\n{filename}")

    labels: list[str] = []

    with filename.open("r", encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 3:
                raise ValueError(
                    f"Niepoprawna linia {line_number} w pliku bazy:\n{raw_line}"
                )

            symmetry_sector = parts[1]
            operator_string = parts[2]
            basis_label = f"{symmetry_sector} {operator_string}"
            labels.append(basis_label)

    if not labels:
        raise ValueError("Plik bazy nie zawiera żadnych elementów.")

    duplicates = sorted(
        label for label in set(labels) if labels.count(label) > 1
    )
    if duplicates:
        raise ValueError(
            "Plik bazy zawiera powtarzające się etykiety:\n"
            + "\n".join(duplicates)
        )

    label_to_index = {
        label: index
        for index, label in enumerate(labels, start=1)
    }

    return labels, label_to_index


# =============================================================================
# WCZYTYWANIE LIOM-ÓW I NONZERO MODES
# =============================================================================

MODE_HEADER_PATTERN = re.compile(
    r"^#\s*(LIOM\s+\d+|NONZERO MODE\s+\d+)\s*$"
)

EIGENVALUE_PATTERN = re.compile(
    r"^#\s*eigenvalue\s*=\s*(.+?)\s*$"
)


def load_liom_modes(filename: Path) -> list[ModeData]:
    """
    Wczytuje wszystkie bloki:
        # LIOM n
        # NONZERO MODE n
    Każda linia danych ma format:
        translated_basis_element <TAB> basis_element <TAB> squared_coefficient
    """
    if not filename.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku LIOM-ów:\n{filename}")

    modes: list[ModeData] = []
    current_mode: ModeData | None = None

    with filename.open("r", encoding="utf-8") as file:
        for line_number, raw_line in enumerate(file, start=1):
            line = raw_line.strip()

            if not line:
                continue

            mode_match = MODE_HEADER_PATTERN.match(line)
            if mode_match:
                current_mode = ModeData(name=mode_match.group(1))
                modes.append(current_mode)
                continue

            eigenvalue_match = EIGENVALUE_PATTERN.match(line)
            if eigenvalue_match and current_mode is not None:
                current_mode.eigenvalue = float(eigenvalue_match.group(1))
                continue

            if line.startswith("#"):
                continue

            if current_mode is None:
                raise ValueError(
                    f"Dane przed nagłówkiem LIOM/NONZERO MODE, linia {line_number}."
                )

            parts = raw_line.rstrip("\n").split("\t")

            if len(parts) < 3:
                raise ValueError(
                    f"Niepoprawna linia {line_number} w pliku LIOM-ów:\n"
                    f"{raw_line}"
                )

            basis_label = parts[-2].strip()
            squared_coefficient = float(parts[-1].strip())

            if basis_label in current_mode.coefficients:
                raise ValueError(
                    f"Etykieta {basis_label!r} występuje dwa razy "
                    f"w bloku {current_mode.name}."
                )

            current_mode.coefficients[basis_label] = squared_coefficient

    if not modes:
        raise ValueError("Nie znaleziono żadnych bloków LIOM ani NONZERO MODE.")

    return modes


# =============================================================================
# ŁĄCZENIE DANYCH Z INDEKSAMI BAZY
# =============================================================================

def mode_to_basis_vector(
    mode: ModeData,
    basis_labels: list[str],
    label_to_index: dict[str, int],
) -> np.ndarray:
    """
    Tworzy wektor długości równej całej bazie.

    Elementy, które nie zostały zapisane w pliku LIOM-ów ze względu
    na coeff_tol, otrzymują wartość 0.
    """
    values = np.zeros(len(basis_labels), dtype=float)
    missing_labels: list[str] = []

    for basis_label, squared_coefficient in mode.coefficients.items():
        basis_index = label_to_index.get(basis_label)

        if basis_index is None:
            missing_labels.append(basis_label)
            continue

        values[basis_index - 1] = squared_coefficient

    if missing_labels:
        raise KeyError(
            f"W bloku {mode.name} znaleziono etykiety nieobecne w pliku bazy:\n"
            + "\n".join(missing_labels)
        )

    return values


# =============================================================================
# WYKRES
# =============================================================================

def safe_filename(text: str) -> str:
    return (
        text.lower()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
    )


def plot_mode(
    mode: ModeData,
    values: np.ndarray,
    output_directory: Path,
) -> Path | None:
    basis_indices = np.arange(1, len(values) + 1)

    # Szerokość rośnie wraz z liczbą elementów bazy.
    figure_width = max(14.0, 0.19 * len(values))

    plt.figure(figsize=(figure_width, 6.5))
    plt.bar(basis_indices, values, width=0.85)

    if use_log_y:
        plt.yscale("log")

    plt.xlim(0.3, len(values) + 0.7)
    plt.xlabel("Indeks stanu bazowego")
    plt.ylabel(r"$|a_i|^2$")

    eigenvalue_text = (
        "brak"
        if mode.eigenvalue is None
        else f"{mode.eigenvalue:.8e}"
    )
    plt.title(
        f"{mode.name}   |   eigenvalue = {eigenvalue_text}\n"
    )

    if show_every_basis_index:
        plt.xticks(
            basis_indices,
            [str(index) for index in basis_indices],
            rotation=90,
            fontsize=6,
        )
    else:
        tick_step = max(1, len(values) // 20)
        shown_indices = basis_indices[::tick_step]
        plt.xticks(shown_indices, shown_indices)

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    output_path: Path | None = None

    if save_plots:
        output_directory.mkdir(parents=True, exist_ok=True)
        source_file_name = lioms_filename.stem
        mode_name = mode.name.replace(" ", "_")
        output_path = output_directory / (f"{mode_name}_{source_file_name}.png")
        plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"Zapisano: {output_path}")

    if show_plots:
        plt.show()
    else:
        plt.close()

    return output_path


# =============================================================================
# PROGRAM GŁÓWNY
# =============================================================================

def main() -> None:
    print("Wczytuję plik bazy:")
    print(basis_filename)

    basis_labels, label_to_index = load_basis(basis_filename)

    print(f"Liczba stanów bazowych: {len(basis_labels)}")

    # Kontrolny przykład z Twojego pliku:
    example_label = "RR 02|30"
    if example_label in label_to_index:
        print(
            f"Kontrola indeksowania: {example_label} "
            f"ma indeks {label_to_index[example_label]}"
        )

    print("\nWczytuję plik LIOM-ów:")
    print(lioms_filename)

    modes = load_liom_modes(lioms_filename)
    print(f"Liczba znalezionych bloków: {len(modes)}")

    available_names = [mode.name for mode in modes]
    print("Dostępne bloki:")
    for name in available_names:
        print(f"  {name}")

    if selected_modes is None:
        modes_to_plot = modes
    else:
        requested = set(selected_modes)
        modes_to_plot = [
            mode for mode in modes
            if mode.name in requested
        ]

        missing_modes = sorted(requested - set(available_names))
        if missing_modes:
            raise ValueError(
                "Nie znaleziono wybranych bloków:\n"
                + "\n".join(missing_modes)
            )

    if not modes_to_plot:
        raise ValueError("Nie wybrano żadnego bloku do narysowania.")

    print("\nTworzę wykresy:")

    for mode in modes_to_plot:
        values = mode_to_basis_vector(
            mode,
            basis_labels,
            label_to_index,
        )
        plot_mode(mode, values, output_directory)

    print("\nGotowe.")


if __name__ == "__main__":
    main()