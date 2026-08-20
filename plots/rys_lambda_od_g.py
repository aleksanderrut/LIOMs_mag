from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA
# =============================================================================

# Może to być dowolny plik z danej serii.
# Numer po "_eig" zostanie automatycznie podmieniony.
filename_template = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_4_20.08.2026\siatka_omega_g\grid_M3_Jp0.0_d0.8_d20.3_Tboth_Pboth_Fyes_Bboth_FIdyes_eig1.txt"
)

# Zakres wartości własnych do narysowania na JEDNYM wykresie
eig_min = 16
eig_max = 17

save_plot = True
show_plot = True

# Folder zapisu
output_directory = Path(r"C:\Users\aleks\Desktop\praca magisterska\meeting_20.08.2026")
output_directory.mkdir(parents=True, exist_ok=True)


# =============================================================================
# SPRAWDZENIE USTAWIEŃ
# =============================================================================

if eig_min < 1:
    raise ValueError("eig_min musi być większe lub równe 1.")

if eig_max < eig_min:
    raise ValueError("eig_max nie może być mniejsze niż eig_min.")

filename_pattern = re.compile(r"_eig\d+\.txt$")

if filename_pattern.search(filename_template.name) is None:
    raise ValueError("Nazwa pliku wzorcowego musi kończyć się np. _eig5.txt")


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


def build_plot_title(filename: Path, eig_min: int, eig_max: int) -> str:
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

    if eig_min == eig_max:
        parts.append(f"eig = {eig_min}")
    else:
        parts.append(f"eig = {eig_min}-{eig_max}")

    return "; ".join(parts)


# =============================================================================
# NAZWA PLIKU WYNIKOWEGO
# =============================================================================

p = extract_plot_parameters(filename_template)

support = p.get("M", "unknown")
j_prime = p.get("Jp", "unknown")
delta = p.get("delta", "unknown")
fermion_identity = p.get("FId", "unknown")

output_filename = output_directory / (
    f"lampda_vs_g_"
    f"M_{support}_"
    f"Jp_{j_prime}_"
    f"d_{delta}_"
    f"FId_{fermion_identity}_"
    f"eig_{eig_min}-{eig_max}.png"
)

print(f"Plik wynikowy:\n{output_filename}")


# =============================================================================
# TWORZENIE JEDNEGO WSPÓLNEGO WYKRESU
# =============================================================================

fig, ax = plt.subplots(figsize=(10, 6))

for eig_index in range(eig_min, eig_max + 1):
    current_filename = filename_template.with_name(
        filename_pattern.sub(f"_eig{eig_index}.txt", filename_template.name)
    )

    if not current_filename.exists():
        raise FileNotFoundError(
            f"Nie znaleziono pliku dla eig = {eig_index}:\n{current_filename}"
        )

    data = np.loadtxt(current_filename)

    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(
            f"Plik dla eig = {eig_index} musi zawierać co najmniej trzy kolumny:\n"
            f"{current_filename}"
        )

    # kolumna 1 -> omega_0
    # kolumna 2 -> g
    # kolumna 3 -> lambda
    omega_0 = data[:, 0]
    g = data[:, 1]
    lambda_values = data[:, 2]

    ax.plot(
        g,
        lambda_values,
        marker="o",
        markersize=3.0,
        linewidth=1.2,
        label=rf"$\lambda_{{{eig_index}}}$",
    )


# =============================================================================
# FORMATOWANIE WYKRESU
# =============================================================================

ax.set_xlabel(r"$g$")
ax.set_ylabel(r"$\lambda$")
ax.set_title(build_plot_title(filename_template, eig_min, eig_max))
ax.legend()
ax.grid(alpha=0.3)

fig.tight_layout(rect=(0.0, 0.10, 1.0, 1.0))

# fig.text(
#     0.5,
#     0.025,
#     r"$S^{z}_{i,1}\left(S^{-}_{i,2}+S^{+}_{i,2}\right)"
#     r"\qquad ; \qquad"
#     r"S^{z}_{i,1}S^{z}_{i+1,1}"
#     r"\qquad\longrightarrow\qquad"
#     r"S^{+}_{i,1}S^{-}_{i+1,1}"
#     r"+S^{-}_{i,1}S^{+}_{i+1,1}$",
#     ha="center",
#     va="bottom",
#     fontsize=12,
#     bbox=dict(
#         boxstyle="round,pad=0.4",
#         facecolor="white",
#         edgecolor="black",
#     ),
# )


# =============================================================================
# ZAPIS I WYŚWIETLENIE
# =============================================================================

if save_plot:
    fig.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Zapisano wykres:\n{output_filename}")

if show_plot:
    plt.show()
else:
    plt.close(fig)