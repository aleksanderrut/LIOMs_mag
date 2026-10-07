from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA

filename_template = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g\grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"
)

eig_min, eig_max = 1, 20

# "fixed_omega" -> stałe omega_0, oś X = g
# "fixed_g"     -> stałe g, oś X = omega_0
plot_mode = "fixed_g"

omega_fixed = 0.4827586206896552
g_fixed = 0.4827586206896552
value_tolerance = 1e-10

save_plot = True
show_plot = True

output_directory = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_05.10.2026"
)
output_directory.mkdir(parents=True, exist_ok=True)

# =============================================================================
# SPRAWDZENIE USTAWIEŃ

if eig_min < 1:
    raise ValueError("eig_min musi być większe lub równe 1.")

if eig_max < eig_min:
    raise ValueError("eig_max nie może być mniejsze niż eig_min.")

if plot_mode not in ["fixed_omega", "fixed_g"]:
    raise ValueError('plot_mode musi być równy "fixed_omega" albo "fixed_g".')

filename_pattern = re.compile(r"_eig\d+\.txt$")

if filename_pattern.search(filename_template.name) is None:
    raise ValueError("Nazwa pliku wzorcowego musi kończyć się np. _eig5.txt")

# =============================================================================
# PARAMETRY Z NAZWY PLIKU

def extract_plot_parameters(filename: Path) -> dict[str, str]:
    patterns = {
        "M": r"_M([^_]+)",
        "Jp": r"_Jp([^_]+)",
        "delta": r"_d([^_]+)",
        "delta2": r"_d2([^_]+)",
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

# =============================================================================
# TYTUŁ

def build_plot_title(filename: Path, eig_min: int, eig_max: int, fixed_value: float) -> str:
    p = extract_plot_parameters(filename)
    parts = []

    if "M" in p:
        parts.append(f"M = {p['M']}")
    if "Jp" in p:
        parts.append(rf"$J' = {p['Jp']}$")
    if "delta" in p:
        parts.append(rf"$\Delta = {p['delta']}$")
    if "delta2" in p:
        parts.append(rf"$\Delta_2 = {p['delta2']}$")
    if "T" in p:
        parts.append(f"T = {p['T']}")
    if "P" in p:
        parts.append(f"P = {p['P']}")
    if "F" in p:
        parts.append(f"F = {p['F']}")
    if "B" in p:
        parts.append(f"B = {p['B']}")
    if "FId" in p:
        parts.append(f"FId = {p['FId']}")

    parts.append(f"eig = {eig_min}" if eig_min == eig_max else f"eig = {eig_min}-{eig_max}")

    if plot_mode == "fixed_omega":
        parts.append(rf"$\omega_0 = {fixed_value:.3g}$")
    else:
        parts.append(rf"$g = {fixed_value:.3g}$")

    return "; ".join(parts)

# =============================================================================
# PARAMETRY DO NAZWY PLIKU

p = extract_plot_parameters(filename_template)

support = p.get("M", "unknown")
j_prime = p.get("Jp", "unknown")
delta = p.get("delta", "unknown")
delta2 = p.get("delta2", "unknown")
fermion_identity = p.get("FId", "unknown")

# =============================================================================
# WYKRES

fig, ax = plt.subplots(figsize=(10, 6))
fixed_value_title = None

for eig_index in range(eig_min, eig_max + 1):

    current_filename = filename_template.with_name(
        filename_pattern.sub(f"_eig{eig_index}.txt", filename_template.name)
    )

    if not current_filename.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku dla eig = {eig_index}:\n{current_filename}")

    data = np.loadtxt(current_filename)

    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(
            f"Plik dla eig = {eig_index} musi zawierać co najmniej trzy kolumny:\n{current_filename}"
        )

    omega_0_all = data[:, 0]
    g_all = data[:, 1]
    lambda_all = data[:, 2]

    # STAŁE OMEGA_0
    if plot_mode == "fixed_omega":

        mask = np.isclose(omega_0_all, omega_fixed, rtol=0.0, atol=value_tolerance)

        if not np.any(mask):
            unique_omega = np.unique(omega_0_all)
            nearest_omega = unique_omega[np.argmin(np.abs(unique_omega - omega_fixed))]
            difference = abs(nearest_omega - omega_fixed)

            raise ValueError(
                f"\nNie znaleziono omega_0 = {omega_fixed} dla eig = {eig_index}.\n"
                f"Najbliższa wartość: omega_0 = {nearest_omega}\n"
                f"różnica = {difference}\n\n"
                f"Ustaw: omega_fixed = {nearest_omega}"
            )

        omega_0 = omega_0_all[mask]
        g = g_all[mask]
        lambda_values = lambda_all[mask]

        sort_indices = np.argsort(g)
        omega_0 = omega_0[sort_indices]
        g = g[sort_indices]
        lambda_values = lambda_values[sort_indices]

        current_fixed_value = omega_0[0]

        if not np.allclose(omega_0, current_fixed_value, rtol=0.0, atol=value_tolerance):
            raise ValueError(
                f"Dla eig={eig_index} omega_0 nie jest stałe: "
                f"{np.min(omega_0)} -> {np.max(omega_0)}"
            )

        if fixed_value_title is None:
            fixed_value_title = current_fixed_value

        elif not np.isclose(
            fixed_value_title, current_fixed_value,
            rtol=0.0, atol=value_tolerance
        ):
            raise ValueError(
                f"Omega_0 dla eig={eig_index} różni się od poprzednich plików: "
                f"{fixed_value_title} != {current_fixed_value}"
            )

        x_values = g
        x_label = r"$g$"

        print(
            f"\neig = {eig_index}\n"
            f"Tryb = stałe omega_0\n"
            f"Żądane omega_0 = {omega_fixed}\n"
            f"Znalezione omega_0 = {current_fixed_value:.16g}\n"
            f"Liczba punktów = {len(g)}\n"
            f"g = {np.min(g)} -> {np.max(g)}"
        )

    # STAŁE g
    elif plot_mode == "fixed_g":

        mask = np.isclose(g_all, g_fixed, rtol=0.0, atol=value_tolerance)

        if not np.any(mask):
            unique_g = np.unique(g_all)
            nearest_g = unique_g[np.argmin(np.abs(unique_g - g_fixed))]
            difference = abs(nearest_g - g_fixed)

            raise ValueError(
                f"\nNie znaleziono g = {g_fixed} dla eig = {eig_index}.\n"
                f"Najbliższa wartość: g = {nearest_g}\n"
                f"różnica = {difference}\n\n"
                f"Ustaw: g_fixed = {nearest_g}"
            )

        omega_0 = omega_0_all[mask]
        g = g_all[mask]
        lambda_values = lambda_all[mask]

        sort_indices = np.argsort(omega_0)
        omega_0 = omega_0[sort_indices]
        g = g[sort_indices]
        lambda_values = lambda_values[sort_indices]

        current_fixed_value = g[0]

        if not np.allclose(g, current_fixed_value, rtol=0.0, atol=value_tolerance):
            raise ValueError(
                f"Dla eig={eig_index} g nie jest stałe: "
                f"{np.min(g)} -> {np.max(g)}"
            )

        if fixed_value_title is None:
            fixed_value_title = current_fixed_value

        elif not np.isclose(
            fixed_value_title, current_fixed_value,
            rtol=0.0, atol=value_tolerance
        ):
            raise ValueError(
                f"g dla eig={eig_index} różni się od poprzednich plików: "
                f"{fixed_value_title} != {current_fixed_value}"
            )

        x_values = omega_0
        x_label = r"$\omega_0$"

        print(
            f"\neig = {eig_index}\n"
            f"Tryb = stałe g\n"
            f"Żądane g = {g_fixed}\n"
            f"Znalezione g = {current_fixed_value:.16g}\n"
            f"Liczba punktów = {len(omega_0)}\n"
            f"omega_0 = {np.min(omega_0)} -> {np.max(omega_0)}"
        )

    ax.plot(
        x_values, lambda_values,
        marker="o", markersize=3.0, linewidth=1.2,
        label=rf"$\lambda_{{{eig_index}}}$"
    )

# =============================================================================
# NAZWA PLIKU

if plot_mode == "fixed_omega":
    plot_name = "lampda_vs_g"
    fixed_tag = f"_omega0_{fixed_value_title:.3g}"

else:
    plot_name = "lampda_vs_omega0"
    fixed_tag = f"_g_{fixed_value_title:.3g}"

output_filename = output_directory / (
    f"{plot_name}_M_{support}_Jp_{j_prime}_d_{delta}_d2_{delta2}_"
    f"FId_{fermion_identity}_eig_{eig_min}-{eig_max}{fixed_tag}.png"
)

print(f"\nPlik wynikowy:\n{output_filename}")

# =============================================================================
# FORMATOWANIE

ax.set_xlabel(x_label)
ax.set_xlim(0, 2)
ax.set_ylabel(r"$\lambda$")
ax.set_title(build_plot_title(filename_template, eig_min, eig_max, fixed_value_title))
ax.legend()
ax.grid(alpha=0.3)

# =============================================================================
# ZAPIS I WYŚWIETLENIE

if save_plot:
    fig.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"\nZapisano wykres:\n{output_filename}")

if show_plot:
    plt.show()
else:
    plt.close(fig)