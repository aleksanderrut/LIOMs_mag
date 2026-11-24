from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# USTAWIENIA
filename = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g_liomsy\lioms_grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3_kompozycja.data"
)

# =============================================================================
# KLASYFIKACJA

# "I_n" -> np. I_0, I_1, I_2
# "P_n" -> np. P_0, P_1, P_2
# "typ" -> m, f, b
klasyfikacja = "I_n"
klasy_do_pokazania = [0, 1, 2, 3]

# =============================================================================
# TRYB WYKRESU

# "fixed_omega" -> stałe omega_0, oś X = g
# "fixed_g"     -> stałe g, oś X = omega_0
plot_mode = "fixed_omega"

# =============================================================================
# STAŁE WARTOŚCI

omega_fixed = 1.5172413793103448
g_fixed = 0.4827586206896552
value_tolerance = 1e-10

# =============================================================================
# ZAPIS / WYŚWIETLANIE

save_plot = True
show_plot = True
output_directory = Path(r"C:\Users\aleks\Desktop\praca magisterska\spotkanie_05.10.2026")
output_directory.mkdir(parents=True, exist_ok=True)

# =============================================================================
# SPRAWDZENIE USTAWIEŃ

if not filename.exists(): raise FileNotFoundError(f"Nie znaleziono pliku:\n{filename}")
if plot_mode not in ["fixed_omega", "fixed_g"]: raise ValueError('plot_mode musi być równy "fixed_omega" albo "fixed_g".')
if klasyfikacja not in ["I_n", "P_n", "typ"]: raise ValueError('klasyfikacja musi być równa "I_n", "P_n" albo "typ".')
if len(klasy_do_pokazania) == 0: raise ValueError("klasy_do_pokazania nie może być puste.")

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
        "eig_index": r"_eig(\d+)",
    }
    parameters = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, filename.stem)
        if match: parameters[key] = match.group(1)
    return parameters

# =============================================================================
# OPIS KLASY

def class_label(grupa) -> str:
    grupa = str(grupa)
    if klasyfikacja == "I_n": return rf"$I_{{{grupa}}}$"
    if klasyfikacja == "P_n": return rf"$P_{{{grupa}}}$"
    return grupa

# =============================================================================
# NAZWA KOLUMNY ODPOWIADAJĄCEJ KLASIE

def find_class_column(header: list[str], grupa) -> str:
    grupa = str(grupa)

    if klasyfikacja == "I_n":
        possible_names = [f"I_n={grupa}"]
    elif klasyfikacja == "P_n":
        possible_names = [f"P_n={grupa}", f"I_0={grupa}"]
    else:
        possible_names = [f"typ={grupa}"]

    for column_name in possible_names:
        if column_name in header: return column_name

    available = [
        name for name in header
        if name.startswith(("I_n=", "P_n=", "I_0=", "typ="))
    ]

    raise ValueError(
        f"\nNie znaleziono klasy {class_label(grupa)}.\n"
        f"Szukano kolumn: {possible_names}\n"
        f"Dostępne kolumny klasyfikacji: {available}"
    )

# =============================================================================
# WCZYTANIE NAGŁÓWKA I DANYCH

with open(filename, "r", encoding="utf-8") as file:
    lines = file.readlines()

header_index = None
header = None

for index, line in enumerate(lines):
    stripped = line.strip()
    if not stripped: continue

    columns = stripped.lstrip("#").strip().split()

    if "omega_0" in columns and "g" in columns:
        header_index = index
        header = columns
        break

if header is None: raise ValueError("Nie znaleziono nagłówka zawierającego kolumny omega_0 oraz g.")

data = np.loadtxt(filename, skiprows=header_index + 1, comments="#")
if data.ndim == 1: data = data.reshape(1, -1)

if data.shape[1] != len(header):
    raise ValueError(
        f"Liczba kolumn danych nie zgadza się z nagłówkiem. "
        f"Nagłówek: {len(header)}, dane: {data.shape[1]}."
    )

# =============================================================================
# INDEKSY PODSTAWOWYCH KOLUMN

omega_column = header.index("omega_0")
g_column = header.index("g")
omega_0_all = data[:, omega_column]
g_all = data[:, g_column]

# =============================================================================
# KOLUMNY WYBRANYCH KLAS

class_columns = {}

for grupa in klasy_do_pokazania:
    column_name = find_class_column(header, grupa)
    class_columns[str(grupa)] = header.index(column_name)

# =============================================================================
# PARAMETRY Z NAZWY PLIKU

p = extract_plot_parameters(filename)

support = p.get("M", "unknown")
j_prime = p.get("Jp", "unknown")
delta = p.get("delta", "unknown")
delta2 = p.get("delta2", "unknown")
fermion_identity = p.get("FId", "unknown")
eig_index = p.get("eig_index", "unknown")

# =============================================================================
# WYBÓR PRZEKROJU

if plot_mode == "fixed_omega":
    mask = np.isclose(omega_0_all, omega_fixed, rtol=0.0, atol=value_tolerance)

    if not np.any(mask):
        unique_omega = np.unique(omega_0_all)
        nearest_omega = unique_omega[np.argmin(np.abs(unique_omega - omega_fixed))]
        difference = abs(nearest_omega - omega_fixed)
        raise ValueError(f"Nie znaleziono omega_0 = {omega_fixed}. Najbliższe: {nearest_omega}, różnica = {difference}. Ustaw omega_fixed = {nearest_omega}")

    omega_0 = omega_0_all[mask]
    g = g_all[mask]

    sort_indices = np.argsort(g)
    omega_0 = omega_0[sort_indices]
    g = g[sort_indices]

    fixed_value_title = omega_0[0]

    if not np.allclose(omega_0, fixed_value_title, rtol=0.0, atol=value_tolerance):
        raise ValueError("Wybrane wartości omega_0 nie są jednakowe.")

    x_values = g
    x_label = r"$g$"

    print(
        f"\neig = {eig_index}\n"
        f"Tryb = stałe omega_0\n"
        f"Żądane omega_0 = {omega_fixed}\n"
        f"Znalezione omega_0 = {fixed_value_title:.16g}\n"
        f"Liczba użytych punktów = {len(g)}\n"
        f"g = {np.min(g)} -> {np.max(g)}"
    )

elif plot_mode == "fixed_g":
    mask = np.isclose(g_all, g_fixed, rtol=0.0, atol=value_tolerance)

    if not np.any(mask):
        unique_g = np.unique(g_all)
        nearest_g = unique_g[np.argmin(np.abs(unique_g - g_fixed))]
        difference = abs(nearest_g - g_fixed)
        raise ValueError(f"Nie znaleziono g = {g_fixed}. Najbliższe: {nearest_g}, różnica = {difference}. Ustaw g_fixed = {nearest_g}")

    omega_0 = omega_0_all[mask]
    g = g_all[mask]

    sort_indices = np.argsort(omega_0)
    omega_0 = omega_0[sort_indices]
    g = g[sort_indices]

    fixed_value_title = g[0]

    if not np.allclose(g, fixed_value_title, rtol=0.0, atol=value_tolerance):
        raise ValueError("Wybrane wartości g nie są jednakowe.")

    x_values = omega_0
    x_label = r"$\omega_0$"

    print(
        f"\neig = {eig_index}\n"
        f"Tryb = stałe g\n"
        f"Żądane g = {g_fixed}\n"
        f"Znalezione g = {fixed_value_title:.16g}\n"
        f"Liczba użytych punktów = {len(omega_0)}\n"
        f"omega_0 = {np.min(omega_0)} -> {np.max(omega_0)}"
    )

# =============================================================================
# TWORZENIE WYKRESU

fig, ax = plt.subplots(figsize=(10, 6))

# =============================================================================
# WYKRES WYBRANYCH KLAS

for grupa in klasy_do_pokazania:
    column_index = class_columns[str(grupa)]
    class_values = data[:, column_index][mask][sort_indices]

    ax.plot(
        x_values, class_values,
        marker="o", markersize=3.0, linewidth=1.2,
        label=class_label(grupa)
    )

# =============================================================================
# TYTUŁ

def build_plot_title(filename: Path, fixed_value: float) -> str:
    p = extract_plot_parameters(filename)
    parts = []

    if "M" in p: parts.append(f"M = {p['M']}")
    if "Jp" in p: parts.append(rf"$J' = {p['Jp']}$")
    if "delta" in p: parts.append(rf"$\Delta = {p['delta']}$")
    if "delta2" in p: parts.append(rf"$\Delta_2 = {p['delta2']}$")
    if "T" in p: parts.append(f"T = {p['T']}")
    if "P" in p: parts.append(f"P = {p['P']}")
    if "F" in p: parts.append(f"F = {p['F']}")
    if "B" in p: parts.append(f"B = {p['B']}")
    if "FId" in p: parts.append(f"FId = {p['FId']}")
    if "eig_index" in p: parts.append(f"eig = {p['eig_index']}")

    if plot_mode == "fixed_omega":
        parts.append(rf"$\omega_0 = {fixed_value:.3g}$")
    else:
        parts.append(rf"$g = {fixed_value:.3g}$")

    labels = [class_label(grupa) for grupa in klasy_do_pokazania]

    if len(labels) == 1:
        parts.append(f"class: {labels[0]}")
    else:
        parts.append("classes: " + ", ".join(labels))

    return "; ".join(parts)

# =============================================================================
# NAZWA PLIKU WYNIKOWEGO

if plot_mode == "fixed_omega":
    plot_name = "classes_vs_g"
    fixed_tag = f"_omega0_{fixed_value_title:.3g}"
else:
    plot_name = "classes_vs_omega0"
    fixed_tag = f"_g_{fixed_value_title:.3g}"

classes_tag = "-".join(str(grupa) for grupa in klasy_do_pokazania)

output_filename = output_directory / (
    f"{plot_name}_M_{support}_Jp_{j_prime}_d_{delta}_d2_{delta2}_"
    f"FId_{fermion_identity}_eig_{eig_index}_{klasyfikacja}_"
    f"{classes_tag}{fixed_tag}.png"
)

print(f"\nKlasyfikacja = {klasyfikacja}")
print("Pokazywane klasy = " + ", ".join(str(x) for x in klasy_do_pokazania))
print(f"\nPlik wynikowy:\n{output_filename}")

# =============================================================================
# FORMATOWANIE

ax.set_xlabel(x_label)
ax.set_xlim(0, 2)
ax.set_ylabel(r"$\sum |c_i|^2$")
ax.set_ylim(0, 1)
ax.set_title(build_plot_title(filename, fixed_value_title))
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