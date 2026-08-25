from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


# =============================================================================
# USTAWIENIA
# =============================================================================

filename_template = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_6_25.08.2026\siatka_omega_g\grid_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig1.txt"
)

eig_min = 1
eig_max = 20


# =============================================================================
# TRYB WYKRESU
# =============================================================================
#
# "fixed_omega"
#     -> ustalone omega_0
#     -> skrypt sam wyszukuje wszystkie odpowiednie wiersze
#     -> oś X = g
#     -> wykres lambda(g)
#
# "fixed_g"
#     -> ustalone g
#     -> skrypt sam wyszukuje wszystkie odpowiednie wiersze
#     -> oś X = omega_0
#     -> wykres lambda(omega_0)
#
# =============================================================================

plot_mode = "fixed_g"


# =============================================================================
# STAŁE WARTOŚCI
# =============================================================================
#
# Używana jest tylko wartość odpowiadająca wybranemu plot_mode.
#
# Przykład:
#
# plot_mode = "fixed_omega"
# omega_fixed = 0.5
#
# albo:
#
# plot_mode = "fixed_g"
# g_fixed = 0.5
#
# =============================================================================

omega_fixed = 0.5063291139240507
g_fixed = 1.4936708860759493

# =============================================================================
# TOLERANCJA PRZY WYSZUKIWANIU WARTOŚCI
# =============================================================================
#
# Używamy np.isclose zamiast zwykłego ==,
# ponieważ liczby w pliku są zapisane jako Float64.
#
# =============================================================================

value_tolerance = 1e-10


# =============================================================================
# ZAPIS / WYŚWIETLANIE
# =============================================================================

save_plot = True
show_plot = True

output_directory = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer"
)

output_directory.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# SPRAWDZENIE USTAWIEŃ
# =============================================================================

if eig_min < 1:
    raise ValueError(
        "eig_min musi być większe lub równe 1."
    )

if eig_max < eig_min:
    raise ValueError(
        "eig_max nie może być mniejsze niż eig_min."
    )

if plot_mode not in [
    "fixed_omega",
    "fixed_g"
]:
    raise ValueError(
        'plot_mode musi być równy "fixed_omega" albo "fixed_g".'
    )


filename_pattern = re.compile(
    r"_eig\d+\.txt$"
)

if filename_pattern.search(
    filename_template.name
) is None:

    raise ValueError(
        "Nazwa pliku wzorcowego musi kończyć się np. _eig5.txt"
    )


# =============================================================================
# PARAMETRY Z NAZWY PLIKU
# =============================================================================

def extract_plot_parameters(
    filename: Path
) -> dict[str, str]:

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

        match = re.search(
            pattern,
            filename.stem
        )

        if match:
            parameters[key] = match.group(1)

    return parameters


# =============================================================================
# TYTUŁ
# =============================================================================

def build_plot_title(
    filename: Path,
    eig_min: int,
    eig_max: int,
    fixed_value: float
) -> str:

    p = extract_plot_parameters(
        filename
    )

    parts = []

    if "M" in p:
        parts.append(
            f"M = {p['M']}"
        )

    if "Jp" in p:
        parts.append(
            rf"$J' = {p['Jp']}$"
        )

    if "delta" in p:
        parts.append(
            rf"$\Delta = {p['delta']}$"
        )

    if "delta2" in p:
        parts.append(
            rf"$\Delta_2 = {p['delta2']}$"
        )

    if "T" in p:
        parts.append(
            f"T = {p['T']}"
        )

    if "P" in p:
        parts.append(
            f"P = {p['P']}"
        )

    if "F" in p:
        parts.append(
            f"F = {p['F']}"
        )

    if "B" in p:
        parts.append(
            f"B = {p['B']}"
        )

    if "FId" in p:
        parts.append(
            f"FId = {p['FId']}"
        )


    if eig_min == eig_max:

        parts.append(
            f"eig = {eig_min}"
        )

    else:

        parts.append(
            f"eig = {eig_min}-{eig_max}"
        )


    # =========================================================================
    # STAŁY PARAMETR
    # 3 CYFRY ZNACZĄCE
    # =========================================================================

    if plot_mode == "fixed_omega":

        parts.append(
            rf"$\omega_0 = {fixed_value:.3g}$"
        )

    elif plot_mode == "fixed_g":

        parts.append(
            rf"$g = {fixed_value:.3g}$"
        )


    return "; ".join(parts)


# =============================================================================
# PARAMETRY DO NAZWY PLIKU WYNIKOWEGO
# =============================================================================

p = extract_plot_parameters(
    filename_template
)


support = p.get(
    "M",
    "unknown"
)

j_prime = p.get(
    "Jp",
    "unknown"
)

delta = p.get(
    "delta",
    "unknown"
)

delta2 = p.get(
    "delta2",
    "unknown"
)

fermion_identity = p.get(
    "FId",
    "unknown"
)


# =============================================================================
# TWORZENIE WYKRESU
# =============================================================================

fig, ax = plt.subplots(
    figsize=(10, 6)
)


# =============================================================================
# STAŁA WARTOŚĆ RZECZYWIŚCIE ZNALEZIONA W PLIKU
# =============================================================================

fixed_value_title = None


# =============================================================================
# PĘTLA PO EIG
# =============================================================================

for eig_index in range(
    eig_min,
    eig_max + 1
):

    current_filename = filename_template.with_name(
        filename_pattern.sub(
            f"_eig{eig_index}.txt",
            filename_template.name
        )
    )


    if not current_filename.exists():

        raise FileNotFoundError(
            f"Nie znaleziono pliku dla eig = {eig_index}:\n"
            f"{current_filename}"
        )


    # =========================================================================
    # WCZYTANIE CAŁEGO PLIKU
    # =========================================================================

    data = np.loadtxt(
        current_filename
    )


    if data.ndim != 2 or data.shape[1] < 3:

        raise ValueError(
            f"Plik dla eig = {eig_index} musi zawierać "
            f"co najmniej trzy kolumny:\n"
            f"{current_filename}"
        )


    # =========================================================================
    # KOLUMNY
    # =========================================================================
    #
    # kolumna 1 -> omega_0
    # kolumna 2 -> g
    # kolumna 3 -> lambda
    #
    # =========================================================================

    omega_0_all = data[:, 0]

    g_all = data[:, 1]

    lambda_all = data[:, 2]


    # =========================================================================
    # TRYB: STAŁE OMEGA_0
    # =========================================================================

    if plot_mode == "fixed_omega":

        # ---------------------------------------------------------------------
        # Znajdujemy wszystkie punkty posiadające zadane omega_0
        # ---------------------------------------------------------------------

        mask = np.isclose(
            omega_0_all,
            omega_fixed,
            rtol=0.0,
            atol=value_tolerance
        )


        # ---------------------------------------------------------------------
        # Jeżeli nie znaleziono wartości omega_0
        # ---------------------------------------------------------------------

        if not np.any(mask):

            unique_omega = np.unique(
                omega_0_all
            )

            nearest_index = np.argmin(
                np.abs(
                    unique_omega - omega_fixed
                )
            )

            nearest_omega = unique_omega[
                nearest_index
            ]

            difference = abs(
                nearest_omega - omega_fixed
            )

            raise ValueError(
                f"\nNie znaleziono omega_0 = {omega_fixed} "
                f"dla eig = {eig_index}.\n"
                f"Najbliższa wartość w pliku to:\n"
                f"omega_0 = {nearest_omega}\n"
                f"różnica = {difference}\n\n"
                f"Jeżeli chcesz jej użyć, ustaw:\n"
                f"omega_fixed = {nearest_omega}"
            )


        # ---------------------------------------------------------------------
        # Wyciągamy wszystkie pasujące punkty
        # ---------------------------------------------------------------------

        omega_0 = omega_0_all[
            mask
        ]

        g = g_all[
            mask
        ]

        lambda_values = lambda_all[
            mask
        ]


        # ---------------------------------------------------------------------
        # Sortowanie po g
        # ---------------------------------------------------------------------

        sort_indices = np.argsort(
            g
        )

        omega_0 = omega_0[
            sort_indices
        ]

        g = g[
            sort_indices
        ]

        lambda_values = lambda_values[
            sort_indices
        ]


        # ---------------------------------------------------------------------
        # Faktyczna wartość omega_0
        # ---------------------------------------------------------------------

        current_fixed_value = omega_0[0]


        # ---------------------------------------------------------------------
        # Sprawdzamy, czy omega_0 jest rzeczywiście stałe
        # ---------------------------------------------------------------------

        if not np.allclose(
            omega_0,
            current_fixed_value,
            rtol=0.0,
            atol=value_tolerance
        ):

            raise ValueError(
                f"Dla eig={eig_index} wybrane wartości omega_0 "
                f"nie są jednakowe.\n"
                f"omega_min = {np.min(omega_0)}\n"
                f"omega_max = {np.max(omega_0)}"
            )


        # ---------------------------------------------------------------------
        # Sprawdzenie zgodności między różnymi eig
        # ---------------------------------------------------------------------

        if fixed_value_title is None:

            fixed_value_title = current_fixed_value

        elif not np.isclose(
            fixed_value_title,
            current_fixed_value,
            rtol=0.0,
            atol=value_tolerance
        ):

            raise ValueError(
                f"Omega_0 dla eig={eig_index} różni się "
                f"od wartości w poprzednich plikach.\n"
                f"Poprzednie omega_0 = {fixed_value_title}\n"
                f"Obecne omega_0 = {current_fixed_value}"
            )


        # ---------------------------------------------------------------------
        # Dane na wykres
        # ---------------------------------------------------------------------

        x_values = g

        x_label = r"$g$"


        # ---------------------------------------------------------------------
        # Informacja w terminalu
        # ---------------------------------------------------------------------

        print()

        print(
            f"eig = {eig_index}"
        )

        print(
            "Tryb = stałe omega_0"
        )

        print(
            f"Żądane omega_0 = {omega_fixed}"
        )

        print(
            f"Znalezione omega_0 = {current_fixed_value:.16g}"
        )

        print(
            f"Liczba użytych punktów = {len(g)}"
        )

        print(
            f"g = {np.min(g)} -> {np.max(g)}"
        )


    # =========================================================================
    # TRYB: STAŁE g
    # =========================================================================

    elif plot_mode == "fixed_g":

        # ---------------------------------------------------------------------
        # Znajdujemy wszystkie punkty posiadające zadane g
        # ---------------------------------------------------------------------

        mask = np.isclose(
            g_all,
            g_fixed,
            rtol=0.0,
            atol=value_tolerance
        )


        # ---------------------------------------------------------------------
        # Jeżeli nie znaleziono wartości g
        # ---------------------------------------------------------------------

        if not np.any(mask):

            unique_g = np.unique(
                g_all
            )

            nearest_index = np.argmin(
                np.abs(
                    unique_g - g_fixed
                )
            )

            nearest_g = unique_g[
                nearest_index
            ]

            difference = abs(
                nearest_g - g_fixed
            )

            raise ValueError(
                f"\nNie znaleziono g = {g_fixed} "
                f"dla eig = {eig_index}.\n"
                f"Najbliższa wartość w pliku to:\n"
                f"g = {nearest_g}\n"
                f"różnica = {difference}\n\n"
                f"Jeżeli chcesz jej użyć, ustaw:\n"
                f"g_fixed = {nearest_g}"
            )


        # ---------------------------------------------------------------------
        # Wyciągamy wszystkie pasujące punkty
        # ---------------------------------------------------------------------

        omega_0 = omega_0_all[
            mask
        ]

        g = g_all[
            mask
        ]

        lambda_values = lambda_all[
            mask
        ]


        # ---------------------------------------------------------------------
        # Sortowanie po omega_0
        # ---------------------------------------------------------------------

        sort_indices = np.argsort(
            omega_0
        )

        omega_0 = omega_0[
            sort_indices
        ]

        g = g[
            sort_indices
        ]

        lambda_values = lambda_values[
            sort_indices
        ]


        # ---------------------------------------------------------------------
        # Faktyczna wartość g
        # ---------------------------------------------------------------------

        current_fixed_value = g[0]


        # ---------------------------------------------------------------------
        # Sprawdzamy, czy g jest rzeczywiście stałe
        # ---------------------------------------------------------------------

        if not np.allclose(
            g,
            current_fixed_value,
            rtol=0.0,
            atol=value_tolerance
        ):

            raise ValueError(
                f"Dla eig={eig_index} wybrane wartości g "
                f"nie są jednakowe.\n"
                f"g_min = {np.min(g)}\n"
                f"g_max = {np.max(g)}"
            )


        # ---------------------------------------------------------------------
        # Sprawdzenie zgodności między różnymi eig
        # ---------------------------------------------------------------------

        if fixed_value_title is None:

            fixed_value_title = current_fixed_value

        elif not np.isclose(
            fixed_value_title,
            current_fixed_value,
            rtol=0.0,
            atol=value_tolerance
        ):

            raise ValueError(
                f"g dla eig={eig_index} różni się "
                f"od wartości w poprzednich plikach.\n"
                f"Poprzednie g = {fixed_value_title}\n"
                f"Obecne g = {current_fixed_value}"
            )


        # ---------------------------------------------------------------------
        # Dane na wykres
        # ---------------------------------------------------------------------

        x_values = omega_0

        x_label = r"$\omega_0$"


        # ---------------------------------------------------------------------
        # Informacja w terminalu
        # ---------------------------------------------------------------------

        print()

        print(
            f"eig = {eig_index}"
        )

        print(
            "Tryb = stałe g"
        )

        print(
            f"Żądane g = {g_fixed}"
        )

        print(
            f"Znalezione g = {current_fixed_value:.16g}"
        )

        print(
            f"Liczba użytych punktów = {len(omega_0)}"
        )

        print(
            f"omega_0 = "
            f"{np.min(omega_0)} -> {np.max(omega_0)}"
        )


    # =========================================================================
    # WYKRES
    # =========================================================================

    ax.plot(
        x_values,
        lambda_values,

        marker="o",
        markersize=3.0,

        linewidth=1.2,

        label=rf"$\lambda_{{{eig_index}}}$",
    )


# =============================================================================
# NAZWA PLIKU WYNIKOWEGO
# =============================================================================

if plot_mode == "fixed_omega":

    plot_name = "lampda_vs_g"

    fixed_tag = (
        f"_omega0_{fixed_value_title:.3g}"
    )


elif plot_mode == "fixed_g":

    plot_name = "lampda_vs_omega0"

    fixed_tag = (
        f"_g_{fixed_value_title:.3g}"
    )


output_filename = output_directory / (

    f"{plot_name}_"

    f"M_{support}_"

    f"Jp_{j_prime}_"

    f"d_{delta}_"

    f"d2_{delta2}_"

    f"FId_{fermion_identity}_"

    f"eig_{eig_min}-{eig_max}"

    f"{fixed_tag}.png"
)


print()

print(
    f"Plik wynikowy:\n{output_filename}"
)


# =============================================================================
# FORMATOWANIE
# =============================================================================

ax.set_xlabel(
    x_label
)

# Zakres osi X od 0 do 2
ax.set_xlim(
    0,
    2
)

ax.set_ylabel(
    r"$\lambda$"
)

ax.set_title(
    build_plot_title(
        filename_template,
        eig_min,
        eig_max,
        fixed_value_title
    )
)

ax.legend()

ax.grid(
    alpha=0.3
)


# =============================================================================
# ZAPIS I WYŚWIETLENIE
# =============================================================================

if save_plot:

    fig.savefig(
        output_filename,
        dpi=300,
        bbox_inches="tight"
    )

    print()

    print(
        f"Zapisano wykres:\n{output_filename}"
    )


if show_plot:

    plt.show()

else:

    plt.close(fig)