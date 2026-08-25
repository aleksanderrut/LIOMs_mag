from pathlib import Path


# =============================================================================
# PLIKI DO PRZESORTOWANIA
# =============================================================================
#
# Możesz podać tutaj zarówno:
#
# 1. plik grid z lambdami:
#    omega_0    g    lambda
#
# 2. plik lioms_grid:
#    omega_0    g    3:...;4:...;...
#
# Skrypt działa na obu formatach.
#
# Możesz podać jeden albo wiele plików.
# =============================================================================

filenames = [

    Path(
        r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_6_25.08.2026\siatka_omega_g\grid_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig1.txt"
    ),

    Path(
        r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_6_25.08.2026\siatka_omega_g_liomsy\lioms_grid_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig4.txt"
    ),
]


# =============================================================================
# SPOSÓB SORTOWANIA
# =============================================================================
#
# True:
#
# najpierw po g,
# potem wewnątrz każdego g po omega_0
#
# czyli:
#
# g = 0.0:
#     omega_0 = 0.0
#     omega_0 = 0.025
#     omega_0 = 0.050
#     ...
#
# g = 0.025:
#     omega_0 = 0.0
#     omega_0 = 0.025
#     omega_0 = 0.050
#     ...
#
# =============================================================================

sort_by_g_then_omega = True


# =============================================================================
# FUNKCJA SORTUJĄCA PLIK
# =============================================================================

def reorder_file(
    filename: Path
):

    # =========================================================================
    # SPRAWDZENIE CZY PLIK ISTNIEJE
    # =========================================================================

    if not filename.exists():

        raise FileNotFoundError(
            f"Nie znaleziono pliku:\n{filename}"
        )


    # =========================================================================
    # WCZYTANIE
    # =========================================================================

    data_lines = []
    other_lines = []


    with open(
        filename,
        "r",
        encoding="utf-8"
    ) as f:

        for line_number, line in enumerate(
            f,
            start=1
        ):

            stripped = line.strip()


            # -----------------------------------------------------------------
            # Pusta linia
            # -----------------------------------------------------------------

            if not stripped:
                continue


            # -----------------------------------------------------------------
            # Podział tylko w celu odczytania omega_0 i g
            #
            # Sama oryginalna linia NIE jest modyfikowana.
            # -----------------------------------------------------------------

            columns = stripped.split()


            if len(columns) < 2:

                other_lines.append(
                    line
                )

                continue


            # -----------------------------------------------------------------
            # Pierwsza kolumna = omega_0
            # Druga kolumna    = g
            # -----------------------------------------------------------------

            try:

                omega_0 = float(
                    columns[0]
                )

                g = float(
                    columns[1]
                )

            except ValueError:

                # np. jakiś nagłówek tekstowy
                other_lines.append(
                    line
                )

                continue


            # -----------------------------------------------------------------
            # Zapamiętujemy:
            #
            # omega_0
            # g
            # numer pierwotnego wiersza
            # CAŁĄ oryginalną linię
            #
            # -----------------------------------------------------------------

            data_lines.append(
                (
                    omega_0,
                    g,
                    line_number,
                    line
                )
            )


    # =========================================================================
    # KONTROLA
    # =========================================================================

    if len(data_lines) == 0:

        raise ValueError(
            f"Nie znaleziono danych numerycznych w pliku:\n"
            f"{filename}"
        )


    # =========================================================================
    # SORTOWANIE
    # =========================================================================
    #
    # Element:
    #
    # (
    #     omega_0,
    #     g,
    #     line_number,
    #     line
    # )
    #
    # Dlatego:
    #
    # x[0] -> omega_0
    # x[1] -> g
    #
    # Chcemy:
    #
    # najpierw g,
    # potem omega_0
    #
    # więc klucz:
    #
    # (x[1], x[0])
    #
    # =========================================================================

    if sort_by_g_then_omega:

        data_lines.sort(
            key=lambda x: (
                x[1],
                x[0]
            )
        )

    else:

        # Opcjonalnie:
        # standardowa kolejność omega_0 -> g

        data_lines.sort(
            key=lambda x: (
                x[0],
                x[1]
            )
        )


    # =========================================================================
    # NAZWA NOWEGO PLIKU
    # =========================================================================
    #
    # np.
    #
    # grid_....txt
    #
    # ->
    #
    # grid_...._reorder.txt
    #
    # =========================================================================

    output_filename = filename.with_name(
        filename.stem
        + "_reorder"
        + filename.suffix
    )


    # =========================================================================
    # ZAPIS
    # =========================================================================

    with open(
        output_filename,
        "w",
        encoding="utf-8"
    ) as f:


        # ---------------------------------------------------------------------
        # Jeżeli w pliku były jakieś linie tekstowe / nagłówki,
        # zachowujemy je na początku.
        # ---------------------------------------------------------------------

        for line in other_lines:

            if line.endswith("\n"):

                f.write(
                    line
                )

            else:

                f.write(
                    line + "\n"
                )


        # ---------------------------------------------------------------------
        # Dane w nowej kolejności
        # ---------------------------------------------------------------------

        for omega_0, g, line_number, line in data_lines:

            if line.endswith("\n"):

                f.write(
                    line
                )

            else:

                f.write(
                    line + "\n"
                )


    # =========================================================================
    # INFORMACJA
    # =========================================================================

    print()
    print("============================================")
    print("PRZESORTOWANO PLIK")
    print("============================================")

    print(
        f"Plik wejściowy:\n{filename}"
    )

    print()

    print(
        f"Liczba punktów: {len(data_lines)}"
    )

    print()

    print(
        "Kolejność: g -> omega_0"
        if sort_by_g_then_omega
        else "Kolejność: omega_0 -> g"
    )

    print()

    print(
        f"Plik wynikowy:\n{output_filename}"
    )


# =============================================================================
# PĘTLA PO PLIKACH
# =============================================================================

for filename in filenames:

    reorder_file(
        filename
    )


print()
print("============================================")
print("GOTOWE")
print("============================================")