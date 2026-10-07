from pathlib import Path

# =============================================================================
# PLIKI DO PRZESORTOWANIA

# Obsługiwane formaty:
# grid:       omega_0    g    lambda
# lioms_grid: omega_0    g    3:...;4:...;...

filenames = [
    # Path(r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_6_25.08.2026\siatka_omega_g\grid_M3_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig1.txt"),

    Path(
        r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\siatka_omega_g_liomsy\lioms_grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"
    ),
]

# =============================================================================
# SPOSÓB SORTOWANIA

# True  -> najpierw g, potem omega_0
# False -> najpierw omega_0, potem g
sort_by_g_then_omega = True

# =============================================================================
# FUNKCJA SORTUJĄCA

def reorder_file(filename: Path):
    if not filename.exists():
        raise FileNotFoundError(f"Nie znaleziono pliku:\n{filename}")

    data_lines = []
    other_lines = []

    with open(filename, "r", encoding="utf-8") as f:
        for line_number, line in enumerate(f, start=1):
            stripped = line.strip()

            if not stripped:
                continue

            columns = stripped.split()

            if len(columns) < 2:
                other_lines.append(line)
                continue

            try:
                omega_0 = float(columns[0])
                g = float(columns[1])
            except ValueError:
                other_lines.append(line)
                continue

            data_lines.append((omega_0, g, line_number, line))

    if not data_lines:
        raise ValueError(f"Nie znaleziono danych numerycznych w pliku:\n{filename}")

    # x[0] = omega_0, x[1] = g
    if sort_by_g_then_omega:
        data_lines.sort(key=lambda x: (x[1], x[0]))
    else:
        data_lines.sort(key=lambda x: (x[0], x[1]))

    output_filename = filename.with_name(filename.stem + "_reorder" + filename.suffix)

    with open(output_filename, "w", encoding="utf-8") as f:
        for line in other_lines:
            f.write(line if line.endswith("\n") else line + "\n")

        for _, _, _, line in data_lines:
            f.write(line if line.endswith("\n") else line + "\n")

    kolejnosc = "g -> omega_0" if sort_by_g_then_omega else "omega_0 -> g"

    print("\n============================================")
    print("PRZESORTOWANO PLIK")
    print("============================================")
    print(f"Plik wejściowy:\n{filename}")
    print(f"\nLiczba punktów: {len(data_lines)}")
    print(f"\nKolejność: {kolejnosc}")
    print(f"\nPlik wynikowy:\n{output_filename}")

# =============================================================================
# PĘTLA PO PLIKACH

for filename in filenames:
    reorder_file(filename)

print("\n============================================")
print("GOTOWE")
print("============================================")