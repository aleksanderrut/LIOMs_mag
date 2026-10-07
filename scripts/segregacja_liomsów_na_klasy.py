from pathlib import Path

FOLDER = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026"
)

PLIK_BAZY = FOLDER / "bazy" / "operators_ladder_mag_M_4_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_0.5_g_0.5_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true_klas.txt"

PLIK_LIOM = FOLDER / "siatka_omega_g_liomsy" / "lioms_grid_M4_Jp0.0_d0.001_d20.0_Tboth_Pboth_Fyes_Bboth_FIdyes_eig3.txt"

PLIK_WYJSCIOWY = PLIK_LIOM.with_name(PLIK_LIOM.stem + "_kompozycja.data")

# ============================================================
# WCZYTANIE KLASYFIKACJI BAZY

klasyfikacja = {}

with open(PLIK_BAZY, "r") as f:
    for linia in f:
        kolumny = linia.split()

        numer = int(kolumny[0])
        I_n = kolumny[3]
        I_0 = kolumny[4]
        typ = kolumny[5]

        klasyfikacja[numer] = {
            "I_n": I_n,
            "I_0": I_0,
            "typ": typ
        }


# ============================================================
# WSZYSTKIE MOŻLIWE KLASY

klasy_I_n = sorted(set(x["I_n"] for x in klasyfikacja.values()), key=int)

klasy_I_0 = sorted(
    set(x["I_0"] for x in klasyfikacja.values() if x["I_0"] != "n"),
    key=int
)

if any(x["I_0"] == "n" for x in klasyfikacja.values()):
    klasy_I_0.append("n")

klasy_typ = ["f", "b", "m"]


# ============================================================
# ANALIZA PLIKU LIOM

wyniki = []

with open(PLIK_LIOM, "r") as f:
    for linia in f:
        kolumny = linia.strip().split()

        omega_0 = float(kolumny[0])
        g = float(kolumny[1])

        suma_I_n = {klasa: 0.0 for klasa in klasy_I_n}
        suma_I_0 = {klasa: 0.0 for klasa in klasy_I_0}
        suma_typ = {klasa: 0.0 for klasa in klasy_typ}

        elementy = kolumny[2].split(";")

        for element in elementy:
            numer, wspolczynnik = element.split(":")

            numer = int(numer)
            wspolczynnik = float(wspolczynnik)

            I_n = klasyfikacja[numer]["I_n"]
            I_0 = klasyfikacja[numer]["I_0"]
            typ = klasyfikacja[numer]["typ"]

            suma_I_n[I_n] += wspolczynnik
            suma_I_0[I_0] += wspolczynnik
            suma_typ[typ] += wspolczynnik

        wyniki.append((omega_0, g, suma_I_n, suma_I_0, suma_typ))


# ============================================================
# ZAPIS

with open(PLIK_WYJSCIOWY, "w") as f:

    naglowek = ["omega_0", "g"]

    naglowek += [f"I_n={x}" for x in klasy_I_n]
    naglowek += [f"I_0={x}" for x in klasy_I_0]
    naglowek += [f"typ={x}" for x in klasy_typ]

    f.write("\t".join(naglowek) + "\n")

    for omega_0, g, suma_I_n, suma_I_0, suma_typ in wyniki:

        wiersz = [
            f"{omega_0:.16e}",
            f"{g:.16e}"
        ]

        wiersz += [f"{suma_I_n[x]:.16e}" for x in klasy_I_n]
        wiersz += [f"{suma_I_0[x]:.16e}" for x in klasy_I_0]
        wiersz += [f"{suma_typ[x]:.16e}" for x in klasy_typ]

        f.write("\t".join(wiersz) + "\n")

print("Gotowe.")
print(f"Wczytano klasyfikację: {PLIK_BAZY}")
print(f"Wczytano LIOM-y:        {PLIK_LIOM}")
print(f"Zapisano:               {PLIK_WYJSCIOWY}")