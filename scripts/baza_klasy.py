from pathlib import Path

FOLDER_BAZY = Path(
    r"C:\Users\aleks\Desktop\praca magisterska\dane_serwer\set_7_03.09.2026\bazy"
)

# Tutaj wpisz nazwę pliku wejściowego
NAZWA_PLIKU = "operators_ladder_mag_M_4_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_0.5_g_0.5_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true.txt"
PLIK_WEJSCIOWY = FOLDER_BAZY / NAZWA_PLIKU

# Dodanie "_klas" przed rozszerzeniem .data
PLIK_WYJSCIOWY = PLIK_WEJSCIOWY.with_name(PLIK_WEJSCIOWY.stem + "_klas" + PLIK_WEJSCIOWY.suffix)


# ============================================================
# I_n - {A_i, |#S^{+}_{i,2} - S^{-}_{i,2}| = n}
def policz_I_n(stan):
    noga_bozonowa = stan.split("|")[1] # po znaku | jest noga bozonowa

    liczba_S_plus = noga_bozonowa.count("1")
    liczba_S_minus = noga_bozonowa.count("3")
    I_n = abs(liczba_S_plus - liczba_S_minus)
    return I_n

# ============================================================
# I_n = 0: |Σ pozycje S^{+}_{i,2} - Σ pozycje S^{-}_{i,2}|
def policz_P_n(stan, I_n):
    if I_n != 0:
        return "n"
    noga_bozonowa = stan.split("|")[1] # po znaku | jest noga bozonowa

    pozycje_S_plus = [i for i, op in enumerate(noga_bozonowa, start=1) if op == "1"]
    pozycje_S_minus = [i for i, op in enumerate(noga_bozonowa, start=1) if op == "3"]
    P_n = abs(sum(pozycje_S_plus) - sum(pozycje_S_minus))
    return P_n

# ============================================================
# typ operatora: b - bozonowy, f - fermionowy, m - mieszany

def policz_typ(stan):
    noga_fermionowa = stan.split("|")[0]
    noga_bozonowa = stan.split("|")[1]

    if all(op == "0" for op in noga_fermionowa):
        return "b"
    elif all(op == "0" for op in noga_bozonowa):
        return "f"
    else:
        return "m"

# ============================================================
# WCZYTANIE BAZY

wyniki = []

with open(PLIK_WEJSCIOWY, "r") as f:

    for numer_wiersza, linia in enumerate(f, start=1):

        linia = linia.strip()
        if not linia:
            continue

        kolumny = linia.split()
        # Format:
        # 1    RR    0000|2000
        sektor = kolumny[1]
        stan = kolumny[2]
        I_n = policz_I_n(stan)
        P_n = policz_P_n(stan, I_n)
        typ = policz_typ(stan)
        wyniki.append((numer_wiersza, sektor, stan, I_n, P_n, typ))


# ============================================================
# ZAPIS

with open(PLIK_WYJSCIOWY, "w") as f:
    for numer, sektor, stan, I_n, P_n, typ in wyniki:
        f.write(f"{numer}\t{sektor}\t{stan}\t{I_n}\t{P_n}\t{typ}\n")


print("Gotowe.")
print(f"Wczytano: {PLIK_WEJSCIOWY}")
print(f"Zapisano:  {PLIK_WYJSCIOWY}")
print(f"Liczba elementów bazy: {len(wyniki)}")