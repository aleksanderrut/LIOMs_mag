
using Pkg
using ArgParse
using SparseArrays
using LinearAlgebra
using DelimitedFiles


# ============================================================
# GENEROWANIE BAZY XXZ + HARD-CORE BOZONY
# ============================================================

function gen_baza_xxz(L::Int; mag=:all)

    rozmiar_pelnej_bazy = 1 << L

    if mag == :all
        liczba_stanow = rozmiar_pelnej_bazy
        sektory = 0:L

    elseif mag isa Real
        N_up = L / 2 + mag
        N_up = Int(N_up)
        liczba_stanow = binomial(L, N_up)
        sektory = N_up:N_up

    else
        error("Nieprawidlowa wartosc mag: $mag")
    end

    # Uwzglednienie drugiej nogi
    liczba_stanow *= rozmiar_pelnej_bazy

    baza_int = Vector{Tuple{Int, Int}}(
        undef, liczba_stanow
    )

    i = 1

    for N_up_sektor in sektory
        for state in 0:(rozmiar_pelnej_bazy - 1)

            if count_ones(state) == N_up_sektor

                for state_bos in 0:(rozmiar_pelnej_bazy - 1)
                    baza_int[i] = (state, state_bos)
                    i += 1
                end

            end
        end
    end

    indeks_odwrotny = Dict{Tuple{Int, Int}, Int}()

    for (i, state) in enumerate(baza_int)
        indeks_odwrotny[state] = i
    end

    println("LICZBA STANOW: ", liczba_stanow)

    baza_bit = [
        (string(state[1], base=2, pad=L),
         string(state[2], base=2, pad=L))
        for state in baza_int
    ]

    return baza_int, baza_bit, indeks_odwrotny
end


# ============================================================
# GENEROWANIE HAMILTONIANU I OPERATOROW
# ============================================================

function gen_ham_xxz(L::Int, J::Real, Delta::Real, PBC::Bool,
    baza_int::Vector{Tuple{Int, Int}},
    indeks_odwrotny::Dict{Tuple{Int, Int}, Int},
    omega0::Real, g::Real)

    ile_stanow = length(baza_int)

    H = zeros(Float64, ile_stanow, ile_stanow)

    Q_1 = zeros(Float64, ile_stanow, ile_stanow)
    Q_2 = zeros(Float64, ile_stanow, ile_stanow)
    Q_3 = zeros(Float64, ile_stanow, ile_stanow)
    Q_4 = zeros(Float64, ile_stanow, ile_stanow)

    # Liczba polaczen
    liczba_polaczen = L - 1

    if PBC == true && L > 2
        liczba_polaczen += 1
    end

    # Petla po stanach
    for i in 1:length(baza_int)

        state, state_bos = baza_int[i]

        # Petla po polaczeniach
        for bit_s in 0:(liczba_polaczen - 1)

            if bit_s == L - 1
                bit_t = 0
            else
                bit_t = bit_s + 1
            end

            # =====================================================
            # HAMILTONIAN XXZ - PIERWSZA NOGA
            # =====================================================

            spin_s = (state >> bit_s) & 1
            spin_t = (state >> bit_t) & 1

            # CZESC DIAGONALNA
            if spin_s == spin_t
                H[i, i] += (J * Delta) / 4
            else
                H[i, i] -= (J * Delta) / 4
            end

            # CZESC POZADIAGONALNA
            if spin_s != spin_t

                new_state = state ⊻ (
                    (1 << bit_s) | (1 << bit_t)
                )

                new_index = indeks_odwrotny[
                    (new_state, state_bos)
                ]

                H[new_index, i] += J / 2
            end

            # =====================================================
            # OPERATOR Q_1
            # 1/2 * (S+_i,2 S-_{i+1,2} + h.c.)
            # =====================================================

            bos_s = (state_bos >> bit_s) & 1
            bos_t = (state_bos >> bit_t) & 1

            if bos_s != bos_t

                new_state_bos = state_bos ⊻ (
                    (1 << bit_s) | (1 << bit_t)
                )

                new_index = indeks_odwrotny[
                    (state, new_state_bos)
                ]

                Q_1[new_index, i] += 0.5
            end

            # =====================================================
            # OPERATOR Q_3
            # Sz_i,2 * Sz_{i+1,2}
            # =====================================================

            if bos_s == bos_t
                Q_3[i, i] += 0.25
            else
                Q_3[i, i] -= 0.25
            end

        end

        # =====================================================
        # CZLON DIAGONALNY BOZONOW
        # omega0 * Sz_l,2
        # =====================================================

        H[i, i] += omega0 * (
            count_ones(state_bos) - L / 2
        )

        # =====================================================
        # SPRZEZENIE SPIN-BOZON
        # g * Sz_l,1 * (S+_l,2 + S-_l,2)
        # =====================================================

        for bit_s in 0:(L - 1)

            spin_s = (state >> bit_s) & 1

            sz_s = spin_s - 0.5

            # Kreacja lub anihilacja bozonu
            new_state_bos = state_bos ⊻ (1 << bit_s)

            new_index = indeks_odwrotny[
                (state, new_state_bos)
            ]

            H[new_index, i] += g * sz_s

            # =================================================
            # OPERATOR Q_2
            # S+_i,2 + S-_i,2
            # =================================================

            Q_2[new_index, i] += 1.0

            # =================================================
            # OPERATOR Q_4
            # 4 * Sz_l,1 * Sz_l,2
            # =================================================

            spin_bos = (state_bos >> bit_s) & 1

            sz_boson = spin_bos - 0.5

            Q_4[i, i] += 4.0 * sz_s * sz_boson

        end

    end

    return H, Q_1, Q_2, Q_3, Q_4
end


# ============================================================
# NORMALIZACJA OPERATORA
# ============================================================

function normalizuj_operator(Q)

    # Wymiar przestrzeni Hilberta
    D = size(Q, 1)

    # Norma Hilberta-Schmidta
    norma_Q = sqrt(sum(abs2, Q) / D)

    if norma_Q == 0
        error("Operator ma zerowa norme")
    end

    # Operator o normie jednostkowej
    return Q / norma_Q

end


# ============================================================
# OBLICZANIE KOMUTATORA
# UWAGA: TA FUNKCJA NIE NORMALIZUJE OPERATORA
# ============================================================

function oblicz_komutator(H, Q)

    # Komutator [H,Q]
    C = H * Q - Q * H

    # Obliczanie Tr(C†C)
    wynik = 0.0

    for i in 1:size(C, 2)
        for j in 1:size(C, 1)

            wynik += abs2(C[j, i])

        end
    end

    return C, wynik / size(H, 1)
end


# ============================================================
# ARGUMENTY PROGRAMU
# ============================================================

function parse_args()

    s = ArgParseSettings(
        description = "Generate XXZ Hamiltonian with hard-core bosons"
    )

    @add_arg_table! s begin

        "-M", "--sites"
            help = "Dlugosc lancucha"
            arg_type = Int
            default = 6
            dest_name = "M"

        "-N", "--n-up"
            help = "Liczba spinow w gore"
            arg_type = Int
            default = 3
            dest_name = "N_up"

        "-J"
            help = "Stala wymiany J"
            arg_type = Float64
            default = 1.0
            dest_name = "J"

        "-d", "--delta"
            help = "Parametr anizotropii Delta"
            arg_type = Float64
            default = 0.3
            dest_name = "Delta"

        "-w", "--omega0"
            help = "Energia bozonow"
            arg_type = Float64
            default = 1.0
            dest_name = "omega0"

        "-P", "--pbc"
            help = "Warunki brzegowe: yes albo no"
            arg_type = String
            required = true
            range_tester = x -> x in ["yes", "no"]
            dest_name = "pbc"

    end

    return ArgParse.parse_args(s)
end


# ============================================================
# ZAPIS DANYCH
# ============================================================

function save_data(data, filename, PBC; kwargs...)

    folder = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"

    mkpath(folder)

    filename_full = "$(filename)_PBC$(PBC)"

    for (key, value) in kwargs
        filename_full *= "_$(key)$(value)"
    end

    filepath = joinpath(folder, filename_full * ".data")

    writedlm(filepath, data)

    println("Zapisano do pliku: ", filepath)

end


# ============================================================
# MAIN
# ============================================================

function main()

    # ========================================================
    # ODCZYT PARAMETROW
    # ========================================================

    args = parse_args()

    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    omega0 = args["omega0"]
    PBC = args["pbc"] == "yes"

    # Zakres sprzezenia g
    G_MIN = 0.0
    G_MAX = 10.0
    G_POINTS = 20

    g_values = range(
        G_MIN, G_MAX, length=G_POINTS
    )

    # ========================================================
    # WSPOLCZYNNIKI Z PAULISTRINGS DLA g = 10
    # ========================================================

    W_3 = 9.9747961865030332e-01
    W_50 = 2.4936031094430042e-03

    # Znak wzgledny skladnikow
    ZNAK_Q4 = +1.0

    a = sqrt(W_3)
    b = ZNAK_Q4 * sqrt(W_50)

    # Magnetyzacja pierwszej nogi
    mag = N_up - M / 2

    # ========================================================
    # GENEROWANIE BAZY
    # ========================================================

    println("Generowanie bazy...")

    baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(
        M; mag=mag
    )

    # ========================================================
    # ZAPIS BAZY
    # ========================================================

    baza_data = Matrix{Any}(
        undef, length(baza_int), 5
    )

    for i in 1:length(baza_int)

        baza_data[i, 1] = i
        baza_data[i, 2] = baza_int[i][1]
        baza_data[i, 3] = baza_int[i][2]
        baza_data[i, 4] = baza_bit[i][1]
        baza_data[i, 5] = baza_bit[i][2]

    end

    save_data(
        baza_data, "Baza", PBC;
        L_x=M, N_up, J, Delta, omega0
    )

    # ========================================================
    # PRZYGOTOWANIE TABELI WYNIKOW
    # ========================================================

    # Kolumny: g, Q_1, Q_2, Q_3, Q_4, Q_slow

    wyniki = Matrix{Any}(undef, G_POINTS + 1, 6)

    wyniki[1, :] = [
        "g", "Q_1", "Q_2", "Q_3", "Q_4", "Q_slow"
    ]

    # ========================================================
    # PETLA PO g
    # ========================================================

    for (k, g) in enumerate(g_values)

        println("========================================")
        println("Punkt: ", k, "/", G_POINTS)
        println("g = ", g)

        # ====================================================
        # GENEROWANIE HAMILTONIANU OD ZERA
        # ====================================================

        H, Q_1, Q_2, Q_3, Q_4 = gen_ham_xxz(
            M, J, Delta, PBC,
            baza_int, indeks_odwrotny,
            omega0, g
        )

        
        # ====================================================
        # SPRAWDZENIE NORM OPERATORÓW PRZED NORMALIZACJĄ
        # ====================================================

        D = size(H, 1)

        println("Tr(Q_1†Q_1)/D = ", sum(abs2, Q_1) / D)
        println("Tr(Q_2†Q_2)/D = ", sum(abs2, Q_2) / D)
        println("Tr(Q_3†Q_3)/D = ", sum(abs2, Q_3) / D)
        println("Tr(Q_4†Q_4)/D = ", sum(abs2, Q_4) / D)


        # ====================================================
        # JEDNORAZOWA NORMALIZACJA OPERATOROW BAZOWYCH
        # ====================================================

        Q_1_norm = normalizuj_operator(Q_1)
        Q_2_norm = normalizuj_operator(Q_2)
        Q_3_norm = normalizuj_operator(Q_3)
        Q_4_norm = normalizuj_operator(Q_4)

        # ====================================================
        # OPERATOR WOLNEGO MODU
        # WSPOLCZYNNIKI Z PAULISTRINGS
        # ====================================================

        Q_slow = a * Q_2_norm + b * Q_4_norm

        # NIE NORMALIZUJEMY PONOWNIE Q_slow

        # ====================================================
        # OBLICZANIE KOMUTATOROW
        # ====================================================

        C_1, wynik_1 = oblicz_komutator(H, Q_1_norm)

        C_2, wynik_2 = oblicz_komutator(H, Q_2_norm)

        C_3, wynik_3 = oblicz_komutator(H, Q_3_norm)

        C_4, wynik_4 = oblicz_komutator(H, Q_4_norm)

        C_slow, wynik_slow = oblicz_komutator(
            H, Q_slow
        )

        # ====================================================
        # ZAPIS WYNIKOW DO TABELI
        # ====================================================

        wyniki[k + 1, :] = [
            g,
            wynik_1,
            wynik_2,
            wynik_3,
            wynik_4,
            wynik_slow
        ]

        println("Q_1 = ", wynik_1)
        println("Q_2 = ", wynik_2)
        println("Q_3 = ", wynik_3)
        println("Q_4 = ", wynik_4)
        println("Q_slow = ", wynik_slow)

    end

    # ========================================================
    # ZAPIS WYNIKOW
    # ========================================================

    save_data(
        wyniki,
        "Komutatory_vs_g_wagi_PS_bez_drugiej_normalizacji",
        PBC;
        L_x=M, N_up, J, Delta, omega0,
        znak_Q4=ZNAK_Q4
    )

    println("========================================")
    println("Zakonczono obliczenia dla wszystkich g.")

end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
