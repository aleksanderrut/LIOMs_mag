using ArgParse
using LinearAlgebra
using DelimitedFiles

# GENEROWANIE BAZY: FERMIONY + HARD-CORE BOZONY
function gen_baza_xxz(L::Int; mag=:all)
    rozmiar_pelnej_bazy = 1 << L
    #noga fermionowa
    if mag == :all
        liczba_stanow = rozmiar_pelnej_bazy
        sektory = 0:L
    elseif mag isa Real
        N_up = Int(L / 2 + mag)
        liczba_stanow = binomial(L, N_up)
        sektory = N_up:N_up
    end
    #noga fermionowa * kombinacje na nodze boznowej
    liczba_stanow = liczba_stanow * rozmiar_pelnej_bazy
    #Tuple na elementy bazy
    baza_int = Vector{Tuple{Int, Int}}(undef, liczba_stanow)

    i = 1 # iterator
    for N_up_sektor in sektory
        for state in 0:(rozmiar_pelnej_bazy - 1)#po wszystkich kobinacjach fermionowych
            if count_ones(state) == N_up_sektor
                for state_bos in 0:(rozmiar_pelnej_bazy - 1)#po wszystkich kombinacjach bozonowych 
                    baza_int[i] = (state, state_bos)
                    i += 1
                end
            end
        end
    end

    indeks_odwrotny = Dict{Tuple{Int, Int}, Int}() #Słownik na indeksy odwrotne
    for (i, state) in enumerate(baza_int)
        indeks_odwrotny[state] = i
    end

    #Zapis w zapisie binarnym bazy
    baza_bit = Vector{Tuple{String, String}}(undef, liczba_stanow)
    for i in 1:liczba_stanow
        state = baza_int[i]
        state_ferm = state[1]
        state_bos = state[2]
        ferm_bit = string(state_ferm, base=2, pad=L)
        bos_bit = string(state_bos, base=2, pad=L)
        baza_bit[i] = (ferm_bit, bos_bit)
    end
    #println("LICZBA STANOW: ", liczba_stanow)

    return baza_int, baza_bit, indeks_odwrotny
end

# GENEROWANIE HAMILTONIANU I OPERATOROW
# Q_3   = RR 000|300
# Q_4   = RR 000|220
# Q_8   = RR 000|130
# Q_12  = RR 000|330
# Q_50  = RR 200|200
# Q_62  = RR 200|230
# Q_118 = RR 020|320
# Q_231 = II 130|320
# Q_237 = II 130|230
# Q_1128 = RR 123|320
# Q_1170 = RR 123|023
# 0 = I, 1 = S+, 2 = Sz, 3 = S-
function gen_ham_xxz(L::Int, J::Real, Delta::Real, PBC::Bool, baza_int::Vector{Tuple{Int, Int}}, indeks_odwrotny::Dict{Tuple{Int, Int}, Int}, omega0::Real, g::Real)
    ile_stanow = length(baza_int)
    liczba_polaczen = L - 1
    if PBC == true && L > 2
        liczba_polaczen += 1
    end
    # ========================================================
    # CZESC I: GENEROWANIE TYLKO HAMILTONIANU H
    # ========================================================
    H = zeros(Float64, ile_stanow, ile_stanow)
    #pętla po stanach
    for i in 1:ile_stanow
        state, state_bos = baza_int[i]
        #pętla po fermionowej nodze
        for bit_s in 0:(liczba_polaczen - 1)
            if bit_s == L - 1
                bit_t = 0
            else
                bit_t = bit_s + 1
            end
            #zczytanie spinów
            spin_s = (state >> bit_s) & 1
            spin_t = (state >> bit_t) & 1
            #diagonala fermiony
            if spin_s == spin_t
                H[i, i] += (J * Delta) / 4
            else
                H[i, i] -= (J * Delta) / 4
            end
            #przeskok fermionowy
            if spin_s != spin_t
                new_state = state ⊻ ((1 << bit_s) | (1 << bit_t))
                new_index = indeks_odwrotny[(new_state, state_bos)]
                H[new_index, i] += J / 2
            end
        end
        #energia bozonowa
        H[i, i] += omega0 * (count_ones(state_bos) - L / 2)
        #sprzezenie za pomocą g
        for bit_s in 0:(L - 1)
            spin_s = (state >> bit_s) & 1
            sz_ferm_s = spin_s - 0.5

            new_state_bos = state_bos ⊻ (1 << bit_s)
            new_index = indeks_odwrotny[(state, new_state_bos)]

            H[new_index, i] += g * sz_ferm_s
        end
    end

    # ========================================================
    # CZESC II: GENEROWANIE TYLKO OPERATOROW Q
    # ========================================================
    if PBC == true
        liczba_trojek = L
    else
        liczba_trojek = L - 2
    end
    Q_3    = zeros(Float64, ile_stanow, ile_stanow)
    Q_4    = zeros(Float64, ile_stanow, ile_stanow)
    Q_8    = zeros(Float64, ile_stanow, ile_stanow)
    Q_12   = zeros(Float64, ile_stanow, ile_stanow)
    Q_50   = zeros(Float64, ile_stanow, ile_stanow)
    Q_62   = zeros(Float64, ile_stanow, ile_stanow)
    Q_118  = zeros(Float64, ile_stanow, ile_stanow)
    Q_231  = zeros(Float64, ile_stanow, ile_stanow)
    Q_237  = zeros(Float64, ile_stanow, ile_stanow)
    Q_1128 = zeros(Float64, ile_stanow, ile_stanow)
    Q_1170 = zeros(Float64, ile_stanow, ile_stanow)

    for i in 1:ile_stanow
        state, state_bos = baza_int[i]
        # 1. OPERATORY JEDNOMIEJSCOWE Q_3 ORAZ Q_50
        for bit_s in 0:(L - 1)
            spin_s = (state >> bit_s) & 1
            bos_s = (state_bos >> bit_s) & 1
            sz_ferm_s = spin_s - 0.5
            sz_bos_s = bos_s - 0.5

            # Q_3 = RR 000|300 ; S+_s,2 + S-_s,2
            new_state_bos = state_bos ⊻ (1 << bit_s)
            new_index = indeks_odwrotny[(state, new_state_bos)]
            Q_3[new_index, i] += 1.0

            # Q_50 = RR 200|200 ; Sz_s,1 * Sz_s,2
            Q_50[i, i] += 2 * sz_ferm_s * sz_bos_s
        end

        # 2. OPERATORY DWUMIEJSCOWE Q_4, Q_8, Q_12, Q_62, Q_118, Q_231, Q_237
        for bit_s in 0:(liczba_polaczen - 1)
            if bit_s == L - 1
                bit_t = 0
            else
                bit_t = bit_s + 1
            end
            spin_s = (state >> bit_s) & 1
            spin_t = (state >> bit_t) & 1
            bos_s = (state_bos >> bit_s) & 1
            bos_t = (state_bos >> bit_t) & 1
            sz_ferm_s = spin_s - 0.5
            sz_ferm_t = spin_t - 0.5
            sz_bos_s = bos_s - 0.5
            sz_bos_t = bos_t - 0.5

            # Q_8 = RR 000|130 ; (S+_s,2*S-_t,2 + S-_s,2*S+_t,2)
            if bos_s != bos_t
                new_state_bos = state_bos ⊻ ((1 << bit_s) | (1 << bit_t))
                new_index = indeks_odwrotny[(state, new_state_bos)]
                Q_8[new_index, i] += 1.0
            end

            # Q_12 = RR 000|330 ; (S-_s,2*S-_t,2 + S+_s,2*S+_t,2)
            if bos_s == bos_t
                new_state_bos = state_bos ⊻ ((1 << bit_s) | (1 << bit_t))
                new_index = indeks_odwrotny[(state, new_state_bos)]
                Q_12[new_index, i] += 1.0
            end

            # Q_4 = RR 000|220 ; Sz_s,2 * Sz_t,2
            Q_4[i, i] += sz_bos_s * sz_bos_t

            # Q_62 = RR 200|230 ; Sz_s,1 * Sz_s,2 * (S+_t,2 + S-_t,2)
            flipped_bos_t = state_bos ⊻ (1 << bit_t)
            indeks_62 = indeks_odwrotny[(state, flipped_bos_t)]
            Q_62[indeks_62, i] += sz_ferm_s * sz_bos_s

            # Q_118 = RR 020|320 ; Sz_t,1 * (S+_s,2 + S-_s,2) * Sz_t,2
            flipped_bos_s = state_bos ⊻ (1 << bit_s)
            indeks_118 = indeks_odwrotny[(state, flipped_bos_s)]
            Q_118[indeks_118, i] += sz_ferm_t * sz_bos_t

            # Q_231 = II 130|320
            # Q_237 = II 130|230
            # Oba wymagaja przeskoku fermionu s <-> t.
            if spin_s != spin_t
                new_state = state ⊻ ((1 << bit_s) | (1 << bit_t))
                znak_ferm = spin_t - spin_s
                znak_bos_s = 2 * bos_s - 1
                znak_bos_t = 2 * bos_t - 1
                indeks_231 = indeks_odwrotny[(new_state, flipped_bos_s)]
                indeks_237 = indeks_odwrotny[(new_state, flipped_bos_t)]

                # Dwa czynniki i daja i*i = -1
                Q_231[indeks_231, i] += (-znak_ferm * znak_bos_s * sz_bos_t)
                Q_237[indeks_237, i] += (-znak_ferm * znak_bos_t * sz_bos_s)
            end
        end

        # 3. OPERATORY TRZYMIEJSCOWE Q_1128 ORAZ Q_1170
        for bit_s in 0:(liczba_trojek - 1)
            bit_t = (bit_s + 1) % L
            bit_u = (bit_s + 2) % L
            spin_s = (state >> bit_s) & 1
            spin_t = (state >> bit_t) & 1
            spin_u = (state >> bit_u) & 1

            if spin_s != spin_u
                new_state = state ⊻ ((1 << bit_s) | (1 << bit_u))
                sz_ferm_t = spin_t - 0.5
                bos_t = (state_bos >> bit_t) & 1
                sz_bos_t = bos_t - 0.5

                wspolczynnik = sz_ferm_t * sz_bos_t

                # Q_1128 = RR 123|320
                new_bos_1128 = state_bos ⊻ (1 << bit_s)
                indeks_1128 = indeks_odwrotny[(new_state, new_bos_1128)]
                Q_1128[indeks_1128, i] += wspolczynnik

                # Q_1170 = RR 123|023
                new_bos_1170 = state_bos ⊻ (1 << bit_u)
                indeks_1170 = indeks_odwrotny[(new_state, new_bos_1170)]
                Q_1170[indeks_1170, i] += wspolczynnik
            end
        end
    end
    return H, Q_3, Q_4, Q_8, Q_12, Q_50, Q_62, Q_118, Q_231, Q_237, Q_1128, Q_1170
end

# NORMALIZACJA OPERATORA
function normalizuj_operator(Q)
    D = size(Q, 1)
    norma_kwadrat = sum(abs2, Q) / D
    return Q / sqrt(norma_kwadrat)
end

# OBLICZANIE KOMUTATORA
function oblicz_komutator(H, Q)
    C = H * Q - Q * H
    D = size(H, 1)
    wynik = sum(abs2, C) / D
    return C, wynik
end

#ARGUMENTY
function parse_args()
    s = ArgParseSettings(description = "XXZ + hard-core bosons: PauliStrings commutators")
    @add_arg_table! s begin
        "-M", "--sites"
            help = "Dlugosc lancucha"
            arg_type = Int
            default = 3
            dest_name = "M"
        "-N", "--n-up"
            help = "Liczba spinow w gore"
            arg_type = Int
            default = 1
            dest_name = "N_up"
        "-J"
            help = "Stala J"
            arg_type = Float64
            default = 1.0
            dest_name = "J"
        "-d", "--delta"
            help = "Anizotropia Delta"
            arg_type = Float64
            default = 0.3
            dest_name = "Delta"
        "-w", "--omega0"
            help = "Czestosc omega0"
            arg_type = Float64
            default = 1.0
            dest_name = "omega0"
        "-P", "--pbc"
            help = "Warunki brzegowe yes lub no"
            arg_type = String
            required = true
            range_tester = x -> x in ["yes", "no"]
            dest_name = "pbc"
    end
    return ArgParse.parse_args(s)
end

#ZAPIS
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

# ZAPIS MACIERZY Tr([H,Q]†[H,Q])/D VS g
function save_wyniki_vs_g(wyniki, kolumny, PBC; kwargs...)
    folder = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
    mkpath(folder)
    filename_full = "Tr_komutatorow_vs_g_PBC$(PBC)"
    for (key, value) in kwargs
        filename_full *= "_$(key)$(value)"
    end
    filepath = joinpath(folder, filename_full * ".data")
    open(filepath, "w") do io
        println(io, join(kolumny, "\t"))
        writedlm(io, wyniki, '\t')
    end
    println("Zapisano macierz komutatorow do pliku: ", filepath)
    return filepath
end

# Losowy stan początkowy
function losuj_q1_gauss(n::Int)
    q1 = randn(n) .+ im .* randn(n)
    return q1 / norm(q1)
end

# Przygotowanie korelacji ED - transformacja operatora wykonywana tylko raz dla kazdego g
function przygotuj_korelacje_ED(wynik_H, A, Psi_0)
    V = wynik_H.vectors
    A_wlasna = adjoint(V) * A * V
    Psi_wlasna = adjoint(V) * Psi_0
    Phi_wlasna = A_wlasna * Psi_wlasna
    W = abs2.(A_wlasna)
    return (A_wlasna, Psi_wlasna, Phi_wlasna, W)
end

# Funkcja korelacji z wykorzystaniem ED
function funkcja_korelacji_ED(dane, fazy)
    A_wlasna, Psi_wlasna, Phi_wlasna, W = dane
    Psi_t = fazy .* Psi_wlasna
    Phi_t = fazy .* Phi_wlasna
    C_t = dot(Psi_t, A_wlasna * Phi_t)
    return real(C_t)
end

# Funkcja korelacji sladowej T -> inf z wykorzystaniem ED
function funkcja_korelacji_trace_ED(dane, fazy)
    A_wlasna, Psi_wlasna, Phi_wlasna, W = dane
    C_t = dot(fazy, W * fazy)
    return real(C_t) / length(fazy)
end

# ZAPIS FUNKCJI KORELACJI ED VS g I t
function save_korelacje_ED(wyniki_korelacji, kolumny_korelacji, PBC; nazwa="Korelacje_ED_vs_g_t", kwargs...)
    folder = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
    mkpath(folder)
    filename_full = "$(nazwa)_PBC$(PBC)"
    for (key, value) in kwargs
        filename_full *= "_$(key)$(value)"
    end
    filepath = joinpath(folder, filename_full * ".data")
    open(filepath, "w") do io
        println(io, join(kolumny_korelacji, "\t"))
        writedlm(io, wyniki_korelacji, '\t')
    end
    println("Zapisano funkcje korelacji ED do pliku: ", filepath)
    return filepath
end

# =========================================================================================================
function main()
    args = parse_args()
    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    omega0 = args["omega0"]
    PBC = args["pbc"] == "yes"

    g_MIN = 0.0
    g_MAX = 10.0
    g_POINTS = 41
    g_values = range(g_MIN, g_MAX, length=g_POINTS)

    # PIERWSZY WOLNY MOD ; INDEKSY 3 I 50
    W_3  = 9.9747961865030332e-01
    W_50 = 2.4936031094430042e-03
    a_3  = sqrt(W_3)
    a_50 = sqrt(W_50)

    # DRUGI WOLNY MOD ; INDEKSY 8, 12, 118 I 62
    W_8   = 4.9786722486762175e-01
    W_12  = 4.9710224964245658e-01
    W_118 = 2.4853980395447483e-03
    W_62  = 2.4853980395447436e-03
    a_8   =  sqrt(W_8)
    a_12  =  sqrt(W_12)
    a_118 =  sqrt(W_118)
    a_62  =  sqrt(W_62)

    # TRZECI WOLNY MOD; INDEKSY 8, 4, 231, 237, 62, 118, 1128, 1170
    W3_8   = 9.8480872701838784e-01
    W3_4   = 1.3452350434232842e-02
    W3_231 = 2.7210764810602230e-04
    W3_237 = 2.7210764810600853e-04
    W3_62  = 1.9838883970923013e-04
    W3_118 = 1.9838883970921172e-04
    W3_1128 = 1.7059543031483211e-04
    W3_1170 = 1.7059543031482842e-04
    a3_8 = sqrt(W3_8)
    a3_4 = sqrt(W3_4)
    a3_231 = sqrt(W3_231)
    a3_237 = sqrt(W3_237)
    a3_62 = sqrt(W3_62)
    a3_118 = sqrt(W3_118)
    a3_1128 = sqrt(W3_1128)
    a3_1170 = sqrt(W3_1170)

    # GENEROWANIE BAZY
    mag = :all
    N_up = "all"
    println("Generowanie pelnej bazy...")
    baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(M; mag=mag)
    println("Wymiar przestrzeni Hilberta: ", length(baza_int))
    baza_data = Matrix{Any}(undef, length(baza_int), 5)
    for i in 1:length(baza_int)
        baza_data[i, 1] = i
        baza_data[i, 2] = baza_int[i][1]
        baza_data[i, 3] = baza_int[i][2]
        baza_data[i, 4] = baza_bit[i][1]
        baza_data[i, 5] = baza_bit[i][2]
    end
    save_data(baza_data, "Baza", PBC; L_x=M, N_up, J, Delta, omega0)

    # Pierwsza kolumna: g. Wszystkie pozostale: Tr([H,Q]†[H,Q])/D.
    kolumny = [
        "g",
        "Q_3", "Q_4", "Q_8", "Q_12", "Q_50",
        "Q_62", "Q_118", "Q_231", "Q_237", "Q_1128", "Q_1170",
        "Q_slow_3_50",
        "Q_slow_8_12_62_118",
        "Q_slow_8_4_231_237_62_118_1128_1170"
    ]
    wyniki = Matrix{Float64}(undef, g_POINTS, length(kolumny))

    # PARAMETRY KORELACJI CZASOWYCH
    t_MIN = 0.0
    t_MAX = 100.0
    t_POINTS = 501
    t_values = range(t_MIN, t_MAX, length=t_POINTS)

    # Pierwsza kolumna: g. Wszystkie pozostałe: C_Q(t) = Re⟨Ψ₀|e^(iHt) Q e^(-iHt) Q|Ψ₀⟩
    kolumny_korelacji = [
        "g",
        "t",
        "C_slow_3_50",
        "C_slow_8_12_62_118",
        "C_slow_8_4_231_237_62_118_1128_1170"
    ]

    # Losowy stan początkowy — ten sam dla wszystkich g
    Psi_0 = losuj_q1_gauss(length(baza_int))

    # Macierz na wyniki funkcji korelacji
    wyniki_korelacji = Matrix{Float64}(undef, g_POINTS * t_POINTS, length(kolumny_korelacji))
    wyniki_korelacji_trace = Matrix{Float64}(undef, g_POINTS * t_POINTS, length(kolumny_korelacji))

    #pętla po g =============================================================================================
    for (k, g) in enumerate(g_values)
        println("Punkt: ", k, "/", g_POINTS)

        #generowane Ham i wszystkich operatorów
        H, Q_3, Q_4, Q_8, Q_12, Q_50, Q_62, Q_118, Q_231, Q_237, Q_1128, Q_1170 =
        gen_ham_xxz(M, J, Delta, PBC, baza_int, indeks_odwrotny, omega0, g)

        if k == 1
            D = size(H, 1)
            println("--- NORMY KWADRATOWE HS ---")
            println("Q_3   = ", sum(abs2, Q_3) / D)
            println("Q_4   = ", sum(abs2, Q_4) / D)
            println("Q_8   = ", sum(abs2, Q_8) / D)
            println("Q_12  = ", sum(abs2, Q_12) / D)
            println("Q_50  = ", sum(abs2, Q_50) / D)
            println("Q_62  = ", sum(abs2, Q_62) / D)
            println("Q_118 = ", sum(abs2, Q_118) / D)
            println("Q_231 = ", sum(abs2, Q_231) / D)
            println("Q_237 = ", sum(abs2, Q_237) / D)
            println("Q_1128 = ", sum(abs2, Q_1128) / D)
            println("Q_1170 = ", sum(abs2, Q_1170) / D)
        end

        Q_3   = normalizuj_operator(Q_3)
        Q_4   = normalizuj_operator(Q_4)
        Q_8   = normalizuj_operator(Q_8)
        Q_12  = normalizuj_operator(Q_12)
        Q_50  = normalizuj_operator(Q_50)
        Q_62  = normalizuj_operator(Q_62)
        Q_118 = normalizuj_operator(Q_118)
        Q_231 = normalizuj_operator(Q_231)
        Q_237 = normalizuj_operator(Q_237)
        Q_1128 = normalizuj_operator(Q_1128)
        Q_1170 = normalizuj_operator(Q_1170)

        Q_slow_3_50 = +a_3 * Q_3 + a_50 * Q_50
        Q_slow_8_12_62_118 = (-a_8 * Q_8 -a_12 * Q_12 -a_118 * Q_118 -a_62 * Q_62)
        Q_slow_8_4_231_237_62_118_1128_1170 = (+a3_8 * Q_8 +a3_4 * Q_4 -a3_231 * Q_231 +a3_237 * Q_237 -a3_62 * Q_62 -a3_118 * Q_118 +a3_1128 * Q_1128 +a3_1170 * Q_1170)

        if k == 1
            D = size(H, 1)
            println("--- NORMY KWADRATOWE HS ---")
            println("Q_3   = ", sum(abs2, Q_3) / D)
            println("Q_4   = ", sum(abs2, Q_4) / D)
            println("Q_8   = ", sum(abs2, Q_8) / D)
            println("Q_12  = ", sum(abs2, Q_12) / D)
            println("Q_50  = ", sum(abs2, Q_50) / D)
            println("Q_62  = ", sum(abs2, Q_62) / D)
            println("Q_118 = ", sum(abs2, Q_118) / D)
            println("Q_231 = ", sum(abs2, Q_231) / D)
            println("Q_237 = ", sum(abs2, Q_237) / D)
            println("Q_1128 = ", sum(abs2, Q_1128) / D)
            println("Q_1170 = ", sum(abs2, Q_1170) / D)

            println("Q_slow_3_50 = ", sum(abs2, Q_slow_3_50) / D)
            println("Q_slow_8_12_62_118 = ", sum(abs2, Q_slow_8_12_62_118) / D)
            println("Q_slow_8_4_231_237_62_118_1128_1170 = ", sum(abs2, Q_slow_8_4_231_237_62_118_1128_1170) / D)

            println("Suma wag 3,50 = ",W_3 + W_50)
            println("Suma wag 8,12,62,118 = ",W_8 + W_12 + W_62 + W_118)
            println("Suma wag 8,4,231,237,62,118,1128,1170 = ",W3_8 + W3_4 + W3_231 + W3_237 + W3_62 + W3_118 + W3_1128 + W3_1170)
        end

        C_3, wynik_3 = oblicz_komutator(H, Q_3)
        C_4, wynik_4 = oblicz_komutator(H, Q_4)
        C_8, wynik_8 = oblicz_komutator(H, Q_8)
        C_12, wynik_12 = oblicz_komutator(H, Q_12)
        C_50, wynik_50 = oblicz_komutator(H, Q_50)
        C_62, wynik_62 = oblicz_komutator(H, Q_62)
        C_118, wynik_118 = oblicz_komutator(H, Q_118)
        C_231, wynik_231 = oblicz_komutator(H, Q_231)
        C_237, wynik_237 = oblicz_komutator(H, Q_237)
        C_1128, wynik_1128 = oblicz_komutator(H, Q_1128)
        C_1170, wynik_1170 = oblicz_komutator(H, Q_1170)
        C_slow_1, wynik_slow_1 = oblicz_komutator(H, Q_slow_3_50)
        C_slow_2, wynik_slow_2 = oblicz_komutator(H, Q_slow_8_12_62_118)
        C_slow_3, wynik_slow_3 = oblicz_komutator(H, Q_slow_8_4_231_237_62_118_1128_1170)

        wyniki[k, :] = [g,wynik_3,wynik_4,wynik_8,wynik_12,wynik_50,wynik_62,wynik_118,wynik_231,wynik_237,wynik_1128,wynik_1170,wynik_slow_1,wynik_slow_2,wynik_slow_3]

        # println("Q_3 = ", wynik_3)
        # println("Q_4 = ", wynik_4)
        # println("Q_8 = ", wynik_8)
        # println("Q_12 = ", wynik_12)
        # println("Q_50 = ", wynik_50)
        # println("Q_62 = ", wynik_62)
        # println("Q_118 = ", wynik_118)
        # println("Q_231 = ", wynik_231)
        # println("Q_237 = ", wynik_237)
        # println("Q_1128 = ", wynik_1128)
        # println("Q_1170 = ", wynik_1170)
        # println("Q_slow_3_50 = ", wynik_slow_1)
        # println("Q_slow_8_12_62_118 = ", wynik_slow_2)
        # println("Q_slow_8_4_231_237_62_118_1128_1170 = ", wynik_slow_3)

        operatory_korelacji = [Q_slow_3_50, Q_slow_8_12_62_118, Q_slow_8_4_231_237_62_118_1128_1170] #WYBOR OPERATOROW
        wynik_H = eigen(Symmetric(H)) # DIAGONALIZACJA ED - RAZ DLA KAZDEGO g

        dane_korelacji = Vector{Any}(undef, length(operatory_korelacji))
        for a in 1:length(operatory_korelacji)
            A = operatory_korelacji[a]
            dane_korelacji[a] = przygotuj_korelacje_ED(wynik_H, A, Psi_0)
        end

        # PETLA PO CZASIE
        for (j, t) in enumerate(t_values)
            fazy = exp.(-im .* wynik_H.values .* t)
            wiersz = (k - 1) * t_POINTS + j
            wyniki_korelacji[wiersz, 1] = g
            wyniki_korelacji[wiersz, 2] = t
            wyniki_korelacji_trace[wiersz, 1] = g
            wyniki_korelacji_trace[wiersz, 2] = t

            for a in 1:length(operatory_korelacji)
                dane = dane_korelacji[a]

                C_t = funkcja_korelacji_ED(dane, fazy)
                wyniki_korelacji[wiersz, a + 2] = C_t

                C_trace_t = funkcja_korelacji_trace_ED(dane, fazy)
                wyniki_korelacji_trace[wiersz, a + 2] = C_trace_t
            end
        end
    end
    save_wyniki_vs_g(wyniki, kolumny, PBC;L_x=M, N_up, J, Delta, omega0,gmin=g_MIN, gmax=g_MAX, ng=g_POINTS)
    save_korelacje_ED(wyniki_korelacji, kolumny_korelacji, PBC; nazwa="Korelacje_ED_vs_g_t", L_x=M, N_up, J, Delta, omega0,gmin=g_MIN, gmax=g_MAX, ng=g_POINTS,tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
    save_korelacje_ED(wyniki_korelacji_trace, kolumny_korelacji, PBC; nazwa="Korelacje_trace_ED_vs_g_t", L_x=M, N_up, J, Delta, omega0, gmin=g_MIN, gmax=g_MAX, ng=g_POINTS, tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
