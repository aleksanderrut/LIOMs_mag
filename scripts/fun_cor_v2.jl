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
# Q_1   = RR 000|200
# Q_3   = RR 000|300
# Q_4   = RR 000|220
# Q_8   = RR 000|130
# Q_12  = RR 000|330
# Q_26  = RR 000|103
# Q_50  = RR 200|200
# Q_60  = RR 200|130
# Q_62  = RR 200|230
# Q_118 = RR 020|320
# Q_231 = II 130|320
# Q_237 = II 130|230
# Q_289 = RR 130|000
# Q_300 = RR 130|130
# Q_304 = RR 130|330
# Q_1128 = RR 123|320
# Q_1170 = RR 123|023
# Q_6    = RR 000|320
# Q_10   = RR 000|230
# Q_235  = II 130|130
# Q_292  = RR 130|300
# Q_298  = RR 130|030
# Q_1121 = RR 123|000
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

            H[new_index, i] += g * (sz_ferm_s + 0.5)
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
    Q_1    = zeros(Float64, ile_stanow, ile_stanow)
    Q_3    = zeros(Float64, ile_stanow, ile_stanow)
    Q_4    = zeros(Float64, ile_stanow, ile_stanow)
    Q_8    = zeros(Float64, ile_stanow, ile_stanow)
    Q_12   = zeros(Float64, ile_stanow, ile_stanow)
    Q_26   = zeros(Float64, ile_stanow, ile_stanow)
    Q_50   = zeros(Float64, ile_stanow, ile_stanow)
    Q_60   = zeros(Float64, ile_stanow, ile_stanow)
    Q_62   = zeros(Float64, ile_stanow, ile_stanow)
    Q_118  = zeros(Float64, ile_stanow, ile_stanow)
    Q_231  = zeros(Float64, ile_stanow, ile_stanow)
    Q_237  = zeros(Float64, ile_stanow, ile_stanow)
    Q_289  = zeros(Float64, ile_stanow, ile_stanow)
    Q_300  = zeros(Float64, ile_stanow, ile_stanow)
    Q_304  = zeros(Float64, ile_stanow, ile_stanow)
    Q_1128 = zeros(Float64, ile_stanow, ile_stanow)
    Q_1170 = zeros(Float64, ile_stanow, ile_stanow)

    # NOWE OPERATORY
    Q_6    = zeros(Float64, ile_stanow, ile_stanow)
    Q_10   = zeros(Float64, ile_stanow, ile_stanow)
    Q_235  = zeros(Float64, ile_stanow, ile_stanow)
    Q_292  = zeros(Float64, ile_stanow, ile_stanow) # RR 130|300
    Q_298  = zeros(Float64, ile_stanow, ile_stanow)
    Q_1121 = zeros(Float64, ile_stanow, ile_stanow)

    for i in 1:ile_stanow
        state, state_bos = baza_int[i]
        # 1. OPERATORY JEDNOMIEJSCOWE Q_1, Q_3 ORAZ Q_50
        for bit_s in 0:(L - 1)
            spin_s = (state >> bit_s) & 1
            bos_s = (state_bos >> bit_s) & 1
            sz_ferm_s = spin_s - 0.5
            sz_bos_s = bos_s - 0.5

            # Q_1 = RR 000|200 ; Sz_s,2
            Q_1[i, i] += sz_bos_s

            # Q_3 = RR 000|300 ; S+_s,2 + S-_s,2
            new_state_bos = state_bos ⊻ (1 << bit_s)
            new_index = indeks_odwrotny[(state, new_state_bos)]
            Q_3[new_index, i] += 1.0

            # Q_50 = RR 200|200 ; Sz_s,1 * Sz_s,2
            Q_50[i, i] += 2 * sz_ferm_s * sz_bos_s
        end

        # 2. OPERATORY DWUMIEJSCOWE Q_4, Q_8, Q_12, Q_60, Q_62, Q_118, Q_231, Q_237, Q_289, Q_300, Q_304
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
            # Q_60 = RR 200|130 ; Sz_s,1 * (S+_s,2*S-_t,2 + S-_s,2*S+_t,2)
            if bos_s != bos_t
                new_state_bos = state_bos ⊻ ((1 << bit_s) | (1 << bit_t))
                new_index = indeks_odwrotny[(state, new_state_bos)]
                Q_8[new_index, i] += 1.0
                Q_60[new_index, i] += sz_ferm_s
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

            # Q_6 = RR 000|320 ; (S-_s,2 + S+_s,2) * Sz_t,2
            # Przeskok bozonu na s, mnozenie przez Sz na t.
            Q_6[indeks_118, i] += sz_bos_t

            # Q_10 = RR 000|230 ; Sz_s,2 * (S-_t,2 + S+_t,2)
            # Przeskok bozonu na t, mnozenie przez Sz na s.
            Q_10[indeks_62, i] += sz_bos_s

            # Q_231 = II 130|320
            # Q_237 = II 130|230
            # Oba wymagaja przeskoku fermionu s <-> t.
            if spin_s != spin_t
                new_state = state ⊻ ((1 << bit_s) | (1 << bit_t))

                # Q_289 = RR 130|000
                # Fermionowy przeskok s <-> t, identycznosc na nodze bozonowej.
                indeks_289 = indeks_odwrotny[(new_state, state_bos)]
                Q_289[indeks_289, i] += 1.0

                znak_ferm = spin_t - spin_s
                znak_bos_s = 2 * bos_s - 1
                znak_bos_t = 2 * bos_t - 1
                indeks_231 = indeks_odwrotny[(new_state, flipped_bos_s)]
                indeks_237 = indeks_odwrotny[(new_state, flipped_bos_t)]

                # Dwa czynniki i daja i*i = -1
                Q_231[indeks_231, i] += (-znak_ferm * znak_bos_s * sz_bos_t)
                Q_237[indeks_237, i] += (-znak_ferm * znak_bos_t * sz_bos_s)

                # Q_235 = II 130|130
                # [i(A_f-A_f^dagger)]*[i(A_b-A_b^dagger)] = -(A_f-A_f^dagger)*(A_b-A_b^dagger)
                # Warunek: obie nogi maja na pozycjach s,t przeciwne spiny.
                if bos_s != bos_t
                    flipped_bos_st = state_bos ⊻ ((1 << bit_s) | (1 << bit_t))
                    indeks_235 = indeks_odwrotny[(new_state, flipped_bos_st)]
                    znak_bos_st = bos_t - bos_s
                    Q_235[indeks_235, i] += -znak_ferm * znak_bos_st

                    # Q_300 = RR 130|130
                    # Przeskok fermionu i bozonu na tej samej parze s,t.
                    Q_300[indeks_235, i] += 1.0
                else
                    # Q_304 = RR 130|330
                    # Przeskok fermionu oraz jednoczesne podniesienie/opuszczenie
                    # obu spinow bozonowych na pozycjach s,t.
                    flipped_bos_st = state_bos ⊻ ((1 << bit_s) | (1 << bit_t))
                    indeks_304 = indeks_odwrotny[(new_state, flipped_bos_st)]
                    Q_304[indeks_304, i] += 1.0
                end

                # Q_292 = RR 130|300
                # (S+_s,1*S-_t,1 + S-_s,1*S+_t,1) * (S-_s,2 + S+_s,2)
                # indeks_231 ma juz wlasciwe przelaczenie nog: ferm s<->t, bos s.
                Q_292[indeks_231, i] += 1.0

                # Q_298 = RR 130|030
                # Fermionowe R: przeskok s<->t; bozonowe R: S-_t,2 + S+_t,2.
                # indeks_237 ma juz wlasciwe przelaczenie nog: ferm s<->t, bos t.
                Q_298[indeks_237, i] += 1.0
            end
        end

        # 3. OPERATORY TRZYMIEJSCOWE Q_26, Q_1128 ORAZ Q_1170
        for bit_s in 0:(liczba_trojek - 1)
            bit_t = (bit_s + 1) % L
            bit_u = (bit_s + 2) % L
            spin_s = (state >> bit_s) & 1
            spin_t = (state >> bit_t) & 1
            spin_u = (state >> bit_u) & 1

            # Q_26 = RR 000|103 ; (S+_s,2*S-_u,2 + S-_s,2*S+_u,2)
            bos_s = (state_bos >> bit_s) & 1
            bos_u = (state_bos >> bit_u) & 1
            if bos_s != bos_u
                new_state_bos = state_bos ⊻ ((1 << bit_s) | (1 << bit_u))
                new_index = indeks_odwrotny[(state, new_state_bos)]
                Q_26[new_index, i] += 1.0
            end

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

                # Q_1121 = RR 123|000
                # (S+_s,1*Sz_t,1*S-_u,1 + S-_s,1*Sz_t,1*S+_u,1) * I_bos
                # Bez zmiany stanow na nodze bozonowej.
                indeks_1121 = indeks_odwrotny[(new_state, state_bos)]
                Q_1121[indeks_1121, i] += sz_ferm_t
            end
        end
    end
    return H, Q_1, Q_3, Q_4, Q_8, Q_12, Q_50, Q_62, Q_118, Q_231, Q_237, Q_289, Q_300, Q_304, Q_1128, Q_1170, Q_26, Q_60, Q_6, Q_10, Q_235, Q_292, Q_298, Q_1121
end

# =============================================================================

# WCZYTYWANIE OPERATORÓW Z ZAPISANEJ BAZY
function wczytaj_baze_operatorow(filename)
    lista_operatorow = []
    for linia in eachline(filename)
        pola = split(strip(linia))
        # Pierwsza kolumna jest ignorowana
        sektor = String(pola[2])
        nogi = split(pola[3], "|")

        # Zamiana znaków na cyfry 0,1,2,3
        fermion_1 = [parse(Int, string(c)) for c in nogi[1]]
        bozon_2  = [parse(Int, string(c)) for c in nogi[2]]
        push!(lista_operatorow, (sektor, fermion_1, bozon_2))
    end
    println("Wczytano $(length(lista_operatorow)) operatorów")
    return lista_operatorow
end

# SPRZĘŻENIE HERMITOWSKIE
function sprzezenie_operatora(lista)
    lista_hc = copy(lista)
    for i in eachindex(lista)
        if lista[i] == 1
            lista_hc[i] = 3       # S+ -> S-
        elseif lista[i] == 3
            lista_hc[i] = 1       # S- -> S+
        else
            lista_hc[i] = lista[i]  # I i Sz bez zmian
        end
    end
    return lista_hc
end

# ROZWINIĘCIE JEDNEJ NOGI W SEKTORZE R LUB I
function rozwin_noge(lista, sektor)
    # Sprzężenie hermitowskie operatora
    lista_hc = sprzezenie_operatora(lista)
    if sektor == 'R'
        wsp_1 = 1.0
        wsp_2 = 1.0
    elseif sektor == 'I'
        wsp_1 =  im
        wsp_2 = -im
    end
    # Każdy element: (operator, współczynnik)
    wyrazy = [(lista, wsp_1), (lista_hc, wsp_2)]
    return wyrazy
end

# ROZWINIĘCIE CAŁEGO OPERATORA NA DWÓCH NOGACH
function rozwin_operator(sektor, fermion_1, bozon_2)
    sektor_fermion = sektor[1]
    fermion_rozwiniety = rozwin_noge(fermion_1, sektor_fermion)
    sektor_bozon = sektor[2]
    bozon_rozwiniety = rozwin_noge(bozon_2, sektor_bozon)
    wyrazy = []

    # Każdy wyraz fermionowy mnożymy przez każdy wyraz bozonowy
    for wyraz_fermion in fermion_rozwiniety
        # Rozdzielenie operatora i jego współczynnika
        operator_fermion = wyraz_fermion[1]
        wspolczynnik_fermion = wyraz_fermion[2]
        for wyraz_bozon in bozon_rozwiniety
            # Rozdzielenie operatora i jego współczynnika
            operator_bozon = wyraz_bozon[1]
            wspolczynnik_bozon = wyraz_bozon[2]
            wspolczynnik = wspolczynnik_fermion * wspolczynnik_bozon # Mnożenie współczynników

            # Zapisanie wyrazu: (operator fermionowy, operator bozonowy, współczynnik)
            push!(wyrazy, (operator_fermion, operator_bozon, wspolczynnik))
        end
    end
    return wyrazy
end

# DZIAŁANIE JEDNEGO STRINGU NA STAN BITOWY
function dzialaj_na_stan(state, operatory, przesuniecie, L)
    new_state = state
    amplituda = 1.0
    for (x, op) in enumerate(operatory)

        op == 0 && continue

        # Numer bitu, z uwzględnieniem translacji
        bit_s = mod(przesuniecie + x - 1, L)

        spin = (new_state >> bit_s) & 1

        if op == 1
            # S+ |0> = |1>
            # S+ |1> = 0
            spin == 1 && return nothing

            new_state = new_state ⊻ (1 << bit_s)

        elseif op == 3
            # S- |1> = |0>
            # S- |0> = 0
            spin == 0 && return nothing

            new_state = new_state ⊻ (1 << bit_s)

        elseif op == 2
            # Sz |0> = -1/2 |0>
            # Sz |1> = +1/2 |1>
            amplituda *= spin - 0.5
        end
    end

    return new_state, amplituda
end


# GENEROWANIE MACIERZY OPERATORA O_i
function generuj_operator_z_bazy(numer, lista_operatorow, L, baza_int, indeks_odwrotny, PBC)
    sektor, gora, dol = lista_operatorow[numer] # Odczyt operatora o podanym numerze
    wyrazy = rozwin_operator(sektor, gora, dol) # Rozwinięcie RR, RI, IR lub II
    zasieg = findlast(x -> gora[x] != 0 || dol[x] != 0, eachindex(gora)) # Faktyczny zasięg operatora

    # PBC: wszystkie L położeń
    # OBC: tylko położenia mieszczące się w łańcuchu
    if PBC
        przesuniecia = 0:L-1
    else
        przesuniecia = 0:L-zasieg
    end
        ile_stanow = length(baza_int)

    Q = zeros(ComplexF64, ile_stanow, ile_stanow)
    # Pętla po stanach początkowych
    for i in 1:ile_stanow
        state, state_bos = baza_int[i]
        # Suma po translacjach
        for przesuniecie in przesuniecia

            # Suma po wyrazach rozwinięcia symetrii
            for (op_gora, op_dol, wsp) in wyrazy

                wynik_f = dzialaj_na_stan(state, op_gora, przesuniecie, L)
                wynik_f === nothing && continue

                wynik_b = dzialaj_na_stan(state_bos, op_dol, przesuniecie, L)
                wynik_b === nothing && continue

                new_state, amp_f = wynik_f
                new_state_bos, amp_b = wynik_b

                new_index = get(indeks_odwrotny, (new_state, new_state_bos), 0) # Szukanie indeksu nowego stanu w bazie
                new_index == 0 && continue # Jeśli stan wypada poza wybrany sektor

                # Wpis do macierzy operatora
                Q[new_index, i] += wsp * amp_f * amp_b
            end
        end
    end
    return Q
end


# =============================================================================

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

# ZAPIS MACIERZY Tr([H,Q]†[H,Q])/D VS g / omega0
function save_wyniki_vs_g(wyniki, kolumny, PBC; nazwa="Tr_komutatorow_vs_g", kwargs...)
    folder = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
    mkpath(folder)
    filename_full = "$(nazwa)_PBC$(PBC)"
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

# STAN POCZATKOWY: JEDEN ELEKTRON NA WYBRANEJ POZYCJI ORAZ LOSOWY GAUSSOWSKI WEKTOR STANU BOZONOWEGO
function losuj_stan_fermion_lokalny(M::Int, pozycja::Int, baza_int, indeks_odwrotny)
    # Jeden elektron na wybranej pozycji
    state_ferm = 1 << (pozycja - 1)
    # Losowy znormalizowany stan bozonowy
    Psi_bos = losuj_q1_gauss(1 << M)
    # Wektor stanu calego ukladu
    Psi_0 = zeros(ComplexF64, length(baza_int))
    for state_bos in 0:((1 << M) - 1)
        indeks = get(indeks_odwrotny,(state_ferm, state_bos),0)
        indeks != 0 || error("Brak stanu w bazie")
        Psi_0[indeks] = Psi_bos[state_bos + 1]
    end
    return Psi_0
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

# ZAPIS FUNKCJI KORELACJI ED VS g / omega0 I t
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

    # g_MIN = 0.0
    # g_MAX = 10.0
    # g_POINTS = 51
    # g_values = range(g_MIN, g_MAX, length=g_POINTS)

    g_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    g_MIN = minimum(g_values)
    g_MAX = maximum(g_values)
    g_POINTS = length(g_values)

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

    # CZWARTY WOLNY MOD: podane WSPOLCZYNNIKI (nie kwadraty wag)
    c4_8   = +9.9654870456974975e-01
    c4_118 = -5.3218725317609773e-02
    c4_62  = -5.3218725317609732e-02
    c4_4   = -3.3646755004379414e-02

    # PIATY WOLNY MOD: podane WSPOLCZYNNIKI (nie kwadraty wag)
    c5_298  = +6.3817744742265958e-01
    c5_292  = +6.3817744742265947e-01
    c5_1121 = -4.2780884225378929e-01
    c5_235  = -3.1705977742643089e-02
    c5_10   = -1.9813416356805943e-02
    c5_6    = -1.9813416356805828e-02

    # SZOSTY WOLNY MOD: podane WSPOLCZYNNIKI (nie kwadraty wag)
    c6_289 = -5.1794706276320912e-01
    c6_50  = +4.7927210307102769e-01
    c6_1   = -4.7651488840914752e-01
    c6_304 = -3.6816385756565068e-01
    c6_300 = -3.6816385756565068e-01
    c6_3   = +4.2433580883231209e-02

    # GENEROWANIE BAZY - cała magentyzacja
    mag = :all
    N_up = "all"
    println("Generowanie pelnej bazy...")
    baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(M; mag=mag)

    # GENEROWANIE BAZY - sektor magnetyzacji
    # mag = N_up - M / 2
    # println("Generowanie bazy...")
    # baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(M; mag=mag)

    #println("Wymiar przestrzeni Hilberta: ", length(baza_int))

    baza_data = Matrix{Any}(undef, length(baza_int), 5)
    for i in 1:length(baza_int)
        baza_data[i, 1] = i
        baza_data[i, 2] = baza_int[i][1]
        baza_data[i, 3] = baza_int[i][2]
        baza_data[i, 4] = baza_bit[i][1]
        baza_data[i, 5] = baza_bit[i][2]
    end
    Psi_0_wspolne = losuj_q1_gauss(length(baza_int))
    pozycja_fermionu=5
    Psi_0_fermion_lokalny = losuj_stan_fermion_lokalny(M, pozycja_fermionu, baza_int, indeks_odwrotny)

    println("Pozycja elektronu: ", pozycja_fermionu)
    println("Norma stanu lokalnego: ", norm(Psi_0_fermion_lokalny))
    println("Wymiar przestrzeni Hilberta: ", length(baza_int))

    # PETLA PO omega_0; tryb pojedynczej wartosci z argumentu -w:
    omega0_values = [omega0]
    # omega0_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    # omega0_values = range(0.0, 2.0, length=11)
    
    skan_omega0 = length(omega0_values) > 1
    skan_g = length(g_values) > 1
    if skan_omega0 && skan_g
        kolumny_parametrow = ["omega0", "g"]
        nazwa_skanu = "omega0_g"
    elseif skan_omega0
        kolumny_parametrow = ["omega0"]
        nazwa_skanu = "omega0"
    else
        kolumny_parametrow = ["g"]
        nazwa_skanu = "g"
    end
    liczba_kolumn_parametrow = length(kolumny_parametrow)
    liczba_punktow = length(omega0_values) * g_POINTS
    parametry_pliku = if skan_omega0
        (L_x=M, N_up=N_up, J=J, Delta=Delta,
         gmin=g_MIN, gmax=g_MAX, ng=g_POINTS,
         omega0min=minimum(omega0_values), omega0max=maximum(omega0_values), nw=length(omega0_values))
    else
        (L_x=M, N_up=N_up, J=J, Delta=Delta, omega0=first(omega0_values),
         gmin=g_MIN, gmax=g_MAX, ng=g_POINTS)
    end

    # Pierwsze kolumny: parametry skanowania. Pozostale: Tr([H,Q]†[H,Q])/D.
    # Zwykle operatory Q_i zapisujemy w kolejnosci rosnacych indeksow.
    kolumny = vcat(kolumny_parametrow, [
        "Q_1", "Q_3", "Q_4", "Q_6", "Q_8", "Q_10", "Q_12", "Q_26", "Q_50",
        "Q_60", "Q_62", "Q_118", "Q_231", "Q_235", "Q_237", "Q_289", "Q_292",
        "Q_298", "Q_300", "Q_304", "Q_1121", "Q_1128", "Q_1170",
        "Q_slow_3_50",
        "Q_slow_8_4_231_237_62_118_1128_1170",
        "Q_slow_8_12_62_118",
        "Q_slow_8_118_62_4",
        "Q_slow_298_292_1121_235_10_6",
        "Q_slow_289_50_1_304_300_3"
    ])
    wyniki = Matrix{Float64}(undef, liczba_punktow, length(kolumny))

    # PARAMETRY KORELACJI CZASOWYCH
    t_MIN = 0.0
    t_MAX = 100.0
    t_POINTS = 501
    t_values = range(t_MIN, t_MAX, length=t_POINTS)

    # Pierwsze kolumny: parametry, potem czas t. Pozostale: C_Q(t) = Re⟨Ψ₀|e^(iHt) Q e^(-iHt) Q|Ψ₀⟩
    kolumny_korelacji = vcat(kolumny_parametrow, [
        "t",
        "Q_1", "Q_3", "Q_6", "Q_8", "Q_10", "Q_26", "Q_50", "Q_60",
        "Q_235", "Q_289", "Q_292", "Q_298", "Q_300", "Q_304", "Q_1121",
        "C_slow_3_50",
        "C_slow_8_4_231_237_62_118_1128_1170",
        "C_slow_8_12_62_118",
        "C_slow_8_118_62_4",
        "C_slow_298_292_1121_235_10_6",
        "C_slow_289_50_1_304_300_3"
    ])

    # Losowy stan poczatkowy — ten sam dla wszystkich omega_0 i g
    Psi_0 = Psi_0_wspolne

    # Macierz na wyniki funkcji korelacji
    wyniki_korelacji = Matrix{Float64}(undef, liczba_punktow * t_POINTS, length(kolumny_korelacji))
    wyniki_korelacji_trace = Matrix{Float64}(undef, liczba_punktow * t_POINTS, length(kolumny_korelacji))
    wyniki_korelacji_lokalny = Matrix{Float64}(undef, liczba_punktow * t_POINTS, length(kolumny_korelacji))

    for (k_omega, omega0) in enumerate(omega0_values)
        println("Obliczenia dla omega_0 = ", omega0)

        save_data(baza_data, "Baza", PBC; L_x=M, N_up, J, Delta, omega0)

        #pętla po g =============================================================================================
        for (k, g) in enumerate(g_values)
            println("Punkt: ", k, "/", g_POINTS)
            wiersz_parametrow = (k_omega - 1) * g_POINTS + k
            parametry_wiersza = skan_omega0 ? (skan_g ? [omega0, g] : [omega0]) : [g]

            #generowane Ham i wszystkich operatorów
            H, Q_1, Q_3, Q_4, Q_8, Q_12, Q_50, Q_62, Q_118, Q_231, Q_237, Q_289, Q_300, Q_304, Q_1128, Q_1170, Q_26, Q_60, Q_6, Q_10, Q_235, Q_292, Q_298, Q_1121 =
            gen_ham_xxz(M, J, Delta, PBC, baza_int, indeks_odwrotny, omega0, g)

            # save_data(H, "Ham", PBC; L_x=M, N_up, J, Delta, omega0)
            # save_data(Q_3, "Q_3", PBC; L_x=M, N_up, J, Delta, omega0)
            # save_data(Q_4, "Q_4", PBC; L_x=M, N_up, J, Delta, omega0)
            # save_data(Q_8, "Q_8", PBC; L_x=M, N_up, J, Delta, omega0)

            if k == 1
                D = size(H, 1)
                println("--- NORMY KWADRATOWE HS ---")
                println("Q_1   = ", sum(abs2, Q_1) / D)
                println("Q_3   = ", sum(abs2, Q_3) / D)
                println("Q_4   = ", sum(abs2, Q_4) / D)
                println("Q_8   = ", sum(abs2, Q_8) / D)
                println("Q_12  = ", sum(abs2, Q_12) / D)
                println("Q_26   = ", sum(abs2, Q_26) / D)
                println("Q_50  = ", sum(abs2, Q_50) / D)
                println("Q_60   = ", sum(abs2, Q_60) / D)
                println("Q_62  = ", sum(abs2, Q_62) / D)
                println("Q_118 = ", sum(abs2, Q_118) / D)
                println("Q_231 = ", sum(abs2, Q_231) / D)
                println("Q_237 = ", sum(abs2, Q_237) / D)
                println("Q_289 = ", sum(abs2, Q_289) / D)
                println("Q_300 = ", sum(abs2, Q_300) / D)
                println("Q_304 = ", sum(abs2, Q_304) / D)
                println("Q_1128 = ", sum(abs2, Q_1128) / D)
                println("Q_1170 = ", sum(abs2, Q_1170) / D)
                println("Q_6 = ", sum(abs2, Q_6) / D)
                println("Q_10 = ", sum(abs2, Q_10) / D)
                println("Q_235 = ", sum(abs2, Q_235) / D)
                println("Q_292 = ", sum(abs2, Q_292) / D)
                println("Q_298 = ", sum(abs2, Q_298) / D)
                println("Q_1121 = ", sum(abs2, Q_1121) / D)
            end

            Q_1   = normalizuj_operator(Q_1)
            Q_3   = normalizuj_operator(Q_3)
            Q_4   = normalizuj_operator(Q_4)
            Q_8   = normalizuj_operator(Q_8)
            Q_12  = normalizuj_operator(Q_12)
            Q_26  = normalizuj_operator(Q_26)
            Q_50  = normalizuj_operator(Q_50)
            Q_60  = normalizuj_operator(Q_60)
            Q_62  = normalizuj_operator(Q_62)
            Q_118 = normalizuj_operator(Q_118)
            Q_231 = normalizuj_operator(Q_231)
            Q_237 = normalizuj_operator(Q_237)
            Q_289 = normalizuj_operator(Q_289)
            Q_300 = normalizuj_operator(Q_300)
            Q_304 = normalizuj_operator(Q_304)
            Q_1128 = normalizuj_operator(Q_1128)
            Q_1170 = normalizuj_operator(Q_1170)
            Q_6 = normalizuj_operator(Q_6)
            Q_10 = normalizuj_operator(Q_10)
            Q_235 = normalizuj_operator(Q_235)
            Q_292 = normalizuj_operator(Q_292)
            Q_298 = normalizuj_operator(Q_298)
            Q_1121 = normalizuj_operator(Q_1121)

            Q_slow_3_50 = +a_3 * Q_3 + a_50 * Q_50
            Q_slow_8_12_62_118 = (-a_8 * Q_8 -a_12 * Q_12 -a_118 * Q_118 -a_62 * Q_62)
            Q_slow_8_4_231_237_62_118_1128_1170 = (+a3_8 * Q_8 +a3_4 * Q_4 -a3_231 * Q_231 +a3_237 * Q_237 -a3_62 * Q_62 -a3_118 * Q_118 +a3_1128 * Q_1128 +a3_1170 * Q_1170)

            # DWA NOWE LIOMSY: wspolczynniki sa podane bezposrednio ze znakami.
            Q_slow_8_118_62_4 = c4_8 * Q_8 + c4_118 * Q_118 + c4_62 * Q_62 + c4_4 * Q_4
            Q_slow_298_292_1121_235_10_6 = c5_298 * Q_298 + c5_292 * Q_292 + c5_1121 * Q_1121 + c5_235 * Q_235 + c5_10 * Q_10 + c5_6 * Q_6
            Q_slow_289_50_1_304_300_3 = c6_289 * Q_289 + c6_50 * Q_50 + c6_1 * Q_1 + c6_304 * Q_304 + c6_300 * Q_300 + c6_3 * Q_3

            if k == 1
                D = size(H, 1)
                println("--- NORMY KWADRATOWE HS ---")
                println("Q_1   = ", sum(abs2, Q_1) / D)
                println("Q_3   = ", sum(abs2, Q_3) / D)
                println("Q_4   = ", sum(abs2, Q_4) / D)
                println("Q_8   = ", sum(abs2, Q_8) / D)
                println("Q_12  = ", sum(abs2, Q_12) / D)
                println("Q_26   = ", sum(abs2, Q_26) / D)
                println("Q_50  = ", sum(abs2, Q_50) / D)
                println("Q_60   = ", sum(abs2, Q_60) / D)
                println("Q_62  = ", sum(abs2, Q_62) / D)
                println("Q_118 = ", sum(abs2, Q_118) / D)
                println("Q_231 = ", sum(abs2, Q_231) / D)
                println("Q_237 = ", sum(abs2, Q_237) / D)
                println("Q_289 = ", sum(abs2, Q_289) / D)
                println("Q_300 = ", sum(abs2, Q_300) / D)
                println("Q_304 = ", sum(abs2, Q_304) / D)
                println("Q_1128 = ", sum(abs2, Q_1128) / D)
                println("Q_1170 = ", sum(abs2, Q_1170) / D)
                println("Q_6 = ", sum(abs2, Q_6) / D)
                println("Q_10 = ", sum(abs2, Q_10) / D)
                println("Q_235 = ", sum(abs2, Q_235) / D)
                println("Q_292 = ", sum(abs2, Q_292) / D)
                println("Q_298 = ", sum(abs2, Q_298) / D)
                println("Q_1121 = ", sum(abs2, Q_1121) / D)

                println("Q_slow_3_50 = ", sum(abs2, Q_slow_3_50) / D)
                println("Q_slow_8_12_62_118 = ", sum(abs2, Q_slow_8_12_62_118) / D)
                println("Q_slow_8_4_231_237_62_118_1128_1170 = ", sum(abs2, Q_slow_8_4_231_237_62_118_1128_1170) / D)
                println("Q_slow_8_118_62_4 = ", sum(abs2, Q_slow_8_118_62_4) / D)
                println("Q_slow_298_292_1121_235_10_6 = ", sum(abs2, Q_slow_298_292_1121_235_10_6) / D)
                println("Q_slow_289_50_1_304_300_3 = ", sum(abs2, Q_slow_289_50_1_304_300_3) / D)

                # println("Suma wag 3,50 = ",W_3 + W_50)
                # println("Suma wag 8,12,62,118 = ",W_8 + W_12 + W_62 + W_118)
                # println("Suma wag 8,4,231,237,62,118,1128,1170 = ",W3_8 + W3_4 + W3_231 + W3_237 + W3_62 + W3_118 + W3_1128 + W3_1170)
            end

            C_1, wynik_1 = oblicz_komutator(H, Q_1)
            C_3, wynik_3 = oblicz_komutator(H, Q_3)
            C_4, wynik_4 = oblicz_komutator(H, Q_4)
            C_8, wynik_8 = oblicz_komutator(H, Q_8)
            C_12, wynik_12 = oblicz_komutator(H, Q_12)
            C_26, wynik_26 = oblicz_komutator(H, Q_26)
            C_50, wynik_50 = oblicz_komutator(H, Q_50)
            C_60, wynik_60 = oblicz_komutator(H, Q_60)
            C_62, wynik_62 = oblicz_komutator(H, Q_62)
            C_118, wynik_118 = oblicz_komutator(H, Q_118)
            C_231, wynik_231 = oblicz_komutator(H, Q_231)
            C_237, wynik_237 = oblicz_komutator(H, Q_237)
            C_289, wynik_289 = oblicz_komutator(H, Q_289)
            C_300, wynik_300 = oblicz_komutator(H, Q_300)
            C_304, wynik_304 = oblicz_komutator(H, Q_304)
            C_1128, wynik_1128 = oblicz_komutator(H, Q_1128)
            C_1170, wynik_1170 = oblicz_komutator(H, Q_1170)
            C_slow_1, wynik_slow_1 = oblicz_komutator(H, Q_slow_3_50)
            C_slow_2, wynik_slow_2 = oblicz_komutator(H, Q_slow_8_12_62_118)
            C_slow_3, wynik_slow_3 = oblicz_komutator(H, Q_slow_8_4_231_237_62_118_1128_1170)
            C_6, wynik_6 = oblicz_komutator(H, Q_6)
            C_10, wynik_10 = oblicz_komutator(H, Q_10)
            C_235, wynik_235 = oblicz_komutator(H, Q_235)
            C_292, wynik_292 = oblicz_komutator(H, Q_292)
            C_298, wynik_298 = oblicz_komutator(H, Q_298)
            C_1121, wynik_1121 = oblicz_komutator(H, Q_1121)
            C_slow_4, wynik_slow_4 = oblicz_komutator(H, Q_slow_8_118_62_4)
            C_slow_5, wynik_slow_5 = oblicz_komutator(H, Q_slow_298_292_1121_235_10_6)
            C_slow_6, wynik_slow_6 = oblicz_komutator(H, Q_slow_289_50_1_304_300_3)

            wyniki[wiersz_parametrow, :] = vcat(parametry_wiersza, [
                wynik_1, wynik_3, wynik_4, wynik_6, wynik_8, wynik_10, wynik_12, wynik_26, wynik_50,
                wynik_60, wynik_62, wynik_118, wynik_231, wynik_235, wynik_237,
                wynik_289, wynik_292, wynik_298, wynik_300, wynik_304, wynik_1121, wynik_1128, wynik_1170,
                wynik_slow_1, wynik_slow_3, wynik_slow_2, wynik_slow_4, wynik_slow_5, wynik_slow_6
            ])

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

            # Matrix[...] utrzymuje zarowno Matrix{Float64}, jak i Matrix{ComplexF64}
            # bez konwersji wszystkich operatorow na zespolone (mniejsze zuzycie pamieci).

            operatory_korelacji = Matrix[
                Q_1, Q_3, Q_6, Q_8, Q_10, Q_26, Q_50, Q_60,
                Q_235, Q_289, Q_292, Q_298, Q_300, Q_304, Q_1121,
                Q_slow_3_50, Q_slow_8_4_231_237_62_118_1128_1170,
                Q_slow_8_12_62_118, Q_slow_8_118_62_4,
                Q_slow_298_292_1121_235_10_6,
                Q_slow_289_50_1_304_300_3,
            ] #WYBOR OPERATOROW
            wynik_H = eigen(Symmetric(H)) # DIAGONALIZACJA ED - RAZ DLA KAZDEGO g

            dane_korelacji = Vector{Any}(undef, length(operatory_korelacji))
            for a in 1:length(operatory_korelacji)
                A = operatory_korelacji[a]
                dane_korelacji[a] = przygotuj_korelacje_ED(wynik_H, A, Psi_0)
            end

            # TRZECIA METODA - ELEKTRON LOKALNY + LOSOWE BOZONY
            Psi_wlasna_lokalny = adjoint(wynik_H.vectors) * Psi_0_fermion_lokalny
            dane_korelacji_lokalny = Vector{Any}(undef, length(operatory_korelacji))
            for a in 1:length(operatory_korelacji)
                # Wykorzystujemy juz obliczona macierz operatora
                # w bazie wlasnej Hamiltonianu
                A_wlasna, _, _, W = dane_korelacji[a]
                Phi_wlasna_lokalny = A_wlasna * Psi_wlasna_lokalny
                dane_korelacji_lokalny[a] = (A_wlasna, Psi_wlasna_lokalny, Phi_wlasna_lokalny, W)
            end

            # PETLA PO CZASIE
            for (j, t) in enumerate(t_values)
                fazy = exp.(-im .* wynik_H.values .* t)
                wiersz = (wiersz_parametrow - 1) * t_POINTS + j
                wyniki_korelacji[wiersz, 1:liczba_kolumn_parametrow] = parametry_wiersza
                wyniki_korelacji[wiersz, liczba_kolumn_parametrow + 1] = t
                wyniki_korelacji_trace[wiersz, 1:liczba_kolumn_parametrow] = parametry_wiersza
                wyniki_korelacji_trace[wiersz, liczba_kolumn_parametrow + 1] = t
                # Parametry i czas dla trzeciej metody
                wyniki_korelacji_lokalny[wiersz, 1:liczba_kolumn_parametrow] = parametry_wiersza
                wyniki_korelacji_lokalny[wiersz, liczba_kolumn_parametrow + 1] = t

                for a in 1:length(operatory_korelacji)
                    dane = dane_korelacji[a]

                    C_t = funkcja_korelacji_ED(dane, fazy)
                    wyniki_korelacji[wiersz, a + liczba_kolumn_parametrow + 1] = C_t

                    C_trace_t = funkcja_korelacji_trace_ED(dane, fazy)
                    wyniki_korelacji_trace[wiersz, a + liczba_kolumn_parametrow + 1] = C_trace_t

                    # Trzecia metoda - lokalny elektron
                    C_lokalny_t = funkcja_korelacji_ED(dane_korelacji_lokalny[a], fazy)
                    wyniki_korelacji_lokalny[wiersz, a + liczba_kolumn_parametrow + 1] = C_lokalny_t
                end
            end
        end
    end # koniec petli po omega_0
    # Zapis calej siatki po zakonczeniu petli: nazwa pliku odpowiada faktycznym kolumnom.
    save_wyniki_vs_g(wyniki, kolumny, PBC;
        nazwa="Tr_komutatorow_vs_$(nazwa_skanu)", parametry_pliku...)
    save_korelacje_ED(wyniki_korelacji, kolumny_korelacji, PBC;
        nazwa="Korelacje_ED_vs_$(nazwa_skanu)_t", parametry_pliku...,
        tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
    save_korelacje_ED(wyniki_korelacji_trace, kolumny_korelacji, PBC;
        nazwa="Korelacje_trace_ED_vs_$(nazwa_skanu)_t", parametry_pliku...,
        tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
    # ZAPIS TRZECIEJ METODY:
    save_korelacje_ED(wyniki_korelacji_lokalny, kolumny_korelacji, PBC;
        nazwa="Korelacje_fermion_lokalny_ED_vs_$(nazwa_skanu)_t", parametry_pliku...,
        x=pozycja_fermionu, tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
