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

# GENEROWANIE HAMILTONIANU
function gen_ham_xxz(L::Int, J::Real, Delta::Real, PBC::Bool, baza_int::Vector{Tuple{Int, Int}}, indeks_odwrotny::Dict{Tuple{Int, Int}, Int}, omega0::Real, g::Real)
    ile_stanow = length(baza_int)
    liczba_polaczen = L - 1
    if PBC == true && L > 2
        liczba_polaczen += 1
    end

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
    return H
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
        # Numer bitu, z uwzględnieniem translacji
        bit_s = mod(przesuniecie + x - 1, L)
        spin = (new_state >> bit_s) & 1

        if op == 0
            # I |0> = |0>
            # I |1> = |1>
            continue

        elseif op == 1
            # S+ |0> = |1>
            # S+ |1> = 0
            if spin == 1
                return nothing
            end
            new_state = new_state ⊻ (1 << bit_s)

        elseif op == 2
            # Sz |0> = -1/2 |0>
            # Sz |1> = +1/2 |1>
            amplituda *= spin - 0.5

        elseif op == 3
            # S- |1> = |0>
            # S- |0> = 0
            if spin == 0
                return nothing
            end
            new_state = new_state ⊻ (1 << bit_s)
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

# OBLICZANIE KOMUTATORA
function oblicz_komutator(H, Q)
    C = H * Q - Q * H
    D = size(H, 1)
    wynik = sum(abs2, C) / D
    return C, wynik
end

# LOSOWY STAN POCZATKOWY Gauss
function losuj_q1_gauss(n::Int)
    q1 = randn(n) .+ im .* randn(n)
    return q1 / norm(q1)
end

# LOSOWY STAN POCZATKOWY - FERMION LOKALNY + BOZON Gauss
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

#ZAPIS - OGOLNY
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
        "-g", "--coupling"
            help = "Sprzezenie g"
            arg_type = Float64
            default = 0.75
            dest_name = "g"
        "-P", "--pbc"
            help = "Warunki brzegowe yes lub no"
            arg_type = String
            required = true
            range_tester = x -> x in ["yes", "no"]
            dest_name = "pbc"
    end
    return ArgParse.parse_args(s)
end

# =========================================================================================================
function main()
    args = parse_args()
    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    omega0 = args["omega0"]
    g = args["g"]
    PBC = args["pbc"] == "yes"

    # =========================================================================================
    # WYBOR SKANOWANIA
    SKAN_OMEGA0 = false
    SKAN_G = true

    # ARGUMENT omega_0
    if SKAN_OMEGA0
        # omega0_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
        omega0_values = range(0.0, 2.0, length=11)
    else
        omega0_values = [omega0]
    end

    # ARGUMENT g
    if SKAN_G
        # g_values = range(0.0, 10.0, length=51)
        g_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    else
        g_values = [g]
    end
    # =========================================================================================

    # =========================================================================================
    # GENEROWANIE BAZY - Sektor magnetyzacji czy czy bez podziału na sektory 
    # cała magentyzacja
    mag = :all
    N_up = "all"
    println("Generowanie pelnej bazy...")
    baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(M; mag=mag)

    # sektor magnetyzacji
    # mag = N_up - M / 2
    # println("Generowanie bazy...")
    # baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(M; mag=mag)

    # WCZYTANIE ZAPISANEJ BAZY OPERATOROW
    plik_bazy_operatorow = joinpath(@__DIR__, "liom_mag_dane", "bazy",
        "operators_ladder_mag_M_3_J_1.0_Jp_0.0_d_0.001_d2_0.0_w_1.0_g_1.0_T_both_P_both_Sz_cons_fermion_yes_Sz_cons_boson_both_FId_true.txt")
    lista_operatorow = wczytaj_baze_operatorow(plik_bazy_operatorow)

    # GENEROWANIE OPERATOROW Q Z BAZY - ich numery
    numery_operatorow = [1, 3, 50, 289, 300, 304]

    operatory_z_bazy = Dict{Int, Matrix{ComplexF64}}()
    for numer in numery_operatorow
        operatory_z_bazy[numer] = generuj_operator_z_bazy(numer, lista_operatorow, M, baza_int, indeks_odwrotny, PBC)
    end

    save_data(operatory_z_bazy[1], "Q_1", PBC; L_x=M, N_up=N_up)
    save_data(operatory_z_bazy[3], "Q_3", PBC; L_x=M, N_up=N_up)

    # WOLNY MOD dla M=3 ten schodzący do zera
    c_289 = -5.1794706276320912e-01
    c_50  = +4.7927210307102769e-01
    c_1   = -4.7651488840914752e-01
    c_304 = -3.6816385756565068e-01
    c_300 = -3.6816385756565068e-01
    c_3   = +4.2433580883231209e-02

    Q_slow_289_50_1_304_300_3 = (
    c_289 * operatory_z_bazy[289] +
    c_50  * operatory_z_bazy[50]  +
    c_1   * operatory_z_bazy[1]   +
    c_304 * operatory_z_bazy[304] +
    c_300 * operatory_z_bazy[300] +
    c_3   * operatory_z_bazy[3]
    )

    baza_data = Matrix{Any}(undef, length(baza_int), 5)
    for i in 1:length(baza_int)
        baza_data[i, 1] = i
        baza_data[i, 2] = baza_int[i][1]
        baza_data[i, 3] = baza_int[i][2]
        baza_data[i, 4] = baza_bit[i][1]
        baza_data[i, 5] = baza_bit[i][2]
    end

    #Losowe stany
    Psi_0 = losuj_q1_gauss(length(baza_int))
    pozycja_fermionu=5
    Psi_0_fermion_lokalny = losuj_stan_fermion_lokalny(M, pozycja_fermionu, baza_int, indeks_odwrotny)
    # println("Pozycja elektronu: ", pozycja_fermionu)
    # println("Norma stanu lokalnego: ", norm(Psi_0_fermion_lokalny))
    # println("Wymiar przestrzeni Hilberta: ", length(baza_int))

    # USTALENIE RODZAJU SKANOWANIA
    skan_omega0 = length(omega0_values) > 1
    skan_g = length(g_values) > 1

    if skan_omega0 && skan_g
        kolumny_parametrow = ["omega0", "g"]
        nazwa_skanu = "omega0_g"
    elseif skan_omega0
        kolumny_parametrow = ["omega0"]
        nazwa_skanu = "omega0"
    elseif skan_g
        kolumny_parametrow = ["g"]
        nazwa_skanu = "g"
    else
        kolumny_parametrow = ["omega0", "g"]
        nazwa_skanu = "punkt"
    end

    liczba_kolumn_parametrow = length(kolumny_parametrow)
    liczba_punktow = length(omega0_values) * length(g_values)

    # PARAMETRY ZAPISYWANE W NAZWIE PLIKU
    if skan_omega0 && skan_g
        parametry_pliku = (
            L_x=M, N_up=N_up, J=J, Delta=Delta,
            omega0min=minimum(omega0_values), omega0max=maximum(omega0_values),
            nw=length(omega0_values),
            gmin=minimum(g_values), gmax=maximum(g_values), ng=length(g_values)
        )
    elseif skan_omega0
        parametry_pliku = (
            L_x=M, N_up=N_up, J=J, Delta=Delta,
            omega0min=minimum(omega0_values), omega0max=maximum(omega0_values),
            nw=length(omega0_values), g=first(g_values)
        )
    elseif skan_g
        parametry_pliku = (
            L_x=M, N_up=N_up, J=J, Delta=Delta,
            omega0=first(omega0_values),
            gmin=minimum(g_values), gmax=maximum(g_values), ng=length(g_values)
        )
    else
        parametry_pliku = (
            L_x=M, N_up=N_up, J=J, Delta=Delta,
            omega0=first(omega0_values), g=first(g_values)
        )
    end

    # Pierwsze kolumny: parametry skanowania. Pozostale: Tr([H,Q]†[H,Q])/D.
    kolumny = vcat(kolumny_parametrow,
        ["Q_$numer" for numer in numery_operatorow],
        ["Q_slow_289_50_1_304_300_3"]
    )
    wyniki = Matrix{Float64}(undef, liczba_punktow, length(kolumny))

    # Pierwsze kolumny: parametry, potem czas t. Pozostale: C_Q(t)
    kolumny_korelacji = vcat(kolumny_parametrow, ["t"],
        ["Q_$numer" for numer in numery_operatorow],
        ["C_slow_289_50_1_304_300_3"]
    )

    # PARAMETRY KORELACJI CZASOWYCH
    t_MIN = 0.0
    t_MAX = 100.0
    t_POINTS = 501
    t_values = range(t_MIN, t_MAX, length=t_POINTS)

    # Macierz na wyniki funkcji korelacji
    wyniki_korelacji = Matrix{Float64}(undef, liczba_punktow * t_POINTS, length(kolumny_korelacji))
    wyniki_korelacji_trace = Matrix{Float64}(undef, liczba_punktow * t_POINTS, length(kolumny_korelacji))
    wyniki_korelacji_lokalny = Matrix{Float64}(undef, liczba_punktow * t_POINTS, length(kolumny_korelacji))

    for (k_omega, omega0) in enumerate(omega0_values)
        println("Obliczenia dla omega_0 = ", omega0)
        for (k_g, g) in enumerate(g_values)
            println("Obliczenia dla g = ", g)

            wiersz_parametrow = (k_omega - 1) * length(g_values) + k_g
            if skan_omega0 && skan_g
                parametry_wiersza = [omega0, g]
            elseif skan_omega0
                parametry_wiersza = [omega0]
            elseif skan_g
                parametry_wiersza = [g]
            else
                parametry_wiersza = [omega0, g]
            end

            # generowanie Hamiltonianu (operatory Q wygenerowano juz z bazy)
            H = gen_ham_xxz(M, J, Delta, PBC, baza_int, indeks_odwrotny, omega0, g)
            save_data(H, "Ham", PBC; L_x=M, N_up, J, Delta, omega0, g)

            #policzenie komutatorów i zapisanie wyników
            wyniki_komutatorow = Float64[]
            for numer in numery_operatorow
                _, wynik = oblicz_komutator(H, operatory_z_bazy[numer])
                push!(wyniki_komutatorow, wynik)
            end
            _, wynik_slow_1 = oblicz_komutator(H, Q_slow_289_50_1_304_300_3)
            wyniki[wiersz_parametrow, :] = vcat(parametry_wiersza, wyniki_komutatorow, [wynik_slow_1])

            #korelacja czasowa ED
            operatory_korelacji = Matrix[operatory_z_bazy[numer] for numer in numery_operatorow] #WYBOR OPERATOROW
            push!(operatory_korelacji, Q_slow_289_50_1_304_300_3)
            
            # DIAGONALIZACJA ED - RAZ DLA KAZDEGO g lub omega_0
            wynik_H = eigen(Symmetric(H)) 

            # LOSOWY STAN POCZATKOWY
            dane_korelacji = Vector{Any}(undef, length(operatory_korelacji))
            for a in 1:length(operatory_korelacji)
                A = operatory_korelacji[a]
                dane_korelacji[a] = przygotuj_korelacje_ED(wynik_H, A, Psi_0)
            end

            # ELEKTRON LOKALNY + LOSOWE BOZONY
            dane_korelacji_lokalny = Vector{Any}(undef, length(operatory_korelacji))
            for a in 1:length(operatory_korelacji)
                A = operatory_korelacji[a]
                dane_korelacji_lokalny[a] = przygotuj_korelacje_ED(wynik_H, A, Psi_0_fermion_lokalny)
            end

            # NORMALIZACJA FUNKCJI KORELACJI: C(t) / C(0)
            fazy_0 = ones(ComplexF64, length(wynik_H.values))
            C_0 = [funkcja_korelacji_ED(dane, fazy_0) for dane in dane_korelacji]
            C_0_trace = [funkcja_korelacji_trace_ED(dane, fazy_0) for dane in dane_korelacji]
            C_0_lokalny = [funkcja_korelacji_ED(dane, fazy_0) for dane in dane_korelacji_lokalny]

            # PETLA PO CZASIE
            for (j, t) in enumerate(t_values)
                fazy = exp.(-im .* wynik_H.values .* t)

                wiersz = (wiersz_parametrow - 1) * t_POINTS + j
                wyniki_korelacji[wiersz, 1:liczba_kolumn_parametrow] = parametry_wiersza
                wyniki_korelacji[wiersz, liczba_kolumn_parametrow + 1] = t
                wyniki_korelacji_trace[wiersz, 1:liczba_kolumn_parametrow] = parametry_wiersza
                wyniki_korelacji_trace[wiersz, liczba_kolumn_parametrow + 1] = t
                wyniki_korelacji_lokalny[wiersz, 1:liczba_kolumn_parametrow] = parametry_wiersza
                wyniki_korelacji_lokalny[wiersz, liczba_kolumn_parametrow + 1] = t

                for a in 1:length(operatory_korelacji)

                    C_t = funkcja_korelacji_ED(dane_korelacji[a], fazy) / C_0[a]
                    wyniki_korelacji[wiersz, a + liczba_kolumn_parametrow + 1] = C_t

                    C_trace_t = funkcja_korelacji_trace_ED(dane_korelacji[a], fazy) / C_0_trace[a]
                    wyniki_korelacji_trace[wiersz, a + liczba_kolumn_parametrow + 1] = C_trace_t

                    C_lokalny_t = funkcja_korelacji_ED(dane_korelacji_lokalny[a], fazy) / C_0_lokalny[a]
                    wyniki_korelacji_lokalny[wiersz, a + liczba_kolumn_parametrow + 1] = C_lokalny_t
                end
            end
        end
    end 
    # Zapis calej siatki po zakonczeniu petli: nazwa pliku odpowiada faktycznym kolumnom.
    save_wyniki_vs_g(wyniki, kolumny, PBC;
        nazwa="Tr_komutatorow_vs_$(nazwa_skanu)", parametry_pliku...)
    save_korelacje_ED(wyniki_korelacji, kolumny_korelacji, PBC;
        nazwa="Korelacje_ED_vs_$(nazwa_skanu)_t", parametry_pliku...,
        tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
    save_korelacje_ED(wyniki_korelacji_trace, kolumny_korelacji, PBC;
        nazwa="Korelacje_trace_ED_vs_$(nazwa_skanu)_t", parametry_pliku...,
        tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
    save_korelacje_ED(wyniki_korelacji_lokalny, kolumny_korelacji, PBC;
        nazwa="Korelacje_fermion_lokalny_ED_vs_$(nazwa_skanu)_t", parametry_pliku...,
        x=pozycja_fermionu, tmin=t_MIN, tmax=t_MAX, nt=t_POINTS)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end