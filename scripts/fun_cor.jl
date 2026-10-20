using Pkg
using ArgParse
using SparseArrays
using LinearAlgebra
using DelimitedFiles

# Generowanie bazy XXZ + hard-core bozony
function gen_baza_xxz(L::Int; mag=:all)
    # Rozmiar bazy w zależności od trybu
    rozmiar_pelnej_bazy = 1 << L # przesówanie bitowe jedynki w lewo o L mijesc
    if mag == :all
        liczba_stanow = rozmiar_pelnej_bazy
        sektory = 0:L #sektor zdefiniowany prz liczbę N_up w nim
    elseif mag isa Real
        N_up = L/2 + mag # N_up = (N_up+N_down)/2 + (N_up - N_down)/2
        N_up = Int(N_up) # tka dla pewność
        liczba_stanow = binomial(L, N_up)
        sektory = N_up:N_up #sektor zdefiniowany prz liczbę N_up w nim
    end

    # NOWE - uwzględnienie drugiej nogi
    liczba_stanow = liczba_stanow * rozmiar_pelnej_bazy

    # Tworzenie baz z uwzgledeniem sektorów magnetyzacji
    baza_int = Vector{Tuple{Int, Int}}(undef, liczba_stanow) # NOWE
    i = 1 # bo Julia
    for N_up_sektor in sektory
        for state in 0:(rozmiar_pelnej_bazy - 1)
            if count_ones(state) == N_up_sektor

                # NOWE - wszystkie konfiguracje bozonów
                for state_bos in 0:(rozmiar_pelnej_bazy - 1)
                    baza_int[i] = (state, state_bos)
                    i = i + 1
                end

            end
        end
    end

    # Odwzorowanie odworotne stan -> indeks w bazie, zapisanw jako słownik
    indeks_odwrotny = Dict{Tuple{Int, Int}, Int}() # NOWE
    for (i, state) in enumerate(baza_int)
        indeks_odwrotny[state] = i
    end
    println("LICZBA STANOW: ")
    println(liczba_stanow)

    
    # Zamiana bazy na zapis binarny
    baza_bit = [
        (string(state[1], base=2, pad=L),
         string(state[2], base=2, pad=L))
        for state in baza_int
    ]

    return baza_int, baza_bit, indeks_odwrotny
end

# Generowanie ham XXZ + hard-core bozony
function gen_ham_xxz(L::Int, J::Real, Delta::Real, PBC::Bool,
    baza_int::Vector{Tuple{Int, Int}},
    indeks_odwrotny::Dict{Tuple{Int, Int}, Int},
    omega0::Real, g::Real)

    #tworzenie macierzy Hamiltonianu
    ile_stanow = length(baza_int)
    H = zeros(Float64, ile_stanow, ile_stanow)

    # Liczba połączeń między sąsiednimi spinami
    liczba_polaczen = L - 1
    if PBC == true && L > 2
        liczba_polaczen += 1
    end

    #pętla po stanch
    for i in 1:length(baza_int)
        state, state_bos = baza_int[i] # NOWE - odczyt obu nóg

        #pętla po możliwych przeskokach
        for bit_s in 0:(liczba_polaczen - 1)
            #struktura poąłczeń (dla PBC ostatni łączy się z pierwszym)
            if bit_s == L - 1
                bit_t = 0
            else
                bit_t = bit_s + 1
            end

            #odczyt spinów z uzyciem bitowego and (bit&1) 0 = down, 1 = up
            spin_s = (state >> bit_s) & 1
            spin_t = (state >> bit_t) & 1

            # CZĘŚĆ DIAGONALNA
            if spin_s == spin_t
                H[i, i] += (J*Delta)/4
            else
                H[i, i] -= (J*Delta)/4
            end

            # CZĘŚĆ POZADIAGONALNA
            if spin_s != spin_t
                # Zamiana dwóch spinów przez XOR - ⊻ a to | - to or
                new_state = state ⊻ ((1 << bit_s) | (1 << bit_t))
                new_index = indeks_odwrotny[(new_state, state_bos)] # NOWE
                H[new_index, i] += J / 2
            end
        end

        # =====================================================
        # NOWE - CZŁON DIAGONALNY BOZONÓW
        # omega0 * n_l,2
        # =====================================================

        H[i, i] += omega0 * count_ones(state_bos)

        # =====================================================
        # NOWE - SPRZĘŻENIE SPIN-BOZON
        # g * n_l,1 * (b†_l,2 + b_l,2)
        # =====================================================

        for bit_s in 0:(L - 1)

            spin_s = (state >> bit_s) & 1

            if spin_s == 1
                # Kreacja lub anihilacja bozonu
                new_state_bos = state_bos ⊻ (1 << bit_s)

                new_index = indeks_odwrotny[(state, new_state_bos)]
                H[new_index, i] += g
            end

        end

    end
    return H
end

function parse_args()
    s = ArgParseSettings(description = "Generate XXZ Hamiltonian with hard-core bosons")
    @add_arg_table! s begin
        "-M", "--sites"
            help = "Długość łańcucha, czyli liczba spinów"
            arg_type = Int
            default = 6
            dest_name = "M"
        "-N", "--n-up"
            help = "Liczba spinów w górę"
            arg_type = Int
            default = 3
            dest_name = "N_up"
        "-J"
            help = "Stała wymiany J"
            arg_type = Float64
            default = 1.0
            dest_name = "J"
        "-d", "--delta"
            help = "Parametr anizotropii Delta"
            arg_type = Float64
            default = 0.3
            dest_name = "Delta"
        "-w", "--omega0"
            help = "Energia bozonów hard-core"
            arg_type = Float64
            default = 1.0
            dest_name = "omega0"
        "-P", "--pbc"
            help = "Periodyczne warunki brzegowe: yes albo no"
            arg_type = String
            required = true
            range_tester = x -> x in ["yes", "no"]
            dest_name = "pbc"
    end
    return ArgParse.parse_args(s)
end

function save_data(data, filename, PBC; kwargs...)

    # Folder zapisu
    folder = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\fun_cor_data"
    mkpath(folder)

    # Nazwa pliku
    filename_full = "$(filename)_PBC$(PBC)"

    # Dodanie parametrów do nazwy
    for (key, value) in kwargs
        filename_full *= "_$(key)$(value)"
    end

    filepath = joinpath(folder, filename_full * ".data")

    # Zapis danych
    writedlm(filepath, data)

    println("Zapisano do pliku: ", filepath)

end

"""
Generating a random normalized vector
"""
function losuj_q1(n::Int)
    q1 = (2 .* rand(n) .- 1) .+ im .* (2 .* rand(n) .- 1)
    q1 = q1 / norm(q1)
    return q1
end

"""
Generating a random normalized vector with Gaussian distribution
"""
function losuj_q1_gauss(n::Int)
    q1 = randn(n) .+ im .* randn(n)
    q1 = q1 / norm(q1)
    return q1
end

"""
Correlation function
"""
function funkcja_korelacji(H, A, delta_t::Real, Psi, Phi)

    # Jednorazowa ewolucja t -> t + delta_t za pomocą funkcji
    Psi_koncowe = dynamika_lanczos(H, Psi, delta_t)
    Phi_koncowe = dynamika_lanczos(H, Phi, delta_t)

    C_t = dot(Psi_koncowe, A * Phi_koncowe) # Funkcja korelacji w nowym punkcie czasu

    return real(C_t), Psi_koncowe, Phi_koncowe
end

"""
The Exact Diagonalization Algorithm for time evolution
"""
function dynamika_ED(H, psi0, delta_t, liczba_krokow)
    n = size(H, 1)
    psi0 = ComplexF64.(psi0)
    stany = zeros(ComplexF64, n, liczba_krokow)
    # Diagonalizacja Hamiltonianu H = VDV†
    wynik_H = eigen(Matrix(H))
    D = Diagonal(wynik_H.values)
    V = wynik_H.vectors
    for krok in 1:liczba_krokow
        psi = V * exp(-im * D * krok * delta_t) * adjoint(V) * psi0
        stany[:, krok] = psi
    end
    return stany
end

"""
Corelation function using Exact Diagonalization
"""
function funkcja_korelacji_ED(H, A, t::Real, Psi_0)
    Psi_0 = ComplexF64.(Psi_0)
    Phi_0 = A * Psi_0
    wynik_H = eigen(Matrix(H)) 
    D = Diagonal(wynik_H.values)
    V = wynik_H.vectors

    Psi_t = V * exp(-im * D * t) * adjoint(V) * Psi_0
    Phi_t = V * exp(-im * D * t) * adjoint(V) * Phi_0
    C_t = dot(Psi_t, A * Phi_t) # C(t) = <Psi(t)| A |Phi(t)>
    return real(C_t)
end

function main()
# Odczyt parametrów
    args = parse_args()

    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    omega0 = args["omega0"]
    PBC = args["pbc"] == "yes"

    g = 1.0

    # Magnetyzacja pierwszej nogi
    mag = N_up - M/2

    # ========================================================
    # GENEROWANIE BAZY
    # ========================================================

    println("Generowanie bazy...")

    baza_int, baza_bit, indeks_odwrotny = gen_baza_xxz(M; mag=mag)

    # ========================================================
    # GENEROWANIE HAMILTONIANU
    # ========================================================

    println("Generowanie Hamiltonianu...")

    @time H = gen_ham_xxz(
        M, J, Delta, PBC,
        baza_int, indeks_odwrotny,
        omega0, g
    )

    println("Rozmiar Hamiltonianu: ", size(H))

    # ========================================================
    # ZAPIS HAMILTONIANU
    # ========================================================

    save_data(H, "Hamiltonian", PBC;
        L_x=M, N_up, J, Delta, omega0, g
    )

    # ========================================================
    # ZAPIS BAZY
    # ========================================================

    baza_data = Matrix{Any}(undef, length(baza_int), 5)

    for i in 1:length(baza_int)

        baza_data[i, 1] = i
        baza_data[i, 2] = baza_int[i][1]
        baza_data[i, 3] = baza_int[i][2]
        baza_data[i, 4] = baza_bit[i][1]
        baza_data[i, 5] = baza_bit[i][2]

    end

    save_data(baza_data, "Baza", PBC;
        L_x=M, N_up, J, Delta, omega0, g
    )


    # # ---- Hamiltonian
    # baza_bin, baza_int, indeks_odwrotny = gen_baza_xxz(M, N_up)
    # # println("Wygenerowano bazę modelu XXZ.")
    # zapisz_baze_do_plikow(M, N_up, baza_bin, baza_int, indeks_odwrotny)
    # # println("Zapisano do plików bazę modelu XXZ.")
    # H, A_1, J_1 = gen_ham_xxz(M, J, Delta, pbc, baza_bin, indeks_odwrotny)
    # # println("Stworzono Hamiltonian modelu XXZ.")    
    # # save_ham(H, M, N_up, J, Delta, pbc_argument)
    # # println("Hamiltonian zapisano do pliku.")

    # # ---- Jaki opertor?
    # if operator_argument == "A_1"
    #     A = A_1
    # elseif operator_argument == "J_1"
    #     A = J_1
    # end
    # println("Wybrany operator: ", operator_argument)

    # # ---- ED ----
    # wynik = eigen(Matrix(H))
    # wartosci_wlasne = wynik.values
    # wektory_wlasne = wynik.vectors
    # save_eigenvalues(wartosci_wlasne, "wartosci_wlasne_ED", M, N_up, J, Delta, pbc_argument)
    # println("Wartości własne zapisano do pliku.")

    # # # ---- dynamika za pomocą diagonalizacji ----
    # # # psi0 = losuj_q1_gauss(size(H, 1))
    # # # delta_t = 0.01
    # # # liczba_krokow = 100
    # # stany_ED = dynamika_ED(H, psi0, delta_t, liczba_krokow)
    # # # zapisanie wyników ewolucji do pliku
    # # output_directory = joinpath(@__DIR__, "Lanczos_XXZ_wyniki")
    # # plik_dynamika_ED = joinpath(output_directory, "dynamika_ED_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument)_dt_$(delta_t)_kroki_$(liczba_krokow).txt")
    # # open(plik_dynamika_ED, "w") do io
    # #     for i in axes(stany_ED, 1)
    # #         println(io, join(stany_ED[i, :], " "))
    # #     end
    # # end

    # # ---- funkcja korealacji -----
    # Psi_0 = losuj_q1_gauss(size(H, 1))
    # # Psi_0 = fill(1 / sqrt(size(H, 1)), size(H, 1))
    # # Psi_0 = zeros(ComplexF64, size(H, 1))
    # # Psi_0[1] = 1.0

    # t_min = 0.0
    # t_max = 50.0
    # liczba_krokow_czasu = 1000
    # delta_t = (t_max - t_min) / liczba_krokow_czasu # Jeden krok czasu = jeden restart Lanczosa

    # println("Przedział czasu: ", t_min, " - ", t_max)
    # println("Liczba kroków czasu: ", liczba_krokow_czasu)
    # println("delta_t = ", delta_t)
    # # Stany początkowe
    # Psi = Psi_0
    # Phi = A * Psi
    # korelacja_cala = zeros(Float64, liczba_krokow_czasu + 1, 2)
    # korelacja_cala[1, 1] = t_min
    # korelacja_cala[1, 2] = korelacja_cala[1, 2]

    # for krok in 1:liczba_krokow_czasu
    #     println("correlation function: ", krok, "/", liczba_krokow_czasu)
    #     C_t_ED = funkcja_korelacji_ED(H, A, t, Psi_0)
    #     korelacja_cala[krok + 1, 2] = C_t_ED
    # end

    # # Zapis funkcji korelacji do pliku
    # output_directory = joinpath(@__DIR__, "model_wyniki")
    # mkpath(output_directory)
    # plik_korelacja = joinpath(output_directory, "korelacja_ED_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument)_operator_$(operator_argument).txt")
    # open(plik_korelacja, "w") do io
    #     for i in axes(korelacja_cala, 1)
    #         println(io, korelacja_cala[i, 1], " ", korelacja_cala[i, 2])
    #     end
    # end
    # println("Funkcja korelacji została zapisana do pliku:")
    # println(plik_korelacja)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end