using Pkg
using ArgParse
using SparseArrays
using LinearAlgebra

"""
Dodaje zera z lewej strony napisu, aż osiągnie długość doc_dlugosc.
Przykład: add_zeros("101", 5) == "00101"
"""
function add_zeros(s::AbstractString, doc_dlugosc::Int)::String
    return lpad(s, doc_dlugosc, '0')
end

"""
Zamienia liczbę dziesiętną na zapis binarny o zadanej minimalnej długości.
Przykład: dec_na_bin(5, 4) == "0101"
"""
function dec_na_bin(dec::Integer, length::Int)::String
    bin = string(dec, base = 2)
    return add_zeros(bin, length)
end

"""
Zamienia napis binarny na liczbę dziesiętną.
Przykład: bin_na_dec("0101") == 5
"""
function bin_na_dec(bin::AbstractString)::Int
    return parse(Int, bin; base = 2)
end

"""
Oblicza symbol Newtona: n nad k.
"""
function dwumian(n::Int, k::Int)::Int
    if k < 0 || k > n
        return 0
    end
    if k == 0 || k == n
        return 1
    end
    # (n nad k) = (n nad n-k)
    k = min(k, n - k)
    result = 1
    for i in 1:k
        result *= n - i + 1
        result ÷= i
    end
    return result
end

"""
Zwraca:
    baza_bin         - stany jako napisy binarne,
    baza_int         - stany jako liczby całkowite,
    indeks_odwrotny  - mapowanie: indeks_odwrotny[state_int + 1] = indeks w bazie
"""
function gen_baza_xxz(M::Int, N_up::Int)
    rozmiar_pelnej_bazy = 2^M
    baza_bin = String[]
    baza_int = Int[]
    # Dla stanu o wartości state_int odwołujemy się przez state_int + 1, ponieważ tablice w Julii są indeksowane od 1.
    indeks_odwrotny = zeros(Int, rozmiar_pelnej_bazy)

    for state_int in 0:(rozmiar_pelnej_bazy - 1)
        if count_ones(state_int) == N_up
            state_bin = dec_na_bin(state_int, M)
            push!(baza_bin, state_bin)
            push!(baza_int, state_int)
            indeks_bazy = length(baza_int)
            indeks_odwrotny[state_int + 1] = indeks_bazy
        end
    end
    return baza_bin, baza_int, indeks_odwrotny
end

"""
Zapisuje trzy tablice bazy do osobnych plików w folderze
Lanczos_XXZ_wyniki znajdującym się w głównym katalogu projektu.
"""
function zapisz_baze_do_plikow(M::Int, N_up::Int, baza_bin::Vector{String}, baza_int::Vector{Int}, indeks_odwrotny::Vector{Int})
    output_directory = joinpath(@__DIR__, "Lanczos_XXZ_wyniki")
    mkpath(output_directory)

    plik_baza_bin = joinpath(output_directory, "baza_bin_M_$(M)_Nup_$(N_up).txt")
    plik_baza_int = joinpath(output_directory, "baza_int_M_$(M)_Nup_$(N_up).txt")
    plik_indeks_odwrotny = joinpath(output_directory, "indeks_odwrotny_M_$(M)_Nup_$(N_up).txt")
    # Zapis bazy binarnej
    open(plik_baza_bin, "w") do io
        for stan_bin in baza_bin
            println(io, stan_bin)
        end
    end
    # Zapis bazy całkowitoliczbowej
    open(plik_baza_int, "w") do io
        for stan_int in baza_int
            println(io, stan_int)
        end
    end
    # Zapis tablicy indeksów odwrotnych
    open(plik_indeks_odwrotny, "w") do io
        for indeks in indeks_odwrotny
            println(io, indeks)
        end
    end
    return nothing
end

"""
Zwraca: H - rzadka macierz Hamiltonianu
"""
function gen_ham_xxz(M::Int, J::Real, Delta::Real, pbc::Bool, baza_bin::Vector{String}, indeks_odwrotny::Vector{Int})
    ile_stanow = length(baza_bin)
    H = spzeros(Float64, ile_stanow, ile_stanow) # Tworzymy rzadką macierz Hamiltonianu wpełniną zerami

    for i in 1:ile_stanow
        stan = collect(baza_bin[i]) # Zamieniamy String na tablicę znaków

        for s in 1:M

            # CZĘŚĆ POZADIAGONALNA
            if stan[s] == '1'

                if s > 1 && stan[s - 1] == '0' # 01 -> 10
                    nowy_stan = copy(stan)
                    nowy_stan[s] = '0'
                    nowy_stan[s - 1] = '1'
                    nowy_stan_bin = join(nowy_stan)
                    nowy_stan_int = bin_na_dec(nowy_stan_bin)
                    nowy_indeks = indeks_odwrotny[nowy_stan_int + 1]
                    H[nowy_indeks, i] += J / 2
                end

                if s < M && stan[s + 1] == '0' # 10 -> 01
                    nowy_stan = copy(stan)
                    nowy_stan[s] = '0'
                    nowy_stan[s + 1] = '1'
                    nowy_stan_bin = join(nowy_stan)
                    nowy_stan_int = bin_na_dec(nowy_stan_bin)
                    nowy_indeks = indeks_odwrotny[nowy_stan_int + 1]
                    H[nowy_indeks, i] += J / 2
                end
                
                # PBC: przeskok z pierwszego węzła na ostatni
                if pbc == 1 && M > 2 
                    if s == 1 && stan[M] == '0'
                        nowy_stan = copy(stan)
                        nowy_stan[1] = '0'
                        nowy_stan[M] = '1'
                        nowy_stan_bin = join(nowy_stan)
                        nowy_stan_int = bin_na_dec(nowy_stan_bin)
                        nowy_indeks = indeks_odwrotny[nowy_stan_int + 1]
                        H[nowy_indeks, i] += J / 2
                    end
                end

                # PBC: przeskok z ostatniego węzła na pierwszy
                if pbc == 1 && M > 2
                    if s == M && stan[1] == '0'
                        nowy_stan = copy(stan)
                        nowy_stan[M] = '0'
                        nowy_stan[1] = '1'
                        nowy_stan_bin = join(nowy_stan)
                        nowy_stan_int = bin_na_dec(nowy_stan_bin)
                        nowy_indeks = indeks_odwrotny[nowy_stan_int + 1]
                        H[nowy_indeks, i] += J / 2
                    end
                end
            end

            # CZĘŚĆ DIAGONALNA
            if s < M
                if stan[s] == stan[s + 1]
                H[i, i] += J * Delta/4 # Pary 11 oraz 00
                else
                H[i, i] -= J * Delta/4 # Pary 10 oraz 01
                end
            end

            # PBC
            if pbc == 1 && M > 2
                if s == M
                    if stan[M] == stan[1]
                    H[i, i] += J * Delta/4 # Pary 11 oraz 00 na brzegu
                    else
                    H[i, i] -= J * Delta/4 # Pary 10 oraz 01 na brzegu
                    end
                end
            end
        end
    end
    return H
end

"""
Zapisuje niezerowe elementy rzadkiego Hamiltonianu XXZ.
"""
function save_ham(H::SparseMatrixCSC, M::Int, N_up::Int, J::Real, Delta::Real, pbc_argument::String)
    output_directory = joinpath(@__DIR__, "Lanczos_XXZ_wyniki")
    mkpath(output_directory)

    plik_hamiltonian = joinpath(output_directory, "hamiltonian_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument).txt")
    wiersze, kolumny, wartosci = findnz(H)     # Pobranie indeksów i wartości wszystkich niezerowych elementów.

    open(plik_hamiltonian, "w") do io
        for k in eachindex(wartosci)
            println(io, wiersze[k], '\t', kolumny[k], '\t', wartosci[k])
        end
    end
    return nothing
end

"""
Zapisuje wartości własne Hamiltonianu do pliku.
Każda wartość własna jest zapisywana w osobnym wierszu.
"""
function save_eigenvalues(wartosci_wlasne::AbstractVector, nazwa_tablicy::String, M::Int, N_up::Int, J::Real, Delta::Real, pbc_argument::String)
    output_directory = joinpath(@__DIR__, "Lanczos_XXZ_wyniki")
    mkpath(output_directory)

    plik_wartosci_wlasne = joinpath(output_directory, "$(nazwa_tablicy)_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument).txt")

    open(plik_wartosci_wlasne, "w") do io
        for wartosc in wartosci_wlasne
            println(io, wartosc)
        end     
    end
    return nothing
end

"""
Specifying model parameters
"""
function parse_args()
  s = ArgParseSettings(description = "Compute LIOMs for two-leg XXZ ladder with leg-sector basis")
  @add_arg_table s begin
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
         "-P", "--pbc"
            help = "Periodyczne warunki brzegowe: yes albo no"
            arg_type = String
            required = true
            range_tester = x -> x in ["yes", "no"]
            dest_name = "pbc"
  end
  return ArgParse.parse_args(s)
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
The Lanczos Algorithm
"""
function lanczos(A, q0, q1, beta1)
    # Tworzymy macierz T tej samej wielkości co A
    n = size(A, 1)
    max_iter = min(n, 1000)
    # Normalizujemy wektor początkowy q₁
    q = ComplexF64.(q1)
    q = q / norm(q)
    # Początkowy poprzedni wektor q₀
    q_previous = ComplexF64.(q0)
    # Początkowa wartość β₁
    beta = zeros(Float64, n)
    beta[1] = Float64(beta1)
    # Macierz wektorów Lanczosa ; Q = [q₁ q₂ ... qₙ]
    Q = zeros(ComplexF64, n, n)
    Q[:, 1] = q
    # Współczynniki macierzy trójdiagonalnej - T
    alpha = zeros(Float64, n)
    energie = zeros(Float64, max_iter, 2)

    for j in 1:max_iter
        # println("\n==============================")
        # println("ITERACJA: ", j)
        # println("==============================")
        # println("q_$(j - 1) = ")
        # display(q_previous)
        # println("q_$j = ")
        # display(q)
        # println("beta_$j = ", beta[j])
        # αⱼ = qⱼᵀ A qⱼ
        alpha[j] = real(dot(q, A * q))
        # println("alpha_$j = ", alpha[j])
        T_j = SymTridiagonal(alpha[1:j], beta[2:j])
        energie[j, 1] = j
        energie[j, 2] = eigmin(T_j)
        # wⱼ₊₁ = A qⱼ - βⱼ qⱼ₋₁ - αⱼ qⱼ
        w = A * q - beta[j] * q_previous - alpha[j] * q
        # println("w_$(j + 1) = ")
        # display(w)
        if j < max_iter 
            # βⱼ₊₁ = ||wⱼ₊₁||
            beta[j + 1] = norm(w)
            # println("beta_$(j + 1) = ", beta[j + 1])
            if beta[j + 1] < 1e-14 # jesil beta robi się zamałe algorytm się zatrzymuje 
                T = SymTridiagonal(alpha[1:j], beta[2:j])
                return T, Q[:, 1:j], energie[1:j, :]
            end
            # qⱼ₊₁ = wⱼ₊₁ / βⱼ₊₁
            q_next = w / beta[j + 1]
            # println("q_$(j + 1) = ")
            # display(q_next)
            # Zapisujemy qⱼ₊₁ w macierzy Q
            Q[:, j + 1] = q_next
            # Przygotowujemy następną iterację
            q_previous = q
            q = q_next
        end
    end
    # Budujemy macierz trójdiagonalną
    T = SymTridiagonal(alpha[1:max_iter], beta[2:max_iter])
    return T, Q[:, 1:max_iter], energie
end

"""
The Lanczos Algorithm for time evolution
"""
function dynamika_lanczos(H, psi0, delta_t, liczba_krokow)
    n = size(H, 1)
    psi = ComplexF64.(psi0)
    stany = zeros(ComplexF64, n, liczba_krokow)

    for krok in 1:liczba_krokow
        q0 = zeros(ComplexF64, n)
        beta_1 = 0.0
        T, Q, energie = lanczos(H, q0, psi, beta_1)
        # Diagonalizacja macierzy T: T = V D V†
        wynik_T = eigen(T)
        D = Diagonal(wynik_T.values)
        V = wynik_T.vectors
        # Stan początkowy zapisany w bazie Kryłowa
        e1 = zeros(ComplexF64, size(T, 1))
        e1[1] = 1.0
        # Ewolucja o jeden krok czasowy
        psi = Q * V * exp(-im * D * delta_t) * adjoint(V) * e1

        stany[:, krok] = psi
    end

    return stany
end

function main()
    args = parse_args()
    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    pbc_argument = args["pbc"]
    pbc = pbc_argument == "yes"

    baza_bin, baza_int, indeks_odwrotny = gen_baza_xxz(M, N_up)
    # println("Wygenerowano bazę modelu XXZ.")
    # zapisz_baze_do_plikow(M, N_up, baza_bin, baza_int, indeks_odwrotny)
    # println("Zapisano do plików bazę modelu XXZ.")

    H = gen_ham_xxz(M, J, Delta, pbc, baza_bin, indeks_odwrotny)
    # println("Stworzono Hamiltonian modelu XXZ.")    
    # save_ham(H, M, N_up, J, Delta, pbc_argument)
    # println("Hamiltonian zapisano do pliku.")

    # ED
    # wynik = eigen(Matrix(H))
    # wartosci_wlasne = wynik.values
    # # wektory_wlasne = wynik.vectors
    # save_eigenvalues(wartosci_wlasne, "wartosci_wlasne_ED", M, N_up, J, Delta, pbc_argument)
    # println("Wartości własne zapisano do pliku.")

    # Lanczos
    # q0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    q0 = zeros(size(H, 1))
    # q1 = ComplexF64[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    q1 = losuj_q1(size(H, 1))
    beta_1 = 0.0
    # println("Wylosowany q1:")
    # display(q1)
    T, Q, energie = lanczos(H, q0, q1, beta_1)
    # println("Macierz T:")
    # display(Matrix(T))
    wynik_LAN = eigen(Matrix(T))
    wartosci_wlasne_LAN = wynik_LAN.values
    # wektory_wlasne = wynik.vectors
    save_eigenvalues(wartosci_wlasne_LAN, "wartosci_wlasne_LAN", M, N_up, J, Delta, pbc_argument)
    println("Wartości własne zapisano do pliku.")

    # Porównanie wyników z ED i Lanczosa
    E_ED = eigmin(Matrix(H))
    roznice_energii = zeros(Float64, size(energie, 1), 2)
    for j in axes(energie, 1)
        roznice_energii[j, 1] = energie[j, 1]
        roznice_energii[j, 2] = abs(energie[j, 2] - E_ED)
    end
    output_directory = joinpath(@__DIR__, "Lanczos_XXZ_wyniki")
    plik_roznice_energii = joinpath(output_directory, "roznica_energii_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument).txt")
    open(plik_roznice_energii, "w") do io
        for j in axes(roznice_energii, 1)
            println(io, roznice_energii[j, 1], " ", roznice_energii[j, 2])
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end