using Pkg
using ArgParse
using SparseArrays
using LinearAlgebra
using ProgressMeter

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
Wykonuje cykliczną translację stanu binarnego o jedno miejsce w lewo.
Przykład: translacja_stanu("0011") == "0110"
"""
function translacja_stanu(stan::String)::String
    return stan[2:end] * stan[1]
end

"""
Wykonuje cykliczną translację stanu binarnego o jedno miejsce w lewo.
Przykład: "0011" -> "0110"
"""
function translacja_stanu(stan::String)::String
    return stan[2:end] * stan[1]
end

"""
Generuje bazę stanów translacyjnych dla układu o M miejscach i N_up spinach w górę.
taki sposób zapisu tyle samo elmentów co stanów w bazie położeniowej
3 6 12 9
5 10
"""
function gen_trans_baza(M::Int, N_up::Int, k::Int)
    baza = Vector{Vector{Int}}()
    wykorzystane = Set{Int}() # Set - bez duplikatów oraz kolejność nie jest ważna
    for state_int in 0:(2^M - 1)
    stan = dec_na_bin(state_int, M)
        if count_ones(state_int) != N_up # Interesuje nas tylko sektor z ustalonym N_up jeśli jest inna pzrechodzimy do naśtepnego obrotu pętli
            continue
        end
        if state_int in wykorzystane # Jeżeli stan został już znaleziony jako translacja wcześniejszego stanu, pzrechodzimy do naśtepnego obrotu pętli
            continue
        end
            element_bazy = Int[] # wektor Int na nowe elementy
            aktualny_stan = stan

            push!(element_bazy, bin_na_dec(aktualny_stan)) 
            push!(wykorzystane, bin_na_dec(aktualny_stan))
            aktualny_stan = translacja_stanu(aktualny_stan)

            while aktualny_stan != stan

                push!(element_bazy, bin_na_dec(aktualny_stan))
                push!(wykorzystane, bin_na_dec(aktualny_stan))
                aktualny_stan = translacja_stanu(aktualny_stan)

            end            
            R_a = length(element_bazy)
            if mod(k, M ÷ R_a) != 0 # warunek na k 
                continue
            end
            push!(baza, element_bazy) # dodanie całego wiersza z translacjami do bazy 
    end
    return baza
end

"""
Generuje reprezentantów orbit translacyjnych dla układu o M miejscach i N_up spinach w górę.
(a - reprezentant, R_a - okres translacji)
"""
function gen_trans_baza_sandvik(M::Int, N_up::Int, k::Int)
    baza = Vector{Tuple{Int,Int}}()
    for state_int in 0:(2^M - 1)
        if count_ones(state_int) != N_up # Interesuje nas tylko sektor z ustalonym N_up jeśli jest inna pzrechodzimy do naśtepnego obrotu pętli
            continue
        end
        stan = dec_na_bin(state_int, M) # analizowany kadydat na reprezentanta
        aktualny_stan = translacja_stanu(stan) # ta funkcja przyjmuje stringa
        reprezentant = state_int # wektor Int na nowe elementy
        okres_translacji = 1 # początkowy okres translacji =1 , bo już została wykonana jedna translacja

        while aktualny_stan != stan
            aktualny_stan_int = bin_na_dec(aktualny_stan)
            if aktualny_stan_int < reprezentant # reprezentami śa najmnijesze inty z danej grupy
                reprezentant = aktualny_stan_int
            end
            # Kolejna translacja stanu
            okres_translacji += 1
            aktualny_stan = translacja_stanu(aktualny_stan)
        end
        if mod(k, M ÷ okres_translacji) != 0 # warunek na k 
            continue
        end
        # Zapisz, gdy sam jest, najmniejszym elementem swojej grupy translacyjnej
        if state_int == reprezentant
            push!(baza, (reprezentant, okres_translacji))
        end
    end
    return baza
end

"""
Dla podanego stanu znajduje reprezentanta oraz liczbę translacji l 
"""
function znajdz_reprezentanta_oraz_l(stan_int::Int, M::Int)
    stan = dec_na_bin(stan_int, M)
    aktualny_stan = stan
    reprezentant = stan_int
    l = 0
    aktualne_l = 0

    aktualny_stan = translacja_stanu(aktualny_stan) # pierwsza translacja
    aktualne_l += 1

    while aktualny_stan != stan
        aktualny_stan_int = bin_na_dec(aktualny_stan)
        if aktualny_stan_int < reprezentant
            reprezentant = aktualny_stan_int
            l = aktualne_l
        end
        aktualny_stan = translacja_stanu(aktualny_stan)
        aktualne_l += 1
    end
    return reprezentant, l # zwracamy dwie rzeczy
end

"""
Dla bazy: [(a_1, R_1), (a_2, R_2), ...]
zwraca słownik: reprezentant => indeks w bazie
Przyklad: baza = [(3, 4), (5, 2), (17, 8)] - > 3  => 1; 5  => 2; 17 => 3
"""
function gen_indeks_odwrotny(baza::Vector{Tuple{Int,Int}})
    indeks_odwrotny = Dict{Int,Int}()
    for i in 1:length(baza)
        reprezentant = baza[i][1]
        indeks_odwrotny[reprezentant] = i
    end
    return indeks_odwrotny
end

"""
Generuje macierz Hamiltonianu XXZ w bazie stanów translacyjnych.
"""
function gen_ham_XXZ_ms(M::Int, J::Real, Delta::Real, k::Int, baza::Vector{Tuple{Int,Int}})
    ile_stanow = length(baza)
    H = spzeros(ComplexF64, ile_stanow, ile_stanow)
    A_1 = spzeros(ComplexF64, ile_stanow, ile_stanow)
    J_1 = spzeros(ComplexF64, ile_stanow, ile_stanow)
    indeks_odwrotny = gen_indeks_odwrotny(baza)     # Tworzenie słownika odwrotnego, aby szybko znaleźć indeks reprezentanta w bazie
    k_fiz = 2π * k / M # fizyczna wartość momentum
    # główna pętla
    for i in 1:ile_stanow
        a = baza[i][1]
        R_a = baza[i][2]
        stan = collect(dec_na_bin(a, M)) # zamienia na chary każdego stringa
        for s in 1:M
            # Numer następnego miejsca w danym wiązaniu, uwzględnienie periodyczności
            if s < M
                s2 = s + 1
            else
                s2 = 1
            end
            # CZĘŚĆ POZADIAGONALNA; 01 -> 10 lub 10 -> 01
            if stan[s] != stan[s2]
                nowy_stan = copy(stan)
                nowy_stan[s] = stan[s2]
                nowy_stan[s2] = stan[s]
                # powrót do stringa i konwersja na int
                nowy_stan_bin = join(nowy_stan)
                nowy_stan_int = bin_na_dec(nowy_stan_bin)
                # Stan po działaniu H nie musi być reprezentantem, znajdujemy jego reprezentanta b oraz l
                reprezentant_b, l = znajdz_reprezentanta_oraz_l(nowy_stan_int, M)
                # Szukamy reprezentanta w aktualnym sektorze momentum
                nowy_indeks = get(indeks_odwrotny, reprezentant_b, 0)
                # Jeżeli reprezentant istnieje w tym sektorze, dodajemy część pozadiagonalną
                if nowy_indeks != 0
                    R_b = baza[nowy_indeks][2]
                    # <b(k)|H_j|a(k)> = h_j(a) exp(-i q l) sqrt(R_a/R_b)
                    H[nowy_indeks, i] += (J / 2) * exp(-im * k_fiz * l) * sqrt(R_a / R_b)
                    # Operator prądu J_1
                    if stan[s] == '1' && stan[s2] == '0'
                        # 10 -> 01
                        J_1[nowy_indeks, i] += (im * J / 2) * exp(-im * k_fiz * l) * sqrt(R_a / R_b)
                    elseif stan[s] == '0' && stan[s2] == '1'
                        # 01 -> 10
                        J_1[nowy_indeks, i] += (-im * J / 2) * exp(-im * k_fiz * l) * sqrt(R_a / R_b)
                    end
                end
            end
            # CZĘŚĆ DIAGONALNA
            if stan[s] == stan[s2]
                H[i, i] += J * Delta / 4
                A_1[i, i] += 1 / 4
            else
                H[i, i] -= J * Delta / 4
                A_1[i, i] -= 1 / 4
            end
        end
    end

    return H, A_1, J_1
end

"""
Zapisuje niezerowe elementy rzadkiego Hamiltonianu w bazie momentum.
"""
function save_ham(H::SparseMatrixCSC, M::Int, N_up::Int, J::Real, Delta::Real, k::Int)
    output_directory = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state"
    mkpath(output_directory)
    plik_hamiltonian = joinpath(output_directory, "hamiltonian_XXZ_ms_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_k_$(k).txt")
    wiersze, kolumny, wartosci = findnz(H)
    open(plik_hamiltonian, "w") do io
        for i in eachindex(wartosci)
            println(io, wiersze[i], '\t', kolumny[i], '\t', wartosci[i])
        end
    end
    println("Hamiltonian został zapisany do pliku:")
    println(plik_hamiltonian)
    return nothing
end

"""
Zapisuje wartości własne Hamiltonianu do pliku.
Każda wartość własna jest zapisywana w osobnym wierszu.
"""
function save_eigenvalues(wartosci_wlasne::AbstractVector, nazwa_tablicy::String, M::Int, N_up::Int, J::Real, Delta::Real, k::Int)
    output_directory = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state"
    mkpath(output_directory)
    plik_wartosci_wlasne = joinpath(output_directory, "$(nazwa_tablicy)_XXZ_ms_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_k_$(k).txt")
    open(plik_wartosci_wlasne, "w") do io
        for wartosc in wartosci_wlasne
            println(io, wartosc)
        end
    end
    println("Wartości własne zostały zapisane do pliku:")
    println(plik_wartosci_wlasne)
    return nothing
end

"""
Zapisuje wszystkie wartości własne ze wszystkich sektorów momentum. Pierwsza kolumna: k_fiz, druga kolumna: wartość własna
"""
function save_all_eigenvalues(wszystkie_wartosci_wlasne::Vector{Tuple{Float64,Float64}},M::Int, N_up::Int, J::Real, Delta::Real)
    output_directory = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state"
    mkpath(output_directory)
    plik_wartosci_wlasne = joinpath(output_directory, "wartosci_wlasne_wszystkie_k_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta).txt")
    open(plik_wartosci_wlasne, "w") do io
        for (k_fiz, wartosc) in wszystkie_wartosci_wlasne
            println(io, k_fiz, '\t', wartosc)
        end
    end
    println()
    println("Wszystkie wartości własne zostały zapisane do pliku:")
    println(plik_wartosci_wlasne)
    return nothing
end

"""
Generates a normalized random state vector in the computational basis, with amplitudes drawn from a Gaussian distribution.
"""
function gen_losowy_stan_gaussowski(M::Int, N_up::Int)
    rozmiar_bazy = binomial(M, N_up)
    psi0 = randn(ComplexF64, rozmiar_bazy)
    norma_przed = norm(psi0)
    psi0 ./= norma_przed
    return psi0
end

"""
The Lanczos Algorithm
"""
function lanczos(A, q0, q1, beta1)
    # Tworzymy macierz T tej samej wielkości co A
    n = size(A, 1)
    max_iter = min(n, 5)
    # Normalizujemy wektor początkowy q₁
    q = ComplexF64.(q1)
    q = q / norm(q)
    # println("norma_na_początku_lanczosa = ", norm(q))
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
The Lanczos Algorithm for one time evolution step
"""
function dynamika_lanczos(H, psi0, delta_t)
    n = size(H, 1)
    psi0 = ComplexF64.(psi0)
    q0 = zeros(ComplexF64, n)
    beta_1 = 0.0
    norma_psi0 = norm(psi0) # Zapamiętujemy normę stanu wejściowego
    T, Q, energie = lanczos(H, q0, psi0, beta_1) # Budujemy przestrzeń Kryłowa dla aktualnego stanu
    # Diagonalizacja macierzy T = V D V†
    wynik_T = wynik_T = eigen(Symmetric(Matrix(T));alg = LinearAlgebra.DivideAndConquer())
    D = Diagonal(wynik_T.values)
    V = wynik_T.vectors
    # Stan początkowy w bazie Kryłowa
    e1 = zeros(ComplexF64, size(T, 1))
    e1[1] = 1.0
    # Jednorazowa ewolucja: t -> t + delta_t. Lanczos pracuje na znormalizowanym q1, dlatego na końcu przywracamy normę stanu wejściowego
    psi_koncowe = norma_psi0 * Q * V * exp(-im * D * delta_t) * adjoint(V) * e1

    return psi_koncowe
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
function dynamika_ED(H, Psi_0, Phi_0, t_min, delta_t, liczba_krokow)
    n = size(H, 1)
    Psi_0 = ComplexF64.(Psi_0)
    Phi_0 = ComplexF64.(Phi_0)
    stany_Psi_ED = zeros(ComplexF64, n, liczba_krokow)
    stany_Phi_ED = zeros(ComplexF64, n, liczba_krokow)
    wynik_H = eigen(Hermitian(Matrix(H))) # Diagonalizacja Hamiltonianu H = V D V†
    E = wynik_H.values
    V = wynik_H.vectors
    D = Diagonal(E)
    # Stany początkowe w bazie własnej Hamiltonianu
    Psi_wlasna = adjoint(V) * Psi_0
    Phi_wlasna = adjoint(V) * Phi_0

    for krok in 1:liczba_krokow
        czas = t_min + krok * delta_t
        ewolucja = exp(-im * D * czas)
        stany_Psi_ED[:, krok] = V * ewolucja * Psi_wlasna
        stany_Phi_ED[:, krok] = V * ewolucja * Phi_wlasna
    end
    return stany_Psi_ED, stany_Phi_ED, E, V
end

"""
Correlation function using Exact Diagonalization
"""
function funkcja_korelacji_ED(A, Psi_t, Phi_t)
    C_t = dot(Psi_t, A * Phi_t)
    return real(C_t)
end

"""
Correlation function using Exact Diagonalization T->inf equation with trace
"""
function funkcja_korelacji_trace_ED(E, A_wlasna, t::Real)
    C_t = 0.0 + 0.0im
    for n in eachindex(E)
        for l in eachindex(E)
            C_t += abs2(A_wlasna[n, l]) * exp(im * (E[n] - E[l]) * t) # C(t) = sum_{n,l} |<n|A|l>|² exp[i(Eₙ - Eₗ)t]
        end
    end
    return real(C_t)
end

"""
Przyjmuje parametry z terminala
"""
function parse_args()
    s = ArgParseSettings(description = "Compute XXZ Hamiltonian in momentum-state basis")
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
            help = "Stała sprzężenia J"
            arg_type = Float64
            default = 1.0
            dest_name = "J"
        "-d", "--delta"
            help = "Anizotropia Delta"
            arg_type = Float64
            default = 0.3
            dest_name = "Delta"
        "-k"
            help = "Numer sektora momentum k"
            arg_type = Int
            default = 0
            dest_name = "k"
        "--all-k"
            help = "Jeśli true, licz dla wszystkich sektorów k = 0,...,M-1; jeśli false, licz tylko dla podanego k"
            arg_type = Bool
            default = false
            dest_name = "all_k"
        "-R", "--realizations"
            help = "Liczba losowych stanów Psi_0 użytych do uśrednienia Lanczosa"
            arg_type = Int
            default = 1
            dest_name = "liczba_realizacji"
    end
    return ArgParse.parse_args(s)
end

function main()
    # args = parse_args()
    # M = args["M"]
    # N_up = args["N_up"]
    # J = args["J"]
    # Delta = args["Delta"]
    # k = args["k"]
    # all_k = args["all_k"]
    # folder_wyniki = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state"
    # mkpath(folder_wyniki)
    # # czy pętla po k czy pojedńcze k 
    # if all_k == true
    #     lista_k = 0:(M - 1)
    # else
    #     lista_k = [k]
    # end

    # ---- Dyspersja ----
    # wszystkie_wartosci_wlasne = Vector{Tuple{Float64,Float64}}()
    # for k in lista_k
    #     czas_k = @elapsed begin
    #     println("===================================")
    #     println("Liczenie sektora k = ", k)
    #     println("===================================")

    #     baza = gen_trans_baza_sandvik(M, N_up, k) # Generowanie bayzy
    #     println("Liczba stanów w bazie: ", length(baza))

    #     nazwa_pliku_baza = joinpath(folder_wyniki, "baza_trans_sandvik_M_$(M)_N_$(N_up)_k_$(k).txt") # Zapis bazy
    #     open(nazwa_pliku_baza, "w") do plik
    #         for element_bazy in baza
    #             println(plik, join(element_bazy, " "))
    #         end
    #     end
    #     println("Baza została zapisana do pliku:")
    #     println(nazwa_pliku_baza)

    #     H, A_1, J_1 = gen_ham_XXZ_ms(M, J, Delta, k, baza) # Generowanie Hamiltonianu
    #     # save_ham(H, M, N_up, J, Delta, k) # Zapis Hamiltonianu

    #     eigen_result = eigen(Hermitian(Matrix(H))) # Diagonalizacja, z powiedzeniem funkcji że jest hermitowska
    #     wartosci_wlasne = eigen_result.values
    #     save_eigenvalues(wartosci_wlasne, "wartosci_wlasne", M, N_up, J, Delta, k) # Zapis wartości własnych

    #     k_fiz = 2π * k / M # fizyczna wartość pędu
    #     for wartosc in wartosci_wlasne
    #         push!(wszystkie_wartosci_wlasne, (k_fiz, wartosc))
    #     end
    # end
    # println("Czas liczenia dla k = $(k): $(czas_k) s")
    # end  # @elapsed begin
    # save_all_eigenvalues(wszystkie_wartosci_wlasne, M, N_up, J, Delta)

    # # ---- Lanczos - funkcja korelacji, ED - funkcja korelacji oraz funkcja korelacji liczona z Tracem ----
    # t_min = 0.0
    # t_max = 50.0
    # liczba_krokow_czasu = 5000
    # delta_t = (t_max - t_min) / liczba_krokow_czasu
    # korelacja_suma = zeros(Float64, liczba_krokow_czasu + 1, 4)
    # # Jeśli liczymy wszystkie k, generujemy jeden znormalizowany stan w całej przestrzeni
    # if all_k == true
    #     Psi_0 = gen_losowy_stan_gaussowski(M, N_up)
    #     pierwszy_indeks = 1
    # end
    # # Pętla po k 
    # for k in lista_k
    #     czas_k = @elapsed begin
    #         println()
    #         println("Liczenie sektora k = ", k)
    #         println("===================================")
    #         baza = gen_trans_baza_sandvik(M, N_up, k)
    #         rozmiar_sektora = length(baza)
    #         println("Rozmiar sektora = ", rozmiar_sektora)
    #         # Stan początkowy, fragment Gaussa znormalizowanego w całej przestrzeni
    #         if all_k == true
    #             ostatni_indeks = pierwszy_indeks + rozmiar_sektora - 1
    #             Psi = copy(Psi_0[pierwszy_indeks:ostatni_indeks])
    #             pierwszy_indeks = ostatni_indeks + 1
    #         # Gaussian tylko w wybranym sektorze k i normalizacja w tym sektorze
    #         else
    #             Psi = randn(ComplexF64, rozmiar_sektora)
    #             Psi ./= norm(Psi)
    #         end
    #         Psi_poczatkowe = copy(Psi)
    #         # println("Norma Psi_k = ", norm(Psi))
    #         H, A_1, J_1 = gen_ham_XXZ_ms(M, J, Delta, k, baza)
    #         Phi = J_1 * Psi
    #         Phi_poczatkowe = copy(Phi)
    #         # println("Norma Phi_k = ", norm(Phi))
            
    #         # ---- Lanczos ----
    #         # t = 0
    #         korelacja_suma[1, 1] = t_min
    #         C_0_k = real(dot(Psi, J_1 * Phi))
    #         korelacja_suma[1, 2] += C_0_k
    #         # Ewolucja czasowa
    #         for krok in 1:liczba_krokow_czasu
    #             czas = t_min + krok * delta_t
    #             C_t_k, Psi, Phi = funkcja_korelacji(H, J_1, delta_t, Psi, Phi)
    #             korelacja_suma[krok + 1, 1] = czas
    #             korelacja_suma[krok + 1, 2] += C_t_k
    #         end
    #     end
    #     println("Czas Lanczosa dla k = $(k): ", round(czas_k, digits = 3), " s")

    #         # ---- ED zgodność z Lanczosem ----
    #         czas_ED_k = @elapsed begin
    #         stany_Psi_ED, stany_Phi_ED, E_ED, V_ED = dynamika_ED(H, Psi_poczatkowe, Phi_poczatkowe, t_min, delta_t, liczba_krokow_czasu)
    #         # t = 0
    #         C_0_ED_k = funkcja_korelacji_ED(J_1, Psi_poczatkowe, Phi_poczatkowe)
    #         korelacja_suma[1, 3] += C_0_ED_k
    #         # Ewolucja czasowa
    #         for krok in 1:liczba_krokow_czasu
    #             Psi_t_ED = stany_Psi_ED[:, krok]
    #             Phi_t_ED = stany_Phi_ED[:, krok]
    #             C_t_ED_k = funkcja_korelacji_ED(J_1, Psi_t_ED, Phi_t_ED)
    #             korelacja_suma[krok + 1, 3] += C_t_ED_k
    #         end
    #     end
    #     println("Czas ED dla k = $(k): ", round(czas_ED_k, digits = 3), " s")

    #         # ---- ED zgodność z Tr T->inf ----
    #         czas_trace_k = @elapsed begin
    #         stany_Psi_ED, stany_Phi_ED, E_ED, V_ED = dynamika_ED(H, Psi_poczatkowe, Phi_poczatkowe, t_min, delta_t, liczba_krokow_czasu)
    #         J_1_wlasna = adjoint(V_ED) * Matrix(J_1) * V_ED
    #         # t = 0
    #         C_0_trace_k = funkcja_korelacji_trace_ED(E_ED, J_1_wlasna, t_min)
    #         korelacja_suma[1, 4] += C_0_trace_k
            
    #         for krok in 1:liczba_krokow_czasu
    #             czas = t_min + krok * delta_t
    #             C_t_trace_k = funkcja_korelacji_trace_ED(E_ED, J_1_wlasna, czas)
    #             korelacja_suma[krok + 1, 4] += C_t_trace_k
    #         end        
    #     end
    #     println("Czas Trace ED dla k = $(k): ", round(czas_trace_k, digits = 3), " s")

    # end
    #     if all_k == true
    #         Z_trace = binomial(M, N_up)
    #     else
    #         Z_trace = size(H, 1)
    #     end
    #     korelacja_suma[:, 4] ./= Z_trace

    # # Zapis funkcji korelacji do pliku
    # output_directory = joinpath(@__DIR__, "momnetu_state")
    # pbc_argument = "yes"
    # operator_argument = "J_1"
    # mkpath(output_directory)
    # plik_korelacja = joinpath(output_directory,"cor_ms_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument)_operator_$(operator_argument).txt")
    # open(plik_korelacja, "w") do io
    #     for i in axes(korelacja_suma, 1)
    #         println(io, korelacja_suma[i, 1], " ", korelacja_suma[i, 2], " ", korelacja_suma[i, 3], " ", korelacja_suma[i, 4])
    #     end
    # end
    # println("Funkcja korelacji została zapisana do pliku:")
    # println(plik_korelacja)

    # ---- Lanczos z uśrednieniem porównie z ED Tr ----
    args = parse_args()
    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    k = args["k"]
    all_k = args["all_k"]
    liczba_realizacji = args["liczba_realizacji"]
    folder_wyniki = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state"
    mkpath(folder_wyniki)
    if all_k == true
        lista_k = 0:(M - 1)
    else
        lista_k = [k]
    end

    # ---- Funkcja korelacji z uśrednieniami ----
    t_min = 0.0
    t_max = 50.0
    liczba_krokow_czasu = 5000
    delta_t = (t_max - t_min) / liczba_krokow_czasu
    korelacja_suma = zeros(Float64, liczba_krokow_czasu + 1, 4)
    for krok in 0:liczba_krokow_czasu
        korelacja_suma[krok + 1, 1] = t_min + krok * delta_t
    end
    # poporawne ustawenie pierwszy_index
    println("Liczba realizacji Psi_0 = ", liczba_realizacji)
    Psi_0_lista = [gen_losowy_stan_gaussowski(M, N_up) for realizacja in 1:liczba_realizacji]
    pierwszy_indeks = 1
    if all_k == false
        for k_tmp in 0:(k - 1)
            pierwszy_indeks += length(gen_trans_baza_sandvik(M, N_up, k_tmp))
        end
    end

    # Pętla po sektorach momentum
    for k in lista_k
        println()
        println("Liczenie sektora k = ", k)
        println("===================================")
        baza = gen_trans_baza_sandvik(M, N_up, k)
        rozmiar_sektora = length(baza)
        println("Rozmiar sektora = ", rozmiar_sektora)
        H, A_1, J_1 = gen_ham_XXZ_ms(M, J, Delta, k, baza)
        ostatni_indeks = pierwszy_indeks + rozmiar_sektora - 1

        # Lanczos - średnia po losowych stanach
        korelacja_lanczos_k = zeros(Float64, liczba_krokow_czasu + 1)
        czas_lanczos_k = @elapsed begin
        for realizacja in 1:liczba_realizacji
            println("Lanczos: realizacja $(realizacja)/$(liczba_realizacji)")
            Psi = copy(Psi_0_lista[realizacja][pierwszy_indeks:ostatni_indeks])
            Phi = J_1 * Psi
            C_0_k = real(dot(Psi, J_1 * Phi))
            korelacja_lanczos_k[1] += C_0_k
            for krok in 1:liczba_krokow_czasu
                C_t_k, Psi, Phi = funkcja_korelacji(H, J_1, delta_t, Psi, Phi)
                korelacja_lanczos_k[krok + 1] += C_t_k
            end
        end
        end
        korelacja_lanczos_k ./= liczba_realizacji
        korelacja_suma[:, 2] .+= korelacja_lanczos_k
        println("Czas Lanczosa dla k = $(k): ", round(czas_lanczos_k, digits = 3), " s")

        # ED Trace T -> inf
        czas_trace_k = @elapsed begin
            wynik_H = eigen(Hermitian(Matrix(H)))
            E_ED = wynik_H.values
            V_ED = wynik_H.vectors
            J_1_wlasna = adjoint(V_ED) * Matrix(J_1) * V_ED
            C_0_trace_k = funkcja_korelacji_trace_ED(E_ED, J_1_wlasna, t_min)
            korelacja_suma[1, 4] += C_0_trace_k
            for krok in 1:liczba_krokow_czasu
                czas = t_min + krok * delta_t
                C_t_trace_k = funkcja_korelacji_trace_ED(E_ED, J_1_wlasna, czas)
                korelacja_suma[krok + 1, 4] += C_t_trace_k
            end
        end
        println("Czas Trace ED dla k = $(k): ", round(czas_trace_k, digits = 3), " s")
        if all_k == true
            pierwszy_indeks = ostatni_indeks + 1
        end
    end
    # Normalizacja wyrażenia z tracem 
    if all_k == true
        Z_trace = binomial(M, N_up)
    else
        Z_trace = size(H, 1)
    end
    korelacja_suma[:, 4] ./= Z_trace


    # Zapis
    output_directory = joinpath(@__DIR__, "momnetu_state")
    pbc_argument = "yes"
    operator_argument = "J_1"
    mkpath(output_directory)
    plik_korelacja = joinpath(output_directory,"cor_ms_XXZ_M_$(M)_Nup_$(N_up)_J_$(J)_Delta_$(Delta)_PBC_$(pbc_argument)_operator_$(operator_argument)_R_$(liczba_realizacji).txt")
    open(plik_korelacja, "w") do io
        for i in axes(korelacja_suma, 1)
            println(io, korelacja_suma[i, 1], " ", korelacja_suma[i, 2], " ", korelacja_suma[i, 3], " ", korelacja_suma[i, 4])
        end
    end
    println("Funkcja korelacji została zapisana do pliku:")
    println(plik_korelacja)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end