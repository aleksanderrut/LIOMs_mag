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
    # Tworzenie słownika odwrotnego, aby szybko znaleźć indeks reprezentanta w bazie
    indeks_odwrotny = gen_indeks_odwrotny(baza)
    k_fiz = 2π*k/M # fizyczna wartość momentum
    # głowna pętla
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
                # powrót to stringa i konwersja na int
                nowy_stan_bin = join(nowy_stan)
                nowy_stan_int = bin_na_dec(nowy_stan_bin)

                # Stan po działaniu H nie musi być reprezentantem, zajdujemy jego reprezentanta b oraz l:
                reprezentant_b, l = znajdz_reprezentanta_oraz_l(nowy_stan_int, M)

                # Jeżeli reprezentanta nie ma w bazie dla danego k, to stan nie jest kompatybilny z tym sektorem momentum.
                nowy_indeks = get(indeks_odwrotny, reprezentant_b, 0)
                if nowy_indeks == 0
                    continue
                end
                R_b = baza[nowy_indeks][2]
                # <b(k)|H_j|a(k)> = h_j(a) exp(-i q l) sqrt(R_a/R_b)
                H[nowy_indeks, i] += (J / 2) * exp(-im * k_fiz * l) * sqrt(R_a / R_b)
            end
            # CZĘŚĆ DIAGONALNA
            if stan[s] == stan[s2]
                H[i, i] += J * Delta / 4 # Pary 11 oraz 00
            else
                H[i, i] -= J * Delta / 4 # Pary 10 oraz 01
            end
        end
    end
    return H
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
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_args()
    M = args["M"]
    N_up = args["N_up"]
    J = args["J"]
    Delta = args["Delta"]
    k = args["k"]
    all_k = args["all_k"]
    folder_wyniki = raw"C:\Users\aleks\Desktop\praca magisterska\LIOMs_mag\scripts\momnetu_state"
    mkpath(folder_wyniki)

    # czy pętla po k czy pojedńcze k 
    if all_k == true
        lista_k = 0:(M - 1)
    else
        lista_k = [k]
    end

    for k in lista_k

        println()
        println("===================================")
        println("Liczenie sektora k = ", k)
        println("===================================")

        # Generowanie bazy
        baza = gen_trans_baza_sandvik(M, N_up, k)

        println("Liczba stanów w bazie: ", length(baza))

        # Zapis bazy
        nazwa_pliku_baza = joinpath(
            folder_wyniki,
            "baza_trans_sandvik_M_$(M)_N_$(N_up)_k_$(k).txt"
        )

        open(nazwa_pliku_baza, "w") do plik
            for element_bazy in baza
                println(plik, join(element_bazy, " "))
            end
        end

        println("Baza została zapisana do pliku:")
        println(nazwa_pliku_baza)

        # Generowanie Hamiltonianu
        H = gen_ham_XXZ_ms(M, J, Delta, k, baza)

        # Zapis Hamiltonianu
        save_ham(H, M, N_up, J, Delta, k)

        # Diagonalizacja
        eigen_result = eigen(Matrix(H))
        wartosci_wlasne = eigen_result.values

        # Zapis wartości własnych
        save_eigenvalues(
            wartosci_wlasne,
            "wartosci_wlasne",
            M,
            N_up,
            J,
            Delta,
            k
        )
    end
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end