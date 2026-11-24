import PauliStrings as ps
using LinearAlgebra
using Printf
using DelimitedFiles
using ArgParse

"""
0 = identity
1 = S+
2 = Sz
3 = S-
"""
symbol_map = Dict(0 => "1", 1 => "S+", 2 => "Sz", 3 => "S-")

@inline hc(o::Int) = (o == 1 ? 3 : (o == 3 ? 1 : o))

function parse_digits_string(s::String)
  return [parse(Int, string(c)) for c in collect(s)]
end

function parse_label(label::String)
  parts = split(label, "|")
  length(parts) != 2 && error("Wrong label format. Expected e.g. 130|222")
  upper = parse_digits_string(parts[1])
  lower = parse_digits_string(parts[2])
  length(upper) != length(lower) && error("Upper and lower labels must have the same length.")
  return upper, lower
end

function sector_to_symbols(sector::String)
  length(sector) != 2 && error("Sector must have length 2, e.g. RI, IR, RR, II")
  upper_T = sector[1] == 'R' ? :even : :odd
  lower_T = sector[2] == 'R' ? :even : :odd
  return upper_T, lower_T
end

"""
Projected one-leg object.
R/even: O + O†
I/odd:  i(O - O†)
"""
function leg_terms(op_list::Vector{Int}, time_reversal::Symbol)
  op_list_hc = hc.(op_list)

  if time_reversal == :even
    return [(op_list, 1.0 + 0.0im), (op_list_hc, 1.0 + 0.0im)]
  elseif time_reversal == :odd
    return [(op_list, 0.0 + 1.0im), (op_list_hc, 0.0 - 1.0im)]
  else
    error("Unknown time-reversal sector: $time_reversal")
  end
end

# ============================================================
# SKOŃCZONY HAMILTONIAN JAKO ps.Operator, BEZ OperatorTS
# ============================================================

"""
Finite two-leg XXZ ladder Hamiltonian.
Tu NIE MA OperatorTS.
Każdy składnik po l jest dodany jawnie do ps.Operator.
"""
function XXZ_ladder_finite(J::Float64, L::Int, Δ::Float64, U_1::Float64, U_2::Float64)
  H = ps.Operator(L * 2)

  for l in 1:L
    lp = mod1(l + 1, L)

    for y in 1:2
      H += (J / 2) * ps.string_2d(("S+", l, y, "S-", lp, y), L, 2)
      H += (J / 2) * ps.string_2d(("S-", l, y, "S+", lp, y), L, 2)
      H += (J * Δ) * ps.string_2d(("Sz", l, y, "Sz", lp, y), L, 2)
    end

    H += (U_1) * ps.string_2d(("Sz", l, 1, "Sz", l, 2), L, 2)

    H += (U_2 / 2) * ps.string_2d(("S+", l, 1, "S-", l, 2), L, 2)
    H += (U_2 / 2) * ps.string_2d(("S+", l, 2, "S-", l, 1), L, 2)
  end

  return H / L
end

"""
Buduje jeden przesunięty string.

Przykład:
upper = [1,3,0]
lower = [2,2,2]

shift = 0:
  x=1 -> l=1
  x=2 -> l=2
  x=3 -> l=3

shift = 1:
  x=1 -> l=2
  x=2 -> l=3
  x=3 -> l=4

Zawijanie periodyczne przez mod1.
"""
function string_2d_from_legs_shifted(upper_list::Vector{Int}, lower_list::Vector{Int}, L::Int, shift::Int; coeff = 1.0 + 0.0im)
  args = Any[]
  M = length(upper_list)

  M > L && error("Observable support M=$M is larger than finite system size L=$L.")

  for x in 1:M
    l = mod1(x + shift, L)

    op = upper_list[x]
    if op != 0
      push!(args, symbol_map[op], l, 1)
    end
  end

  for x in 1:M
    l = mod1(x + shift, L)

    op = lower_list[x]
    if op != 0
      push!(args, symbol_map[op], l, 2)
    end
  end

  length(args) == 0 && error("Tried to build identity-only string from observable label.")

  return coeff * ps.string_2d(Tuple(args), L, 2)
end

"""
Skończona translacyjna suma:

O = sum_{shift=0}^{L-1} O_shift

To zastępuje OperatorTS, ale jawnie dodaje każdy przesunięty string.
"""
function translated_string_from_legs(upper_list::Vector{Int}, lower_list::Vector{Int}, L::Int; coeff = 1.0 + 0.0im)
  O = ps.Operator(L * 2)

  for shift in 0:L-1
    O += string_2d_from_legs_shifted(upper_list, lower_list, L, shift; coeff = coeff)
  end

  return O
end

"""
Skończony odpowiednik build_ladder_operator z kodu LIOM.

Dla np. RI oraz 130|222:
- górę projektujemy jako R/even,
- dół jako I/odd,
- mnożymy projekcje między nogami,
- jawnie sumujemy po przesunięciach na skończonym L.
"""
function basis_operator_finite(upper_list::Vector{Int}, upper_T::Symbol, lower_list::Vector{Int}, lower_T::Symbol, L::Int)
  O = ps.Operator(L * 2)

  upper_terms = leg_terms(upper_list, upper_T)
  lower_terms = leg_terms(lower_list, lower_T)

  for (u_list, cu) in upper_terms
    for (l_list, cl) in lower_terms
      O += translated_string_from_legs(u_list, l_list, L; coeff = cu * cl)
    end
  end

  return O
end

function basis_operator_finite(sector::String, label::String, L::Int)
  upper_list, lower_list = parse_label(label)
  upper_T, lower_T = sector_to_symbols(sector)
  return basis_operator_finite(upper_list, upper_T, lower_list, lower_T, L)
end

"""
Obserwabla jako kombinacja kilku elementów bazy.

Format:
observable_terms = [
  (wspolczynnik, "RI", "130|222"),
  (wspolczynnik, "IR", "132|000"),
]
"""
function build_observable_finite(observable_terms, L::Int)
  A = ps.Operator(L * 2)

  for (coeff, sector, label) in observable_terms
    A += coeff * basis_operator_finite(sector, label, L)
  end

  return A
end

function dense_matrix_from_operator(O)
  return Matrix(O)
end

function commutator_norm_operator(H, A)
  C = ps.com(H, A)
  numerator = ps.trace_product(ps.dagger(C), C)
  denominator = ps.trace_product(ps.dagger(A), A)

  abs(denominator) < 1e-14 && error("Observable has zero norm.")

  return real(numerator / denominator)
end

function autocorrelation_from_ed(H, A, times::Vector{Float64})
  println("Converting H to matrix...")
  Hmat = dense_matrix_from_operator(H)

  println("Converting A to matrix...")
  Amat = dense_matrix_from_operator(A)

  println("Diagonalizing Hamiltonian...")
  F = eigen(Hermitian(Hmat))

  E = F.values
  V = F.vectors

  println("Transforming observable to energy basis...")
  A_E = V' * Amat * V

  weights = abs2.(A_E)
  denominator = sum(weights)

  abs(denominator) < 1e-14 && error("Observable has zero matrix norm.")

  C = zeros(Float64, length(times))

  println("Computing autocorrelation...")
  for (k, t) in enumerate(times)
    phase = exp.(im .* ((E .- E') .* t))
    value = sum(phase .* weights) / denominator
    C[k] = real(value)
  end

  return C
end

function parse_args()
  s = ArgParseSettings(description = "Finite ED autocorrelation checker for LIOM observables")
  @add_arg_table s begin
    "--delta", "-d"
    help = "Anisotropy parameter Δ"
    arg_type = Float64
    default = 0.3
    "--U_1"
    help = "Rung Sz-Sz coupling U_1"
    arg_type = Float64
    default = 0.75
    "--U_2"
    help = "Rung XY coupling U_2"
    arg_type = Float64
    default = 0.75
    "--L", "-L"
    help = "Finite ladder length"
    arg_type = Int
    default = 4
    "--tmax"
    help = "Maximal time"
    arg_type = Float64
    default = 100.0
    "--nt"
    help = "Number of time points"
    arg_type = Int
    default = 501
    "--output-dir"
    help = "Output directory"
    arg_type = String
    default = "ed_autocorr_dane"
  end
  return ArgParse.parse_args(s)
end

function main()
  args = parse_args()

  L = args["L"]

  J = 1.0
  Delta = args["delta"]
  U_1 = args["U_1"]
  U_2 = args["U_2"]

  tmax = args["tmax"]
  nt = args["nt"]
  output_dir = args["output-dir"]

  # ==========================================================
  # TUTAJ  OBSERWABLĘ
  # Format jednego składnika:
  # (wspolczynnik, "sektor", "gorna|dolna")
  # Przykłady:
  # (1.0, "RI", "130|222")
  # (-0.43, "IR", "132|000")
  # (0.12, "RI", "222|130")
  # Czyli obserwabla:
  # A = suma_i coeff_i * O_i
  # Na start możesz wpisać jeden operator z bazy,
  # a potem ręcznie przepisać kilka największych współczynników
  # z wektora własnego LIOM-a.
  # ==========================================================

  observable_terms = [
    (1.0, "RI", "130|222"),
  ]

  # ==========================================================

  println("Parameters:")
  println("L = ", L)
  println("N spins = ", 2L)
  println("Hilbert space dimension = ", 2^(2L))
  println("J = ", J)
  println("Delta = ", Delta)
  println("U_1 = ", U_1)
  println("U_2 = ", U_2)
  println("tmax = ", tmax)
  println("nt = ", nt)

  println()
  println("Observable terms:")
  for item in observable_terms
    println(item)
  end

  println()
  println("Building finite Hamiltonian as ps.Operator...")
  H = XXZ_ladder_finite(J, L, Delta, U_1, U_2)

  println("Building observable A as finite ps.Operator...")
  A = build_observable_finite(observable_terms, L)

  # Jeśli obserwabla ma być hermitowska, zostawiamy tę linię.
  # Dla sektorów R/I zwykle powinna wyjść hermitowska,
  # ale ta symetryzacja usuwa drobne problemy numeryczne.
  A = (A + ps.dagger(A)) / 2

  println()
  println("Commutator norm Tr([H,A]†[H,A]) / Tr(A†A):")
  lambda_A = commutator_norm_operator(H, A)
  println(lambda_A)

  times = collect(range(0.0, tmax, length = nt))

  C = autocorrelation_from_ed(H, A, times)

  out_path = joinpath(@__DIR__, output_dir)
  mkpath(out_path)

  file_tag = "ED_autocorr_L_$(L)_J_$(J)_d_$(Delta)_U_1_$(U_1)_U_2_$(U_2)"
  output_file = joinpath(out_path, "$(file_tag).txt")

  open(output_file, "w") do io
    println(io, "# t\tC_A(t)")
    println(io, "# L = ", L)
    println(io, "# J = ", J)
    println(io, "# Delta = ", Delta)
    println(io, "# U_1 = ", U_1)
    println(io, "# U_2 = ", U_2)
    println(io, "# commutator_norm = ", lambda_A)
    println(io, "# observable_terms:")
    for item in observable_terms
      println(io, "# ", item)
    end
    println(io, "#")
    for i in eachindex(times)
      println(io, times[i], "\t", C[i])
    end
  end
  println()
  println("Saved autocorrelation to:")
  println(output_file)
  println()
  println("First values:")
  for i in 1:min(10, length(times))
    @printf("%10.5f  % .12f\n", times[i], C[i])
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end