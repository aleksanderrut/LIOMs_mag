import PauliStrings as ps
using LinearAlgebra
using Printf
using ProgressMeter
using Base.Threads
using DelimitedFiles
using ArgParse

"""
One-leg XXZ ladder Hamiltonian with periodic boundary conditions along the leg and rung phonon coupling.
Hamiltonian written in the S+, S-, Sz basis.
"""
function XXZ_ladder(J::Float64, L::Int, Δ::Float64, ω_0::Float64, g::Float64, J_prime::Float64)
  H = ps.Operator(L * 2)

  for l in 1:L
    lp = mod1(l + 1, L)

    H += (J / 2) * ps.string_2d(("S+", l, 1, "S-", lp, 1), L, 2)
    H += (J / 2) * ps.string_2d(("S-", l, 1, "S+", lp, 1), L, 2)
    H += (J * Δ) * ps.string_2d(("Sz", l, 1, "Sz", lp, 1), L, 2)

    H += ω_0 * ps.string_2d(("Sz", l, 2), L, 2)
    H += g * ps.string_2d(("Sz", l, 1, "S+", l, 2), L, 2)
    H += g * ps.string_2d(("Sz", l, 1, "S-", l, 2), L, 2)

    H += (J_prime / 2) * ps.string_2d(("S+", l, 2, "S-", lp, 2), L, 2)
    H += (J_prime / 2) * ps.string_2d(("S-", l, 2, "S+", lp, 2), L, 2)
    #H += (J_prime * Δ) * ps.string_2d(("Sz", l, 2, "Sz", lp, 2), L, 2)
  end
  
  return ps.OperatorTS{(L, 2), (true, false)}(H) / L
end

"""
hs_product(op1, op2)
Hilbert-Schmidt inner product normalized by 2^N.
"""
@inline hs_product(op1::ps.OperatorTS, op2::ps.OperatorTS) = ps.trace_product(ps.dagger(op1), op2) / UInt128(2)^ps.qubitlength(op1)

"""
hs_norm(op) -> Float64
Hilbert-Schmidt norm.
"""
@inline hs_norm(op::ps.OperatorTS) = sqrt(real(hs_product(op, op)))

"""
0 = identity
1 = S+
2 = Sz
3 = S-
"""
@inline hc(o::Int) = (o == 1 ? 3 : (o == 3 ? 1 : o))

"""
undigit(olist::Vector{Int}) -> Int
Convert a base-4 digit list to its integer representation.
"""
@inline undigit(olist) = sum([olist[i] * 4^(i - 1) for i in eachindex(olist)])

symbol_map = Dict(0 => "1", 1 => "S+", 2 => "Sz", 3 => "S-")
superscript_map = Dict("S+" => "P", "Sz" => "Z", "S-" => "M", "1" => "1")
@inline hc(o::String) = (o == "P" ? "M" : (o == "M" ? "P" : o))

function digits_string(op_list::Vector{Int})
  return join(string.(op_list), "")
end

"""
is_odd_under_parity(op_list, time_reversal) -> Bool
Determine if an operator list is odd under parity assuming given time-reversal symmetry.
"""
function is_odd_under_parity(op_list::Vector{Int}, time_reversal::Symbol)
  cnt_z = count(x -> x == 2, op_list)
  sign = (-1)^cnt_z
  if time_reversal == :odd
    return sign == 1
  end
  return sign == -1
end

"""
is_even_under_parity(op_list, time_reversal) -> Bool
Determine if an operator list is even under parity assuming given time-reversal symmetry.
"""
function is_even_under_parity(op_list::Vector{Int}, time_reversal::Symbol)
  return !is_odd_under_parity(op_list, time_reversal)
end

function parity_sector(op_list::Vector{Int}, time_reversal::Symbol)
  if is_odd_under_parity(op_list, time_reversal)
    return :odd
  else
    return :even
  end
end

"""
Compose sectors from upper and lower leg.
even * even = even; odd * odd = even; even * odd = odd; odd * even = odd
"""
function compose_sector(a::Symbol, b::Symbol)
  return a == b ? :even : :odd
end

"""
Saved basis convention.
odd(imag) -> I, even(real) -> R
"""
function sector_char(s::Symbol)
  s == :odd && return "I"
  s == :even && return "R"
  error("Unknown sector: $s")
end

"""
Build one local ladder string using ps.string_2d.
The support starts at l = 1. Translation along legs is handled by OperatorTS.
"""
function string_2d_from_legs(upper_list::Vector{Int}, lower_list::Vector{Int}, L::Int; coeff = 1 + 0im)
  args = Any[]
  M = length(upper_list)

  for x in 1:M
    op = upper_list[x]
    if op != 0
      push!(args, symbol_map[op], x, 1)
    end
  end

  for x in 1:M
    op = lower_list[x]
    if op != 0
      push!(args, symbol_map[op], x, 2)
    end
  end

  return coeff * ps.string_2d(Tuple(args), L, 2)
end

"""
Projected one-leg object. T even: O + O† / T odd: i(O - O†)
"""
function leg_terms(op_list::Vector{Int}, time_reversal::Symbol)
  op_list_hc = hc.(op_list)
  if time_reversal == :even
    return [(op_list, 1 + 0im), (op_list_hc, 1 + 0im)]
  elseif time_reversal == :odd
    return [(op_list, 0 + 1im), (op_list_hc, 0 - 1im)]
  end
end

"""
Build full ladder operator from projected upper and lower legs.
"""
function build_ladder_operator(upper_list::Vector{Int}, upper_T::Symbol, lower_list::Vector{Int}, lower_T::Symbol, L::Int)
  op = ps.Operator(L * 2)
  upper_terms = leg_terms(upper_list, upper_T)
  lower_terms = leg_terms(lower_list, lower_T)
  #takie trochę mnożenie kazdego elmentu z każdym
  """
  przykład       góra:  +i 130    oraz   -i 310;        dół:   +i 132    oraz   -i 312
  (+i)(+i) 130|132  = - 130|132
  (+i)(-i) 130|312  = + 130|312
  (-i)(+i) 310|132  = + 310|132
  (-i)(-i) 310|312  = - 310|312
  wynik: - 130|132 + 130|312 + 310|132 - 310|312 
  """ 
  for (u_list, cu) in upper_terms
    for (l_list, cl) in lower_terms
      op += string_2d_from_legs(u_list, l_list, L; coeff = cu * cl)
    end
  end
  return ps.OperatorTS{(L, 2), (true, false)}(op)
end

"""
Generate all one-leg blocks and remember their sectors.
Each element is: (op_list, T, P, Sz_check)
"""
function generate_leg_blocks(M::Int, conserve_Sz::Symbol)
  legs = []
  for i in 0:4^M-1
    op_list = digits(i, base = 4, pad = M)
    op_list_hc = hc.(op_list)
    op_hc = undigit(op_list_hc)
    op_hc > i && continue

    Sz_check = sum(op_list .== 1) - sum(op_list .== 3) # już rostrzygające czy Sz jest zachowane czy nie
    (conserve_Sz == :yes) && Sz_check != 0 && continue
    (conserve_Sz == :no) && Sz_check == 0 && continue

    if sum(in.(op_list, Ref([1, 3]))) != 0 # czy mozna zbudowac T odd (imag) jesli idziemy do drugiego warunku 
      T = :odd
      P = parity_sector(op_list, T) # już rostrzygające czy P - ecen czy P - odd
      push!(legs, (op_list, T, P, Sz_check))
    end

    T = :even
    P = parity_sector(op_list, T) # już rostrzygające czy P - ecen czy P - odd
    push!(legs, (op_list, T, P, Sz_check))
  end
  return legs
end

"""
Generate ladder basis. M means support on one leg.
"""
function all_M_ladder_leg_sym(L::Int, M::Int, time_reversal::Symbol, parity::Symbol, conserve_Sz_fermion::Symbol, conserve_Sz_boson::Symbol, include_fermion_identity::Bool)
  ops = ps.OperatorTS[]
  ops_list = []
  ops_list_rows = String[]

  upper_legs = generate_leg_blocks(M, conserve_Sz_fermion)
  lower_legs = generate_leg_blocks(M, conserve_Sz_boson)

  #pęta po wszystkich kombinacjach opertorów na górenj i dolnej nodze
  for upper in upper_legs
    upper_list, upper_T, upper_P, upper_Sz = upper
    for lower in lower_legs
      lower_list, lower_T, lower_P, lower_Sz = lower

      #wyrzucenie 000|000
      all(upper_list .== 0) && all(lower_list .== 0) && continue

      #opcjonalne wyrzucenie elementów działających wyłącznie na nodze bozonowej
      !include_fermion_identity && all(upper_list .== 0) && continue

      upper_list[1] == 0 && lower_list[1] == 0 && continue

      #policzenie współych smetrii opertora już na dwóch nogach
      total_T = compose_sector(upper_T, lower_T)
      total_P = compose_sector(upper_P, lower_P)

      (time_reversal == :even) && total_T != :even && continue
      (time_reversal == :odd) && total_T != :odd && continue

      (parity == :even) && total_P != :even && continue
      (parity == :odd) && total_P != :odd && continue

      op = build_ladder_operator(upper_list, upper_T, lower_list, lower_T, L)
      #println("Operator basis element:")
      #println(sector_char(upper_T) * sector_char(lower_T), "\t", digits_string(upper_list), "|", digits_string(lower_list))
      #println(op)
      #println()
      nrm = hs_norm(op)
      abs(nrm) < 1e-12 && continue

      push!(ops, op)
      push!(ops_list, (upper_list, lower_list, upper_T, lower_T))
      push!(ops_list_rows, "\t" * sector_char(upper_T) * sector_char(lower_T) * "\t" * digits_string(upper_list) * "|" * digits_string(lower_list))
    end
  end

  return ops, ops_list, ops_list_rows
end

"""
compute_lioms(H, L, max_supp, time_reversal, parity, conserve_Sz_fermion, conserve_Sz_boson, include_fermion_identity)
"""
function compute_lioms(H::ps.OperatorTS, L::Int, max_supp::Int, time_reversal::Symbol, parity::Symbol, conserve_Sz_fermion::Symbol, conserve_Sz_boson::Symbol, include_fermion_identity::Bool)
  println("Generating operators...")
  ops, ops_list, ops_list_rows = all_M_ladder_leg_sym(L, max_supp, time_reversal, parity, conserve_Sz_fermion, conserve_Sz_boson, include_fermion_identity)
  println("Generated $(length(ops)) operators")
  println("Computing norms...")
  if length(ops) == 0
    error("Empty operator basis after symmetry filtering.")
  end

  norms = hs_norm.(ops)
  good = findall(x -> abs(x) > 1e-12, norms)
  ops = ops[good]
  ops_list = ops_list[good]
  ops_list_rows = ops_list_rows[good]
  norms = norms[good]
  if length(ops) == 0
    error("Empty operator basis after norm filtering.")
  end
  ops ./= norms

  println("Precomputing commutators...")
  comms = Vector{ps.OperatorTS}(undef, length(ops))
  @showprogress for i in eachindex(ops)
    comms[i] = im * ps.com(H, ops[i])
  end

  corr_mat = zeros(Float64, length(ops), length(ops))
  total = (length(ops) * (length(ops) + 1)) ÷ 2
  p = Progress(total; desc = "Computing correlation matrix...", showspeed = true)

  @threads :greedy for i in eachindex(ops)
    comm_i = comms[i]
    for j in i:length(ops)
      comm_j = comms[j]
      corr_mat[i, j] = real(hs_product(comm_i, comm_j))
      i != j && (corr_mat[j, i] = corr_mat[i, j])
      next!(p)
    end
  end
  finish!(p)

  F = eigen(Symmetric(corr_mat))
  return F.values, F.vectors, ops_list, ops_list_rows
end

"""
Save basis in requested short format.
"""
function save_operator_labels(filename, ops_list)
  open(filename, "w") do io
    for item in ops_list
      upper_list, lower_list, upper_T, lower_T = item
      println(io, "1", "\t", sector_char(upper_T) * sector_char(lower_T), "\t", digits_string(upper_list), "|", digits_string(lower_list))
    end
  end
end

"""
Save all found LIOMs to one file, Output format:
# LIOM k, eigenvalue = ...
basis_element    coefficient
"""
function site_suffix(pos::Int)
  pos == 1 && return "i"
  return "i+$(pos - 1)"
end

function spin_string_from_list(op_list::Vector{Int})
  parts = String[]

  for (pos, op) in enumerate(op_list)
    op_name = symbol_map[op]
    site = site_suffix(pos)

    if op_name == "1"
      push!(parts, "1_{$site}")
    elseif op_name == "Sz"
      push!(parts, "S^z_{$site}")
    elseif op_name == "S+"
      push!(parts, "S^+_{$site}")
    elseif op_name == "S-"
      push!(parts, "S^-_{$site}")
    end
  end

  return join(parts, " ")
end

"""
Save exact LIOMs and the first few nonzero modes to one file.

Exact LIOMs are eigenvectors with eigenvalue abs(lambda) < eval_tol.
Nonzero modes are the first n_nonzero eigenvectors with abs(lambda) >= eval_tol.

Output format:
# LIOM k
# eigenvalue = ...
translated_basis_element    basis_element    coefficient
"""
function save_lioms(filename, evals, evecs, ops_list, ops_list_rows; eval_tol = 1e-10, coeff_tol = 1e-8, n_nonzero = 3)
  liom_indices = findall(x -> abs(x) < eval_tol, evals)
  nonzero_indices = findall(x -> abs(x) >= eval_tol, evals)
  nonzero_indices = nonzero_indices[1:min(n_nonzero, length(nonzero_indices))]

  open(filename, "w") do io
    println(io, "# Number of exact LIOMs found = ", length(liom_indices))
    println(io, "# Number of saved nonzero modes = ", length(nonzero_indices))
    println(io, "# eval_tol = ", eval_tol)
    println(io, "# coeff_tol = ", coeff_tol)
    println(io)

    println(io, "# =========================")
    println(io, "# EXACT LIOMs")
    println(io, "# =========================")
    println(io)

    for (liom_number, j) in enumerate(liom_indices)
      squared_coeffs = abs2.(evecs[:, j])
      sorted_indices = sortperm(squared_coeffs, rev = true)

      println(io, "# LIOM ", liom_number)
      println(io, "# eigenvalue = ", evals[j])
      println(io, "# sum of squared coefficients = ", sum(squared_coeffs))
      println(io, "# translated_basis_element\tbasis_element\tsquared_coefficient")

      for i in sorted_indices
        squared_coeff = squared_coeffs[i]
        if squared_coeff > coeff_tol^2
          basis_label = strip(ops_list_rows[i])
          basis_label = replace(basis_label, "\t" => " ")
          upper_list, lower_list, upper_T, lower_T = ops_list[i]
          translated_upper = spin_string_from_list(upper_list)
          translated_lower = spin_string_from_list(lower_list)
          translated_label = translated_upper * " | " * translated_lower
          @printf(
            io,
            "%s\t%s\t%.16e\n",
            translated_label,
            basis_label,
            squared_coeff
          )
        end
      end
      println(io)
    end
    println(io)
    println(io, "# =========================")
    println(io, "# FIRST NONZERO MODES")
    println(io, "# =========================")
    println(io)
    for (mode_number, j) in enumerate(nonzero_indices)
      squared_coeffs = abs2.(evecs[:, j])
      sorted_indices = sortperm(squared_coeffs, rev = true)
      println(io, "# NONZERO MODE ", mode_number)
      println(io, "# eigenvalue = ", evals[j])
      println(io, "# sum of squared coefficients = ", sum(squared_coeffs))
      println(io, "# translated_basis_element\tbasis_element\tsquared_coefficient")
      for i in sorted_indices
        squared_coeff = squared_coeffs[i]
        if squared_coeff > coeff_tol^2
          basis_label = strip(ops_list_rows[i])
          basis_label = replace(basis_label, "\t" => " ")
          upper_list, lower_list, upper_T, lower_T = ops_list[i]
          translated_upper = spin_string_from_list(upper_list)
          translated_lower = spin_string_from_list(lower_list)
          translated_label = translated_upper * " | " * translated_lower
          @printf(
            io,
            "%s\t%s\t%.16e\n",
            translated_label,
            basis_label,
            squared_coeff
          )
        end
      end
      println(io)
    end
  end
end

function parse_args()
  s = ArgParseSettings(description = "Compute LIOMs for two-leg XXZ ladder with leg-sector basis")
  @add_arg_table s begin
    "--delta", "-d"
    help = "Anisotropy parameter Δ"
    arg_type = Float64
    default = 0.3
    "--omega", "-w"
    help = "Omega_0 parameter"
    arg_type = Float64
    default = 1.0
    "--coupling", "-g"
    help = "interaction strength between phonon and fermions"
    arg_type = Float64
    default = 1.0
    "--J-prime"
    help = "Coupling strength along the second leg"
    arg_type = Float64
    default = 0.0
    "--max-supp", "-M"
    help = "Support M on one leg"
    arg_type = Int
    default = 2
    "--time-reversal", "-T"
    help = "Total ladder time-reversal sector: even, odd, both"
    arg_type = String
    default = "odd"
    "--parity", "-P"
    help = "Total ladder parity sector: even, odd, both"
    arg_type = String
    default = "even"
    "--conserve-Sz-fermion", "-F"
    help = "Conserve Sz on the fermionic leg? yes, no, both"
    arg_type = String
    default = "yes"
    "--conserve-Sz-boson", "-B"
    help = "Conserve Sz on the bosonic leg? yes, no, both"
    arg_type = String
    default = "both"
    "--include-fermion-identity"
    help = "Include basis elements containing only identity operators on the fermionic leg? yes or no"
    arg_type = String
    default = "yes"
  end
  return ArgParse.parse_args(s)
end

function main()
  args = parse_args()

  allowed_time_rev = ["even", "odd", "both"]
  !(args["time-reversal"] in allowed_time_rev) && error("Invalid --type: $(args["time-reversal"]). Must be one of $(allowed_time_rev)")
  allowed_parity = ["even", "odd", "both"]
  !(args["parity"] in allowed_parity) && error("Invalid --parity: $(args["parity"]). Must be one of $(allowed_parity)")
  allowed_conserve_Sz = ["yes", "no", "both"]
  !(args["conserve-Sz-fermion"] in allowed_conserve_Sz) &&
  error("Invalid --conserve-Sz-fermion: $(args["conserve-Sz-fermion"]). Must be one of $(allowed_conserve_Sz)")
  !(args["conserve-Sz-boson"] in allowed_conserve_Sz) &&
  error("Invalid --conserve-Sz-boson: $(args["conserve-Sz-boson"]). Must be one of $(allowed_conserve_Sz)")
  allowed_identity_options = ["yes", "no"]
  !(args["include-fermion-identity"] in allowed_identity_options) &&
  error("Invalid --include-fermion-identity: $(args["include-fermion-identity"]). Must be one of $(allowed_identity_options)")
  time_reversal = Symbol(args["time-reversal"])
  conserve_Sz_fermion = Symbol(args["conserve-Sz-fermion"])
  conserve_Sz_boson = Symbol(args["conserve-Sz-boson"])
  include_fermion_identity = args["include-fermion-identity"] == "yes"
  parity = Symbol(args["parity"])

  # Parameters
  J = 1.0
  J_prime = args["J-prime"]
  Delta = args["delta"]
  omega_0 = args["omega"]
  max_supp = args["max-supp"]
  g = args["coupling"]
  L = 2 * max_supp + 2

  println("Parameters:")
  println("J = ", J)
  println("J_prime = ", J_prime)
  println("Δ = ", Delta)
  println("ω_0 = ", omega_0)
  println("g = ", g)
  println("max_supp = ", max_supp)
  println("L = ", L)
  println("Behavior under time-reversal symmetry = ", time_reversal)
  println("Behavior under parity = ", parity)
  println("Conserve Sz on fermionic leg? = ", conserve_Sz_fermion)
  println("Conserve Sz on bosonic leg? = ", conserve_Sz_boson)
  println("Include fermionic identity-only leg? = ", include_fermion_identity)

  H = XXZ_ladder(J, L, Delta, omega_0, g, J_prime)
  #println("Hamiltonian:")
  #println(H)

  tic = time_ns()
  evals, evecs, ops_list, ops_list_rows = compute_lioms(H, L, max_supp, time_reversal, parity, conserve_Sz_fermion, conserve_Sz_boson, include_fermion_identity)
  toc = time_ns()

base_path = "$(@__DIR__)/liom_mag_dane"
evals_path = "$(base_path)/wartosci_wlasne"
evecs_path = "$(base_path)/wektory_wlasne"
basis_path = "$(base_path)/bazy"
logs_path = "$(base_path)/logi"
lioms_path = "$(base_path)/LIOMs_found"
mkpath(evals_path)
mkpath(evecs_path)
mkpath(basis_path)
mkpath(logs_path)
mkpath(lioms_path)
fid_tag = include_fermion_identity ? "yes" : "no"
file_tag = "M$(max_supp)_Jp$(J_prime)_d$(Delta)_w$(omega_0)_g$(g)_T$(time_reversal)_P$(parity)_F$(conserve_Sz_fermion)_B$(conserve_Sz_boson)_FId$(fid_tag)"
filename_evals = "$(evals_path)/eigenvalues_$(file_tag).txt"
filename_evecs = "$(evecs_path)/eigenvectors_$(file_tag).txt"
filename_operators = "$(basis_path)/operators_$(file_tag).txt"
filename_log = "$(logs_path)/log_$(file_tag).txt"
filename_lioms = "$(lioms_path)/LIOMs_$(file_tag).txt"
writedlm(filename_evals, evals)
writedlm(filename_evecs, evecs)
save_operator_labels(filename_operators, ops_list)
save_lioms(filename_lioms, evals, evecs, ops_list, ops_list_rows)
  open(filename_log, "w") do io
    println(io, "J = ", J)
    println(io, "J_prime = ", J_prime)
    println(io, "Δ = ", Delta)
    println(io, "ω_0 = ", omega_0)
    println(io, "g = ", g)
    println(io, "max_supp = ", max_supp)
    println(io, "L = ", L)
    println(io, "Behavior under time-reversal symmetry = ", time_reversal)
    println(io, "Behavior under parity = ", parity)
    println(io, "Conserve Sz on fermionic leg? = ", conserve_Sz_fermion)
    println(io, "Conserve Sz on bosonic leg? = ", conserve_Sz_boson)
    println(io, "Include fermionic identity-only leg? = ", include_fermion_identity)
    println(io, "Operators basis size = ", length(ops_list))
    println(io, "Number of LIOMs found = ", count(x -> abs(x) < 1e-10, evals))
    println(io, "Elapsed time = ", (toc - tic) / 1e9, " s")
  end

  max_evals = min(10, length(evals))
  println("$(max_evals) smallest eigenvalues:")
  println(evals[1:max_evals])
  #List of operator basis elements with their coefficients in the LIOMs corresponding to the smallest eigenvalues.
  
  # println("$(max_evals) eigenvectors corrsponding to smallest eigenvalues:")
  # max_label_width = maximum(length.(ops_list_rows))

  # for i in eachindex(ops_list_rows)
  #   op_str = ops_list_rows[i]
  #   print("∑_l ", rpad(op_str, max_label_width), "   ")
  #   for j in 1:max_evals
  #     num_str = @sprintf("% .5f", evecs[i, j])
  #     num_str = replace(num_str, " -" => "-")
  #     print(" ", num_str)
  #   end
  #   println()
  # end
  
  println("Size of basis for M=$max_supp, " * "time reversal = $time_reversal, " * "parity = $parity, " * "fermionic Sz = $conserve_Sz_fermion, " * "bosonic Sz = $conserve_Sz_boson, " * "fermionic identity included = $include_fermion_identity: ", length(ops_list))
  println("Number of LIOMs found: ", count(x -> abs(x) < 1e-10, evals))
  println("Elapsed time: ", (toc - tic) / 1e9, " s")
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
