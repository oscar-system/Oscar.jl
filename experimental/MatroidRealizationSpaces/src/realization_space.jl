
@attributes mutable struct MatroidRealizationSpace{BaseRingType, RingType} <: AbsAffineScheme{BaseRingType, RingType}
  defining_ideal::Union{Ideal,NumFieldOrderIdeal}
  inequations::Vector{RingElem}
  ambient_ring::Ring
  realization_matrix::Union{MatElem,Nothing}
  char::Union{Int,Nothing}
  q::Union{Int,Nothing}
  ground_ring::Ring
  one_realization::Bool

  # Fields for caching
  underlying_scheme::AbsAffineScheme{BaseRingType, RingType}

  function MatroidRealizationSpace(
    I::Union{Ideal,NumFieldOrderIdeal},
    ineqs::Vector{<:RingElem},
    R::Ring,
    mat::Union{MatElem,Nothing},
    char::Union{Int,Nothing},
    q::Union{Int,Nothing},
    ground_ring::Ring)
    BaseRingType = typeof(ground_ring)
    RingType = typeof(R)
    if R isa MPolyRing
      PolyRingType = typeof(R)
      MultSetType = MPolyPowersOfElement{BaseRingType, elem_type(BaseRingType), PolyRingType, elem_type(PolyRingType)}
      RingType = MPolyQuoLocRing{BaseRingType, elem_type(BaseRingType), PolyRingType, elem_type(PolyRingType), MultSetType}
    end
    return new{BaseRingType, RingType}(I, ineqs, R, mat, char, q, ground_ring, false)
  end
end

function Base.show(io::IO, ::MIME"text/plain", RS::MatroidRealizationSpace)
  if has_attribute(RS, :is_realizable) && !is_realizable(RS)
    if RS.char === nothing && RS.q === nothing
      print(io, "The matroid is not realizable.")
    else
      print(io, "The matroid is not realizable over the specified field or characteristic.")
    end
  else
    io = pretty(io)
    if RS.one_realization
      println(io, "One realization is given by")
    elseif has_attribute(RS, :is_realizable) && is_realizable(RS)
      println(io, "The realizations are parametrized by")
    elseif !has_attribute(RS, :is_realizable)
      println(io, "The realization space is")
    end
    print(io, Indent())
    show(io, MIME("text/plain"), RS.realization_matrix)
    print(io, "\n", Dedent(), "in the ", Lowercase(), RS.ambient_ring)
    I = RS.defining_ideal
    if !iszero(I)
      print(io, "\nwithin the vanishing set of the ideal\n", I)
    end
    if length(RS.inequations) > 0
      print(io, "\navoiding the zero loci of the polynomials\n", RS.inequations)
    end
  end
end

@doc raw"""
    is_realizable(M; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)

* If char = nothing, then this function determines whether the matroid is realizable over some field.

* If `char == 0`, then this function determines whether the matroid is realizable over some field of
    characteristic 0.

* If char = p is prime, this function determines whether the matroid is realizable
    over the finite field ``GF(p)``.

* If `char == p` and `q` is a power of `p`, this function determines whether the matroid is realizable over the
    finite field ``GF(q)``.
"""
function is_realizable(
  M::Matroid; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing
)::Bool
  RS = realization_space(M; char=char, q=q)
  return is_realizable(RS)
end

@attr Bool function is_realizable(RS::MatroidRealizationSpace)
  if !(RS.ambient_ring isa MPolyRing)
    return true
  end
  for p in minimal_primes(RS.defining_ideal)
    component_non_trivial = all(!in(p), RS.inequations)
    if component_non_trivial
      return true
    end
  end
  return false
end

@doc raw"""
    defining_ideal(RS::MatroidRealizationSpace)

The ideal of the matroid realization space `RS`.
"""
defining_ideal(RS::MatroidRealizationSpace) = RS.defining_ideal

@doc raw"""
    inequations(RS::MatroidRealizationSpace)

Generators of the localizing semigroup of `RS`.
These are the polynomials that need to be nonzero in any realization.
"""
inequations(RS::MatroidRealizationSpace) = RS.inequations

@doc raw"""
    ambient_ring(RS::MatroidRealizationSpace)

The polynomial ring containing the ideal `defining_ideal(RS)` and the polynomials in `inequations(RS)`.
"""
ambient_ring(RS::MatroidRealizationSpace) = RS.ambient_ring

### The following method uses types which are declared only later.
# Hence hit has been moved to
#
#   src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Methods.jl.
#=
function underlying_scheme(RS::MatroidRealizationSpace{<:Field, <:MPolyQuoLocRing})
  isdefined(RS, :underlying_scheme) && return RS.underlying_scheme

  P = ambient_ring(RS)::MPolyRing
  I = defining_ideal(RS)::MPolyIdeal
  U = powers_of_element(inequations)::MPolyPowersOfElement
  RS.underlying_scheme = spec(R, I, U)
  return RS.underlying_scheme
end
=#

function underlying_scheme(RS::MatroidRealizationSpace)
  error("method not implemented for realization spaces of type $(typeof(RS))")
end
@doc raw"""
    realization_matrix(RS::MatroidRealizationSpace)

A matrix with entries in ambient_ring(RS) whose columns, when filled in with values satisfying equalities
from `defining_ideal(RS)` and inequations from `inequations(RS)`, form a realization for the matroid.
"""
realization_matrix(RS::MatroidRealizationSpace) = RS.realization_matrix

function realization_space_matrix(M::Matroid, B::Vector{Int}, F::Ring)
  # prepare the combinatorial data

  circs = fundamental_circuits_of_basis(M, B)

  nonIdCols = setdiff(matroid_groundset(M), B)
  circs = [setdiff(c, nonIdCols) for c in circs]

  rk = rank(M)
  n = length(M)

  # we start by computing the number of variables:
  numVars = 0
  unUsedRowsForOnes = collect(2:rk)
  for col in 1:(n - rk), row in 1:rk
    circ = circs[col]
    if !(B[row] == minimum(circ)) && B[row] in circ
      if row in unUsedRowsForOnes
        unUsedRowsForOnes = setdiff(unUsedRowsForOnes, [row])
      else
        numVars += 1
      end
    end
  end

  if numVars > 0
    R, x = polynomial_ring(F, numVars)
  else
    R = F
    x = Vector{MPolyRingElem}()
  end

  unUsedRowsForOnes = collect(2:rk)

  # create the matrix and fill it with entries

  mat = zero_matrix(R, rk, n)

  for i in 1:rk
    mat[i, B[i]] = R(1)
  end

  varCounter = 1

  for col in 1:(n - rk), row in 1:rk
    circ = circs[col]
    c = nonIdCols[col]

    if B[row] == minimum(circ)
      mat[row, c] = R(1)
    elseif B[row] in circ
      if row in unUsedRowsForOnes
        mat[row, c] = R(1)
        unUsedRowsForOnes = setdiff(unUsedRowsForOnes, [row])
      else
        mat[row, c] = x[varCounter]
        varCounter = varCounter + 1
      end
    else
      mat[row, c] = R(0)
    end
  end
  return (R, mat)
end

# Given a matroid M with a basis B, this functions computes for a i in groundset(M)\B a circuit in B \cup i.
# This is a function used in the construction of the matrix underlying the realization space computation.
function fundamental_circuits_of_basis(M::Matroid, B::Vector{Int})
  remaining_elts = setdiff(matroid_groundset(M), B)
  fund_circs = Vector{Vector{Int}}()
  for i in remaining_elts
    push!(fund_circs, fundamental_circuit(M, B, i))
  end
  return fund_circs
end

@doc raw"""
    realization_space(
      M::Matroid;
      B::Union{GroundsetType,Nothing}=nothing,
      saturate::Bool=false,
      simplify::Bool=true,
      char::Union{Int,Nothing}=nothing,
      q::Union{Int,Nothing}=nothing,
      ground_ring::Ring=ZZ
    )::MatroidRealizationSpace

This function returns the data for the coordinate ring of the matroid realization space of the matroid `M`
as a `MatroidRealizationSpace`. This function has several optional parameters.

* `B` is a basis of M that specifies which columns of `realization_matrix(M)` form an identity matrix.
  The default is `nothing`, in which case the basis is chosen for you.

* `saturate` determines whether `defining_ideal(M)` should be saturated with respect to the semigroup
  generated by `inequations(M)`. The default is `false`. The saturation can be rather slow for large instances.

* `simplify` determines whether a reduced realization space is returned which means that the equations
  are used to eliminate variables as far as possible. The default is `true`.

* `char` specifies the characteristic of the coefficient ring. The returned realization space is then the space of all
  realizations over fields of characteristic `char`. The default is `nothing`.

* `q` is an integer and assumed to be a prime power `q=p^k`. The returned realization space is then the space of all
  realizations over the field ``GF(p^k)``. The default is `nothing`.

* `ground_ring` is a ring and specifies the ground_ring over which one wants to consider the realization space,
  e.g. `QQ` or `GF(p)`. The groud_ring `ZZ` means that we compute the space of realizations over all fields.
  The default is `ZZ`.

# Examples
```jldoctest
julia> M = fano_matroid();

julia> RS = realization_space(M)
The realization space is
  [0   1   1   1   1   0   0]
  [1   0   1   1   0   1   0]
  [1   0   1   0   1   0   1]
in the integer ring
within the vanishing set of the ideal
2ZZ

julia> realization_space(non_fano_matroid())
The realization space is
  [1   1   0   0   1   1   0]
  [0   1   1   1   1   0   0]
  [0   1   1   0   0   1   1]
in the integer ring
avoiding the zero loci of the polynomials
RingElem[2]

julia> realization_space(pappus_matroid(), char=0)
The realization space is
  [1   0   1   0   x2   x2                 x2^2    1    0]
  [0   1   1   0    1    1   -x1*x2 + x1 + x2^2    1    1]
  [0   0   0   1   x2   x1                x1*x2   x1   x2]
in the multivariate polynomial ring in 2 variables over QQ
avoiding the zero loci of the polynomials
RingElem[x1 - x2, x2, x1, x2 - 1, x1 + x2^2 - x2, x1 - 1, x1*x2 - x1 - x2^2]

julia> realization_space(uniform_matroid(3,6))
The realization space is
  [1   0   0   1    1    1]
  [0   1   0   1   x1   x3]
  [0   0   1   1   x2   x4]
in the multivariate polynomial ring in 4 variables over ZZ
avoiding the zero loci of the polynomials
RingElem[x1*x4 - x2*x3, x2 - x4, x1 - x3, x1*x4 - x1 - x2*x3 + x2 + x3 - x4, x3 - x4, x4 - 1, x3 - 1, x3, x4, x1 - x2, x2 - 1, x1 - 1, x1, x2]
```
"""
function realization_space(
  M::Matroid;
  B::Union{GroundsetType,Nothing}=nothing,
  saturate::Bool=false,
  simplify::Bool=true,
  char::Union{Int,Nothing}=nothing,
  q::Union{Int,Nothing}=nothing,
  ground_ring::Ring=ZZ
)::MatroidRealizationSpace
  if char != nothing && !is_prime(char) && char != 0
    error("The characteristic has to be 0 or a prime number.")
  end

  #Construct the base ring as F_p if q=p^k
  if q != nothing
    isprimepower, k, p = is_prime_power_with_data(q)
    if !isprimepower
      error("The given q has to be a prime power.")
    end
    if char != nothing && char != p
      error("The given characteristic doesn't match q.")
    else
      char = p
    end
  end

  if char == 0
    ground_ring = QQ
  elseif char != nothing
    ground_ring = GF(char)
  end

  rk = rank(M)
  n = length(M)

  goodM = isomorphic_matroid(M, [i for i in 1:n])

  Bs = bases(goodM)

  if !isnothing(B)
    goodB = sort!(Int.([M.gs2num[j] for j in B]))
  else
    goodB = find_good_basis_heuristically(goodM)
  end
  polyR, mat = realization_space_matrix(goodM, goodB, ground_ring)

  eqs = Vector{RingElem}()
  ineqs = Vector{RingElem}()

  #need to catch the corner-case if there are no variables at all
  if !(polyR isa MPolyRing)
    RS = MatroidRealizationSpace(ideal(polyR, [0]), ineqs, polyR, mat, char, q, ground_ring)
    set_attribute!(RS, :is_realizable, :true)
    return RS
  end

  for col in subsets(Vector(1:n), rk)
    col_det = det(mat[:, col])

    if total_degree(col_det) <= 0
      if col_det != 0 && col in Bs
        if is_unit(col_det)
          continue
        end
      elseif col_det != 0 # and col is not a basis
        # determinant nonzero but set not a basis
        push!(eqs, col_det)
      elseif col in Bs
        push!(ineqs, col_det)
        #determinant zero but set is a basis, i.e. M is not realizable
        RS = MatroidRealizationSpace(ideal(polyR, eqs), ineqs, polyR, nothing, char, q, ground_ring)
        set_attribute!(RS, :is_realizable, :false)
        return RS
      else
        continue
      end
    end

    if col in Bs
      push!(ineqs, col_det)
    else
      push!(eqs, col_det)
    end
  end

  def_ideal = ideal(polyR, eqs)
  def_ideal = ideal(groebner_basis(def_ideal))
  if isone(def_ideal)
    RS = MatroidRealizationSpace(def_ideal, ineqs, polyR, nothing, char, q, ground_ring)
    set_attribute!(RS, :is_realizable, :false)
    return RS
  end

  ineqs = gens_2_prime_divisors(ineqs)

  RS = MatroidRealizationSpace(def_ideal, ineqs, polyR, mat, char, q, ground_ring)

  if simplify
    RS = reduce_realization_space(RS)
  end

  if q != nothing && RS.ambient_ring isa MPolyRing
    I = RS.defining_ideal
    R = RS.ambient_ring
    eqs = Vector{RingElem}()
    for x in gens(R)
      push!(eqs, x^q - x)
    end
    I = I + ideal(R, eqs)
    RS.defining_ideal = I
    if isone(RS.defining_ideal)
      set_attribute!(RS, :is_realizable, :false)
      return RS
    end
  end

  if saturate
    RS.defining_ideal = stepwise_saturation(RS.defining_ideal, RS.inequations)
    if isone(RS.defining_ideal)
      set_attribute!(RS, :is_realizable, :false)
      return RS
    end
  end

  return RS
end

# A heuristic function that tries to find a sensible basis for the moduli space computation for which the defining ideal is not too complicated
function find_good_basis_heuristically(M::Matroid)
  bs = bases(M)
  cs = circuits(M)
  min_num_vars = length(cs) * rank(M)
  min_basis = bs[1]
  for b in bs
    current_num_vars = 0
    for c in cs, e in c
      if e in b
        current_num_vars += 1
      end
    end
    if current_num_vars < min_num_vars
      min_num_vars = current_num_vars
      min_basis = b
    end
  end
  return min_basis
end

# Return the prime divisors of f.
function poly_2_prime_divisors(f::RingElem)
  return map(first, factor(f))
end

# Return the unique prime divisors of the elements of Sgen, again no exponents.
function gens_2_prime_divisors(Sgens::Vector{<:RingElem})
  return unique!(vcat([poly_2_prime_divisors(f) for f in Sgens]...))
end

function stepwise_saturation(I::MPolyIdeal, Sgens::Vector{<:RingElem})
  for f in Sgens
    I = saturation(I, ideal([f]))
  end
  return I
end

@doc raw"""
    realization(M::Matroid; B::Union{GroundsetType,Nothing} = nothing,
      saturate::Bool=false,
      char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing
    )::MatroidRealizationSpace

This function tries to find one realization in the matroid realization space of
the matroid `M`. The output is again a `MatroidRealizationSpace`.

If the matroid is only realizable over an extension of the prime field the
extension field is specified as a splitting field of an irreducible polynomial.
Every root of this polynomial gives an equivalent realization of the matroid.

This function has several optional parameters. Note that one must input either
the characteristic or a specific field of definition for the realization.

* `B` is a basis of M that specifies which columns of `realization_matrix(M)`
  form the identity matrix. The default is `nothing`, in which case the basis is
  chosen for you.

* `char` specifies the characteristic of the coefficient ring, and is used to determine if the matroid
  is realizable over a field of this characteristic. The default is `nothing`.

* `q` is an integer, and when char = p, this input is used to determine whether the matroid
  is realizable over the finite field ``GF(p^{q})``. The default is `nothing`.

* `reduce` determines whether a reduced realization space is returned which means that the equations
  are used to eliminate variables as far as possible. The default is `true`.

* `saturate` determines whether `defining_ideal(M)` should be saturated with respect to the semigroup
  generated by `inequations(M)`. The default is `false`. This can be rather slow for large instances.

# Examples
```jldoctest
julia> realization(pappus_matroid(), char=0)
One realization is given by
  [1   0   1   0   2   2   4   1   0]
  [0   1   1   0   1   1   1   1   1]
  [0   0   0   1   2   3   6   3   2]
in the rational field

julia> realization(pappus_matroid(), q=4)
One realization is given by
  [1   0   1   0   x1 + 1   x1 + 1   x1    1        0]
  [0   1   1   0        1        1    1    1        1]
  [0   0   0   1   x1 + 1       x1    1   x1   x1 + 1]
in the multivariate polynomial ring in 1 variable over GF(2)
within the vanishing set of the ideal
Ideal (x1^2 + x1 + 1)

julia> realization(uniform_matroid(3,6), char=5)
One realization is given by
  [1   0   0   1   1   1]
  [0   1   0   1   4   3]
  [0   0   1   1   3   2]
in the prime field of characteristic 5
```
"""
function realization(
  M::Matroid;
  B::Union{GroundsetType,Nothing}=nothing,
  saturate::Bool=false,
  simplify::Bool=true,
  char::Union{Int,Nothing}=nothing,
  q::Union{Int,Nothing}=nothing,
)
  RS = realization_space(M; B=B, saturate=saturate, simplify=simplify, char=char, q=q)

  return realization(RS)
end

@doc raw"""
    realization(RS::MatroidRealizationSpace)

This function tries to find one realization in the matroid realization `RS`.
The output is again a `MatroidRealizationSpace`.
"""
function realization(RS::MatroidRealizationSpace)
  # If the matroid is not realizable we stop
  if RS.char == nothing && RS.q == nothing
    error("A field or characteristic must be specified")
  end

  if !is_realizable(RS)
    return RS
  end

  # If the ambient ring is not a polynomial ring we can't reduce and we stop
  R = RS.ambient_ring

  !(R isa MPolyRing) && return RS
  Inew = RS.defining_ideal
  eqs = copy(gens(Inew))

  if dim(Inew) == 0
    for p in minimal_primes(Inew)
      if !any(in(p), RS.inequations)
        Inew = p
        break
      end
    end
    RSnew = MatroidRealizationSpace(Inew, Vector{RingElem}(), R, RS.realization_matrix, RS.char, RS.q, RS.ground_ring)
    RSnew = reduce_realization_space(RSnew)
    RSnew.one_realization = true
    return RSnew
  end

  d = min(dim(Inew), nvars(R))
  ineqsnew = RS.inequations

  counter = 0
  base = 7
  if RS.char != nothing && RS.char > 0
    base = RS.char
  end
  upperbound = min(base^d, 10^3)
  while counter < upperbound
    eqsnew = copy(eqs)
    values = digits(counter; base=base, pad=d)
    counter += 1

    for i in 1:d
      push!(eqsnew, R[i] - values[i])
    end
    Inew = ideal(groebner_basis(ideal(R, eqsnew)))
    isone(Inew) && continue
    ineqsnew = copy(RS.inequations)
    for i in 1:length(ineqsnew)
      ineqsnew[i] = reduce(ineqsnew[i], gens(Inew))
    end
    (R(0) in ineqsnew) && continue
    Inew = stepwise_saturation(Inew, ineqsnew)
    isone(Inew) && continue
    break
  end
  counter == upperbound && d != 0 && return RS

  RSnew = MatroidRealizationSpace(Inew, ineqsnew, R, RS.realization_matrix, RS.char, RS.q, RS.ground_ring)
  RSnew = reduce_realization_space(RSnew)
  ineqsnew = RSnew.inequations
  if length(ineqsnew) > 0
    Inew = RSnew.defining_ideal
    Rnew = RSnew.ambient_ring
    ineqsnew = filter(p -> !isone(ideal(groebner_basis(Inew + ideal(Rnew, p)))), ineqsnew)
    RSnew = MatroidRealizationSpace(
      Inew, ineqsnew, Rnew, RSnew.realization_matrix, RSnew.char, RSnew.q, RSnew.ground_ring
    )
  end

  RSnew.one_realization = true

  return RSnew
end

#####################
# full reduction    #
#####################

function coefficient_v(v::RingElem, f::RingElem)
  isone(degree(f, v)) || return (false, parent(f)(0))
  return (true, coeff(f, [v], [1]))
end

function find_solution_v(
  v::RingElem, Igens::Vector{<:RingElem}, Sgens::Vector{<:RingElem}, R::MPolyRing
)
  with_v_deg_1 = [g for g in Igens if isone(degree(g, v))]
  length(with_v_deg_1) != 0 || return (false, R(0))

  for f in with_v_deg_1
    (a, den) = coefficient_v(v, f)
    a || return (false, R(0))
    fac_den = poly_2_prime_divisors(den)
    !issubset(fac_den, Sgens) && continue
    no_v = coeff(f, [v], [0])
    iszero(length(no_v)) && continue
    h = R(-1) * no_v
    return (true, h//den)
  end
  return (false, R(0))
end

# v is replaced by t in f
function sub_map(v::RingElem, t::RingElem, R::MPolyRing, xs::Vector{<:RingElem})
  xs_v = map(x -> x == v ? t : x, xs)
  return hom(R, fraction_field(R), identity, xs_v)
end

# replace v by t in f, only return the numerator.
function sub_v(v::RingElem, t::RingElem, f::RingElem, R::Ring, xs::Vector{<:RingElem})
  m = sub_map(v, t, R, xs)
  new_f = numerator(m(f))
  return new_f
end

# removes factors that are in the semigroup generated by Sgens
function clean(f::RingElem, R::MPolyRing, Sgens::Vector{<:RingElem})
  is_zero(f) && return f # TODO: move to the right place
  fFactors = factor(f)
  cleanf_arr = [k^e for (k, e) in fFactors if !(k in Sgens) || is_unit(k)]
  length(cleanf_arr) > 0 ? prod(cleanf_arr) : unit(fFactors)
end

# variables in ideal
function ideal_vars(Igens::Vector{<:RingElem})
  return unique!(vcat([vars(gen) for gen in Igens]...))
end

function n_new_Sgens(
  x::RingElem, t::RingElem, Sgens::Vector{<:RingElem}, R::Ring, xs::Vector{<:RingElem}
)
  preSgens = unique!([sub_v(x, t, f, R, xs) for f in Sgens])
  if R(0) in preSgens
    return [R(0)]
  end
  return gens_2_prime_divisors(preSgens)
end

function n_new_Igens(
  x::RingElem,
  t::RingElem,
  Igens::Vector{<:RingElem},
  Sgens::Vector{<:RingElem},
  R::Ring,
  xs::Vector{<:RingElem},
)
  preIgens = unique!([clean(sub_v(x, t, f, R, xs), R, Sgens) for f in Igens])
  return filter(!iszero, preIgens)
end

function matrix_clear_den_in_col(X::Oscar.MatElem, c::Int)
  Xc = [denominator(f) for f in X[:, c]]
  t = lcm(Xc)
  result = multiply_column!(X, t, c)
  return result
end

function matrix_clear_den(X::Oscar.MatElem)
  rs, cs = size(X)
  for c in 1:cs
    X = matrix_clear_den_in_col(X, c)
  end
  return X
end

function reduce_ideal_one_step(
  MRS::MatroidRealizationSpace, elim::Vector{<:RingElem}, fullyReduced::Bool
)
  Igens = gens(MRS.defining_ideal)
  Sgens = MRS.inequations
  R = MRS.ambient_ring
  FR = fraction_field(R)
  xs = gens(R)
  X = MRS.realization_matrix
  nr, nc = size(X)

  Ivars = ideal_vars(Igens)

  t = R(0)
  for x in Ivars
    a, t = find_solution_v(x, Igens, Sgens, R)
    a || continue

    phi = sub_map(x, t, R, xs)
    Sgens_new = n_new_Sgens(x, t, Sgens, R, xs)
    if length(Sgens_new) == 0
      Sgens_new = Vector{RingElem}()
    elseif R(0) in Sgens_new
      # 0 is in the inequations, hence M is not realizable.
      MRS.inequations = Sgens_new
      return (MRS, elim, true)
    end

    Igens_new = n_new_Igens(x, t, Igens, Sgens_new, R, xs)
    push!(elim, x)

    phiX = matrix(FR, [phi(X[i, j]) for i in 1:nr, j in 1:nc])
    nX_FR = matrix_clear_den(phiX)
    nX = matrix(R, [numerator(nX_FR[i, j]) for i in 1:nr, j in 1:nc])

    GBnew = collect(groebner_basis(ideal(R, Igens_new)))

    MRS_new = MatroidRealizationSpace(ideal(R, GBnew), Sgens_new, R, nX, MRS.char, MRS.q, MRS.ground_ring)

    return (MRS_new, elim, fullyReduced)
  end

  return (MRS, elim, true)
end

function reduce_realization_space(
  MRS::MatroidRealizationSpace,
  elim::Vector{RingElem}=Vector{RingElem}(),
  fullyReduced::Bool=false,
)

  #If there are no variables left, we don't reduce anything
  if !(MRS.ambient_ring isa MPolyRing)
    return MRS
  end

  (MRS, elim, fullyReduced) = reduce_ideal_one_step(MRS, elim, fullyReduced)

  !fullyReduced && return reduce_realization_space(MRS, elim, fullyReduced)

  R = MRS.ambient_ring
  xs = gens(R)
  cR = coefficient_ring(R)
  X = MRS.realization_matrix
  nr, nc = size(X)
  Igens = gens(MRS.defining_ideal)
  Sgens = MRS.inequations

  # 0 is in the inequations, thus it is not realizable
  if R(0) in Sgens
    MRS.realization_matrix = nothing
    set_attribute!(MRS, :is_realizable, :false)
    return MRS
  end
  xnew_str = ["x$i" for i in 1:length(xs) if !(xs[i] in elim)]

  if length(xnew_str) == 0
    phi = hom(R, cR, [cR(0) for i in 1:length(xs)])
    ambR = codomain(phi)

    if length(Igens) == 0
      Inew = ideal(ambR, [ambR(0)])
    else
      Inew = ideal(ambR, phi.(Igens))
    end
    normal_Sgens = phi.(Sgens)
    # This is a hack as phi.(Sgens) returns a list without types on the empty list
    if normal_Sgens == Any[]
      normal_Sgens = Vector{RingElem}()
    end
  else
    Rnew, xnew = polynomial_ring(coefficient_ring(R), length(xnew_str))

    zero_elim_var = elem_type(Rnew)[]
    j = 1
    for i in 1:length(xs)
      if xs[i] in elim
        push!(zero_elim_var, Rnew(0))
      else
        push!(zero_elim_var, xnew[j])
        j += 1
      end
    end

    phi = hom(R, Rnew, zero_elim_var)

    ambR = codomain(phi)
    if length(Igens) == 0
      Inew = ideal(ambR, ambR(0))
    else
      Inew = ideal(ambR, phi.(Igens))
    end

    if length(Sgens) == 0
      normal_Sgens = Vector{RingElem}()
    else
      Sgens_new = phi.(Sgens)
      normal_Sgens = [normal_form(g, Inew) for g in Sgens_new]
      if !(ambR(0) in normal_Sgens)
        normal_Sgens = gens_2_prime_divisors(Sgens_new)
      end
    end
  end

  if isone(Inew) || ambR(0) in normal_Sgens
    MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, nothing, MRS.char, MRS.q, MRS.ground_ring)
    set_attribute!(MRS_new, :is_realizable, :false)
    return MRS_new
  end

  Xnew = matrix(ambR, [phi(X[i, j]) for i in 1:nr, j in 1:nc])

  #Try to reduce the matrix one last time using the ideal and the inequations
  m, n = size(Xnew)
  A = base_ring(R)
  if length(gens(Inew)) > 0 && gens(Inew)[1] != 0
    for i in 1:m, j in 1:n
      if A == ZZ
        Xnew[i, j] = mod(Xnew[i, j], gens(Inew)[1])
      else
        Xnew[i, j] = reduce(Xnew[i, j], gens(Inew))
      end
    end
  end

  for j in 1:n
    g = gcd(Xnew[:, j]...)
    prime_divisors = poly_2_prime_divisors(g)
    for f in prime_divisors
      if f in normal_Sgens
        for i in 1:m
          Xnew[i, j] = Xnew[i, j] / f
        end
      end
    end
  end
  MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, Xnew, MRS.char, MRS.q, MRS.ground_ring)
  return MRS_new
end
