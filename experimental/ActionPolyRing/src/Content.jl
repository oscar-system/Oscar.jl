###############################################################################
#
#  Construction 
#
###############################################################################

##### Algebras #####

### Difference ###
@doc raw"""
    difference_polynomial_ring(R::Ring, elementary_symbols::Union{Vector{Symbol}, Int}, n_action_maps::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Construct the difference polynomial ring over the base ring `R` with the given elementary symbols and 
`n_action_maps` commuting endomorphisms. 

- If `elementary_symbols` is a vector of symbols, those names are used.  
- If it is an integer `m`, the symbols `u1, …, um` are generated automatically.  

In both cases, the jet variables that are initially available are those with jet `[0,…,0]`, one for each elementary symbol.

This method returns a tuple `(dpr, gens)` where `dpr` is the resulting difference polynomial ring and `gens` is the 
vector of initial jet variables.  

This constructor also accepts all keyword arguments of [`set_ranking!`](@ref) to control the ranking.

# Examples

```jldoctest
julia> R, variablesR = difference_polynomial_ring(QQ, 3, 4)
(Difference polynomial ring in 3 elementary symbols over QQ, DifferencePolyRingElem{QQFieldElem}[u1[0,0,0,0], u2[0,0,0,0], u3[0,0,0,0]])

julia> R
Difference polynomial ring in 3 elementary symbols u1, u2, u3
with 4 commuting endomorphisms
  over rational field

julia> variablesR
3-element Vector{DifferencePolyRingElem{QQFieldElem}}:
 u1[0,0,0,0]
 u2[0,0,0,0]
 u3[0,0,0,0]

julia> S, variablesS = difference_polynomial_ring(QQ, [:a, :b, :c], 4)
(Difference polynomial ring in 3 elementary symbols over QQ, DifferencePolyRingElem{QQFieldElem}[a[0,0,0,0], b[0,0,0,0], c[0,0,0,0]])

julia> S
Difference polynomial ring in 3 elementary symbols a, b, c
with 4 commuting endomorphisms
  over rational field

julia> variablesS
3-element Vector{DifferencePolyRingElem{QQFieldElem}}:
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]
```
"""
function difference_polynomial_ring(R::Ring, n_elementary_symbols::Int, n_action_maps::Int; kwargs...)
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, n_elementary_symbols, n_action_maps)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, [(i, zeros(Int, n_action_maps)) for i in 1:n_elementary_symbols])))
end

function difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, n_action_maps::Int; kwargs...)
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, elementary_symbols, n_action_maps)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, [(i, zeros(Int, n_action_maps)) for i in 1:length(elementary_symbols)])))
end

### Differential ###
@doc raw"""
    differential_polynomial_ring(R::Ring, elementary_symbols::Union{Vector{Symbol}, Int}, n_action_maps::Int) -> Tuple{DifferentialPolyRing, Vector{DifferentialPolyRingElem}}

Construct the differential polynomial ring over the base ring `R` with the given elementary symbols and 
`n_action_maps` commuting derivations. 

- If `elementary_symbols` is a vector of symbols, those names are used.  
- If it is an integer `m`, the symbols `u1, …, um` are generated automatically.  

In both cases, the jet variables that are initially available are those with jet `[0,…,0]`, one for each elementary symbol.

This method returns a tuple `(dpr, gens)` where `dpr` is the resulting differential polynomial ring and `gens` is the 
vector of initial jet variables.  

This constructor also accepts all keyword arguments of [`set_ranking!`](@ref) to control the ranking.

# Examples

```jldoctest
julia> R, variablesR = differential_polynomial_ring(QQ, 3, 4)
(Differential polynomial ring in 3 elementary symbols over QQ, DifferentialPolyRingElem{QQFieldElem}[u1[0,0,0,0], u2[0,0,0,0], u3[0,0,0,0]])

julia> R
Differential polynomial ring in 3 elementary symbols u1, u2, u3
with 4 commuting derivations
  over rational field

julia> variablesR
3-element Vector{DifferentialPolyRingElem{QQFieldElem}}:
 u1[0,0,0,0]
 u2[0,0,0,0]
 u3[0,0,0,0]

julia> S, variablesS = differential_polynomial_ring(QQ, [:a, :b, :c], 4)
(Differential polynomial ring in 3 elementary symbols over QQ, DifferentialPolyRingElem{QQFieldElem}[a[0,0,0,0], b[0,0,0,0], c[0,0,0,0]])

julia> S
Differential polynomial ring in 3 elementary symbols a, b, c
with 4 commuting derivations
  over rational field

julia> variablesS
3-element Vector{DifferentialPolyRingElem{QQFieldElem}}:
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]
```
"""
function differential_polynomial_ring(R::Ring, n_elementary_symbols::Int, n_action_maps::Int; kwargs...)
  dpr = DifferentialPolyRing{elem_type(typeof(R))}(R, n_elementary_symbols, n_action_maps)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, [(i, zeros(Int, n_action_maps)) for i in 1:n_elementary_symbols])))
end

function differential_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, n_action_maps::Int; kwargs...)
  dpr = DifferentialPolyRing{elem_type(typeof(R))}(R, elementary_symbols, n_action_maps)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, [(i, zeros(Int, n_action_maps)) for i in 1:length(elementary_symbols)])))
end

##### Elements #####

### Union ###
(apr::ActionPolyRing)() = elem_type(apr)(apr)
(apr::ActionPolyRing)(upre::AbstractAlgebra.Generic.UniversalPolyRingElem) = elem_type(apr)(apr, upre)
(apr::ActionPolyRing)(mpre::MPolyRingElem) = elem_type(apr)(apr, mpre)
(apr::ActionPolyRing)(a::T) where {T<:RingElement} = apr(base_ring(apr)(a))

### Difference ###
function (dpr::DifferencePolyRing{T})(a::DifferencePolyRingElem{T}) where {T}
  @req parent(a) === dpr "Wrong parent"
  return a
end

### Differential ###
function (dpr::DifferentialPolyRing{T})(a::DifferentialPolyRingElem{T}) where {T}
  @req parent(a) === dpr "Wrong parent"
  return a
end

###############################################################################
#
#  Basic ring functionality 
#
###############################################################################

### Union ###
function Base.deepcopy_internal(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem}, dict::IdDict)
    # Avoid deepcopying the parent as it may refer back to it in one of its dictionaries 
    pp = deepcopy_internal(data(dpre), dict)
    return typeof(dpre)(parent(dpre), pp)
end

### Difference ###
elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}

### Differential ###
elem_type(::Type{DifferentialPolyRing{T}}) where {T} = DifferentialPolyRingElem{T}

parent_type(::Type{DifferentialPolyRingElem{T}}) where {T} = DifferentialPolyRing{T}

### generic ###
@doc raw"""
    zero(A::ActionPolyRing)

Return the zero element of the action polynomial ring `A`.
"""
zero(apr::ActionPolyRing) = apr()

@doc raw"""
    one(A::ActionPolyRing)

Return the multiplicitive identity of the action polynomial ring `A`.
"""
one(apr::ActionPolyRing) = apr(one(base_ring(apr)))

base_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = AbstractAlgebra.Generic.UniversalPolyRing{T}

coefficient_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = parent_type(T)

is_square(apre::ActionPolyRingElem) = is_square(data(apre))

Base.sqrt(apre::ActionPolyRingElem; check::Bool = true) = parent(apre)(sqrt(data(apre); check = check))

is_irreducible(apre::ActionPolyRingElem) = is_irreducible(data(apre))

is_unit(apre::ActionPolyRingElem) = is_unit(data(apre))

characteristic(apr::ActionPolyRing) = characteristic(base_ring(apr))

###############################################################################
#
#  Factorisation
#
###############################################################################

function __wrap_factorization_apr(f::Fac{<:UniversalPolyRingElem}, apr::ActionPolyRing)
   res = Fac{elem_type(apr)}()
   res.unit = apr(f.unit)
   for (fact, expo) in f
      AbstractAlgebra.mulpow!(res, apr(fact), expo)
   end
   return res
end

factor_squarefree(apr::ActionPolyRingElem) = __wrap_factorization_apr(factor_squarefree(data(apr)), parent(apr))

factor(apr::ActionPolyRingElem) = __wrap_factorization_apr(factor(data(apr)), parent(apr))

###############################################################################
#
#  Field access and related stuff
#
###############################################################################

##### Rings #####

coefficient_ring(apr::ActionPolyRing) = coefficient_ring(base_ring(apr))

@doc raw"""
    n_elementary_symbols(A::ActionPolyRing) -> Int

Return the number of elementary symbols of the action polynomial ring `A`.
"""
n_elementary_symbols(apr::ActionPolyRing) = length(elementary_symbols(apr))

@doc raw"""
    elementary_symbols(A::ActionPolyRing) -> Vector{Symbol}

Return the elementary symbols of the action polynomial ring `A` as a vector.
"""
elementary_symbols(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.elementary_symbols

@doc raw"""
    n_action_maps(A::ActionPolyRing) -> Int

Return the number of elementary symbols of the action polynomial ring `A`.
"""
n_action_maps(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.n_action_maps

##### Elements #####

parent(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem}) = dpre.parent

###############################################################################
#
#  Basic polynomial functionality 
#
###############################################################################

@doc raw"""
    coeff(p::ActionPolyRingElem, i::Int)

Return coefficient of the `i`-th term of `p`.
"""
coeff(apre::ActionPolyRingElem, i::Int) = coeff(data(apre), __perm_for_sort_poly(apre)[i])

@doc raw"""
    exponent_vector(p::ActionPolyRingElem, i::Int)

Return the exponent vector of the `i`-th term of `p`.
"""
function exponent_vector(apre::ActionPolyRingElem, i::Int)
  v = exponent_vector(data(apre), __perm_for_sort_poly(apre)[i])
  perm = __perm_for_sort(parent(apre))
  L = length(perm)

  return (length(v) == L) ? v[perm] : vcat(v, zeros(Int, L - length(v)))[perm]
end

@doc raw"""
    monomial(p::ActionPolyRingElem, i::Int)

Return the `i`-th monomial of `p`.
"""
monomial(apre::ActionPolyRingElem, i::Int) = parent(apre)(monomial(data(apre), __perm_for_sort_poly(apre)[i]))

@doc raw"""
    term(p::ActionPolyRingElem, i::Int)

Return the `i`-th term of `p`.
"""
term(apre::ActionPolyRingElem, i::Int) = parent(apre)(term(data(apre), __perm_for_sort_poly(apre)[i]))

@doc raw"""
    is_monomial(p::ActionPolyRingElem)

Return `true` if `p` is a monomial and `false` otherwise.
"""
is_monomial(apre::ActionPolyRingElem) = is_monomial(data(apre))

@doc raw"""
    is_term(p::ActionPolyRingElem)

Return `true` if `p` is a term, i.e. a non-zero multiple of a monomial, and `false` otherwise.
"""
is_term(apre::ActionPolyRingElem) = is_term(data(apre))

@doc raw"""
    length(p::ActionPolyRingElem) -> Int

Return the length of `p`, i.e. the number of terms of `p`.
"""
length(apre::ActionPolyRingElem) = length(data(apre))

@doc raw"""
    total_degree(p::ActionPolyRingElem) -> Int

Return the total degree of `p`.
"""
total_degree(apre::ActionPolyRingElem) = total_degree(data(apre))

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int, jet::Vector{Int}) -> Int

Return the degree of the polynomial `p` in the jet variable specified by `i` and `jet`. If this jet variable
is valid but still untracked, return $0$. This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables).
"""
function degree(apre::ActionPolyRingElem, i::Int, jet::Vector{Int})
  apr = parent(apre)
  upr = base_ring(apr)
  @req __is_valid_jet(apr, i, jet) "Invalid jet variable"
  if is_zero(apre)
    return -1
  end
  jtv = __jtv(apr)
  if haskey(jtv, (i,jet))
    idx = findfirst(var -> var == data(jtv[(i,jet)]), gens(upr))
    return degree(data(apre), gen(upr, idx))
  end
  return 0
end

function degree(apre::ActionPolyRingElem{T}, var::ActionPolyRingElem{T}) where {T}
  check_parent(apre, var)
  @req is_gen(var) "Not a jet variable"
  return degree(apre, __vtj(parent(apre))[var])
end

degree(apre::ActionPolyRingElem, jet_idx::Tuple{Int, Vector{Int}}) = degree(apre, jet_idx...)

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int)

Return the degree of the polynomial `p` in the, among the currently tracked jet variables, 'i'-th largest one.
The index of the jet variable may also be passed as a tuple. Alternatively, the jet variable may be passed right
away instead of its index.
"""
degree(apre::ActionPolyRingElem, i::Int) = degree(apre, __vtj(parent(apre))[gen(parent(apre), i)])

@doc raw"""
    degrees(p::ActionPolyRingElem)

Return an array of the degrees of the polynomial `p` in terms of each jet variable. The jet variables are sorted with
respect to the ranking of the action polynomial ring containing `p`, leading with the largest variable.
"""
degrees(apre::ActionPolyRingElem) = [degree(apre, i) for i in 1:ngens(parent(apre))]

@doc raw"""
    is_constant(p::ActionPolyRingElem)

Return `true` if `p` is a degree zero polynomial or the zero polynomial, i.e. a constant polynomial. 
"""
is_constant(apre::ActionPolyRingElem) = is_constant(data(apre))

@doc raw"""
    vars(p::ActionPolyRingElem)

Return the jet variables actually occuring in `p` as a vector. The jet variables are sorted with respect to the
ranking of the action polynomial ring containing `p`, leading with the largest jet variable.
"""
vars(apre::ActionPolyRingElem) = sort!(parent(apre).(vars(data(apre))); rev = true)

@doc raw"""
    is_gen(p::ActionPolyRingElem)

Return true if `p` is a jet variable in an action polynomial ring.
"""
is_gen(apre::ActionPolyRingElem) = is_gen(data(apre))

@doc raw"""
    gen(A::ActionPolyRing, i::Int, jet::Vector{Int})

Return the jet variable of the action polynomial ring `A` specified by `i` and `jet`. If this jet variable was untracked,
it is tracked afterwards. The jet variable may also be specified by a tuple, see [Specifying jet variables](@ref specifying_jet_variables).
Additionally, the jet variable may also be specified by an integer but the corresponding method [`gen(A::ActionPolyRing, i::Int)`](@ref)
slightly differs in its functionality, as it cannot create new jet variables.

# Examples

```jldoctests
julia> dpr = differential_polynomial_ring(ZZ, [:a, :b, :c], 4)[1]; gen(dpr, 1, [3,1,0,0])
a[3,1,0,0]

julia> gen(dpr, (1, [3,2,0,0]))
a[3,2,0,0]

julia> gens(dpr)
5-element Vector{DifferentialPolyRingElem{ZZRingElem}}:
 a[3,2,0,0]
 a[3,1,0,0]
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]
```
"""
function gen(apr::ActionPolyRing, i::Int, jet::Vector{Int})
  @req __is_valid_jet(apr, i, jet) "invalid jet variable"
  jtv = __jtv(apr)
  if haskey(jtv, (i, jet))
    return deepcopy(jtv[(i, jet)])
  end
  __add_new_jetvar!(apr, i, jet)
  return deepcopy(jtv[(i, jet)])
end

gen(apr::ActionPolyRing, jet_idx::Tuple{Int, Vector{Int}}) = gen(apr, jet_idx...)

@doc raw"""
    gen(A::ActionPolyRing, i::Int)

Among the currently tracked jet variables of `A`, return the `i`-th largest one.

# Examples

```jldoctests
julia> dpr = difference_polynomial_ring(ZZ, [:a, :b, :c], 4)[1]; gen(dpr, 2)
b[0,0,0,0]

julia> set_ranking!(dpr; partition = [[0,1,1],[1,0,0]]); gen(dpr, 2)
c[0,0,0,0]
```
"""
gen(apr::ActionPolyRing, i::Int) = gens(apr)[i]

@doc raw"""
    gens(A::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})

Return the jet variables of the action polynomial ring `A` specified by the entries of `jet_idxs` as
a vector and track all new jet variables.

# Examples

```jldoctests
julia> dpr = differential_polynomial_ring(ZZ, [:a, :b, :c], 4)[1]; gens(dpr, [(1, [3,1,0,0]), (1, [3,2,0,0])])
2-element Vector{DifferentialPolyRingElem{ZZRingElem}}:
 a[3,1,0,0]
 a[3,2,0,0]

julia> gens(dpr)
5-element Vector{DifferentialPolyRingElem{ZZRingElem}}:
 a[3,2,0,0]
 a[3,1,0,0]
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]
```
"""
gens(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}}) = [gen(apr, jet_idx) for jet_idx in jet_idxs]

@doc raw"""
    gens(A::ActionPolyRing)

Return the currently tracked jet variables of the action polynomial ring `A` as a vector. The jet variables are
sorted with respect to the ranking of `A`, leading with the largest jet variable.

# Examples

```jldoctests
julia> dpr = difference_polynomial_ring(ZZ, [:a, :b, :c], 4)[1]; gens(dpr)
3-element Vector{DifferencePolyRingElem{ZZRingElem}}:
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]

julia> set_ranking!(dpr; partition = [[0,1,1],[1,0,0]]); gens(dpr)
3-element Vector{DifferencePolyRingElem{ZZRingElem}}:
 b[0,0,0,0]
 c[0,0,0,0]
 a[0,0,0,0]

julia> gens(dpr, [(1, [1,1,1,1]), (2, [1,1,1,1])]);

julia> gens(dpr)
5-element Vector{DifferencePolyRingElem{ZZRingElem}}:
 b[1,1,1,1]
 b[0,0,0,0]
 c[0,0,0,0]
 a[1,1,1,1]
 a[0,0,0,0]
```
"""
gens(apr::ActionPolyRing) = sort!(collect(values(deepcopy(__jtv(apr)))); rev = true)

number_of_generators(apr::ActionPolyRing) = number_of_generators(base_ring(apr))

number_of_variables(apr::ActionPolyRing) = number_of_variables(base_ring(apr))

@doc raw"""
    getindex(A::ActionPolyRing, i::Int, jet::Vector{Int})

Alias for [`gen(A, i, jet)`](@ref gen(A::ActionPolyRing, i::Int, jet::Vector{Int})).
"""
getindex(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = gen(apr, i, jet)

getindex(apr::ActionPolyRing, jet_idx::Tuple{Int, Vector{Int}}) = gen(apr, jet_idx...)

@doc raw"""
    var_index(p::ActionPolyRingElem)

Return the integer `i` such that `p` is the `i`-th largest currently tracked jet variable. If `p` is not a jet variable an exception is raised. 
"""
function var_index(x::ActionPolyRingElem)
  @req is_gen(x) "Not a jet variable in var_index"
  return findfirst(==(x), gens(parent(x)))
end

@doc raw"""
    constant_coefficient(p::ActionPolyRingElem{T}) -> T

Return the constant coefficient of `p`. Does not throw an error for the zero polynomial like [`trailing_coefficient`](@ref trailing_coefficient(apre::ActionPolyRingElem)).
"""
constant_coefficient(apre::ActionPolyRingElem) = constant_coefficient(data(apre))

@doc raw"""
    leading_coefficient(p::ActionPolyRingElem{T}) -> T 

Return the leading coefficient of the polynomial `p`, i.e. the coefficient of the first (with respect to the ranking of the action polynomial ring containing it) nonzero term.
"""
function leading_coefficient(apre::ActionPolyRingElem{T}) where {T}
  @req length(apre) > 0 "Zero polynomial does not have a leading coefficient"
  return coeff(apre, 1) 
end

@doc raw"""
    leading_monomial(p::ActionPolyRingElem)

Return the leading monomial of the polynomial `p` with respect to the ranking of the action polynomial ring containing it.
"""
function leading_monomial(apre::ActionPolyRingElem) 
  @req length(apre) > 0 "Zero polynomial does not have a leading monomial"
  return monomial(apre, 1)
end

@doc raw"""
    leading_term(p::ActionPolyRingElem)

Return the leading term of the polynomial `p` with respect to the ranking of the action polynomial ring containing it.
"""
function leading_term(apre::ActionPolyRingElem) 
  @req length(apre) > 0 "Zero polynomial does not have a leading term"
  return term(apre, 1)
end

@doc raw"""
    trailing_coefficient(p::ActionPolyRingElem{T}) -> T

Return the trailing coefficient of the polynomial `p`, i.e. the coefficient of the last (with respect to the ranking of the action polynomial ring containing it) nonzero term, or zero if the polynomial is zero.
"""
function trailing_coefficient(apre::ActionPolyRingElem{T}) where {T}
  len = length(apre)
  @req len > 0 "Zero polynomial does not have a trailing coefficient"
  return coeff(apre, len) 
end

@doc raw"""
    trailing_monomial(p::ActionPolyRingElem)

Return the trailing monomial of the polynomial `p` with respect to the ranking of the action polynomial ring containing it.
"""
function trailing_monomial(apre::ActionPolyRingElem) 
  len = length(apre)
  @req len > 0 "Zero polynomial does not have a trailing monomial"
  return monomial(apre, len)
end

@doc raw"""
    trailing_term(p::ActionPolyRingElem)

Return the leading term of the polynomial `p` with respect to the ranking of the action polynomial ring containing it.
"""
function trailing_term(apre::ActionPolyRingElem) 
  len = length(apre)
  @req len > 0 "Zero polynomial does not have a trailing term"
  return term(apre, len)
end

@doc raw"""
    tail(p::ActionPolyRingElem)

Return the tail of `p`, i.e. `p`, i.e. return `p` without its leading term with respect to the ranking of the action polynomial ring containing it.
"""
tail(apre::ActionPolyRingElem) = apre - leading_term(apre)

function derivative(apre::ActionPolyRingElem{T}, var::ActionPolyRingElem{T}) where {T}
  check_parent(apre, var)
  @req is_gen(var) "Not a jet variable"
  return derivative(apre, __vtj(parent(var))[var])
end

@doc raw"""
    derivative(p::ActionPolyRing, i::Int, jet::Vector{Int})

Return the derivative of `p` with respect to the jet variable specified by `i` and `jet`.
This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables).
"""
function derivative(apre::ActionPolyRingElem, i::Int, jet::Vector{Int})
  apr = parent(apre)
  @req __is_valid_jet(apr, i, jet) "Invalid jet variable"
  jtv = __jtv(apr)
  if haskey(jtv, (i, jet))
    return apr(derivative(data(apre), data(jtv[(i, jet)])))
  end
  return zero(apr)
end

derivative(apre::ActionPolyRingElem, jet_idx::Tuple{Int, Vector{Int}}) = derivative(apre, jet_idx...)

derivative(apre::ActionPolyRingElem, i::Int) = derivative(apre, gen(parent(apre), i))

###############################################################################
#
#  Action polynomial functionality 
#
###############################################################################

@doc raw"""
    diff_action(p::DifferencePolyRingElem, i::Int)

Apply the `i`-th endomorphism to the polynomial `p`.

# Examples

```jldoctests
julia> dpr, (a,b,c) = difference_polynomial_ring(ZZ, [:a, :b, :c], 2); f = -2*a*b + 3*a*b^2;

julia> diff_action(3*a, 1)
3*a[1,0]

julia> diff_action(3*a, 2)
3*a[0,1]

julia> diff_action(f, 1)
(3*b[1,0]^2 - 2*b[1,0])*a[1,0]

julia> diff_action(f, 2)
(3*b[0,1]^2 - 2*b[0,1])*a[0,1]
```
"""
function diff_action(dpre::DifferencePolyRingElem{T}, i::Int) where {T}
  dpr = parent(dpre)
  @req i in 1:n_action_maps(dpr) "index out of range"
  d = fill(0, n_action_maps(parent(dpre)))
  d[i] = 1
  return diff_action(dpre, d)
end

@doc raw"""
    diff_action(p::DifferentialPolyRingElem, i::Int)

Apply the `i`-th derivation to the polynomial `p`.

# Examples

```jldoctests
julia> dpr, (a,b,c) = differential_polynomial_ring(ZZ, [:a, :b, :c], 2); f = -2*a*b + 3*a*b^2;

julia> diff_action(3*a, 1)
3*a[1,0]

julia> diff_action(3*a, 2)
3*a[0,1]

julia> diff_action(f, 1)
(3*b[0,0]^2 - 2*b[0,0])*a[1,0] + 6*b[1,0]*a[0,0]*b[0,0] - 2*b[1,0]*a[0,0]

julia> diff_action(f, 2)
(3*b[0,0]^2 - 2*b[0,0])*a[0,1] + 6*b[0,1]*a[0,0]*b[0,0] - 2*b[0,1]*a[0,0]
```
"""
function diff_action(dpre::DifferentialPolyRingElem{T}, i::Int) where {T}
  dpr = parent(dpre)
  @req i in 1:n_action_maps(dpr) "index out of range"

  # Remove constant term: d_i(1) = 0
  dpre -= constant_coefficient(dpre)

  if is_zero(dpre)
    return dpre
  end
 
  d = fill(0, n_action_maps(dpr))
  d[i] = 1
  dpre_vars_idxs = map(var -> __vtj(dpr)[var], vars(dpre))
  add_new_vars_idx = filter(new_idx -> !haskey(__jtu_idx(dpr), new_idx), map(idx -> (idx[1], idx[2] + d), dpre_vars_idxs))
  if !is_empty(add_new_vars_idx)
    __add_new_jetvar!(dpr, add_new_vars_idx)
  end
  jtu = __jtu_idx(dpr)
  old_to_new_pos = Dict{Int, Int}() #Dictionary that links the positions of old with new variables
  for i in 1:length(dpre_vars_idxs)
    (i, idx) = dpre_vars_idxs[i]
    new_idx = jtu[(i, idx + d)]
    old_to_new_pos[jtu[(i, idx)]] = new_idx
  end 

  upr = base_ring(dpr)
  upre = data(dpre)
  C = MPolyBuildCtx(upr)

  for term in terms(upre)
    coeff_t = coeff(term, 1)
    vars_t = vars(term)
    ev = append!(exponent_vector(term, 1), fill(0, ngens(upr) - length(exponent_vector(term, 1))))

    for var in vars_t
      var_pos = findfirst(==(var), gens(upr))
      exp = ev[var_pos]
      shifted_var_pos = old_to_new_pos[var_pos]

      new_exp_vec = copy(ev)
      new_exp_vec[var_pos] -= 1
      new_exp_vec[shifted_var_pos] += 1

      push_term!(C, coeff_t * exp, new_exp_vec)
    end
  end
  return dpr(finish(C))
end

@doc raw"""
    diff_action(p::ActionPolyRingElem, d::Vector{Int}) -> ActionPolyRingElem

Successively apply the `i`-th diff-action `d[i]`-times to the polynomial `p`, where $i = 1, \ldots, \mathrm{length}(d)$. 
"""
function diff_action(dpre::DifferencePolyRingElem{T}, d::Vector{Int}) where {T}
  dpr = parent(dpre)
  @req length(d) == n_action_maps(dpr) && all(>=(0), d) "Invalid vector of multiplicities"
  if is_zero(d)
    return dpre
  end
  dpre_vars_idxs = map(var -> __vtj(dpr)[var], vars(dpre))
  add_new_vars_idx = filter(new_idx -> !haskey(__jtu_idx(dpr), new_idx), map(idx -> (idx[1], idx[2] + d), dpre_vars_idxs))
  if !is_empty(add_new_vars_idx)
    __add_new_jetvar!(dpr, add_new_vars_idx)
  end

  jtu = __jtu_idx(dpr)
  old_to_new_pos = Dict{Int, Int}() #Dictionary that links the positions of old with new variables
  for i in 1:length(dpre_vars_idxs)
    (i, idx) = dpre_vars_idxs[i]
    new_idx = jtu[(i, idx + d)]
    old_to_new_pos[jtu[(i, idx)]] = new_idx
  end 
  
  upr = base_ring(dpr)
  upre = data(dpre)
  C = MPolyBuildCtx(upr)

  for term in terms(upre)
    coeff_t = coeff(term, 1)
    vars_t = vars(term)
    ev = append!(exponent_vector(term, 1), fill(0, ngens(upr) - length(exponent_vector(term, 1))))
    new_exp_vec = copy(ev)

    for var in vars_t
      var_pos = findfirst(==(var), gens(upr))
      shifted_var_pos = old_to_new_pos[var_pos]
      
      new_exp_vec[shifted_var_pos] = ev[var_pos]
      new_exp_vec[var_pos] = 0 
    end
    
    push_term!(C, coeff_t , new_exp_vec)
  end
  return dpr(finish(C))
end

function diff_action(dpre::DifferentialPolyRingElem{T}, d::Vector{Int}) where {T}
  len = length(d)
  @req len == n_action_maps(parent(dpre)) && all(>=(0), d) "Invalid vector of diff multiplicities"
  res = dpre
  for i in 1:len
    for j in 1:d[i]
      res = diff_action(res, i)
    end
  end
  return res
end

@doc"""
    initial(p::ActionPolyRingElem)

Return the initial of the polynomial `p`, i.e. the leading coefficient of `p` regarded as a univariate polynomial in its leader.
"""
function initial(apre::ActionPolyRingElem)
  if is_constant(apre)
    return apre
  end
  res = parent(apre)()
  ld = leader(apre)
  ld_ind = var_index(ld)
  d = degree(apre, ld_ind)
  for (t,e) in zip(terms(apre), exponents(apre))
    if e[ld_ind] < d
      break
    end
    res += remove(t, ld)[2] 
  end
  return res
end

@doc raw"""
    leader(p::ActionPolyRingElem)

Return the leader of the polynomial `p`, that is the largest jet variable with respect to the ranking of `parent(p)`. If `p` is constant, an error is raised.
"""
function leader(apre::ActionPolyRingElem)
  @req !is_constant(apre) "A constant polynomial has no leader"
  return maximum(vars(apre))
end

###############################################################################
#
#  Discriminant and resultant
#
###############################################################################

function resultant(r1::ActionPolyRingElem, r2::ActionPolyRingElem, var::ActionPolyRingElem)
  check_parent(r1, r2) && check_parent(r1, var)
  @req is_gen(var) "Not a jet variable"
  return resultant(r1, r2, __vtj(parent(var))[var])
end

@doc raw"""
    resultant(f::ActionPolyRingElem, g::ActionPolyRingElem, i::Int, jet::Vector{Int})

Return the resultant of `f` and `g` regarded as univariate polynomials in the jet variable specified by `i` and
`jet`. This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables).
"""
function resultant(r1::ActionPolyRingElem, r2::ActionPolyRingElem, i::Int, jet::Vector{Int})
  check_parent(r1, r2)
  S = parent(r1)
  @req __is_valid_jet(S, i, jet) "Invalid jet variable"
  jtv = __jtv(S)
  if haskey(jtv, (i, jet))
    return S(resultant(data(r1), data(r2), data(jtv[(i, jet)])))
  end
  if is_zero(r1) || is_zero(r2)
    return zero(S)
  end
  return one(S)
end

resultant(r1::ActionPolyRingElem, r2::ActionPolyRingElem, jet_idx::Tuple{Int,Vector{Int}}) = resultant(r1, r2, jet_idx...)
resultant(r1::ActionPolyRingElem, r2::ActionPolyRingElem, i::Int) = resultant(r1, r2, gen(parent(r1), i))

@doc raw"""
    discriminant(p::ActionPolyRingElem)

Return the discriminant of `p`.
"""
function discriminant(p::ActionPolyRingElem)
  is_constant(p) && return zero(parent(p))

  ld = leader(p)
  
  if degree(p, ld) % 4 in (0,1)
    return divexact(resultant(p, derivative(p, ld), ld), initial(p))
  end
  
  return -divexact(resultant(p, derivative(p, ld), ld), initial(p))
end

###############################################################################
#
#  Univariate functionality 
#
###############################################################################

@doc raw"""
    is_univariate(p::ActionPolyRingElem)

Return `true` if `p` is a polynomial in a single jet variable and `false` otherwise.
"""
is_univariate(apre::ActionPolyRingElem) = is_univariate(data(apre))

function is_univariate_with_data(apre::ActionPolyRingElem)
  flag, gen_idx = is_univariate_with_data(data(apre))
  is_zero(gen_idx) && return flag, gen_idx
  return flag, findfirst(==(gen_idx), __perm_for_sort(parent(apre)))
end

@doc raw"""
    is_univariate(p::ActionPolyRing)

Return `false`, since an action polynomial ring cannot be univariate.
"""
is_univariate(apr::ActionPolyRing) = false

@doc raw"""
    to_univariate(R::PolyRing{T}, p::ActionPolyRingElem{T}) where {T <: RingElement}

Assuming the polynomial `p` is actually a univariate polynomial, convert the polynomial to a univariate polynomial in the
given univariate polynomial ring `R`. An exception is raised if the polynomial `p` involves more than one jet variable.
"""
to_univariate(R::PolyRing{T}, apre::ActionPolyRingElem{T}) where {T <: RingElement} = to_univariate(R, data(apre))

@doc raw"""
    to_univariate(p::ActionPolyRingElem)

Assuming the polynomial `p` is actually a univariate polynomial in the jet variable `x`, convert the polynomial to a univariate
polynomial in a univariate polynomial ring over the same coefficient ring in the variable `x`. If `p` is constant, it is considered
to be a polynomial in the largest tracked jet variable of its parent. An exception is raised if the polynomial `p` involves
more than one jet variable.
"""
function to_univariate(apre::ActionPolyRingElem)
  flag, var = is_univariate_with_data(apre)
  flag || error("Polynomial is not univariate.")
  if is_zero(var)
    var = 1
  end
  apr = parent(apre)
  x = symbols(apr)[var]
  R, _ = coefficient_ring(apr)[x]
  return to_univariate(R, apre)
end

function univariate_coefficients(r::ActionPolyRingElem, var::ActionPolyRingElem)
  check_parent(r, var)
  @req is_gen(var) "Not a jet variable"
  return univariate_coefficients(r, __vtj(parent(var))[var])
end
  
@doc raw"""
    univariate_coefficients(p::ActionPolyRingElem, i::Int, jet::Vector{Int}) 

Return the coefficient vector of `p` regarded as a univariate polynomial in the jet variable specified by `i` and
`jet`, leading with the constant coefficient. This method allows all versions described in
[Specifying jet variables](@ref specifying_jet_variables).
"""
function univariate_coefficients(r::ActionPolyRingElem, i::Int, jet::Vector{Int})
  d = degree(r, i, jet)
  d == 0 && return [r]
  res = [zero(r) for _ in 1:d+1]
  var = __jtv(parent(r))[(i, jet)]
  v_idx = var_index(var)
  for (t, e) in zip(terms(r), exponents(r)) 
    @inbounds res[e[v_idx] + 1] += remove(t, var)[2]
  end
  return res
end

univariate_coefficients(r::ActionPolyRingElem, jet_idx::Tuple{Int, Vector{Int}}) = univariate_coefficients(r, jet_idx...)
univariate_coefficients(r::ActionPolyRingElem, i::Int) = univariate_coefficients(r, gen(parent(r), i))

###############################################################################
#
#  Evaluation
#
###############################################################################

# We let UnivPoly handle the calculations, so we need to permute argument vectors
function __permute_vals(S::ActionPolyRing, A::Vector{T}) where {T}
  perm = invperm(__perm_for_sort(S))
  n = nvars(S)
  m = length(A) 
  m > n && error("Too many values")
  return vcat(A, [zero(T) for _ in 1:(n - m)])[perm]
end

# Full substitutions with zero padding
evaluate(a::ActionPolyRingElem{T}, A::Vector{T}) where {T <: RingElement} = evaluate(data(a), __permute_vals(parent(a), A))

@doc raw"""
    evaluate(a::ActionPolyRingElem{T}, vals::Vector{V}) where {T <: RingElement, V <: Ringelement}

Evaluate the polynomial expression by substituting in the supplied values in the array `vals` for each of the
tracked jet variables. The evaluation will succeed if multiplication is defined between elements of the
coefficient ring of `a` and elements of `vals`.
"""
evaluate(a::ActionPolyRingElem{T}, vals::Vector{V}) where {T <: RingElement, V <: RingElement} = evaluate(data(a), __permute_vals(parent(a), vals))

(a::ActionPolyRingElem{T})() where {T <: RingElement} = evaluate(a, T[])
(a::ActionPolyRingElem{T})(vals::T...) where {T <: RingElement} = evaluate(a, collect(vals))
(a::ActionPolyRingElem{T})(val::V, vals::V...) where {T <: RingElement, V <: Union{Integer, Rational, AbstractFloat}} = evaluate(a, [val, vals...])
(a::ActionPolyRingElem{T})(vals::Union{NCRingElem, RingElement}...) where {T <: RingElement} = evaluate(a, collect(vals))

# Partial substitutions
@doc raw"""
    evaluate(a::ActionPolyRingElem{T}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, V <: Ringelement}

Evaluate the polynomial expression by substituting in the supplied values in the array `vals` for
the corresponding jet variables specified by the indices given by the array `vars`; see
[Specifying jet variables](@ref specifying_jet_variables). The evaluation will succeed if
multiplication is defined between elements of the coefficient ring of `a` and elements of `vals`.
"""
function evaluate(a::ActionPolyRingElem{T}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, V <: RingElement}
    S = parent(a)
  per = __perm_for_sort(S)
  return S(evaluate(data(a), map(x -> per[x], vars), vals))
end

@doc raw"""
    evaluate(a::PolyT, vars::Vector{PolyT}, vals::Vector{V}) where {PolyT <: ActionPolyRingElem, V <: Ringelement}

Evaluate the polynomial expression by substituting in the supplied values in the array `vals` for
the corresponding jet variables from the vector `vars`; see
[Specifying jet variables](@ref specifying_jet_variables). The evaluation will succeed if
multiplication is defined between elements of the coefficient ring of `a` and elements of `vals`.
"""
evaluate(a::PolyT, vars::Vector{PolyT}, vals::Vector{V}) where {PolyT <: ActionPolyRingElem, V <: RingElement} = parent(a)(evaluate(a, [var_index(x) for x in vars], vals))

###############################################################################
#
#  Iterators
#
###############################################################################

@doc raw"""
    coefficients(p::ActionPolyRingElem)

Return an iterator for the coefficients of `p` with respect to the ranking of the parent of `p`.

# Examples

```jldoctests
julia> dpr, (a,b,c) = difference_polynomial_ring(ZZ, [:a, :b, :c], 4; partition = [[0,1,1],[1,0,0]]); f = -2*a*b + a*c + 3*b^2;

julia> cf = coefficients(f)
Coefficients iterator of 3*b[0,0,0,0]^2 - 2*a[0,0,0,0]*b[0,0,0,0] + c[0,0,0,0]*a[0,0,0,0]

julia> collect(cf)
3-element Vector{ZZRingElem}:
 3
 -2
 1
```
"""
coefficients(a::ActionPolyRingElem{T}) where {T} = ActionPolyCoeffs{typeof(a)}(a)

@doc raw"""
    exponents(p::ActionPolyRingElem)

Return an iterator for the exponents of `p` with respect to the ranking of the parent of `p`.

# Examples

```jldoctests
julia> dpr, (a,b,c) = difference_polynomial_ring(ZZ, [:a, :b, :c], 4; partition = [[0,1,1],[1,0,0]]); f = -2*a*b + a*c + 3*b^2;

julia> ef = exponents(f)
Exponents iterator of 3*b[0,0,0,0]^2 - 2*a[0,0,0,0]*b[0,0,0,0] + c[0,0,0,0]*a[0,0,0,0]

julia> collect(ef)
3-element Vector{Vector{Int64}}:
 [2, 0, 0]
 [1, 0, 1]
 [0, 1, 1]
```
"""
exponents(a::ActionPolyRingElem{T}) where {T} = ActionPolyExponentVectors{typeof(a)}(a)

@doc raw"""
    monomials(p::ActionPolyRingElem)

Return an iterator for the monomials of `p` with respect to the ranking of the parent of `p`.

# Examples

```jldoctests
julia> dpr, (a,b,c) = difference_polynomial_ring(ZZ, [:a, :b, :c], 4; partition = [[0,1,1],[1,0,0]]); f = -2*a*b + a*c + 3*b^2;

julia> mf = monomials(f)
Monomials iterator of 3*b[0,0,0,0]^2 - 2*a[0,0,0,0]*b[0,0,0,0] + c[0,0,0,0]*a[0,0,0,0]

julia> collect(mf)
3-element Vector{DifferencePolyRingElem{ZZRingElem}}:
 b[0,0,0,0]^2
 a[0,0,0,0]*b[0,0,0,0]
 a[0,0,0,0]*c[0,0,0,0]
```
"""
monomials(a::ActionPolyRingElem{T}) where {T} = ActionPolyMonomials{typeof(a)}(a)

@doc raw"""
    terms(p::ActionPolyRingElem)

Return an iterator for the terms of `p` with respect to the ranking of the parent of `p`.

# Examples

```jldoctests
julia> dpr, (a,b,c) = difference_polynomial_ring(ZZ, [:a, :b, :c], 4; partition = [[0,1,1],[1,0,0]]); f = -2*a*b + a*c + 3*b^2;

julia> tf = terms(f)
Terms iterator of 3*b[0,0,0,0]^2 - 2*a[0,0,0,0]*b[0,0,0,0] + c[0,0,0,0]*a[0,0,0,0]

julia> collect(tf)
3-element Vector{DifferencePolyRingElem{ZZRingElem}}:
 3*b[0,0,0,0]^2
 -2*a[0,0,0,0]*b[0,0,0,0]
 a[0,0,0,0]*c[0,0,0,0]
```
"""
terms(a::ActionPolyRingElem{T}) where {T} = ActionPolyTerms{typeof(a)}(a)

__iter_helper(f, poly, state) = state < length(poly) ? (f(poly, state + 1), state + 1) : nothing

Base.iterate(x::ActionPolyCoeffs, state=0) = __iter_helper(coeff, x.poly, state)
Base.iterate(x::ActionPolyExponentVectors, state=0) = __iter_helper(exponent_vector, x.poly, state)
Base.iterate(x::ActionPolyMonomials, state=0) = __iter_helper(monomial, x.poly, state)
Base.iterate(x::ActionPolyTerms, state=0) = __iter_helper(term, x.poly, state)

Base.length(x::Union{ActionPolyCoeffs, ActionPolyExponentVectors, ActionPolyMonomials, ActionPolyTerms}) = length(x.poly)

Base.eltype(::Type{ActionPolyCoeffs{PolyT}}) where {PolyT<:ActionPolyRingElem} = elem_type(coefficient_ring_type(PolyT))
Base.eltype(::Type{ActionPolyExponentVectors{PolyT}}) where {PolyT<:ActionPolyRingElem} = Vector{Int}
Base.eltype(::Type{ActionPolyMonomials{PolyT}}) where {PolyT<:ActionPolyRingElem} = PolyT
Base.eltype(::Type{ActionPolyTerms{PolyT}}) where {PolyT<:ActionPolyRingElem} = PolyT

###############################################################################
#
#  Promotion rules
#
###############################################################################

function AbstractAlgebra.promote_rule(::Type{PolyT}, ::Type{V}) where
        {T<:RingElement, V<:RingElement, PolyT<:ActionPolyRingElem{T}}
    AbstractAlgebra.promote_rule(T, V) == T ? PolyT : Union{}
end

AbstractAlgebra.promote_rule(::Type{PolyT}, ::Type{PolyT}) where
        {T<:RingElement, PolyT<:ActionPolyRingElem{T}} = PolyT

###############################################################################
#
#  Random generation
#
###############################################################################

rand(rng::AbstractRNG, apr::ActionPolyRing, term_range::AbstractUnitRange{Int},
     exp_bound::AbstractUnitRange{Int}, v...) = apr(rand(rng, base_ring(apr), term_range, exp_bound, v...))

rand(apr::ActionPolyRing, term_range, exp_bound, v...) = rand(Random.default_rng(), apr, term_range, exp_bound, v...)

###############################################################################
#
#  Conformance test element generation
#
###############################################################################

ConformanceTests.generate_element(R::ActionPolyRing{ZZRingElem}) = rand(R, 0:4, 0:10, -10:10)
ConformanceTests.generate_element(R::ActionPolyRing{ZZModRingElem}) = rand(R, 0:4, 0:10, -10:10)

###############################################################################
#
#  Misc 
#
###############################################################################

symbols(apr::ActionPolyRing) = symbols(base_ring(apr))[__perm_for_sort(apr)]

canonical_unit(apre::ActionPolyRingElem) = canonical_unit(data(apre))

###############################################################################
#
#  Aux action polynomial rings 
#
###############################################################################

# =========================================
# Getters for ring elements and rings
# =========================================

# Element getters
data(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem}) = dpre.p

__is_perm_up_to_date(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem}) = dpre.is_perm_up_to_date

# Ring getters
base_ring(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.upoly_ring
__jtv(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.jet_to_var::Dict{Tuple{Int, Vector{Int}}, elem_type(dpr)}
__vtj(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.var_to_jet::Dict{elem_type(dpr), Tuple{Int, Vector{Int}}}
__jtu_idx(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.jet_to_upoly_idx
__are_perms_up_to_date(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) = dpr.are_perms_up_to_date

# =========================================
# Permutation getters
# =========================================

function __perm_for_sort(dpr::Union{DifferencePolyRing, DifferentialPolyRing})
  if !__are_perms_up_to_date(dpr) || !isdefined(dpr, :permutation)
    __set_perm_for_sort!(dpr)
  end
  return dpr.permutation
end

function __perm_for_sort_poly(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem})
  dpr = parent(dpre)
  if !__are_perms_up_to_date(dpr) || !__is_perm_up_to_date(dpre)
    __set_perm_for_sort!(dpr)
    __set_perm_for_sort_poly!(dpre)
  end
  return dpre.permutation
end

# =========================================
# Setters for internal state
# =========================================

function __set_are_perms_up_to_date!(dpr::Union{DifferencePolyRing, DifferentialPolyRing}, update::Bool)
    dpr.are_perms_up_to_date = update
end

function __set_is_perm_up_to_date!(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem}, update::Bool)
    dpre.is_perm_up_to_date = update
end

function __set_perm_for_sort!(dpr::Union{DifferencePolyRing, DifferentialPolyRing})
    dpr.permutation = sortperm(dpr.(gens(base_ring(dpr))); rev = true)
    __set_are_perms_up_to_date!(dpr, true)
end

# Assumes are_perms_up_to_date == true
function __set_perm_for_sort_poly!(dpre::Union{DifferencePolyRingElem, DifferentialPolyRingElem})
  exps = collect(exponents(data(dpre)))
  n = length(exps)
  if n <= 1
    dpre.permutation = collect(1:n)
    return __set_is_perm_up_to_date!(dpre, true)
  end

  perm = (parent(dpre).permutation)[1:min(end, length(exps[1]))] #trim unused indices (avoids padding with zeros)

  dpre.permutation = sortperm(exps; lt=__my_lt_for_vec(perm), rev=true)
  __set_is_perm_up_to_date!(dpre, true)
end

function __my_lt_for_vec(perm::Vector{Int})
  return function ___my_lt_for_vec(ei::Vector{Int}, ej::Vector{Int})
    @inbounds for k in perm
      vi, vj = ei[k], ej[k]
      vi != vj && return vi < vj 
    end
    return false
  end
end

#Check if the jet_to_var dictionary of apr could contain the key (i,jet).
__is_valid_jet(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = i in 1:n_elementary_symbols(apr) && length(jet) == n_action_maps(apr) && all(>=(0), jet)

function __add_new_jetvar!(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})
  __set_are_perms_up_to_date!(apr, false)
  s_vec = map(jet_idx -> string(elementary_symbols(apr)[jet_idx[1]]) * "[" * join(jet_idx[2], ",") * "]", jet_idxs)::Vector{String}
  upr = base_ring(apr)
  ng = ngens(upr)
  new_vars = apr.(gens(upr, s_vec))
  for k in 1:length(new_vars)
    var = new_vars[k]
    jetvaridx = jet_idxs[k]
    __jtv(apr)[jetvaridx] = var
    __vtj(apr)[var] = jetvaridx
    __jtu_idx(apr)[jetvaridx] = ng + k
  end
  return new_vars
end

__add_new_jetvar!(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = __add_new_jetvar!(apr, [(i, jet)])

###############################################################################
#
#  Construction / Getter 
#
###############################################################################

@doc raw"""
    ranking(dpr::Union{DifferencePolyRing, DifferentialPolyRing}) -> ActionPolyRingRanking

Return the ranking of the jet variables of the difference or differential polynomial ring `dpr`.
"""
ranking(dpr::DifferencePolyRing{T}) where {T} = dpr.ranking::ActionPolyRingRanking{DifferencePolyRing{T}}
ranking(dpr::DifferentialPolyRing{T}) where {T} = dpr.ranking::ActionPolyRingRanking{DifferentialPolyRing{T}}

@doc raw"""
    set_ranking!(A::ActionPolyRing;
                 partition_name::Symbol = :default,
                 index_ordering_name::Symbol = :default,
                 partition::Vector{Vector{Int}} = Vector{Int}[],
                 index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)) 

This method configures the ranking of the action polynomial ring `A`, using an ordered partition of the elementary symbols and a monomial ordering on the indices. The ranking can be specified either by choosing predefined naming options or by explicitly providing a custom configuration.

# Keyword Arguments
- `partition_name`: Determines the partition of the elementary symbols of `dpr`. Supported values are:
  - `:top`: groups all elementary symbols into a single block,
  - `:pot`: separates each elementary symbol into its own block,
  - `:default`: uses `:top` unless a custom partition is specified.

- `index_ordering_name`: Specifies the ordering on the multiindices. Supported values are:
  - `:lex`: lexicographic ordering,
  - `:deglex`: degree lexicographic ordering,
  - `:invlex`: inverse lexicographic ordering,
  - `:deginvlex`: degree inverse lexicographic ordering,
  - `:degrevlex`: degree reverse lexicographic ordering,
  - `:default`: uses `:lex` unless a custom matrix is specified.

- `partition`: A custom partition of the elementary symbols, represented as a vector of characteristic vectors. The elementary symbols corresponding to the first characteristic vectors are considered largest and so on.

- `index_ordering_matrix`: A custom matrix representing a monomial ordering on the indices. Its number of columns must equal `n_action_maps(A)`.

# Examples

```jldoctests
julia> dpr = differential_polynomial_ring(ZZ, [:a, :b, :c], 4; partition_name=:pot, index_ordering_name = :degrevlex)[1]; ranking(dpr)
Ranking of differential polynomial ring in 3 elementary symbols over ZZ
with elementary symbols partitioned by
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
and ordering of the indices defined by
  [1    1    1    1]
  [0    0    0   -1]
  [0    0   -1    0]
  [0   -1    0    0]

julia> set_ranking!(dpr; partition = [[0,1,1],[1,0,0]], index_ordering_matrix = identity_matrix(ZZ, 4))
Ranking of differential polynomial ring in 3 elementary symbols over ZZ
with elementary symbols partitioned by
  [[0, 1, 1], [1, 0, 0]]
and ordering of the indices defined by
  [1   0   0   0]
  [0   1   0   0]
  [0   0   1   0]
  [0   0   0   1]

julia> set_ranking!(dpr)
Ranking of differential polynomial ring in 3 elementary symbols over ZZ
with elementary symbols partitioned by
  [[1, 1, 1]]
and ordering of the indices defined by
  [1   0   0   0]
  [0   1   0   0]
  [0   0   1   0]
  [0   0   0   1]
```
"""
function set_ranking!(dpr::PolyT;
    partition_name::Symbol = :default,
    index_ordering_name::Symbol = :default,
    partition::Vector{Vector{Int}} = Vector{Int}[],
    index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)
  ) where {PolyT <: ActionPolyRing}
  
  __set_are_perms_up_to_date!(dpr, false)
  m, n = n_elementary_symbols(dpr), n_action_maps(dpr)
  dpr.ranking = ActionPolyRingRanking{PolyT}(dpr, __compute_ranking_params(m, n, partition_name, index_ordering_name, partition, index_ordering_matrix)...)
  return ranking(dpr)
end

function __compute_ranking_params(m::Int, n::Int,
    partition_name,
    index_ordering_name,
    partition,
    index_ordering_matrix
  )
  
  @req partition_name in [:top, :pot, :default] "Invalid name of partition"
  if partition_name == :default
    if is_empty(partition)
      partition = [fill(1, m)] #Use :top by default
    else
      # Otherwise the input is used. Check its validity:
      @req __is_valid_partition(partition, m) "Not a partition of the number of elementary symbols"
    end
  else
    if is_empty(partition)
      if partition_name == :top
        partition = [fill(1, m)]
      else
        partition = [[i == j ? 1 : 0 for j in 1:m] for i in 1:m]
      end
    else # This case is only accessed if both a partition and a name are provided. Then a consistency check is required.
      @req __is_valid_partition(partition, m) "Not a partition of the number of elementary symbols"
      if partition_name == :top
        @req partition == [fill(1, m)] "The partition provided does not match its name"
      else
        @req partition == [[i == j ? 1 : 0 for j in 1:m] for i in 1:m] "The partition provided does not match its name"
      end
    end
  end

  R, _ = polynomial_ring(QQ, n; cached = false) # We use a dummy polynomial ring to extract canonical matrices of monomial orderings
  if index_ordering_name == :default
    if is_empty(index_ordering_matrix)
      index_ordering_matrix = canonical_matrix(lex(R)) #Use :lex by default
    else
      # Otherwise the input is used. Check its validity:
      @req ncols(index_ordering_matrix) == n "The number of columns of the matrix provided must equal $n" 
    end
  else
    @req is_empty(index_ordering_matrix) "Providing both a name and a matrix is not supported. Please just choose one."
    index_ordering_matrix = canonical_matrix(monomial_ordering(R, index_ordering_name))
  end
  return (partition, index_ordering_matrix)
end

###############################################################################
#
#  Field access and related stuff (rankings)
#
###############################################################################

### Difference ###

@doc raw"""
    parent(r::ActionPolyRingRanking)

Return the action polynomial ring `A` with `r = ranking(A)`.
"""
parent(ran::ActionPolyRingRanking) = ran.ring

@doc raw"""
    partition(r::ActionPolyRingRanking) -> Vector{Vector{Int}}

Return the partition of the elementary symbols defined by the ranking `r`
of the action polynomial ring `A`, where `r = ranking(A)`.
"""
partition(ran::ActionPolyRingRanking) = ran.partition

@doc raw"""
    index_ordering_matrix(r::ActionPolyRingRanking) -> ZZMatrix

Return the matrix inducing the monomial ordering of the multiindices defined
by the ranking `r` of the action polynomial ring `A`, where `r = ranking(A)`.
"""
index_ordering_matrix(ran::ActionPolyRingRanking) = ran.index_ordering_matrix

@doc raw"""
    riquier_matrix(r::ActionPolyRingRanking) -> ZZMatrix

Return a Riquier matrix that induces the ranking `r` of the action polynomial ring `A`, where `r = ranking(A)`.
"""
function riquier_matrix(ran::ActionPolyRingRanking)
  if !isdefined(ran, :riquier_matrix)
    par = partition(ran)
    apr = parent(ran)
    upper_part = block_diagonal_matrix([matrix(ZZ, length(par)-1, n_elementary_symbols(apr), vcat(par[1:end-1]...)), index_ordering_matrix(ran)])
    lower_part = block_diagonal_matrix([__in_block_tie_breaking_matrix(par), zero_matrix(ZZ, 0, n_action_maps(apr))])
    ran.riquier_matrix = vcat(upper_part, lower_part)
  end
  return ran.riquier_matrix
end

###############################################################################
#
#  Aux rankings 
#
###############################################################################

__is_valid_partition(par::Vector{Vector{Int}}, m::Int) = all(par_elt -> !is_zero(par_elt), par) && sum(par) == fill(1, m) && all(par_elt -> all(j -> par_elt[j] == 0 || par_elt[j] == 1, 1:m), par)

function __in_block_tie_breaking_matrix(partition::Vector{Vector{Int}})
  m = length(partition[1])
  max_ones = maximum(count(isone, vec) for vec in partition)
  result = zero_matrix(ZZ, max_ones - 1, m)
  for vec in partition
    pos = 1
    ones_seen = 0
    total_ones = count(isone, vec)
    for j in 1:m
      if isone(vec[j])
        ones_seen += 1
        if ones_seen == total_ones
          break
        end
        result[pos, j] = ZZ(1)
        pos += 1
      end
    end
  end
  return result
end

