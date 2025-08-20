###############################################################################
#
#  Construction 
#
###############################################################################

##### Algebras #####

### Difference ###
@doc raw"""
    difference_polynomial_ring(R::Ring, n_elementary_symbols::Int, ndiffs::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Return a tuple consisting of the difference polynomial ring over the ring `R` with the specified number of elementary variables and commuting endomorphisms, and the vector of
these elementary variables. This methods allows for all the keywords of the `set_ranking!` method. See there for further information.

# Examples

```jldoctests
julia> S, variables = difference_polynomial_ring(QQ, 3, 4)
(Difference polynomial ring in 3 elementary symbols over QQ, DifferencePolyRingElem{QQFieldElem}[u1[0,0,0,0], u2[0,0,0,0], u3[0,0,0,0]])

julia> S
Difference polynomial ring in 3 elementary symbols u1, u2, u3
with 4 commuting endomorphisms
  over Rational field

julia> variables
3-element Vector{DifferencePolyRingElem{QQFieldElem}}:
 u1[0,0,0,0]
 u2[0,0,0,0]
 u3[0,0,0,0]
```
"""
function difference_polynomial_ring(R::Ring, n_elementary_symbols::Int, ndiffs::Int; kwargs...)
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, n_elementary_symbols, ndiffs)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, map(i -> (i, fill(0, ndiffs)), 1:n_elementary_symbols))))
end

@doc raw"""
    difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) -> Tuple{DifferencePolyRing, Vector{DifferencePolyRingElem}}

Return a tuple consisting of the difference polynomial ring over the ring `R` with the specified elementary variables and number of commuting endomorphisms, and the vector of
these elementary variables. Note that the multiindex [0..0] of length 'ndiffs' is appended to the variable names provided. This methods allows for all the keywords of the `set_ranking!`
method. See there for further information.

# Examples

```jldoctest
julia> S, variables = difference_polynomial_ring(QQ, [:a, :b, :c], 4)
(Difference polynomial ring in 3 elementary symbols over QQ, DifferencePolyRingElem{QQFieldElem}[a[0,0,0,0], b[0,0,0,0], c[0,0,0,0]])

julia> S
Difference polynomial ring in 3 elementary symbols a, b, c
with 4 commuting endomorphisms
  over Rational field

julia> variables
3-element Vector{DifferencePolyRingElem{QQFieldElem}}:
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]
```
"""
function difference_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int; kwargs...)
  dpr = DifferencePolyRing{elem_type(typeof(R))}(R, elementary_symbols, ndiffs)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, map(i -> (i, fill(0, ndiffs)), 1:length(elementary_symbols)))))
end

### Differential ###
@doc raw"""
    differential_polynomial_ring(R::Ring, n_elementary_symbols::Int, ndiffs::Int) -> Tuple{DifferentialPolyRing, Vector{DifferentialPolyRingElem}}

Return a tuple consisting of the differential polynomial ring over the ring `R` with the specified number of elementary variables and commuting endomorphisms, and the vector of
these elementary variables. This methods allows for all the keywords of the `set_ranking!` method. See there for further information.

# Examples

```jldoctests
julia> S, variables = differential_polynomial_ring(QQ, [:a, :b, :c], 4)
(Differential polynomial ring in 3 elementary symbols over QQ, DifferentialPolyRingElem{QQFieldElem}[a[0,0,0,0], b[0,0,0,0], c[0,0,0,0]])

julia> S
Differential polynomial ring in 3 elementary symbols a, b, c
with 4 commuting derivations
  over Rational field

julia> variables
3-element Vector{DifferentialPolyRingElem{QQFieldElem}}:
 a[0,0,0,0]
 b[0,0,0,0]
 c[0,0,0,0]
```
"""
function differential_polynomial_ring(R::Ring, n_elementary_symbols::Int, ndiffs::Int; kwargs...)
  dpr = DifferentialPolyRing{elem_type(typeof(R))}(R, n_elementary_symbols, ndiffs)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, map(i -> (i, fill(0, ndiffs)), 1:n_elementary_symbols))))
end

@doc raw"""
    differential_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) -> Tuple{DifferentialPolyRing, Vector{DifferentialPolyRingElem}}

Return a tuple consisting of the differential polynomial ring over the ring `R` with the specified elementary variables and number of commuting endomorphisms, and the vector of
these elementary variables. Note that the multiindex [0..0] of length 'ndiffs' is appended to the variable names provided. This methods allows for all the keywords of the `set_ranking!`
method. See there for further information.

# Examples

```jldoctests
julia> S, variables = differential_polynomial_ring(QQ, 3, 4)
(Differential polynomial ring in 3 elementary symbols over QQ, DifferentialPolyRingElem{QQFieldElem}[u1[0,0,0,0], u2[0,0,0,0], u3[0,0,0,0]])

julia> S
Differential polynomial ring in 3 elementary symbols u1, u2, u3
with 4 commuting derivations
  over Rational field

julia> variables
3-element Vector{DifferentialPolyRingElem{QQFieldElem}}:
 u1[0,0,0,0]
 u2[0,0,0,0]
 u3[0,0,0,0]
```
"""
function differential_polynomial_ring(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int; kwargs...)
  dpr = DifferentialPolyRing{elem_type(typeof(R))}(R, elementary_symbols, ndiffs)
  set_ranking!(dpr; kwargs...)
  return (dpr, deepcopy.(__add_new_jetvar!(dpr, map(i -> (i, fill(0, ndiffs)), 1:length(elementary_symbols)))))
end

##### Elements #####

### Difference ###
(dpr::DifferencePolyRing)() = DifferencePolyRingElem{elem_type(base_ring(dpr))}(dpr)

(dpr::DifferencePolyRing)(upre::AbstractAlgebra.Generic.UniversalPolyRingElem) = DifferencePolyRingElem{elem_type(base_ring(dpr))}(dpr, upre)

function (dpr::DifferencePolyRing{T})(a::DifferencePolyRingElem{T}) where {T}
  @req parent(a) === dpr "Wrong parent"
  return a
end

### Differential ###
(dpr::DifferentialPolyRing)() = DifferentialPolyRingElem{elem_type(base_ring(dpr))}(dpr)

(dpr::DifferentialPolyRing)(upre::AbstractAlgebra.Generic.UniversalPolyRingElem) = DifferentialPolyRingElem{elem_type(base_ring(dpr))}(dpr, upre)

function (dpr::DifferentialPolyRing{T})(a::DifferentialPolyRingElem{T}) where {T}
  @req parent(a) === dpr "Wrong parent"
  return a
end

### generic ###
(apr::ActionPolyRing)(a::T) where {T<:RingElement} = apr(__upr(apr)(a))

###############################################################################
#
#  Basic ring functionality 
#
###############################################################################

### Difference ###
elem_type(::Type{DifferencePolyRing{T}}) where {T} = DifferencePolyRingElem{T}

parent_type(::Type{DifferencePolyRingElem{T}}) where {T} = DifferencePolyRing{T}

function Base.deepcopy_internal(dpre::DifferencePolyRingElem{T}, dict::IdDict) where {T}
    # Avoid deepcopying the parent as it may refer back to it in one of its dictionaries 
    pp = deepcopy_internal(data(dpre), dict)
    return DifferencePolyRingElem{T}(parent(dpre), pp)
end

zero(dpr::DifferencePolyRing) = dpr()

one(dpr::DifferencePolyRing) = dpr(one(__upr(dpr)))

### Differential ###
elem_type(::Type{DifferentialPolyRing{T}}) where {T} = DifferentialPolyRingElem{T}

parent_type(::Type{DifferentialPolyRingElem{T}}) where {T} = DifferentialPolyRing{T}

function Base.deepcopy_internal(dpre::DifferentialPolyRingElem{T}, dict::IdDict) where {T}
    # Avoid deepcopying the parent as it may refer back to it in one of its dictionaries 
    pp = deepcopy_internal(data(dpre), dict)
    return DifferentialPolyRingElem{T}(parent(dpre), pp)
end

zero(dpr::DifferentialPolyRing) = dpr()

one(dpr::DifferentialPolyRing) = dpr(one(__upr(dpr)))

### generic ###
base_ring_type(::Type{<:ActionPolyRing{T}}) where {T} = parent_type(T)

is_square(apre::ActionPolyRingElem) = is_square(data(apre))

Base.sqrt(apre::ActionPolyRingElem; check::Bool = true) = parent(apre)(sqrt(data(apre); check = check))

is_irreducible(apre::ActionPolyRingElem) = is_irreducible(data(apre))

is_unit(apre::ActionPolyRingElem) = is_unit(data(apre))

characteristic(apr::ActionPolyRing) = characteristic(__upr(apr))

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

##### Algebras #####

### Difference ###
ndiffs(dpr::DifferencePolyRing) = dpr.ndiffs

elementary_symbols(dpr::DifferencePolyRing) = dpr.elementary_symbols

### Differential ###
ndiffs(dpr::DifferentialPolyRing) = dpr.ndiffs

elementary_symbols(dpr::DifferentialPolyRing) = dpr.elementary_symbols

##### Elements #####

parent(dpre::DifferencePolyRingElem) = dpre.parent

parent(dpre::DifferentialPolyRingElem) = dpre.parent

##### Related stuff #####
base_ring(apr::ActionPolyRing) = base_ring(__upr(apr))

n_elementary_symbols(apr::ActionPolyRing) = length(elementary_symbols(apr))

###############################################################################
#
#  Basic polynomial functionality 
#
###############################################################################

coeff(apre::ActionPolyRingElem, i::Int) = coeff(data(apre), __perm_for_sort_poly(apre)[i])

exponent_vector(apre::ActionPolyRingElem, i::Int) = exponent_vector(data(apre), __perm_for_sort_poly(apre)[i])[__perm_for_sort(parent(apre))]

monomial(apre::ActionPolyRingElem, i::Int) = parent(apre)(monomial(data(apre), __perm_for_sort_poly(apre)[i]))

term(apre::ActionPolyRingElem, i::Int) = parent(apre)(term(data(apre), __perm_for_sort_poly(apre)[i]))

is_monomial(apre::ActionPolyRingElem) = is_monomial(data(apre))

is_term(apre::ActionPolyRingElem) = is_monomial(data(apre))

is_univariate(apre::ActionPolyRingElem) = is_univariate(data(apre))

is_univariate(apr::ActionPolyRing) = false

length(apre::ActionPolyRingElem) = length(data(apre))

@doc raw"""
    total_degree(p::ActionPolyRingElem) -> Int

Return the total degree of `p`.
"""
total_degree(apre::ActionPolyRingElem) = total_degree(data(apre))

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int, jet::Vector{Int}) -> Int

Return the degree of the polynomial `p` in the `i`-th elementary variable with
multiindex `jet`. If this jet variable is valid but still untracked, return $0$. Alternatively, the variable may be passed right away instead of its index.
"""
function degree(apre::ActionPolyRingElem, i::Int, jet::Vector{Int})
  apr = parent(apre)
  upr = __upr(apr)
  @req __is_valid_jet(apr, i, jet) "Invalid jet variable"
  jtv = __jtv(apr)
  if haskey(jtv, (i,jet))
    idx = findfirst(var -> var == data(jtv[(i,jet)]), gens(upr))
    return degree(data(apre), gen(upr, idx))
  end
  return 0
end

function degree(apre::ActionPolyRingElem{T}, var::ActionPolyRingElem{T}) where {T}
  check_parent(apre, var)
  @req is_gen(var) "Not a variable"
  return degree(apre, __vtj(parent(apre))[var])
end

degree(apre::ActionPolyRingElem, jet_idx::Tuple{Int, Vector{Int}}) = degree(apre, jet_idx...)

@doc raw"""
    degree(p::ActionPolyRingElem, i::Int)

Return the degree of the polynomial `p` in the, among the currently tracked variables, 'i'-th largest one. The index of the jet variable may also be passed as a tuple. Alternatively, the variable may be passed right away instead of its index.
"""
degree(apre::ActionPolyRingElem, i::Int) = degree(apre, __vtj(parent(apre))[gen(parent(apre), i)])

@doc raw"""
    degrees(p::ActionPolyRingElem)

Return an array of the degrees of the polynomial `p` in terms of each variable. The variables are sorted with respect to the ranking of the action polynomial ring containing `p`, leading with the largest variable.
"""
degrees(apre::ActionPolyRingElem) = map(i -> degree(apre, i), 1:ngens(parent(apre)))

@doc raw"""
    is_constant(p::ActionPolyRingElem)

Return `true` if `p` is a degree zero polynomial or the zero polynomial, i.e. a constant polynomial. 
"""
is_constant(apre::ActionPolyRingElem) = is_constant(data(apre))

@doc raw"""
    vars(p::ActionPolyRingElem)

Return the variables actually occuring in `p`. The variables are sorted with respect to the ranking of the action polynomial ring containing `p`, leading with the largest variable.
"""
vars(apre::ActionPolyRingElem) = sort(parent(apre).(vars(data(apre))); rev = true)

is_gen(apre::ActionPolyRingElem) = is_gen(data(apre))

@doc raw"""
    gen(apr::ActionPolyRing, i::Int, jet::Vector{Int})

Return the `i`-th elementary variable with multiindex `jet` in the action polynomial
`apr`. If this jet variable was untracked, it is tracked afterwards. The index of the jet variable may also be passed as a tuple.

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
function gen(apr::ActionPolyRing, i::Int, jet::Vector)
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
    gen(apr::ActionPolyRing, i::Int)

Return the, among the currently tracked variables of `apr`, 'i'-th largest one.

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
    gens(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})

Return the jet variables of the action polynomial ring `apr` specified by the entries of `jet_idxs` as
a vector and track all new variables.

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
function gens(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})
  jtv = __jtv(apr)
  new_jet_idxs = filter(jet_idx -> !haskey(jtv, jet_idx), jet_idxs)
  if !is_empty(new_jet_idxs)
    __add_new_jetvar!(apr, new_jet_idxs)
  end
  return map(jet_idx -> jtv[jet_idx], jet_idxs)
end

@doc raw"""
    gens(apr::ActionPolyRing)

Return the currently tracked variables of the action polynomial ring `apr` as a vector. The variables are sorted with respect to the ranking of `apr`, leading with the largest variable.

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
function gens(apr::ActionPolyRing)
  return sort(collect(deepcopy(values(__jtv(apr)))); rev = true)
end

number_of_generators(apr::ActionPolyRing) = number_of_generators(__upr(apr))

number_of_variables(apr::ActionPolyRing) = number_of_variables(__upr(apr))

@doc raw"""
    getindex(apr::ActionPolyRing, i::Int, jet::Vector{Int})

Alias for `gen(apr, i, jet)`.
"""
getindex(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = gen(apr, i, jet)

getindex(apr::ActionPolyRing, jet_idx::Tuple{Int, Vector{Int}}) = gen(apr, jet_idx...)

@doc raw"""
    var_index(x::ActionPolyRingElem)

Return the integer `i` such that `x` is the `i`-th largest currently tracked variable. If `x` is not a variable an exception is raised. 
"""
function var_index(x::ActionPolyRingElem)
  @req is_gen(x) "Not a variable in var_index"
  return findfirst(==(x), gens(parent(x)))
end

constant_coefficient(apre::ActionPolyRingElem) = constant_coefficient(data(apre))

@doc raw"""
    leading_coefficient(p::ActionPolyRingElem{T}) -> T

Return the leading coefficient of the polynomial `p`, i.e. the coefficient of the first (with respect to the ranking of the action polynomial ring containing it) nonzero term, or zero if the polynomial is zero.
"""
function leading_coefficient(apre::ActionPolyRingElem{T}) where {T}
  if length(apre) == 0
    return zero(T)
  end
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
  if len == 0
    return zero(T)
  end
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

Return the tail of `p` with respect to the ranking of the action polynomial ring containing it.
"""
tail(apre::ActionPolyRingElem) = apre - leading_term(apre)

###############################################################################
#
#  Action polynomial functionality 
#
###############################################################################

@doc raw"""
    diff_action(p::ActionPolyRingElem, i::Int)

Apply the `i`-th diff-action to the polynomial `p`. If `p` is a difference polynomial then this action is multiplicative on the variables. If `p` is a differential polynomial then this action adheres to the Leibniz rule.
"""
function diff_action(dpre::DifferencePolyRingElem{T}, i::Int) where {T}
  d = fill(0, ndiffs(parent(dpre)))
  d[i] = 1
  return diff_action(dpre, d)
end

function diff_action(dpre::DifferentialPolyRingElem{T}, i::Int) where {T}
  dpr = parent(dpre)
  @req i in 1:ndiffs(dpr) "index out of range"

  # Remove constant term: d_i(1) = 0
  dpre -= constant_coefficient(dpre)

  if is_zero(dpre)
    return dpre
  end
 
  d = fill(0, ndiffs(dpr))
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

  upr = __upr(dpr)
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

Successively apply the `i`-th diff-action `d[i]`-times to the polynomial `p`, where $i = 1, \ldots, length(d)$. 
"""
function diff_action(dpre::DifferencePolyRingElem{T}, d::Vector{Int}) where {T}
  dpr = parent(dpre)
  @req length(d) == ndiffs(dpr) && all(j -> j >= 0, d) "Invalid vector of multiplicities"
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
  
  upr = __upr(dpr)
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
  @req len == ndiffs(parent(dpre)) && all(j -> j >= 0, d) "Invalid vector of diff multiplicities"
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
  return divexact(leading_term(apre), leader(apre)^degree(apre, leader(apre)); check = true)
end

@doc raw"""
    leader(p::ActionPolyRingElem)

Return the leader of the polynomial `p`, that is the largest variable with respect to the ranking of `parent(p)`. If `p` is constant, an error is raised.
"""
function leader(apre::ActionPolyRingElem)
  @req !is_constant(apre) "A constant polynomial has no leader"
  if is_univariate(apre)
    return vars(apre)[1]
  end
  return max(vars(apre)...)
end
###############################################################################
#
#   Iterators
#
###############################################################################

function coefficients(a::DifferencePolyRingElem)
   return DifferencePolyCoeffs(a)
end

function coefficients(a::DifferentialPolyRingElem)
   return DifferentialPolyCoeffs(a)
end

function AbstractAlgebra.exponent_vectors(a::DifferencePolyRingElem)
   return DifferencePolyExponentVectors(a)
end

function AbstractAlgebra.exponent_vectors(a::DifferentialPolyRingElem)
   return DifferentialPolyExponentVectors(a)
end

function monomials(a::DifferencePolyRingElem)
   return DifferencePolyMonomials(a)
end

function monomials(a::DifferentialPolyRingElem)
   return DifferentialPolyMonomials(a)
end

function terms(a::DifferencePolyRingElem)
   return DifferencePolyTerms(a)
end

function terms(a::DifferentialPolyRingElem)
   return DifferentialPolyTerms(a)
end

function Base.iterate(x::Union{DifferencePolyCoeffs, DifferentialPolyCoeffs})
   if length(x.poly) >= 1
      return coeff(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyCoeffs, DifferentialPolyCoeffs}, state)
   state += 1
   if length(x.poly) >= state
      return coeff(x.poly, state), state
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyExponentVectors, DifferentialPolyExponentVectors})
   if length(x.poly) >= 1
      return exponent_vector(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyExponentVectors, DifferentialPolyExponentVectors}, state)
   state += 1
   if length(x.poly) >= state
      return exponent_vector(x.poly, state), state
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyTerms, DifferentialPolyTerms})
   if length(x.poly) >= 1
      return term(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyTerms, DifferentialPolyTerms}, state)
   state += 1
   if length(x.poly) >= state
      return term(x.poly, state), state
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyMonomials, DifferentialPolyMonomials})
   if length(x.poly) >= 1
      return monomial(x.poly, 1), 1
   else
      return nothing
   end
end

function Base.iterate(x::Union{DifferencePolyMonomials, DifferentialPolyMonomials}, state)
   state += 1
   if length(x.poly) >= state
      return monomial(x.poly, state), state
   else
      return nothing
   end
end

function Base.length(x::Union{DifferencePolyCoeffs, DifferencePolyExponentVectors, DifferencePolyTerms, DifferencePolyMonomials, DifferentialPolyCoeffs, DifferentialPolyExponentVectors, DifferentialPolyTerms, DifferentialPolyMonomials})
   return length(x.poly)
end

function Base.eltype(::Type{DifferencePolyCoeffs{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return S
end

function Base.eltype(::Type{DifferencePolyExponentVectors{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return Vector{Int}
end

function Base.eltype(::Type{DifferencePolyMonomials{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return T
end

function Base.eltype(::Type{DifferencePolyTerms{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return T
end

function Base.eltype(::Type{DifferentialPolyCoeffs{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return S
end

function Base.eltype(::Type{DifferentialPolyExponentVectors{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return Vector{Int}
end

function Base.eltype(::Type{DifferentialPolyMonomials{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return T
end

function Base.eltype(::Type{DifferentialPolyTerms{T}}) where T <: ActionPolyRingElem{S} where S <: RingElement
   return T
end

###############################################################################
#
#   Promotion rules
#
###############################################################################

AbstractAlgebra.promote_rule(::Type{DifferencePolyRingElem{T}}, ::Type{DifferencePolyRingElem{T}}) where {T <: RingElement} = DifferencePolyRingElem{T}

function AbstractAlgebra.promote_rule(::Type{DifferencePolyRingElem{T}}, ::Type{V}) where {T <: RingElement, V <: RingElement}
   AbstractAlgebra.promote_rule(T, V) == T ? DifferencePolyRingElem{T} : Union{}
end

AbstractAlgebra.promote_rule(::Type{DifferentialPolyRingElem{T}}, ::Type{DifferentialPolyRingElem{T}}) where {T <: RingElement} = DifferentialPolyRingElem{T}

function AbstractAlgebra.promote_rule(::Type{DifferentialPolyRingElem{T}}, ::Type{V}) where {T <: RingElement, V <: RingElement}
   AbstractAlgebra.promote_rule(T, V) == T ? DifferentialPolyRingElem{T} : Union{}
end

###############################################################################
#
#   Random generation
#
###############################################################################

rand(rng::AbstractRNG, apr::ActionPolyRing, term_range::AbstractUnitRange{Int},
     exp_bound::AbstractUnitRange{Int}, v...) = apr(rand(rng, __upr(apr), term_range, exp_bound, v...))

rand(apr::ActionPolyRing, term_range, exp_bound, v...) = rand(Random.default_rng(), apr, term_range, exp_bound, v...)


###############################################################################
#
#   Conformance test element generation
#
###############################################################################

function ConformanceTests.generate_element(R::ActionPolyRing{ZZRingElem})
    return rand(R, 0:4, 0:10, -10:10)
end

function ConformanceTests.generate_element(R::ActionPolyRing{ZZModRingElem})
    return rand(R, 0:4, 0:10, -10:10)
end

###############################################################################
#
#  Misc 
#
###############################################################################

symbols(apr::ActionPolyRing) = symbols(__upr(apr))[__perm_for_sort(apr)]

canonical_unit(apre::ActionPolyRingElem) = canonical_unit(data(apre))

###############################################################################
#
#  Update polynomials 
#
###############################################################################

# The multivariate polynomial ring parent(data(data(dpre))) will always equal what has been parent(dpre).upoly_ring.mpoly_ring at the time of its creation, so
# we need a way update this, e.g. after adding variables. This is mandatory for some methods to work proper and also just good practise.

function __update_internals!(dpre::DifferencePolyRingElem)
  upr = __upr(parent(dpre))
  p = data(dpre)
  if upr.mpoly_ring !== parent(data(p))
    dpre.p = upr(p)
  end
end  

function __update_internals!(dpre::DifferentialPolyRingElem)
  upr = __upr(parent(dpre))
  p = data(dpre)
  if upr.mpoly_ring !== parent(data(p))
    dpre.p = upr(p)
  end
end  

###############################################################################
#
#  Aux action polynomial rings 
#
###############################################################################

### Difference ###

#Getters for internals
data(dpre::DifferencePolyRingElem) = dpre.p

__upr(dpr::DifferencePolyRing) = dpr.upoly_ring

__jtv(dpr::DifferencePolyRing) = dpr.jet_to_var

__vtj(dpr::DifferencePolyRing) = dpr.var_to_jet

__jtu_idx(dpr::DifferencePolyRing) = dpr.jet_to_upoly_idx

__are_perms_up_to_date(dpr::DifferencePolyRing) = dpr.are_perms_up_to_date

function __perm_for_sort(dpr::DifferencePolyRing)
  if __are_perms_up_to_date(dpr) && isdefined(dpr, :permutation)
    return dpr.permutation
  end
  __set_perm_for_sort!(dpr)
  return dpr.permutation
end

function __perm_for_sort_poly(dpre::DifferencePolyRingElem)
  dpr = parent(dpre)
  if __are_perms_up_to_date(dpr) && isdefined(dpre, :permutation)
    return dpre.permutation
  end
  __set_perm_for_sort!(dpr)
  __set_perm_for_sort_poly!(dpre)
  return dpre.permutation
end

#Setters for internals
function __set_are_perms_up_to_date!(dpr::DifferencePolyRing, update::Bool)
  dpr.are_perms_up_to_date = update
end

function __set_perm_for_sort!(dpr::DifferencePolyRing)
  dpr.permutation = sortperm(dpr.(gens(__upr(dpr))); rev = true)
  __set_are_perms_up_to_date!(dpr, true)
end

#If this is called, we assume that are_perms_up_to_date == true
function __set_perm_for_sort_poly!(dpre::DifferencePolyRingElem)
  __update_internals!(dpre)
  dpre.permutation = sortperm(map(expo -> expo[parent(dpre).permutation], collect(exponents(data(dpre)))); rev = true)
end

### Differential ###
data(dpre::DifferentialPolyRingElem) = dpre.p

__upr(dpr::DifferentialPolyRing) = dpr.upoly_ring

__jtv(dpr::DifferentialPolyRing) = dpr.jet_to_var

__vtj(dpr::DifferentialPolyRing) = dpr.var_to_jet

__jtu_idx(dpr::DifferentialPolyRing) = dpr.jet_to_upoly_idx

__are_perms_up_to_date(dpr::DifferentialPolyRing) = dpr.are_perms_up_to_date

function __perm_for_sort(dpr::DifferentialPolyRing)
  if __are_perms_up_to_date(dpr) && isdefined(dpr, :permutation)
    return dpr.permutation
  end
  __set_perm_for_sort!(dpr)
  return dpr.permutation
end

function __perm_for_sort_poly(dpre::DifferentialPolyRingElem)
  dpr = parent(dpre)
  if __are_perms_up_to_date(dpr) && isdefined(dpre, :permutation)
    return dpre.permutation
  end
  __set_perm_for_sort!(dpr)
  __set_perm_for_sort_poly!(dpre)
  return dpre.permutation
end

#Setters for internals
function __set_are_perms_up_to_date!(dpr::DifferentialPolyRing, update::Bool)
  dpr.are_perms_up_to_date = update
end

function __set_perm_for_sort!(dpr::DifferentialPolyRing)
  dpr.permutation = sortperm(dpr.(gens(__upr(dpr))); rev = true)
  __set_are_perms_up_to_date!(dpr, true)
end

#If this is called, we assume that are_perms_up_to_date == true
function __set_perm_for_sort_poly!(dpre::DifferentialPolyRingElem)
  __update_internals!(dpre)
  dpre.permutation = sortperm(map(expo -> expo[parent(dpre).permutation], collect(exponents(data(dpre)))); rev = true)
end

#Check if the jet_to_var dictionary of apr could contain the key (i,jet).
__is_valid_jet(apr::ActionPolyRing, i::Int, jet::Vector{Int}) = i in 1:n_elementary_symbols(apr) && length(jet) == ndiffs(apr) && all(j -> j >= 0, jet)

function __add_new_jetvar!(apr::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})
  __set_are_perms_up_to_date!(apr, false)
  s_vec = map(jet_idx -> string(elementary_symbols(apr)[jet_idx[1]]) * "[" * join(jet_idx[2], ",") * "]", jet_idxs)::Vector{String}
  upr = __upr(apr)
  ng = ngens(upr)
  new_vars = apr.(gens(upr, s_vec))
  for k in 1:length(new_vars)
    var = new_vars[k]
    jetvaridx = jet_idxs[k]
    __jtv(apr)[jetvaridx] = var
    __vtj(apr)[var] = jetvaridx
    __jtu_idx(apr)[jetvaridx] = ng + k
  end
  for key in keys(__vtj(apr))
    __update_internals!(key)
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
    ranking(dpr::DifferencePolyRing) -> ActionPolyRingRanking

Return the ranking of the jet variables of the difference polynomial ring `dpr`.
"""
ranking(dpr::DifferencePolyRing) = dpr.ranking
  
@doc raw"""
    ranking(dpr::DifferentialPolyRing) -> ActionPolyRingRanking

Return the ranking of the jet variables of the differential polynomial ring `dpr`
"""
ranking(dpr::DifferentialPolyRing) = dpr.ranking

@doc raw"""
    set_ranking!(apr::ActionPolyRing;
                 partition_name::Symbol = :default,
                 index_ordering_name::Symbol = :default,
                 partition::Vector{Vector{Int}} = Vector{Int}[],
                 index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)) 

This method configures the ranking of the action polynomial ring `apr`, using an ordered partition of the elementary symbols and a monomial ordering on the indices. The ranking can be specified either by choosing predefined naming options or by explicitly providing a custom configuration.

# Keyword Arguments
- `partition_name`: Determines the partition of the elementary symbols of `dpr`. Supported values are:
  - `:top`: groups all variables into a single block,
  - `:pot`: separates each variable into its own block,
  - `:default`: uses `:top` unless a custom partition is specified.

- `index_ordering_name`: Specifies the ordering on the multiindices. Supported values are:
  - `:lex`: lexicographic ordering,
  - `:deglex`: degree lexicographic ordering,
  - `:invlex`: inverse lexicographic ordering,
  - `:deginvlex`: degree inverse lexicographic ordering,
  - `:degrevlex`: degree reverse lexicographic ordering,
  - `:default`: uses `:lex` unless a custom matrix is specified.

- `partition`: A custom partition of the elementary symbols, represented as a vector of characteristic vectors. The elementary symbols corresponding to the first characteristic vectors are considered largest and so on.

- `index_ordering_matrix`: A custom matrix representing a monomial ordering on the indices. Its number of columns must equal `ndiffs(dpr)`.

# Examples

```jldoctests
julia> dpr = differential_polynomial_ring(ZZ, [:a, :b, :c], 4; partition_name=:pot, index_ordering_name = :degrevlex)[1]; ranking(dpr)
Ranking of Differential polynomial ring in 3 elementary symbols over ZZ
with elementary symbols partitioned by
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
and ordering of the indices defined by
  [1    1    1    1]
  [0    0    0   -1]
  [0    0   -1    0]
  [0   -1    0    0]

julia> set_ranking!(dpr; partition = [[0,1,1],[1,0,0]], index_ordering_matrix = identity_matrix(ZZ, 4))
Ranking of Differential polynomial ring in 3 elementary symbols over ZZ
with elementary symbols partitioned by
  [[0, 1, 1], [1, 0, 0]]
and ordering of the indices defined by
  [1   0   0   0]
  [0   1   0   0]
  [0   0   1   0]
  [0   0   0   1]

julia> set_ranking!(dpr)
Ranking of Differential polynomial ring in 3 elementary symbols over ZZ
with elementary symbols partitioned by
  [[1, 1, 1]]
and ordering of the indices defined by
  [1   0   0   0]
  [0   1   0   0]
  [0   0   1   0]
  [0   0   0   1]
```
"""
function set_ranking!(dpr::DifferencePolyRing;
    partition_name::Symbol = :default,
    index_ordering_name::Symbol = :default,
    partition::Vector{Vector{Int}} = Vector{Int}[],
    index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)
  )
  
  __set_are_perms_up_to_date!(dpr, false)
  m, n = n_elementary_symbols(dpr), ndiffs(dpr)
  dpr.ranking = ActionPolyRingRanking{typeof(dpr)}(dpr, __compute_ranking_params(m, n, partition_name, index_ordering_name, partition, index_ordering_matrix)...)
end

function set_ranking!(dpr::DifferentialPolyRing;
    partition_name::Symbol = :default,
    index_ordering_name::Symbol = :default,
    partition::Vector{Vector{Int}} = Vector{Int}[],
    index_ordering_matrix::ZZMatrix = zero_matrix(ZZ, 0, 0)
  )
  
  __set_are_perms_up_to_date!(dpr, false)
  m, n = n_elementary_symbols(dpr), ndiffs(dpr)
  dpr.ranking = ActionPolyRingRanking{typeof(dpr)}(dpr, __compute_ranking_params(m, n, partition_name, index_ordering_name, partition, index_ordering_matrix)...)
end

function __compute_ranking_params(m::Int, n::Int,
    partition_name,
    index_ordering_name,
    partition,
    index_ordering_matrix
  )::Tuple{Vector{Vector{Int64}}, ZZMatrix}
  
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

  @req index_ordering_name in [:default, :lex, :deglex, :invlex, :deginvlex, :degrevlex] "Invalid name of index ordering"
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
    #if is_empty(index_ordering_matrix)
      if index_ordering_name == :lex
        index_ordering_matrix = canonical_matrix(lex(R))
      elseif index_ordering_name == :deglex
        index_ordering_matrix = canonical_matrix(deglex(R))
      elseif index_ordering_name == :invlex
        index_ordering_matrix = canonical_matrix(invlex(R))
      elseif index_ordering_name == :deginvlex
        index_ordering_matrix = canonical_matrix(deginvlex(R))
      elseif index_ordering_name == :degrevlex
        index_ordering_matrix = canonical_matrix(degrevlex(R))
      end
    #end
  end
  return (partition, index_ordering_matrix)
end

###############################################################################
#
#  Field access and related stuff (rankings)
#
###############################################################################

### Difference ###
base_ring(ran::ActionPolyRingRanking) = ran.ring

@doc raw"""
    partition(r::ActionPolyRingRanking) -> Vector{Vector{Int}}

Return the partition of the elementary symbols defined by the ranking `r`
of the action polynomial ring `apr`, where `r = ranking(apr)`.
"""
partition(ran::ActionPolyRingRanking) = ran.partition

@doc raw"""
    index_ordering_matrix(r::ActionPolyRingRanking) -> ZZMatrix

Return the matrix inducing the monomial ordering of the multiindices defined
by the ranking `r` of the action polynomial ring `apr`, where `r = ranking(apr)`.
"""
index_ordering_matrix(ran::ActionPolyRingRanking) = ran.index_ordering_matrix

@doc raw"""
    riquier_matrix(r::ActionPolyRingRanking) -> ZZMatrix

Return a Riquier matrix that induces the ranking `r` of the action polynomial ring `apr`, where `r = ranking(apr)`.
"""
function riquier_matrix(ran::ActionPolyRingRanking)
  if !isdefined(ran, :riquier_matrix)
    par = partition(ran)
    dpr = base_ring(ran)
    upper_part = block_diagonal_matrix([matrix(ZZ, length(par)-1, n_elementary_symbols(dpr), vcat(par[1:end-1]...)), index_ordering_matrix(ran)])
    lower_part = block_diagonal_matrix([__in_block_tie_breaking_matrix(par), zero_matrix(ZZ, 0, ndiffs(dpr))])
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

