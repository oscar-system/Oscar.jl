###############################################################################
#
#  Pseudo-division and Ritt ordering
#
###############################################################################

function pseudorem(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  check_parent(p, var)
  return __core_pseudorem(p, q, var, false)[1]
end

@doc raw"""
    pseudorem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem} -> PolyT

Return the algebraic pseudo-remainder of `p` divided by `q` with respect to the jet variable specified by `i` and
`jet`. If no jet variable is specified then division is performed with respect to the leader of `q`, even allowing
`q` to be a nonzero constant. This method performs division by using a lazy pre-multiplication by the initial of `q`
at each step, only multiplying the remainder when necessary.

This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables); see the online documentation.

# Examples

```jldoctest pseudorem_example
julia> dpr, (x, y) = differential_polynomial_ring(QQ, [:x, :y], 1); p, q = (x^2 + y, y*x + 1)
(x[0]^2 + y[0], y[0]*x[0] + 1)

julia> pseudorem(p, q, x)
y[0]^3 + 1

julia> pseudorem(p, q, y)
x[0]^3 - 1
```
Note that in this example, the leader of `q` is `x`, so it need not be specified for pseudo-division with respect to `x`:
```jldoctest pseudorem_example
julia> leader(q) == x
true

julia> pseudorem(p, q)
y[0]^3 + 1
```
Finally, the pseudoremainder of any polynomial with respect to a non-zero constant is always zero.
```jldoctest pseudorem_example
julia> pseudorem(p, dpr(1))
0
```
"""
function pseudorem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem}
  deg_q = degree(q, i, jet)
  deg_q < 0 && throw(DivideError()) # By convention, (only) the zero polynomial has degree -1 in each jet variable
  @req deg_q > 0 "Cannot pseudo-divide by a polynomial with degree 0 in the specified division variable"
  
  # positive degree ensures existence of the key
  return pseudorem(p, q, __jtv(parent(q))[(i, jet)])
end

pseudorem(p::PolyT, q::PolyT, jet_idx::Tuple{Int, Vector{Int}}) where {PolyT <: ActionPolyRingElem} = pseudorem(p, q, jet_idx...)
pseudorem(p::PolyT, q::PolyT, i::Int) where {PolyT <: ActionPolyRingElem} = pseudorem(p, q, gen(parent(p), i))

function pseudorem(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  is_zero(q) && throw(DivideError())
  is_constant(q) && return zero(q)
  return pseudorem(p, q, leader(q))
end

### with div ###

function pseudodivrem(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  check_parent(p, var)
  return __core_pseudodivrem(p, q, var, false)[1:2]
end

@doc raw"""
    pseudodivrem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem} -> Tuple{PolyT, PolyT}

Return the pair `(s, r)` where `s` is the pseudo-quotient and `r` is the pseudo-remainder of `p`
by `q` with respect to the jet variable specified by `i` and `jet`. If no jet variable is specified
then division is performed with respect to the leader of `q`, even allowing `q` to be a nonzero constant.
The number of pre-multiplications by the leading coefficient of `q` in this jet variable is minimised,
i.e. we have `lc(q)^k * p = s * q + r` where the integer `k >= 0` is minimal and `lc(q)` is the above
mentioned leading coefficient.

This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables); see the online documentation.

# Examples

```jldoctest pseudodivrem_example
julia> dpr, (x, y) = differential_polynomial_ring(QQ, [:x, :y], 1); p, q = (x^2 + y, y*x + 1)
(x[0]^2 + y[0], y[0]*x[0] + 1)

julia> pseudodivrem(p, q, x)
(y[0]*x[0] - 1, y[0]^3 + 1)

julia> y^2 * p == (y*x - 1) * q + (y^3 + 1)
true

julia> pseudodivrem(p, q, y)
(1, x[0]^3 - 1)

julia> x * p == 1 * q + (x^3 - 1)
true
```
Note that in this example, the leader of `q` is `x`, so it need not be specified for pseudo-division with respect to `x`:
```jldoctest pseudodivrem_example
julia> leader(q) == x
true

julia> pseudodivrem(p, q)
(y[0]*x[0] - 1, y[0]^3 + 1)
```
Finally, if the second argument is a non-zero constant, then it depends on the first argument being divisible by said constant,
whether the pseudo-quotient is the first argument divided by the second one or just the first argument. The pseudo-remainder is always
equal to zero in this case.
```jldoctest pseudodivrem_example
julia> pseudodivrem(p, dpr(2))
(1//2*x[0]^2 + 1//2*y[0], 0)
```
"""
function pseudodivrem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem}
  deg_q = degree(q, i, jet)
  deg_q < 0 && throw(DivideError()) # By convention, (only) the zero polynomial has degree -1 in each jet variable
  @req deg_q > 0 "Cannot pseudo-divide by a polynomial with degree 0 in the specified division variable"
  
  # positive degree ensures existence of the key
  return pseudodivrem(p, q, __jtv(parent(q))[(i, jet)])
end

pseudodivrem(p::PolyT, q::PolyT, jet_idx::Tuple{Int, Vector{Int}}) where {PolyT <: ActionPolyRingElem} = pseudodivrem(p, q, jet_idx...)
pseudodivrem(p::PolyT, q::PolyT, i::Int) where {PolyT <: ActionPolyRingElem} = pseudodivrem(p, q, gen(parent(p), i))

function pseudodivrem(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  is_zero(q) && throw(DivideError())
  
  if is_constant(q)
    flag, quo = divides(p, q)
    if flag
      return (quo, zero(p))
    else
      return (p, zero(p))
    end
  end
  
  return pseudodivrem(p, q, leader(q))
end

### Ritt ordering ###

@doc raw"""
    ritt_is_less(p::ActionPolyRingElem, q::ActionPolyRingElem) -> Bool

Return `true` if `p` is smaller than `q` with respect to the Ritt ordering (associated to the
ranking) on the action polynomial ring containing `p` and `q`), otherwise return `false`.
"""
function ritt_is_less(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  if is_constant(p)
    if is_zero(p)
      is_zero(q) && return false
      return true
    end
    is_constant(q) && return false
    return true
  end
  
  is_constant(q) && return false

  ld_p = leader(p)
  ld_q = leader(q)

  ld_p < ld_q && return true
  ld_q < ld_p && return false

  return degree(p, ld_p) < degree(q, ld_p)
end

###############################################################################
#
#  Reduction methods
#
###############################################################################

@doc raw"""
    partially_reduce(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the partial remainder of `p` with respect to the non-zero polynomial `q`, by performing successive pseudo-divisions
of `p` by proper derivatives (or shifts) of `q`. This means that the returned polynomial has strictly smaller degree in each jet
variable that is a proper derivative (or shift) of the leader of `q`. If `q` is a non-zero constant then the zero polynomial is
returned.
"""
partially_reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem} = __core_partially_reduce(p, q, false)[1]

@doc raw"""
    partially_reduce(p::ActionPolyRingElem, S::Vector{ActionPolyRingElem})

Partially reduce the action polynomial `p` with respect to the vector `S`. This is done by pre-sorting `S`
with respect to Ritt ordering, filtering out zero polynomials and then performing top-down partial
reductions of `p` by the remaining elements of `S` until no further reductions are possible.
"""
partially_reduce(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem} = __core_partially_reduce(p, S, false)[1]

@doc raw"""
    reduce(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the full remainder of `p` with respect to the nonzero polynomial `q`. This remainder is
reduced with respect to `q` in the sense that the degree of `p` in each derivative (or shift) of the
leader of `q` is strictly smaller than the degree of the respective derivative (or shift) of `q` in that
same jet variable. If `q` is a nonzero constant then the zero polynomial is returned.
"""
reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem} = __core_reduce(p, q, false)[1]

@doc raw"""
    reduce(p::ActionPolyRingElem, S::Vector{ActionPolyRingElem})

Reduce the action polynomial `p` with respect to the vector `S`. This is done by pre-sorting `S`
with respect to Ritt ordering and then performing top-down reductions of `p` by the elements of
`S` until no further reductions are possible.
"""
reduce(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem} = __core_reduce(p, S, false)[1]

@doc raw"""
    autoreduce(S::Vector{ActionPolyRingElem})

Compute an autoreduced set from the vector of action polynomials `S`. If at any point in the reduction
process a nonzero constant is discovered, the vector containing just this constant is returned.
"""
function autoreduce(S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
  # S will serve as the working list throughout this algorithm
  S = filter(!is_zero, S)
  sort!(S, lt=ritt_is_less)
  
  A = PolyT[] # The "successively" built up autoreduced set 
  while !isempty(S)
    if !isempty(A) && is_constant(A[1])
      return [A[1]]
    elseif !isempty(S) && is_constant(S[1])
      return [S[1]]
    end
    
    p = popfirst!(S)
    original_p = p
    
    rank_dropped = false
    for q in Iterators.reverse(A)
      p = reduce(p, q) 
      
      if is_zero(p)
        break
      end
      
      # Reduce by coefficients in the base ring 
      c_gcd = foldl(gcd, coefficients(p))
      if !isone(c_gcd)
        p = divexact(p, c_gcd)
      end
      
      # Detect rank drops of p
      if ritt_is_less(p, original_p)
        rank_dropped = true
        break 
      end
    end
    
    # Nothing to do here
    if is_zero(p)
      continue 
    end
    
    # Reinsert both p and all element of A greater than p back into the working list S. 
    if rank_dropped
      # inserting elts of A
      while !isempty(A) && ritt_is_less(p, A[end])
        q_popped = pop!(A)
        idx = searchsortedfirst(S, q_popped, lt=ritt_is_less)
        insert!(S, idx, q_popped)
      end
      
      # inserting p
      idx = searchsortedfirst(S, p, lt=ritt_is_less)
      insert!(S, idx, p)
      
      continue
    end
    
    # No rank drop ever occurred, so appending p to A will give an autoreduced set
    push!(A, p)
  end
  
  return A
end

###############################################################################
#
#  Reduction verifiers
#
###############################################################################

#=
@doc raw"""
    leader_shift_for_partial_reduction(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the jet difference vector `k` such that `apply_action(q, k)` has the same leader
as the highest-ranking unreduced jet variable in `p`. If `p` is already partially reduced 
with respect to `q`, return `nothing`. Note that the notion of an 'unreduced' jet variable
depends on the type of `p` (and `q`):

- If the inputs are of type `DifferentialPolyRingElem`, then the occurrence of the jet variable in
  question is sufficient.
- If the inputs are of type `DifferencePolyRingElem`, then the degree of `p` in this variable must be
  greater or equal to the degree of `q` in its leader.
"""=#
function __leader_shift_for_partial_reduction(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  # parent checks are always performed before this method is called from an exported method
  @req !is_constant(q) "Cannot compute leader shift with respect to a constant polynomial"
  is_constant(p) && return nothing
  
  apr = parent(q)
  ld_q = leader(q)
  ld_i, ld_idx = __vtj(apr)[ld_q]
  
  for var in vars(p)
    var_i, var_idx = __vtj(apr)[var]
    
    if var_i == ld_i
      is_derivative = true
      for k in 1:length(ld_idx)
        if var_idx[k] < ld_idx[k]
          is_derivative = false
          break
        end
      end 
    
      if is_derivative && var_idx != ld_idx
        if __is_proper_shift_reducible(p, q, var)
          return var_idx .- ld_idx 
        end
      end 
    end 
  end 
  
  return nothing
end

__is_proper_shift_reducible(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: DifferencePolyRingElem} = degree(p, var) >= degree(q, leader(q))
__is_proper_shift_reducible(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: DifferentialPolyRingElem} = true

@doc raw"""
    is_partially_reduced(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return `true` if the action polynomial `p` is partially reduced with respect to the nonzero action polynomial `q`.
"""
function is_partially_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  @req !is_zero(q) "Cannot partially reduce with respect to the zero polynomial"
  is_constant(q) && return is_zero(p)
  return isnothing(__leader_shift_for_partial_reduction(p, q))
end

@doc raw"""
    is_partially_reduced(p::ActionPolyRingElem, S::Vector{<:ActionPolyRingElem})

Return `true` if `p` is partially reduced with respect to every element in the set `S`.
"""
is_partially_reduced(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem} = all(q -> is_partially_reduced(p, q), S)

@doc raw"""
    is_reduced(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return `true` if the polynomial `p` is fully reduced with respect to the non-constant
polynomial `q`. This means `p` is partially reduced with respect to `q`, and the degree 
of `p` in the leader of `q` is strictly less than the degree of `q` in its leader.
"""
function is_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_zero(q) "Cannot reduce with respect to the zero polynomial"
  is_constant(q) && return is_zero(p)
  
  is_constant(p) && return true
  is_partially_reduced(p, q) || return false
  
  ld_q = leader(q)
  return degree(p, ld_q) < degree(q, ld_q)
end

@doc raw"""
    is_reduced(p::ActionPolyRingElem, S::Vector{<:ActionPolyRingElem})

Return `true` if `p` is fully reduced with respect to every element in the set `S`.
"""
is_reduced(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem} = all(q -> is_reduced(p, q), S)

@doc raw"""
    is_autoreduced(S::Vector{<:ActionPolyRingElem})

Return `true` if `S` is an autoreduced set. A set is autoreduced if every polynomial 
in the set is fully reduced with respect to all other polynomials in the set, and 
the set is sorted by Ritt ordering.
"""
function is_autoreduced(S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
  any(is_constant, S) && return (length(S) == 1 && !is_zero(S[1]) && is_constant(S[1]))
  !issorted(S, lt=ritt_is_less) && return false
  
  for i in 1:length(S)
    for j in 1:length(S)
      i == j && continue
      !is_reduced(S[i], S[j]) && return false
    end
  end
  
  return true
end

###############################################################################
#
#  Internal Thomas decomposition API
#
###############################################################################

__pseudorem_with_factors(p::P, q::P, v::P) where {P <: Union{MPolyRingElem, ActionPolyRingElem}} = __core_pseudorem(p, q, v, true)
__pseudodivrem_with_factors(p::P, q::P, v::P) where {P <: Union{MPolyRingElem, ActionPolyRingElem}} = __core_pseudodivrem(p, q, v, true)
__partially_reduce_with_factors(p::P, q::P) where {P <: Union{MPolyRingElem, ActionPolyRingElem}} = __core_partially_reduce(p, q, true)
__partially_reduce_with_factors(p::P, S::Vector{P}) where {P <: Union{MPolyRingElem, ActionPolyRingElem}} = __core_partially_reduce(p, S, true)
__reduce_with_factors(p::P, q::P) where {P <: Union{MPolyRingElem, ActionPolyRingElem}} = __core_reduce(p, q, true)
__reduce_with_factors(p::P, S::Vector{P}) where {P <: Union{MPolyRingElem, ActionPolyRingElem}} = __core_reduce(p, S, true)

###############################################################################
#
#  Core methods
#
###############################################################################

# For the Thomas decomposition it is necessary to store certain factors that occur during computation of
# pseudo-remainders (and propagate them through further methods). Moreover, the Thomas decomposition of
# a differential system (represented by vectors of ActionPolyRingElems) uses the Thomas decomposition of
# an algebraic system (represented by vectors of MPolyRingElems) as a subroutine. Thus, we need many
# subroutines that work for both polynomial types. These methods get the prefix '__core_'. Of course, a
# lot of these methods are of interest in their own.

function __core_pseudorem(p::P, q::P, v::P, track_factors::Bool) where {P <: Union{MPolyRingElem, ActionPolyRingElem}}
  deg_q = degree(q, v)
  deg_q < 0 && throw(DivideError())
  @req deg_q > 0 "Cannot pseudo-divide by a polynomial with degree 0 in the specified division variable"
  
  factors = P[]
  degree(p, v) < deg_q && return (p, factors)

  lc_q = __univariate_leading_coefficient(q, v)
  rem = deepcopy(p)

  while !is_zero(rem) && (deg_rem = degree(rem, v)) >= deg_q
    lc_rem = __univariate_leading_coefficient(rem, v)
    
    flag, c = divides(lc_rem, lc_q)
    
    if flag
      sub!(rem, rem, c * (v^(deg_rem - deg_q)) * q)
    else
      mul!(rem, lc_q, rem)
      sub!(rem, rem, lc_rem * (v^(deg_rem - deg_q)) * q)
      if track_factors
        push!(factors, lc_q)
      end
    end
  end

  return (rem, factors)
end

function __core_pseudodivrem(p::P, q::P, v::P, track_factors::Bool) where {P <: Union{MPolyRingElem, ActionPolyRingElem}}
  deg_q = degree(q, v)
  deg_q < 0 && throw(DivideError())
  @req deg_q > 0 "Cannot pseudo-divide by a polynomial with degree 0 in the specified division variable"
    
  factors = P[]
  quo = zero(parent(p))
  degree(p, v) < deg_q && return (quo, p, factors)

  lc_q = __univariate_leading_coefficient(q, v)
  rem = deepcopy(p)

  while !is_zero(rem) && (deg_rem = degree(rem, v)) >= deg_q
    lc_rem = __univariate_leading_coefficient(rem, v)
        
    flag, c = divides(lc_rem, lc_q)

    if flag
      quo_term = c * (v^(deg_rem - deg_q))
      add!(quo, quo, quo_term)
      sub!(rem, rem, quo_term * q)
    else
      mul!(rem, lc_q, rem)
      mul!(quo, lc_q, quo)
      quo_term = lc_rem * (v^(deg_rem - deg_q))
      add!(quo, quo, quo_term)
      sub!(rem, rem, quo_term * q)
      if track_factors
        push!(factors, lc_q)
      end
    end
  end
    
  return (quo, rem, factors)
end

# MPolys are always partially reduced
__core_partially_reduce(p::P, q::P, track_factors::Bool) where {P <: MPolyRingElem} = (p, P[])
__core_partially_reduce(p::P, S::Vector{P}, track_factors::Bool) where {P <: MPolyRingElem} = (p, P[])

function __core_partially_reduce(p::P, q::P, track_factors::Bool) where {P <: ActionPolyRingElem}
  @req !is_zero(q) "Cannot partially reduce with respect to the zero polynomial"
  factors = P[]
  is_constant(q) && return (zero(p), factors)  

  p_red = deepcopy(p)
  shift = __leader_shift_for_partial_reduction(p_red, q)
    
  while shift !== nothing
    q_shifted = apply_action(q, shift)
    ld_q_shifted = leader(q_shifted)
    p_red, step_factors = __core_pseudorem(p_red, q_shifted, ld_q_shifted, track_factors)
    
    if track_factors
      append!(factors, step_factors)
    end
    
    shift = __leader_shift_for_partial_reduction(p_red, q)
  end
    
  return (p_red, factors)
end

function __core_partially_reduce(p::P, S::Vector{P}, track_factors::Bool) where {P <: ActionPolyRingElem}
  has_const = false
  for q in S
    @req !is_zero(q) "Cannot partially reduce with respect to a set containing the zero polynomial"
    if !has_const && is_constant(q)
      has_const = true
    end
  end
  
  factors = P[]
  has_const && return (zero(p), factors)

  sorted_S = sort(S, lt=ritt_is_less)

  res = deepcopy(p)
  changed = true
  while changed && !is_zero(res)
    changed = false
    for q in Iterators.reverse(sorted_S)
      original_res = deepcopy(res)
      res, step_factors = __core_partially_reduce(res, q, track_factors)
      
      if track_factors
        append!(factors, step_factors)
      end
      if res != original_res
        changed = true
      end
    end
  end

  return (res, factors)
end

function __core_reduce(p::P, q::P, track_factors::Bool) where {P <: Union{MPolyRingElem, ActionPolyRingElem}}
  @req !is_zero(q) "Cannot reduce with respect to the zero polynomial"
  factors = P[]
  is_constant(q) && return (zero(p), factors)
  
  p_red, p_factors = __core_partially_reduce(p, q, track_factors)
  if track_factors
    append!(factors, p_factors)
  end
  
  res, b_factors = __core_pseudorem(p_red, q, __leader(q), track_factors)
  if track_factors
    append!(factors, b_factors)
  end
  
  return (res, factors)
end

function __core_reduce(p::P, S::Vector{P}, track_factors::Bool) where {P <: Union{MPolyRingElem, ActionPolyRingElem}}
  has_const = false
  for q in S
    @req !is_zero(q) "Cannot reduce with respect to a set containing the zero polynomial"
    if !has_const && is_constant(q)
      has_const = true
    end
  end
  
  factors = P[]
  has_const && return (zero(p), factors)
  
  sorted_S = sort(S, lt=__ritt_is_less)

  res = deepcopy(p)
  changed = true
  
  while changed && !is_zero(res)
    changed = false
    for q in Iterators.reverse(sorted_S)
      original_res = deepcopy(res)
      
      # 1. Partial Reduction
      res, p_factors = __core_partially_reduce(res, q, track_factors)
      if track_factors
        append!(factors, p_factors)
      end
      
      # 2. Algebraic Reduction
      ld_q = __leader(q)
      if degree(res, ld_q) >= degree(q, ld_q) 
        res, b_step = __core_pseudorem(res, q, ld_q, track_factors)
        if track_factors
          append!(factors, b_step)
        end
      end
      
      if res != original_res
        changed = true
        is_zero(res) && break
      end
    end
  end

  return (res, factors)
end

###############################################################################
#
#  Unifying helpers
#
###############################################################################

# This section extends standard functionality from action polynomial rings to mpoly rings
# Since these are not intended to be exported, for each method we also provide a wrapper for
# the associated action polynomial method.

# Return the highest-ranked variable present in `p`. 
# Variables are assumed to be sorted `x_1 > x_2 > ... > x_n` matching `gens(parent(p))`.
function __leader(p::MPolyRingElem)
  @req !is_zero(p) "The zero polynomial has no leader"
  is_constant(p) && return one(parent(p))
    
  R = parent(p)
  for i in 1:nvars(R)
    v = gen(R, i)
    degree(p, v) > 0 && return v
  end
  return one(R)
end
__leader(p::ActionPolyRingElem) = leader(p)

# Return the leading coefficient of `p` regarded as a univariate polynomial in the variable `v`.
function __univariate_leading_coefficient(p::P, v::P) where {P <: MPolyRingElem}
  @req is_gen(v) "Not a variable"
  d = degree(p, v)
  @req d > -1 "The zero polynomial has no leading coefficient"
  d == 0 && return p
    
  res = zero(p)
  v_idx = var_index(v)
    
  for (t, e) in zip(terms(p), exponents(p))
    if @inbounds e[v_idx] == d
      res += remove(t, v)[2]
    end
  end
  return res
end
__univariate_leading_coefficient(p::P, v::P) where {P <: ActionPolyRingElem} = univariate_leading_coefficient(p, v)

# Return the initial of the polynomial `p`.
function __initial(p::MPolyRingElem)
  @req !is_zero(p) "The zero polynomial has no initial"
  is_constant(p) && return p
  return __univariate_leading_coefficient(p, __leader(p))
end
__initial(p::ActionPolyRingElem) = initial(p)

# Return the discriminant of `p`.
function __discriminant(p::MPolyRingElem)
  is_constant(p) && return is_zero(p) ? p : one(p)
  
  ld_p = __leader(p)
  i = var_index(ld_p)
  
  if degree(p, ld_p) % 4 in (0, 1)
    return divexact(resultant(p, derivative(p, i), i), __initial(p))
  else
    return -divexact(resultant(p, derivative(p, i), i), __initial(p))
  end
end
__discriminant(p::ActionPolyRingElem) = discriminant(p)

function __ritt_is_less(p::P, q::P) where {P <: MPolyRingElem}
  if is_constant(p)
    if is_zero(p)
      is_zero(q) && return false
      return true
    end
    is_constant(q) && return false
    return true
  end
  
  is_constant(q) && return false

  ld_p = __leader(p)
  ld_q = __leader(q)
  idx_p = var_index(ld_p)
  idx_q = var_index(ld_q)

  idx_p > idx_q && return true
  idx_q > idx_p && return false

  return degree(p, ld_p) < degree(q, ld_q)
end
__ritt_is_less(p::P, q::P) where {P <: ActionPolyRingElem} = ritt_is_less(p, q)

