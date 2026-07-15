###############################################################################
#
#  Unary Operators 
#
###############################################################################

Base.:-(apre::ActionPolyRingElem) = parent(apre)(-data(apre))

###############################################################################
#
#  Binary Operators 
#
###############################################################################

function Base.:+(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) + data(apre2))
end

function Base.:-(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) - data(apre2))
end

function Base.:*(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) * data(apre2))
end

function Base.:/(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(data(apre1) / data(apre2))
end

function divexact(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}; check::Bool = true) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(divexact(data(apre1), data(apre2)))
end

function divides(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T}
  check_parent(apre1, apre2)
  flag, res = divides(data(apre1), data(apre2))
  return flag, parent(apre1)(res)
end

function Base.div(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T}
  check_parent(apre1, apre2) 
  return parent(apre1)(div(data(apre1), data(apre2)))
end

function Base.divrem(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T}
  check_parent(apre1, apre2) 
  return parent(apre1).(divrem(data(apre1), data(apre2)))
end

function mod(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2)
  return parent(apre1)(mod(data(apre1), data(apre2)))
end

#----

function pseudorem(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  check_parent(p, var)
  @req is_gen(var) "Not a jet variable"
  return pseudorem(p, q, __vtj(parent(var))[var]...)
end

@doc raw"""
    pseudorem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem} -> PolyT

Return the algebraic pseudo-remainder of `p` divided by `q` with respect to the jet variable specified by `i` and
`jet`. If no jet variable is specified then division is performed with respect to the leader of `q`.
This method performs division by using a lazy pre-multiplication by the initial of `q` at each step, only
multiplying the remainder when necessary.

This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables); see the online documentation.

# Examples

```jldoctest pseudorem_example
julia> dpr, (x,y) = differential_polynomial_ring(QQ, [:x,:y], 1); p, q = (x^2 + y, y*x + 1)
(x[0]^2 + y[0], y[0]*x[0] + 1)

julia> pseudorem(p,q,x)
y[0]^3 + 1

julia> pseudorem(p,q,y)
x[0]^3 - 1
```
Note that in this example, the leader of `q` is `x`, so it need not be specified for pseudo-division with respect to `x`:
```jldoctest pseudorem_example
julia> leader(q) == x
true

julia> pseudorem(p,q)
y[0]^3 + 1
```
"""
function pseudorem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  
  deg_q = degree(q, i, jet)
  
  deg_q < 0 && throw(DivideError()) # By convention, (only) the zero polynomial has degree -1 in each jet variable
  @req deg_q > 0 "Cannot pseudo-divide by a polynomial with degree 0 in the provided jet variable"
  degree(p, i, jet) < deg_q && return p # p already algebraically reduced wrt q

  var = __jtv(parent(q))[(i, jet)]
  lc_q = univariate_leading_coefficient(q, i, jet)
  rem = deepcopy(p)

  while !is_zero(rem) && (deg_rem = degree(rem, i, jet)) >= deg_q
    lc_rem = univariate_leading_coefficient(rem, i, jet)
    
    # Check if pre-multiplying by lc_q is necessary
    flag, c = divides(lc_rem, lc_q)
    
    if flag
      sub!(rem, rem, c * (var^(deg_rem - deg_q)) * q)
    else
      mul!(rem, lc_q, rem)
      sub!(rem, rem, lc_rem * (var^(deg_rem - deg_q)) * q)
    end
  end

  return rem
end

pseudorem(p::PolyT, q::PolyT, jet_idx::Tuple{Int, Vector{Int}}) where {PolyT <: ActionPolyRingElem} = pseudorem(p, q, jet_idx...)
pseudorem(p::PolyT, q::PolyT, i::Int) where {PolyT <: ActionPolyRingElem} = pseudorem(p, q, gen(parent(p), i))

function pseudorem(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  is_zero(q) && throw(DivideError())
  return pseudorem(p, q, leader(q))
end

#----

function pseudodivrem(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  check_parent(p, var)
  @req is_gen(var) "Not a jet variable"
  return pseudodivrem(p, q, __vtj(parent(var))[var]...)
end

@doc raw"""
    pseudodivrem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem} -> Tuple{PolyT, PolyT}

Return the pair `(s,r)` where `s` is the pseudo-quotient and `r` is the pseudo-remainder of `p`
by `q` with respect to the jet variable specified by `i` and `jet`. If no jet variable is specified
then division is performed with respect to the leader of `q`. The number of pre-multiplications
by the leading coefficient of `q` in this jet variable is minimised, i.e. we have
`lc(q)^k * p = s * q + r` where the integer `k >= 0` is minimal and `lc(q)` is the above mentioned
leading coefficient.

This method allows all versions described in [Specifying jet variables](@ref specifying_jet_variables); see the online documentation.

# Examples

```jldoctest pseudodivrem_example
julia> dpr, (x,y) = differential_polynomial_ring(QQ, [:x,:y], 1); p, q = (x^2 + y, y*x + 1)
(x[0]^2 + y[0], y[0]*x[0] + 1)

julia> pseudodivrem(p,q,x)
(y[0]*x[0] - 1, y[0]^3 + 1)

julia> y^2 * p == (y*x - 1) * q + (y^3 + 1)
true

julia> pseudodivrem(p,q,y)
(1, x[0]^3 - 1)

julia> x * p == 1 * q + (x^3 - 1)
true
```
Note that in this example, the leader of `q` is `x`, so it need not be specified for pseudo-division with respect to `x`:
```jldoctest pseudodivrem_example
julia> leader(q) == x
true

julia> pseudodivrem(p,q)
(y[0]*x[0] - 1, y[0]^3 + 1)
```
"""
function pseudodivrem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  
  deg_q = degree(q, i, jet)
  
  deg_q < 0 && throw(DivideError()) # By convention, (only) the zero polynomial has degree -1 in all jet variables
  @req deg_q > 0 "Cannot pseudo-divide by a polynomial with degree 0 in the division variable"
  degree(p, i, jet) < deg_q && return (zero(parent(p)), p)

  var = __jtv(parent(q))[(i, jet)]
  
  lc_q = univariate_leading_coefficient(q, i, jet)
  rem = deepcopy(p)
  quo = zero(parent(p))

  while !is_zero(rem) && (deg_rem = degree(rem, i, jet)) >= deg_q
    lc_rem = univariate_leading_coefficient(rem, i, jet)
    
    # Check if pre-multiplying by lc_q is necessary
    flag, c = divides(lc_rem, lc_q)

    if flag
      quo_term = c * (var^(deg_rem - deg_q))

      add!(quo, quo, quo_term)
      sub!(rem, rem, quo_term * q)
    else
      mul!(rem, lc_q, rem)
      mul!(quo, lc_q, quo)

      quo_term = lc_rem * (var^(deg_rem - deg_q))

      add!(quo, quo, quo_term)
      sub!(rem, rem, quo_term * q)
    end
  end
  
  return (quo, rem)
end

pseudodivrem(p::PolyT, q::PolyT, jet_idx::Tuple{Int, Vector{Int}}) where {PolyT <: ActionPolyRingElem} = pseudodivrem(p, q, jet_idx...)
pseudodivrem(p::PolyT, q::PolyT, i::Int) where {PolyT <: ActionPolyRingElem} = pseudodivrem(p, q, gen(parent(p), i))

function pseudodivrem(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  is_zero(q) && throw(DivideError())
  return pseudodivrem(p, q, leader(q))
end

###############################################################################
#
#  gcd and lcm 
#
###############################################################################

function gcd(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
  check_parent(apre1, apre2) 
  return parent(apre1)(gcd(data(apre1), data(apre2)))
end

function lcm(apre1::ActionPolyRingElem{T}, apre2::ActionPolyRingElem{T}) where {T<:RingElement}
   check_parent(apre1, apre2)
   return parent(apre1)(lcm(data(apre1), data(apre2)))
end

###############################################################################
#
#  Remove and valuation 
#
###############################################################################

function remove(z::ActionPolyRingElem{T}, p::ActionPolyRingElem{T}) where {T}
  check_parent(z, p)
  val, q = remove(data(z), data(p))
  return val, parent(z)(q)
end

valuation(z::ActionPolyRingElem{T}, p::ActionPolyRingElem{T}) where {T} = remove(z, p)[1]

###############################################################################
#
#  Comparison 
#
###############################################################################

### Difference ###
Base.:(==)(dpre1::DifferencePolyRingElem{T}, dpre2::DifferencePolyRingElem{T}) where {T<:RingElement} = parent(dpre1) === parent(dpre2) && data(dpre1) == data(dpre2)

function Base.hash(dpre::DifferencePolyRingElem, h::UInt)
  b = 0x475b3fa701aa3148 % UInt
  h = hash(parent(dpre), h)
  h = hash(data(dpre), h)
  return xor(h, b)
end

function Base.isless(dpre1::DifferencePolyRingElem{T}, dpre2::DifferencePolyRingElem{T}) where {T<:RingElement}
  check_parent(dpre1, dpre2)
  dpr = parent(dpre1)
  vtj = __vtj(dpr)
  @req haskey(vtj, dpre1) && haskey(vtj, dpre2) "Not jet variables in comparison"
  m = n_elementary_symbols(dpr)
  ind1, ind2 = vtj[dpre1], vtj[dpre2]
  v1, v2 = vcat([ind1[1] == j ? one(ZZ) : zero(ZZ) for j in 1:m], ind1[2]), vcat([ind2[1] == j ? one(ZZ) : zero(ZZ) for j in 1:m], ind2[2])
  M = riquier_matrix(ranking(dpr))
  return isless(M * v1, M * v2)
end

### Differential ###
Base.:(==)(dpre1::DifferentialPolyRingElem{T}, dpre2::DifferentialPolyRingElem{T}) where {T<:RingElement} = data(dpre1) == data(dpre2) && parent(dpre1) === parent(dpre2)

function Base.hash(dpre::DifferentialPolyRingElem, h::UInt)
  b = 0x5c93cee72de560dd % UInt
  h = hash(parent(dpre), h)
  h = hash(data(dpre), h)
  return xor(h, b)
end

function Base.isless(dpre1::DifferentialPolyRingElem{T}, dpre2::DifferentialPolyRingElem{T}) where {T<:RingElement}
  check_parent(dpre1, dpre2)
  dpr = parent(dpre1)
  vtj = __vtj(dpr)
  @req haskey(vtj, dpre1) && haskey(vtj, dpre2) "Not jet variables in comparison"
  m = n_elementary_symbols(dpr)
  ind1, ind2 = vtj[dpre1], vtj[dpre2]
  v1, v2 = vcat([ind1[1] == j ? one(ZZ) : zero(ZZ) for j in 1:m], ind1[2]), vcat([ind2[1] == j ? one(ZZ) : zero(ZZ) for j in 1:m], ind2[2])
  M = riquier_matrix(ranking(dpr))
  return isless(M * v1, M * v2)
end

###############################################################################
#
#  Unsafe functions 
#
###############################################################################

### Difference ###
function zero!(a::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = zero!(a.p)
  __set_is_perm_up_to_date!(a, false)
  return a
end

function one!(a::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = one!(a.p)
  __set_is_perm_up_to_date!(a, false)
  return a
end

function neg!(z::PolyT, a::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  z.p = neg!(z.p, a.p)
  __set_is_perm_up_to_date!(z, false)
  return z
end

function fit!(a::PolyT, n::Int) where {PolyT <: ActionPolyRingElem}
  fit!(data(a), n)
  __set_is_perm_up_to_date!(a, false)
end

function add!(a::PolyT, b::PolyT, c::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = add!(data(a), data(b), data(c))
  __set_is_perm_up_to_date!(a, false)
  return a
end

function add!(a::PolyT, b::PolyT, c::RingElement) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = add!(data(a), data(b), c)
  __set_is_perm_up_to_date!(a, false)
  return a
end

add!(a::PolyT, b::RingElement, c::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}} = add!(a, c, b)

function sub!(a::PolyT, b::PolyT, c::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = sub!(data(a), data(b), data(c))
  __set_is_perm_up_to_date!(a, false)
  return a
end

function sub!(a::PolyT, b::PolyT, c::RingElement) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = sub!(data(a), data(b), c)
  __set_is_perm_up_to_date!(a, false)
  return a
end

function sub!(a::PolyT, b::RingElement, c::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = sub!(data(a), b, data(c))
  __set_is_perm_up_to_date!(a, false)
  return a
end

function mul!(a::PolyT, b::PolyT, c::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = mul!(data(a), data(b), data(c))
  __set_is_perm_up_to_date!(a, false)
  return a
end

function mul!(a::PolyT, b::PolyT, c::RingElement) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}}
  a.p = mul!(data(a), data(b), c)
  __set_is_perm_up_to_date!(a, false)
  return a
end

mul!(a::PolyT, b::RingElement, c::PolyT) where {T <: RingElement, PolyT <: ActionPolyRingElem{T}} = mul!(a, c, b)

