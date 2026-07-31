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

function Base.:+(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  return parent(p)(data(p) + data(q))
end

function Base.:-(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  return parent(p)(data(p) - data(q))
end

function Base.:*(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  return parent(p)(data(p) * data(q))
end

function Base.:/(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  return parent(p)(data(p) / data(q))
end

function divexact(p::PolyT, q::PolyT; check::Bool = true) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  return parent(p)(divexact(data(p), data(q)))
end

function divides(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  flag, res = divides(data(p), data(q))
  return flag, parent(p)(res)
end

function Base.div(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q) 
  return parent(p)(div(data(p), data(q)))
end

function Base.divrem(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q) 
  return parent(p).(divrem(data(p), data(q)))
end

function mod(p::PolyT, q::PolyT) where {PolyT<:ActionPolyRingElem}
  check_parent(p, q)
  return parent(p)(mod(data(p), data(q)))
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
`jet`. If no jet variable is specified then division is performed with respect to the leader of `q`, even allowing
`q` to be a nonzero constant. This method performs division by using a lazy pre-multiplication by the initial of `q`
at each step, only multiplying the remainder when necessary.

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
Finally, the pseudoremainder of any polynomial with respect to a non-zero constant is always zero.
```jldoctest pseudorem_example
julia> pseudorem(p,dpr(1))
0
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
  is_constant(q) && return zero(q)
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
then division is performed with respect to the leader of `q`, even allowing `q` to be a nonzero constant.
The number of pre-multiplications by the leading coefficient of `q` in this jet variable is minimised,
i.e. we have `lc(q)^k * p = s * q + r` where the integer `k >= 0` is minimal and `lc(q)` is the above
mentioned leading coefficient.

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
Finally, if the second argument is a non-zero constant, then it depends on the first argument being divisible by said constant,
whether the pseudo-quotient is the first argument divided by the second one or just the first argument. The pseudo-remainder is always
equal to zero in this case.
```jldoctest pseudodivrem_example
julia> pseudodivrem(p,dpr(2))
(1//2*x[0]^2 + 1//2*y[0], 0)
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

###############################################################################
#
#  gcd and lcm 
#
###############################################################################

function gcd(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q) 
  return parent(p)(gcd(data(p), data(q)))
end

function lcm(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
   check_parent(p, q)
   return parent(p)(lcm(data(p), data(q)))
end

###############################################################################
#
#  Remove and valuation 
#
###############################################################################

function remove(z::PolyT, p::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(z, p)
  val, q = remove(data(z), data(p))
  return val, parent(z)(q)
end

valuation(z::PolyT, p::PolyT) where {PolyT <: ActionPolyRingElem} = remove(z, p)[1]

###############################################################################
#
#  Comparison 
#
###############################################################################

Base.:(==)(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem} = parent(p) === parent(q) && data(p) == data(q)

@doc raw"""
    isless(var1::ActionPolyRingElem, var2::ActionPolyRingElem) -> Bool

Return `true`, if `var1` is less than `var2` with respect to the Riquier ranking on
the action polynomial ring that contains the jet variables `var1` and `var2`. For comparing
of action polynomials with respect to Ritt ordering, use [`ritt_is_less`](@ref ritt_is_less).
"""
function Base.isless(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  apr = parent(p)
  vtj = __vtj(apr)
  
  @req haskey(vtj, p) && haskey(vtj, q) "Not jet variables in comparison"
  
  ind1 = vtj[p]
  ind2 = vtj[q]
  i1, J1 = ind1[1], ind1[2]
  i2, J2 = ind2[1], ind2[2]
  
  M = riquier_matrix(ranking(apr))
  m = n_action_indeterminates(apr)
  
  for k in 1:size(M, 1)
    val1 = M[k, i1]
    for l in 1:length(J1)
      val1 += M[k, m + l] * J1[l]
    end
    val2 = M[k, i2]
    for l in 1:length(J2)
      val2 += M[k, m + l] * J2[l]
    end
    val1 != val2 && return val1 < val2
  end
  
  return false
end

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

### Hashes ###
function Base.hash(p::DifferencePolyRingElem, h::UInt)
  b = 0x475b3fa701aa3148 % UInt
  h = hash(parent(p), h)
  h = hash(data(p), h)
  return xor(h, b)
end

function Base.hash(p::DifferentialPolyRingElem, h::UInt)
  b = 0x5c93cee72de560dd % UInt
  h = hash(parent(p), h)
  h = hash(data(p), h)
  return xor(h, b)
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

