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

