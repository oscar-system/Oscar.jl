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

