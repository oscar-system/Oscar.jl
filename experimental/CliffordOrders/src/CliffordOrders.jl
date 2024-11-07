export
  CliffordOrder,
  CliffordOrderElem,
  clifford_order,
  ambient_algebra,
  coefficient_ideals

#############################################################
#
#  Datatypes
#
#############################################################

mutable struct CliffordOrder{T, C} <: Hecke.AbstractAssociativeAlgebra{T}

  base_ring::Ring
  ambient_algebra::C
  rank::Int
  lattice::Any
  gram::MatElem
  coefficient_ideals::Vector{NumFieldOrderFractionalIdeal}
  centroid::Any
  disq::T
  center::Any

  function CliffordOrder{T, C}(ls::ZZLat) where {T, C}
    @req is_integral(1//2 * norm(ls)) "The given lattice is not even!"
    qs = rational_span(ls)
    coeff_ids = fill(fractional_ideal(ZZ, 1), 2^rank(ls))
    
    return new{T, C}(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs), coeff_ids)
  end

  function CliffordOrder{T, C}(ls::QuadLat{S, M}) where {T, C, S<:Field, M<:MatElem}
    @req is_integral(fractional_ideal(base_ring(ls), base_field(ls)(1//2)) * norm(ls)) "The given lattice is not even!"
    qs = rational_span(ls)
    coeff_ids = _set_coefficient_ideals(ls)
    return new{T, C}(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs), coeff_ids)
  end

end


### Elements ###
mutable struct CliffordOrderElem{T, C, P} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::P
  coeffs::Vector{S} where S<:FieldElem
  even_coeffs::Vector{S} where S<:FieldElem
  odd_coeffs::Vector{S} where S<:FieldElem

  #Returns the 0-element of the Clifford order C
  function CliffordOrderElem{T, C, P}(CO::CliffordOrder{T, C}) where {T, C, P}
    ambalg = CO.ambient_algebra
    newelt = new{T, C, P}(CO, fill(ambalg.base_ring(), CO.rank))
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  CliffordOrderElem(CO::CliffordOrder) = CliffordOrderElem{elem_type(CO.base_ring), typeof(CO.ambient_algebra), typeof(CO)}(CO)

  #Returns the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function CliffordOrderElem{T, C, P}(CO::CliffordOrder{T, C}, coeff::Vector{S}) where {T, C, P, S<:FieldElem}
    @req length(coeff) == CO.rank "invalid length of coefficient vector"
    
    if CO.base_ring == ZZ
      for i in 1:CO.rank
        @req is_integer(coeff[i]) "The element does not lie in the Clifford order."
      end
    else
      for i in 1:CO.rank
        @req coeff[i] in CO.coefficient_ideals[i] "The element does not lie in the Clifford order."
      end
    end

    newelt = new{T, C, CliffordOrder{T, C}}(CO, coeff)
    _set_even_odd_coeff!(newelt)
    return newelt
  end

  function CliffordOrderElem(CO::CliffordOrder{T, C}, coeff::Vector{S}) where {T, C, S}
    K = CO.ambient_algebra.base_ring
    @req _can_convert_coeff(coeff, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordOrderElem{elem_type(CO.base_ring), typeof(CO.ambient_algebra), typeof(CO)}(CO, K.(coeff))
  end

end

elem_type(::Type{CliffordOrder{T, C}}) where {T, C} = CliffordOrderElem{T, C, CliffordOrder{T, C}}

parent_type(::Type{CliffordOrderElem{T, C, P}}) where {T, C, P} = P

base_ring_type(CO::CliffordOrder{T, C}) where {T, C} = typeof(base_ring(CO))

#############################################################
#
#  Construction
#
#############################################################

### Order ###
@doc raw"""
    clifford_order(ls::Union{ZZLat, QuadLat}) -> CliffordOrder

Return the Clifford order of the even quadratic lattice $ls$.
"""
clifford_order(ls::ZZLat) = CliffordOrder{ZZRingElem, CliffordAlgebra{QQFieldElem, QQMatrix}}(ls)

clifford_order(ls::QuadLat) = CliffordOrder{elem_type(base_ring(ls)), typeof(clifford_algebra(rational_span(ls)))}(ls)

### Elements ###
(C::CliffordOrder)() = CliffordOrderElem(C)

function (C::CliffordOrder)(a::S) where {S<:RingElem}
  res = fill(zero(a), rank(C))
  res[1] = a
  return CliffordOrderElem(C, res)
end

function (C::CliffordOrder)(a::S) where {S<:Number}
  res = fill(zero(a), rank(C))
  res[1] = a
  return CliffordOrderElem(C, res)
end

(C::CliffordOrder)(coeff::Vector{S}) where {S} = CliffordOrderElem(C, coeff)

################################################################################
#
#  Basic field access
#
################################################################################

### Order  ###
@doc raw"""
    base_ring(C::CliffordOrder) -> Ring

Return the base ring of the Clifford order $C$.
"""
base_ring(C::CliffordOrder) = C.base_ring

@doc raw"""
    ambient_algebra(C::CliffordOrder) -> CliffordAlgebra

Return the Clifford algebra of the ambient space of the underlying lattice of $C$.
"""
ambient_algebra(C::CliffordOrder) = C.ambient_algebra

@doc raw"""
    rank(C::CliffordOrder) -> Int

Return the rank of the Clifford order $C$ over its base ring.
"""

rank(C::CliffordOrder) = C.rank
@doc raw"""
    lattice(C::CliffordOrder) -> Union{ZZLat, QuadLat}

Return the underlying even quadratic lattice of the Clifford order $C$.
"""
lattice(C::CliffordOrder) = C.lattice

@doc raw"""
    gram_matrix(C::CliffordOrder) -> MatElem

Return the Gram matrix with respect to the fixed pseudo-basis of the Clifford order $C$.
"""
gram_matrix(C::CliffordOrder) = C.gram

@doc raw"""
    coefficient_ideals(C::CliffordOrder) -> Vector{NumFieldOrderFractionalIdeal}

Return the vector of coefficient ideals of the canonical pseudo-basis of $C$.
"""
coefficient_ideals(C::CliffordOrder) = C.coefficient_ideals

### Elements ###
@doc raw"""
    parent(x::CliffordOrderElem) -> CliffordOrder

Return the Clifford order containing $x$.
"""
parent(x::CliffordOrderElem) = x.parent

@doc raw"""
  coeff(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ with respect to the
canonical pseudo-basis of its parent Clifford order.
"""
coeff(x::CliffordOrderElem) = x.coeffs

@doc raw"""
    even_coeff(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with odd grading are set to zero. This also updates
the field `x.even_coeffs`.
"""
function even_coeff(x::CliffordOrderElem)
  if isdefined(x, :even_coeffs)
    return x.even_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.even_coeffs
end

@doc raw"""
    odd_coeff(x::CliffordOrderElem) -> Vector

Return the coefficient vector of $x$ wrt the
canonical pseudo-basis of its parent Clifford order,
but all coefficients corresponding to basis elements
with even grading are set to zero. This also updates
the field `x.odd_coeffs`.
"""
function odd_coeff(x::CliffordOrderElem)
  if isdefined(x, :odd_coeffs)
    return x.odd_coeffs
  end
  _set_even_odd_coeff!(x)
  return x.odd_coeffs
end

#############################################################
#
#  String I/O
#
#############################################################

### Algebra ###
function Base.show(io::IO, ::MIME"text/plain", C::CliffordOrder)
  io = pretty(io)
  print(io, "Clifford order of even lattice over $(base_ring(C)) with Gram matrix\n")
  print(io, Indent())
  show(io, "text/plain", gram_matrix(C))
  print(io, Dedent(), "\nand coefficient ideals\n")
  print(io, Indent())
  show(io, "text/plain", _coefficient_ideals_of_lattice(lattice(C)))
  print(io, Dedent())
end

Base.show(io::IO, C::CliffordOrder) = print(io, "Clifford order over $(base_ring(C))")

### Elements ###
function Base.show(io::IO, x::CliffordOrderElem)
  print(io, "(")
  foreach(y -> print(io, "$y "), coeff(x)[1:(end - 1)])
  print(io, "$(coeff(x)[end]))")
end

################################################################################
#
#  basic functionality
#
################################################################################

@doc raw"""
    zero(C::CliffordOrder) -> CliffordOrderElem

Return the additive identity of the Clifford order $C$.
"""
zero(C::CliffordOrder) = C()

@doc raw"""
    one(C::CliffordOrder) -> CliffordOrderElem

Return the multiplicative identity of the Clifford order $C$.
"""
function one(C::CliffordOrder)
  v = fill(base_ring(ambient_algebra(C))(), rank(C))
  v[1] = base_ring(ambient_algebra(C))(1)
  return C(v)
end

@doc raw"""
    ambient_basis(C::CliffordOrder, i::Int) -> CliffordOrderElem

Return the $i$-th canonical basis element of the ambient algebra of $C$
as an element of $C$.
"""
function ambient_basis(C::CliffordOrder, i::Int)
  res = C()
  res.coeffs[i] = base_ring(ambient_algebra(C))(1)
  _set_even_odd_coeff!(res)
  return (res, coefficient_ideals(C)[i])
end

@doc raw"""
    gen(C::CliffordOrder, i::Int) -> CliffordOrderElem

Return the $i$-th canonical multiplicative generator of the Clifford
order $C$. This is just the $i$-th basis vector of ambient quadratic
space.
"""
function gen(C::CliffordOrder, i::Int)
  res = C()
  if i <= 0
    res.coeffs[i] #Throws a BoundsError instead of a DomainError for consistency
  end
  res.coeffs[2^(i - 1) + 1] = base_ring(ambient_algebra(C))(1)
  _set_even_odd_coeff!(res)
  return res
end

@doc raw"""
    ambient_basis(C::CliffordOrder) -> Vector{CliffordOrderElem}

Returns the canonical basis of the Clifford order $C$.
"""
basis(C::CliffordOrder) = map(i -> basis(C, i), 1:rank(C))

@doc raw"""
    gens(C::CliffordOrder) -> Vector{CliffordOrderElem}

Returns the vector of canonical multiplicative generators of the Clifford order $C$.
"""
gens(C::CliffordOrder) = map(i -> gen(C, i), 1:rank(lattice(C)))

@doc raw"""
    is_commutative(C::CliffordOrder) -> Bool

Returns `true` if $C$ is commutative and `false` otherwise.
"""
is_commutative(C::CliffordOrder) = rank(C) == 1 || rank(C) == 2

#############################################################
#
#  Auxillary functions
#
#############################################################

# Computes the ideals present in the canonical pseudo-basis of the Clifford order
# of the lattice with ideals 'lat_ids' in its fixed pseudo-basis.
function _set_coefficient_ideals(ls::QuadLat)
  coeff_ids = coefficient_ideals(ls)
  res = [fractional_ideal(base_ring(ls), base_ring(ls)(1))]
  for i in 1:rank(ls)
    tmpres = map(x -> x*coeff_ids[i], res)
    res = vcat(res, tmpres)
  end
  return res
end

function _set_even_odd_coeff!(x::CliffordOrderElem)
  x.even_coeffs = map(
                      y -> if sum(digits(y - 1; base=2, pad=rank(x.parent.lattice))) % 2 == 0
      x.coeffs[y]
    else
      x.parent.ambient_algebra.base_ring()
    end, 1:(x.parent.rank)
  )
  x.odd_coeffs = x.coeffs - x.even_coeffs
  return x
end

function _can_convert_coeff(coeff::Vector{S}, K::Field) where {S}
  if length(coeff) == 0
    return true
  end
  try
    K.(coeff)  # Try converting the coefficient vector to K
      true
  catch
      false
  end
end

#Returns the coefficient_ideals of the given lattice 'ls'
_coefficient_ideals_of_lattice(ls::QuadLat) = coefficient_ideals(ls)::Vector{<:NumFieldOrderFractionalIdeal}

_coefficient_ideals_of_lattice(ls::ZZLat) = fill(fractional_ideal(ZZ,1), rank(ls))::Vector{<:NumFieldOrderFractionalIdeal}
