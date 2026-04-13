
################################################################################
#
#  Datatypes
#
################################################################################

### Algebra ###
# Data structure for Clifford algebras. The type variable 'T' represents the element type
# of the base ring, i.e. it is usually QQFieldElem or AbsSimpleNumFieldElem.
# The type variable 'S' represents the type of the Gram matrix of the underlying quadratic space,
# i.e. it is usually QQMatrix or AbstactAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}.
mutable struct CliffordAlgebra{T, S} <: Hecke.AbstractAssociativeAlgebra{T}
  base_ring::Ring
  space::Hecke.QuadSpace{K, S} where {K} # K = parent_type(T), e.g. T is of type AbsSimpleNumFieldElem and K is of type AbsSimpleNumField
  gram::S
  dim::Int
  orth_elt::Any #elem_type{C}
  disq::T

  #Return the Clifford algebra of the quadratic space 'qs' 
  function CliffordAlgebra{T,S}(qs::Hecke.QuadSpace{K,S}) where {T,S,K}
    gram = gram_matrix(qs)
    return new{T,S}(base_ring(qs), qs, gram, 2^ncols(gram))
  end
end

### Elements ###
# Data structure for the elements of a Clifford algebra. The type variables serve the
# same purpose as they do for Clifford algebras.
mutable struct CliffordAlgebraElem{T,S} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::CliffordAlgebra{T, S}
  coeffs::Vector{T}

  #Return the 0-element of the Clifford algebra C
  CliffordAlgebraElem{T,S}(C::CliffordAlgebra{T,S}) where {T,S} = new{T,S}(C, [C.base_ring() for _ in 1:C.dim])

  CliffordAlgebraElem(C::CliffordAlgebra) =
    CliffordAlgebraElem{elem_type(C.base_ring), typeof(C.gram)}(C)

  #Return the element in the Clifford algebra C with coefficient vector coeff wrt. the canonical basis
  function CliffordAlgebraElem{T,S}(C::CliffordAlgebra{T,S}, coeffs::Vector{R}) where {T,S,R<:FieldElem}
    @req length(coeffs) == C.dim "invalid length of coefficient vector"
    return new{T,S}(C, coeffs)
  end

  CliffordAlgebraElem(C::CliffordAlgebra{T,S}, coeffs::Vector{T}) where {T,S} =
    CliffordAlgebraElem{elem_type(C.base_ring), typeof(C.gram)}(C, coeffs)
  
  function CliffordAlgebraElem(C::CliffordAlgebra{T,S}, coeff::Vector{R}) where {T,S,R}
    K = C.base_ring
    return CliffordAlgebraElem{elem_type(K), typeof(C.gram)}(C, K.(coeff))
  end
end

# Data structure for Clifford orders over rings distinct from the integers. The type variable 'T' represents the element type
# of the base ring, i.e. it is usually AbsSimpleNumFieldOrderElem. The type variable 'C' represents the type of its ambient algebra,
# i.e. it is usually CliffordAlgebra{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}}.
mutable struct CliffordOrder{T, C} <: Hecke.AbstractAssociativeAlgebra{T}

  base_ring::Ring
  ambient_algebra::C
  rank::Int
  lattice::QuadLat
  gram::MatElem

  # In the 3 lines below let CO be an instance of CliffordOrder
  coefficient_ideals::Any # Vector{Hecke.fractional_ideal_type(base_ring_type(C))}
  disq::Any # Tuple{Hecke.fractional_ideal_type(base_ring_type(CO)), elem_type(base_ring(ambient_algebra(CO)))}
  _raw_orth_data::Any #Always of type Tuple{elem_type(ambient_algebra(CO)), Hecke.fractional_ideal_type(base_ring_type(CO))}

  function CliffordOrder{T, C}(ls::QuadLat{S, M}) where {T, C, S<:NumField, M<:MatElem}
    !is_zero(rank(ls)) && @req is_integral(fractional_ideal(base_ring(ls), base_field(ls)(1//2)) * norm(ls)) "The given lattice is not even!"
    qs = rational_span(ls)
    coeff_ids = _coefficient_ideals_of_CO(ls)
    return new{T, C}(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs), coeff_ids)
  end
end

mutable struct ZZCliffordOrder <: Hecke.AbstractAssociativeAlgebra{ZZRingElem}

  base_ring::ZZRing
  ambient_algebra::CliffordAlgebra{QQFieldElem, QQMatrix}
  rank::Int
  lattice::ZZLat
  gram::QQMatrix
  disq::ZZRingElem
  max_orth_elt::Any #Always of type ZZCliffordOrderElem

  function ZZCliffordOrder(ls::ZZLat)
    !is_zero(rank(ls)) && @req is_even(ls) "The given lattice is not even!"
    qs = rational_span(ls) 
    return new(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs))
  end

end

##### Elements #####
# Data structure for elements of Clifford orders over rings distinct from the integers. The type
# variables serve the same purpose as they do for Clifford orders.
mutable struct CliffordOrderElem{T, C, S} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::CliffordOrder{T, C}
  coeffs::Vector{S}

  #Return the 0-element of the Clifford order C
  CliffordOrderElem{T, C, S}(CO::CliffordOrder{T, C}) where {T, C, S} =
    new{T, C, S}(CO, [CO.ambient_algebra.base_ring() for _ in 1:CO.rank])

  function CliffordOrderElem(CO::CliffordOrder)
    CT = typeof(CO.ambient_algebra)
    return CliffordOrderElem{elem_type(CO.base_ring), CT, elem_type(base_ring_type(CT))}(CO)
  end

  #Return the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function CliffordOrderElem{T, C, S}(CO::CliffordOrder{T, C}, coeffs::Vector{S}) where {T, C, S}
    @req length(coeffs) == CO.rank "invalid length of coefficient vector"
    
    for i in 1:CO.rank
      ci = coeffs[i]
      is_zero(ci) || @req ci in coefficient_ideals(CO)[i] "The element does not lie in the Clifford order."
    end
    
    return new{T, C, S}(CO, coeffs)
  end

  function CliffordOrderElem(CO::CliffordOrder{T, C}, coeffs::Vector{S}) where {T, C, S}
    ambalg = CO.ambient_algebra
    K = base_ring(ambalg)
    return CliffordOrderElem{elem_type(CO.base_ring), typeof(ambalg), elem_type(K)}(CO, K.(coeffs))
  end

end

##### Elements #####
mutable struct ZZCliffordOrderElem <: Hecke.AbstractAssociativeAlgebraElem{ZZRingElem}
  parent::ZZCliffordOrder
  coeffs::Vector{QQFieldElem}

  #Return the 0-element of the Clifford order CO
  ZZCliffordOrderElem(CO::ZZCliffordOrder) = new(CO, [QQ() for _ in 1:CO.rank])

  #Return the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeffs::Vector{QQFieldElem})
    @req length(coeffs) == CO.rank "invalid length of coefficient vector"
    for i in 1:CO.rank
      @req is_integer(coeffs[i]) "The element does not lie in the Clifford order."
    end
    return new(CO, coeffs)
  end

  ZZCliffordOrderElem(CO::ZZCliffordOrder, coeffs::Vector{S}) where {S} = ZZCliffordOrderElem(CO, QQ.(coeffs))

end
