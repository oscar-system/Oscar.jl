
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
mutable struct CliffordAlgebra{T,S} <: Hecke.AbstractAssociativeAlgebra{T}
  base_ring::Ring
  space::Hecke.QuadSpace{K,S} where {K} # K = parent_type(T), e.g. T is of type AbsSimpleNumFieldElem and K is of type AbsSimpleNumField
  gram::S
  dim::Int
  basis_of_centroid::Any # Vector{elem_type(C)}, with C an instance of CliffordAlgebra
  disq::T
  basis_of_center::Any # Vector{elem_type(C)}, with C an instance of CliffordAlgebra

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
  CliffordAlgebraElem{T,S}(C::CliffordAlgebra{T,S}) where {T,S} = new{T,S}(C, fill(C.base_ring(), C.dim))

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
    @req __can_convert_coefficients(coeff, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordAlgebraElem{elem_type(C.base_ring), typeof(C.gram)}(C, K.(coeff))
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

  # In the 4 lines below let CO be an instance of CliffordOrder
  coefficient_ideals::Any # Vector{Hecke.fractional_ideal_type(base_ring_type(C))}
  pseudo_basis_of_centroid::Any # Vector{Tuple{elem_type(CO), Hecke.fractional_ideal_type(base_ring_type(CO))}}
  disq::Any # Tuple{Hecke.fractional_ideal_type(base_ring_type(CO)), elem_type(base_ring(ambient_algebra(CO)))}
  pseudo_basis_of_center::Any # Vector{Tuple{elem_type(CO), Hecke.fractional_ideal_type(base_ring_type(CO))}}

  function CliffordOrder{T, C}(ls::QuadLat{S, M}) where {T, C, S<:NumField, M<:MatElem}
    if !is_zero(rank(ls))
      @req is_integral(fractional_ideal(base_ring(ls), base_field(ls)(1//2)) * norm(ls)) "The given lattice is not even!"
    end
    qs = rational_span(ls)
    coeff_ids = _set_coefficient_ideals!(ls)
    return new{T, C}(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs), coeff_ids)
  end
end

mutable struct ZZCliffordOrder <: Hecke.AbstractAssociativeAlgebra{ZZRingElem}

  base_ring::ZZRing
  ambient_algebra::CliffordAlgebra{QQFieldElem, QQMatrix}
  rank::Int
  lattice::ZZLat
  gram::QQMatrix
  basis_of_centroid::Any #Always of type Vector{ZZCliffordOrderElem} 
  disq::ZZRingElem
  basis_of_center::Any #Always of type Vector{ZZCliffordOrderElem}

  function ZZCliffordOrder(ls::ZZLat)
    if !is_zero(rank(ls))
      @req is_even(ls) "The given lattice is not even!"
    end
    qs = rational_span(ls) 
    return new(base_ring(ls), clifford_algebra(qs), 2^rank(ls), ls, gram_matrix(qs))
  end

end

##### Elements #####
# Data structure for elements of Clifford orders over rings distinct from the integers. The type
# variables serve the same purpose as they do for Clifford orders.
mutable struct CliffordOrderElem{T, C} <: Hecke.AbstractAssociativeAlgebraElem{T}
  parent::CliffordOrder{T, C}
  coeffs::Any 

  #Return the 0-element of the Clifford order C
  CliffordOrderElem{T, C}(CO::CliffordOrder{T, C}) where {T, C} =
    new{T, C}(CO, fill(CO.ambient_algebra.base_ring(), CO.rank))

  CliffordOrderElem(CO::CliffordOrder) = CliffordOrderElem{elem_type(CO.base_ring), typeof(CO.ambient_algebra)}(CO)

  #Return the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function CliffordOrderElem{T, C}(CO::CliffordOrder{T, C}, coeffs::Vector{S}) where {T, C, S<:NumFieldElem}
    @req length(coeffs) == CO.rank "invalid length of coefficient vector"
    
    for i in 1:CO.rank
      ci = coeffs[i]
      is_zero(ci) || @req ci in coefficient_ideals(CO)[i] "The element does not lie in the Clifford order."
    end
    
    return new{T, C}(CO, coeffs)
  end

  function CliffordOrderElem(CO::CliffordOrder{T, C}, coeffs::Vector{S}) where {T, C, S}
    K = base_ring(ambient_algebra(CO))
    @req _can_convert_coefficients(coeffs, K) "entries of coefficient vector are not contained in $(K)"
    return CliffordOrderElem{elem_type(base_ring(CO)), typeof(ambient_algebra(CO))}(CO, K.(coeffs))
  end

end

##### Elements #####
mutable struct ZZCliffordOrderElem <: Hecke.AbstractAssociativeAlgebraElem{ZZRingElem}
  parent::ZZCliffordOrder
  coeffs::Vector{QQFieldElem}

  #Return the 0-element of the Clifford order CO
  ZZCliffordOrderElem(CO::ZZCliffordOrder) = new(CO, fill(QQ(), CO.rank))

  #Return the element in the Clifford order CO with coefficient vector coeff with respect to the canonical basis
  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeffs::Vector{QQFieldElem})
    @req length(coeffs) == CO.rank "invalid length of coefficient vector"
    for i in 1:CO.rank
      @req is_integer(coeffs[i]) "The element does not lie in the Clifford order."
    end
    return new(CO, coeffs)
  end

  function ZZCliffordOrderElem(CO::ZZCliffordOrder, coeffs::Vector{S}) where {S}
    @req _can_convert_coefficients(coeffs, QQ) "entries of coefficient vector are not contained in $(QQField)"
    return ZZCliffordOrderElem(CO, QQ.(coeffs))
  end

end
