export AffineCurve
export underlying_scheme
export ambient_space

mutable struct AffineCurve{BRT, BRET, SpecType<:Spec, PolyType<:MPolyElem} <: AbsSpec{BRT, BRET}
  C::SpecType # the underlying scheme
  A::SpecType # the ambient space
  f::Vector{PolyType} # the defining equations

  # Construct the affine curve in Spec(R) which is the zero locus of fⱼ in the hypersurface 
  # complement of Πᵢ gᵢ.
  function AffineCurve(R::MPolyRing, f::Vector{PolyType}, g::Vector{PolyType}) where {PolyType<:MPolyElem}
    all(x->(parent(x)==R), f) || error("polynomials do not belong to the correct ring")
    all(x->(parent(x)==R), g) || error("polynomials do not belong to the correct ring")
    A = Spec(R)
    I = ideal(R, f)
    U = MPolyPowersOfElement(R, g)
    L = MPolyQuoLocalizedRing(R, I, U)
    C = Spec(L)
    dim(C) == 1 || error("the given generators do not define a subscheme of dimension 1")
    return new{typeof(coefficient_ring(R)), elem_type(coefficient_ring(R)), typeof(C), PolyType}(C, A, f)
  end
end

### overwriting the required getters
function OO(C::AffineCurve) 
  return OO(C.C)
end

### overwriting the required type-getters
ring_type(::Type{CurveType}) where {BRT, BRET, SpecType, PolyType, CurveType<:AffineCurve{BRT, BRET, SpecType, PolyType}} = ring_type(SpecType)

### additional getters
underlying_scheme(C::AffineCurve) = C.C
ambient_space(C::AffineCurve) = C.A

### overwriting further generic methods 
function hypersurface_complement(C::AffineCurve, f::MPolyElem) 
  parent(f) == base_ring(OO(C)) || error("polynomial does not belong to the correct ring")
  return AffineCurve(base_ring(OO(C)), gens(modulus(OO(C))), push!(denominators(inverted_set(OO(C))), f))
end

