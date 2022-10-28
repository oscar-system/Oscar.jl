export OO, ambient_ring, base_ring, defining_ideal, dim, codim

export defining_ideal, strict_modulus

export ring_type, base_ring_type, base_ring_elem_type, poly_type, ring_type



########################################################################
# (1) Attributes of AbsSpec
########################################################################

# Here is the inferface for AbsSpec

@Markdown.doc """
    OO(X::AbsSpec) 

On an affine scheme ``X = Spec(R)`` this returns the ring `R`.
"""
function OO(X::AbsSpec{BRT, RT}) where {BRT, RT} 
  OO(underlying_scheme(X))::RT
end


@Markdown.doc """
    ambient_ring(X::AbsSpec)

On an affine scheme ``X = Spec(R)`` over ``𝕜`` this returns a 
polynomial ring ``P = 𝕜[x₁,…,xₙ]`` with natural coercion 
``P → R`` and the property that for every other (commutative) 
ring ``S`` and any homomorphism ``φ : R → S`` there is a morphism 
``ψ : P → S`` factoring through ``φ`` and such that ``φ`` 
is uniquely determined by ``ψ``.
"""
function ambient_ring(X::AbsSpec)
  return ambient_ring(underlying_scheme(X))::MPolyRing
end


@Markdown.doc """
    base_ring(X::AbsSpec)

On an affine scheme ``X/𝕜`` over ``𝕜`` this returns the ring `𝕜`.
"""
function base_ring(X::AbsSpec{BRT, RT}) where {BRT, RT}
  return base_ring(underlying_scheme(X))::BRT
end


@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  return dim(saturated_ideal(modulus(OO(X))))
end


@attr function dim(X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  # the following line is supposed to refer the problem to the
  # algebra side of the problem
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))]))
end


@attr function dim(X::AbsSpec{<:Ring, <:MPolyRing})
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))]))
end


@attr function dim(X::AbsSpec{<:Ring, <:MPolyQuo})
  return dim(modulus(OO(X)))
end


@attr function codim(X::Spec)
  return dim(ideal(ambient_ring(X), [zero(ambient_ring(X))])) - dim(X)
end


@attr String function name(X::Spec)
  return "unnamed affine variety"
end



########################################################################
# (2) Further attributes for AbsSpec
########################################################################

# TODO: Needed?

@attr defining_ideal(X::AbsSpec{<:Any, <:MPolyRing}) = ideal(OO(X), [zero(OO(X))])
defining_ideal(X::AbsSpec{<:Any, <:MPolyQuo}) = modulus(OO(X))

@attr defining_ideal(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = ideal(OO(X), [zero(OO(X))])
defining_ideal(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = modulus(OO(X))

strict_modulus(X::AbsSpec) = saturated_ideal(modulus(OO(X)))



########################################################################
# (3) Further attributes for Spec
########################################################################

strict_modulus(X::Spec) = saturated_ideal(modulus(OO(X)))



########################################################################
# (4) Implementation of the AbsSpec interface for the basic Spec           #
########################################################################

OO(X::Spec) = X.OO
base_ring(X::Spec) = X.kk
ambient_ring(X::Spec{<:Any, <:MPolyRing}) = OO(X)
ambient_ring(X::Spec{<:Any, <:MPolyQuo}) = base_ring(OO(X))
ambient_ring(X::Spec{<:Any, <:MPolyLocalizedRing}) = base_ring(OO(X))
ambient_ring(X::Spec{<:Any, <:MPolyQuoLocalizedRing}) = base_ring(OO(X))
ambient_ring(X::Spec{T, T}) where {T<:Field} = base_ring(X)



########################################################################
# (5) Type getters
########################################################################

# TODO: Needed?

ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = RT
ring_type(X::AbsSpec) = ring_type(typeof(X))

base_ring_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = BRT
base_ring_type(X::AbsSpec) = base_ring_type(typeof(X))
base_ring_elem_type(::Type{SpecType}) where {BRT, RT, SpecType<:AbsSpec{BRT, RT}} = elem_type(BRT)
base_ring_elem_type(X::AbsSpec) = base_ring_elem_type(typeof(X))

poly_type(::Type{SpecType}) where {BRT, RT<:MPolyRing, SpecType<:AbsSpec{BRT, RT}} = elem_type(RT)
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyQuo{T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyLocalizedRing{<:Any, <:Any, <:Any, T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(::Type{SpecType}) where {BRT, T, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, T}, SpecType<:AbsSpec{BRT, RT}} = T
poly_type(X::AbsSpec) = poly_type(typeof(X))

ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = RT
ring_type(X::Spec) = ring_type(typeof(X))
base_ring_type(::Type{Spec{BRT, RT}}) where {BRT, RT} = BRT
base_ring_type(X::Spec) = base_ring_type(typeof(X))
