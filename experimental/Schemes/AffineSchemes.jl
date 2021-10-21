import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose, AbstractAlgebra.Generic.Frac
import Base: ∘
import Oscar: base_ring
import AbstractAlgebra: FPModule, FPModuleElem
import Base.copy

using Oscar.Multiindices
using Oscar.Misc

export AffineScheme, Spec, SpecPrincipalOpen
export base_ring, ambient_ring, defining_ideal, imgs_frac, pullback, pullback_from_parent, pullback_from_root, inclusion_in_parent, inclusion_in_root, set_name!, inverted_element, identity_morphism, denoms, inverses, divide_by_units
export affine_space, localize, subscheme

export AffSchMorphism
export domain, codomain, pullback

@Markdown.doc """
    Scheme{ BaseRingType <: Ring }

A scheme over a ring ``k`` of type `BaseRingType`.
"""
abstract type Scheme{ S <: Ring } end

@doc Markdown.doc"""
    AffineScheme{BaseRingType, RingType <: MPolyRing, RingElemType <: MPolyElem} <: Scheme{BaseRingType}

An affine scheme over a base ring ``k`` of type `BaseRingType`, given by a ring ``R/I`` with 
``R`` a polynomial ring of type `RingType` and elements of type `RingElemType`.
"""
abstract type AffineScheme{BaseRingType, RingType <: MPolyRing, RingElemType <: MPolyElem} <: Scheme{BaseRingType} end

@Markdown.doc """
    SchemeMorphism{BaseRingType <: Ring}

Morphism of Schemes over a ring ``k`` of type `BaseRingType`.
"""
abstract type SchemeMorphism{BaseRingType <: Ring} end

@doc Markdown.doc"""
    Spec{BaseRingType, RingType, RingElemType} <: AffineScheme{BaseRingType, RingType, RingElemType}

An affine scheme ``X = Spec R/I`` with ``R = k[x₁,…,xₙ]`` a free 
polynomial algebra of type `RingType` over a base ring ``k`` of type 
`BaseRingType` and ``I ⊂ R`` a finitely generated ideal 
with elements of type ``RingElemType``.
"""
mutable struct Spec{BaseRingType, RingType, RingElemType} <: AffineScheme{BaseRingType, RingType, RingElemType}
  # the basic fields 
  k::BaseRingType		# the base ring (usually a field) of definition for the scheme
  R::RingType 			# the ambient polynomial ring to model this affine scheme
  I::MPolyIdeal{RingElemType}	# The ideal in R defining the scheme

  # fields for caching
  name::String # the name of this scheme for printing

  function Spec(k::BaseRingType, R::RingType, I::MPolyIdeal{RingElemType} ) where{
			BaseRingType <: Ring, RingType <:MPolyRing , RingElemType <: MPolyElem}
    k == coefficient_ring(R) || error("Base ring of the affine scheme does not coincide with the base ring of the associated algebra")
    base_ring(I) == R || error("Ideal does not belong to the ring")
    return new{BaseRingType, RingType, RingElemType}(k, R, I )
  end
end

###################################################################################
# Getter functions
#
# No caching is needed in this case, since all these variables need to be assigned 
# at instantiation
#
@Markdown.doc """
    base_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}

Returns the base ring ``k`` over which the affine scheme ``A`` is defined.
"""
function base_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return A.k
end

@Markdown.doc """
    ambient_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}

Returns the "ambient ring" ``R`` of the affine scheme ``A``, where ``A = Spec R/I``.
"""
function ambient_ring(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return A.R
end

@Markdown.doc """
     defining_ideal(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
Returns the "defining ideal" ``I`` of the affine scheme ``A``, where ``A = Spec R/I``.
"""
function defining_ideal(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return A.I
end

@Markdown.doc """
    denoms(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}

This function simply extends the `denoms(D::SpecPrincipalOpen)` so that it 
can conveniently be used on all of `AffineScheme`. It returns an empty list.
"""
function denoms(A::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType}
  return typeof(A.R)[]
end


@Markdown.doc """
    Spec( X::Spec{BaseRingType, RingType, RingElemType} ) where {BaseRingType, RingType, RingElemType}

Returns a copy of the affine scheme ``X``.
"""
Spec(X::Spec{BaseRingType, RingType, RingElemType}) where {BaseRingType, RingType, RingElemType} = Spec(base_ring(X), ambient_ring(X), deepcopy(defining_ideal(X)))

# the routine for deepcopies internally used by julia. 
# Note that rings don't need to be deepcopied. Ideals on 
# the other hand, posess information that can be altered. 
Base.deepcopy_internal(X::Spec{BaseRingType, RingType, RingElemType}, dict::IdDict) where {BaseRingType, RingType, RingElemType} = Spec(X)

@Markdown.doc """
    subscheme(X::Spec{BaseRingType, RingType, RingElemType}, J::MPolyIdeal{RingElemType}) where {BaseRingType, RingType, RingElemType}

Returns the subscheme of ``X = Spec R/I`` defined by the ideal ``I + J`` in the ring ``R``.
"""
function subscheme(X::Spec{BaseRingType, RingType, RingElemType}, J::MPolyIdeal{RingElemType}) where {BaseRingType, RingType, RingElemType}
  base_ring(J) == ambient_ring(X) || error( "Ideal does not belong to the ambient ring of the variety" )
  return Spec(base_ring(X), ambient_ring(X), defining_ideal(X) + J)
end

@Markdown.doc """
    subscheme(X::Spec{BaseRingType, RingType, RingElemType}, f::RingElemType) where {BaseRingType, RingType, RingElemType}

Returns the subscheme of ``X = Spec R/I`` defined by the ideal`` I + ⟨f⟩`` in the ring ``R``.
"""
function subscheme(X::Spec{BaseRingType, RingType, RingElemType}, f::RingElemType) where {BaseRingType, RingType, RingElemType}
  return subscheme(X, ideal(ambient_ring(X), f))
end

