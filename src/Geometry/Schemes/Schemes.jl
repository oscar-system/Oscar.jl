import AbstractAlgebra.Ring, Oscar.AlgHom, Oscar.compose, AbstractAlgebra.Generic.Frac
import Base: ∘
import Oscar: base_ring
import AbstractAlgebra: FPModule, FPModuleElem
import Base.copy

include("./Misc.jl")
using Oscar.Misc

include("./Multiindices.jl")
using Oscar.Multiindices

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

include("./AffineSchemes.jl")
