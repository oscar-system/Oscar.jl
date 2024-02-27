@doc raw"""
    AbsAffineAlgebraicSet <: AbsAffineScheme

An affine, geometrically reduced subscheme of an affine space over a field.
"""
abstract type AbsAffineAlgebraicSet{BaseField<:Field, RingType} <:AbsAffineScheme{BaseField, RingType} end

@doc raw"""
    AbsProjectiveAlgebraicSet <: AbsProjectiveScheme

A projective, geometrically reduced scheme of finite type over a field.
"""
abstract type AbsProjectiveAlgebraicSet{BaseField<:Field, RingType} <:AbsProjectiveScheme{BaseField, RingType} end

################################################################################
#
# Abstract types for varieties
#
################################################################################

@doc raw"""
    AbsAffineVariety <: AbsAffineAlgebraicSet

An affine, geometrically integral subscheme of an affine space over a field.
"""
abstract type AbsAffineVariety{BaseField<:Field, RingType} <:AbsAffineAlgebraicSet{BaseField, RingType} end

@doc raw"""
    AbsProjectiveVariety <: AbsProjectiveAlgebraicSet

A geometrically integral subscheme of a projective space over a field.
"""
abstract type AbsProjectiveVariety{BaseField<:Field, GradedRingType<:Ring} <: AbsProjectiveAlgebraicSet{BaseField, GradedRingType} end


@doc raw"""
    AbsCoveredVariety <: Scheme

A separated, geometrically integral scheme of finite type over a field.
"""
abstract type AbsCoveredVariety{BaseField<:Field} <: AbsCoveredScheme{BaseField} end

################################################################################
#
# Abstract types for curves
#
################################################################################

@doc raw"""
    AbsAffineCurve <: AbsAffineAlgebraicSet

A curve in affine space.

An affine curve is an affine algebraic set of dimension one.
"""
abstract type AbsAffineCurve{BaseField<:Field, RingType<:Ring} <: AbsAffineAlgebraicSet{BaseField, RingType} end

@doc raw"""
    AbsProjectiveCurve <: AbsProjectiveVariety

A projective curve embedded in an ambient projective space.
"""
abstract type AbsProjectiveCurve{BaseField<:Field, RingType<:Ring} <: AbsProjectiveAlgebraicSet{BaseField, RingType} end

@doc raw"""
    AbsCoveredCurve

A curve represented in terms of a covering.

A curve is a reduced scheme of dimension 1 of finite type over a field.
"""
abstract type AbsCoveredCurve{BaseField} <: Scheme{BaseField} end


################################################################################
#
# Abstract types for surfaces
#
################################################################################

@doc raw"""
    AbsAffineSurface <: AbsAffineVariety

A surface in affine space.

An affine surface is an affine variety of dimension two.
"""
abstract type AbsAffineSurface{BaseField<:Field, RingType<:Ring} <: AbsAffineVariety{BaseField, RingType} end

@doc raw"""
    AbsProjectiveSurface <: AbsProjectiveVariety

A projective surface embedded in an ambient projective space.
"""
abstract type AbsProjectiveSurface{BaseField<:Field, RingType<:Ring} <: AbsProjectiveVariety{BaseField, RingType} end

@doc raw"""
    AbsCoveredSurface

A geometrically integral scheme of dimension 2 of finite type over a field
represented in terms of a covering.
"""
abstract type AbsCoveredSurface{BaseField<:Field} <: AbsCoveredVariety{BaseField} end

################################################################################
#
# Abstract types for rational points
#
################################################################################

@doc raw"""
    AbsRationalPointSet{P<:AbsAffineScheme,T<:Scheme}

Abstract type for the set of rational points of a scheme.

Let ``k`` be a ring, ``L`` a ``k`-algebra and ``X`` a scheme over ``k``. The set ``X(L)`` of ``L``-rational points of ``X`` is defined as the set of homomorphisms over ``Spec k`` from ``Spec L`` to ``X``, i.e. as ``X(L) = Hom_{Spec k}(Spec L, X)``.

We refer to ``L`` as the `coefficient_ring`, to
``Spec L`` as the `domain` and to ``X`` as the `codomain` of this
homset. For traditional reasons we keep the name
`ambient_scheme` for the codomain.
"""
abstract type AbsRationalPointSet{P<:AbsAffineScheme,T<:Scheme} end

abstract type AbsRationalPoint{RingElemType, ParentType} end

domain(p::AbsRationalPoint) = domain(parent(p))
codomain(p::AbsRationalPoint) = codomain(parent(p))
ambient_scheme(p::AbsRationalPoint) = codomain(p)
ambient_space(p::AbsRationalPoint) = ambient_space(ambient_scheme(p))
coefficient_ring(P::AbsRationalPoint) = coefficient_ring(parent(P))

@doc raw"""
    AbsAffineRationalPoint{CoefficientType, ParentType}

A rational point ``P`` of an affine scheme ``X``.
We refer to ``X`` as the parent of ``P``.

Let ``X \subseteq \mathbb{A}^n_k`` be an algebraic set or more generally a subscheme
defined by the ideal ``I = (f_1, \dots f_r) \subseteq k[x_1,\dots x_n]``.
A rational point ``p`` of ``X`` is a tuple ``p = (p_1, \dots , p_n) \in k^n`` such that
``f_1(p) = \dots = f_n(p) = 0``.
"""
abstract type AbsAffineRationalPoint{RingElemType, ParentType} <: AbsRationalPoint{RingElemType, ParentType} end

@doc raw"""
    AbsProjectiveRationalPoint

A rational point ``P`` of a projective scheme ``X``.
We refer to ``X`` as the parent of ``P``.

Let ``k`` be a field. A rational point is an element of ``\mathbb{P}^n(k) = k^{n+1} \setminus \{0\} / k^*`` where
two vectors ``v,w`` in ``k^{n+1} \setminus \{0\}`` are identified if ``v = \alpha w`` for a non-zero scalar ``\alpha \in k^*``.

Let ``X \subseteq \mathbb{P}^n_k`` be an algebraic set or more generally a closed subscheme defined by the homogeneous ideal ``I = (f_1, \dots f_r)``. Then a rational point of ``X`` is ``p \in \mathbb{P}^n(k)`` such that ``f_1(p) = \dots = f_n(p) = 0``.

This type includes points in weighted projective space.
"""
abstract type AbsProjectiveRationalPoint{RingElemType, ParentType} <: AbsRationalPoint{RingElemType, ParentType} end

abstract type AbsCoveredRationalPoint{RingElemType, ParentType} <: AbsRationalPoint{RingElemType, ParentType} end

