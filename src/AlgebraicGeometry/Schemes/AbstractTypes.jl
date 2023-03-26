@doc Markdown.doc"""
    AbsAffineAlgebraicSet <: AbsSpec

An affine, geometrically reduced subscheme of an affine space over a field.
"""
abstract type AbsAffineAlgebraicSet{BaseField<:Field, RingType} <:AbsSpec{BaseField, RingType} end

@doc Markdown.doc"""
    AbsProjectiveAlgebraicSet <: AbsProjectiveScheme

A projective, geometrically reduced scheme of finite type over a field.
"""
abstract type AbsProjectiveAlgebraicSet{BaseField<:Field, RingType} <:AbsProjectiveScheme{BaseField, RingType} end

################################################################################
#
# Abstract types for varieties
#
################################################################################

@doc Markdown.doc"""
    AbsAffineVariety <: AbsAffineAlgebraicSet

An affine, geometrically integral subscheme of an affine space over a field.
"""
abstract type AbsAffineVariety{BaseField<:Field, RingType} <:AbsAffineAlgebraicSet{BaseField, RingType} end

@doc Markdown.doc"""
    AbsProjectiveVariety <: AbsProjectiveAlgebraicSet

A projective variety over a field.

That is a geometrically integral subscheme of a projective space over a field.
"""
abstract type AbsProjectiveVariety{BaseField<:Field, GradedRingType<:Ring} <: AbsProjectiveAlgebraicSet{BaseField, GradedRingType} end

@doc Markdown.doc"""
    AbsCoveredVariety <: Scheme

A separated, geometrically integral scheme of finite type over a field.
"""
abstract type AbsCoveredVariety{BaseField<:Field} <: AbsCoveredScheme{BaseField} end

################################################################################
#
# Abstract types for curves
#
################################################################################

@doc Markdown.doc"""
    AbsAffineCurve <: AbsAffineVariety

A curve in affine space.

An affine curve is an affine variety of dimension one.
"""
abstract type AbsAffineCurve{BaseField<:Field, RingType<:Ring} <: AbsAffineVariety{BaseField, RingType} end

@doc Markdown.doc"""
    AbsProjectiveCurve <: AbsProjectiveVariety

A projective curve embedded in an ambient projective space.
"""
abstract type AbsProjectiveCurve{BaseField<:Field, RingType<:Ring} <: AbsProjectiveVariety{BaseField, RingType} end

@doc Markdown.doc"""
    AbsCoveredCurve

A curve represented in terms of a covering.

A curve is a geometrically integral scheme of dimension 1 of finite type over a field.
"""
abstract type AbsCoveredCurve{BaseField} <: Scheme{BaseField} end

