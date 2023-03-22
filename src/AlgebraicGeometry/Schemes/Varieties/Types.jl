################################################################################
#
# Abstract types for varieties
#
################################################################################

@doc Markdown.doc"""
    AbsAffineVariety <: AbsSpec{<:Field}

An affine, geometrically integral scheme of finite type over a field.
"""
abstract type AbsAffineVariety{BaseField<:Field, RingType} <:AbsSpec{BaseField, RingType}# where {BaseField <:Field, RingType}
end

@doc Markdown.doc"""
    AbsProjectiveVariety <: AbsProjectiveScheme{<:Field}

A projective variety over a field.

That is a projective, geometrically integral scheme over a field.
"""
abstract type AbsProjectiveVariety{BaseField<:Field, GradedRingType<:Ring} <: AbsProjectiveScheme{BaseField, GradedRingType} end

@doc Markdown.doc"""
    AbsCoveredVariety <: Scheme{<:Field}

A separated, geometrically integral scheme of finite type over a field.

Note that we allow reducible varieties.
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

A curve is a geometrically integral scheme of finite type over a field which is
of dimension one.
"""
abstract type AbsCoveredCurve{BaseField} <: Scheme{BaseField} end

#=
################################################################################
# Plane curves
################################################################################
@doc Markdown.doc"""
    ProjectivePlaneCurve

A plane projective curve embedded in an ambient projective space.
"""
mutable struct PlaneProjectiveCurve{BaseField<:Field, GradedRingType<:Ring} <: AbsProjectiveCurve{BaseField, GradedRingType}
  X::ProjectiveScheme{BaseField, GradedRingType}
end
=#
