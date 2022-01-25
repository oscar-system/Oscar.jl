mutable struct CoveredScheme{CoveringType, BRT, BRET} <: Scheme{BRT, BRET}
  coverings::Vector{CoveringType}
end

mutable struct EffectiveWeilDivisor{SpecType<:Spec, CoveredSchemeType<:CoveredScheme, CoveringType<:Covering}
  X::CoveredSchemeType # parent
  C::CoveringType # the reference covering
  D::Dict{SpecType, SpecType} # keys: patches of X in C; values: the subscheme given by D in that patch

  function EffectiveWeilDivisor(X::CoveredSchemeType, C::CoveringType, D::Dict{SpecType, SpecType}; check::Bool=true) where {SpecType<:Spec, CoveredSchemeType<:CoveredScheme, CoveringType<:Covering}
    C in coverings(X) || error("covering can not be found")
    for X in keys(D)
      X in patches(C) || error("patch can not be found in the given covering")
      if check
        Z = C[X]
        is_closed_embedding(X, Z) || error("the divisor on the patch is not a closed subscheme")
      end
    end
    return new{SpecType, CoveredSchemeType, CoveringType}(X, D, C)
  end
end

parent(D::EffectiveWeilDivisor) = D.X
covering(D::EffectiveWeilDivisor) = D.C
getindex(D::EffectiveWeilDivisor, X::Spec) = (D.D)[X]

parent_type(D::EffectiveWeilDivisor{S, T, U}) where{S, T, U} = T
parent_type(::Type{EffectiveWeilDivisor{S, T, U}}) where{S, T, U} = T
covering_type(D::EffectiveWeilDivisor{S, T, U}) where{S, T, U} = U
covering_type(::Type{EffectiveWeilDivisor{S, T, U}}) where{S, T, U} = U
affine_patch_type(D::EffectiveWeilDivisor{S, T, U}) where{S, T, U} = S
affine_patch_type(::Type{EffectiveWeilDivisor{S, T, U}}) where{S, T, U} = S

# TODO: Implement a reliable hash-function that allows an effective 
# Weil divisor to be used in a dictionary.

mutable struct WeilDivisor{CoveredSchemeType<:CoveredScheme, EffectiveWeilDivisorType<:EffectiveWeilDivisor, CoefficientRingElemType<:AbstractAlgebra.Ring}
  X::CoveredSchemeType # the parent
  coefficients::Dict{EffectiveWeilDivisorType, CoefficientRingElemType} # the formal linear combination

  function(X::CoveredSchemeType, coefficients::Dict{EffectiveWeilDivisorType, CoefficientRingElemType}) where {CoveredSchemeType, EffectiveWeilDivisorType, CoefficientRingElemType}
    # TODO: Do we want to require that the different effective divisors 
    # have the same underlying covering? Probably not.
    for D in keys(coefficients)
      parent(D) == X || error("component of divisor does not lay in the given scheme")
    end
    return new{CoveredSchemeType, EffectiveWeilDivisorType, CoefficientRingElemType}(X, D)
  end
end

parent_type(D::WeilDivisor{S, T, U}) where{S, T, U} = T
parent_type(::Type{EffectiveWeilDivisor{S, T, U}}) where{S, T, U} = T

parent(D::WeilDivisor) = D.X
getindex(D::WeilDivisor, E::EffectiveWeilDivisor) = (D.coefficients)[E]
components(D::WeilDivisor) = keys(D.coefficients)



