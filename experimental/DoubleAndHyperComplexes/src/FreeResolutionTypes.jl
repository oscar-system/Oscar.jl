########################################################################
# FreeResolution -- type declarations
########################################################################
#
@doc raw"""
    FreeResolution{T}

Data structure for free resolutions.
"""
mutable struct FreeResolution{ChainType, MorphismType} <: AbsHyperComplex{ChainType, MorphismType}
    C::Hecke.ComplexOfMorphisms
    underlying_complex::SimpleComplexWrapper{ChainType, MorphismType}

    function FreeResolution(C::Hecke.ComplexOfMorphisms{T}) where {T}
        FR = new{T, morphism_type(T)}()
        FR.C = C

        return FR
    end
end

Base.getindex(FR::FreeResolution, i::Int) = FR.C[i]

function Base.show(io::IO, FR::FreeResolution)
    C = FR.C
    show(io, C)
end

function underlying_complex(res::FreeResolution)
  if !isdefined(res, :underlying_complex)
    res.underlying_complex = SimpleComplexWrapper(res.C)
  end
  return res.underlying_complex
end

