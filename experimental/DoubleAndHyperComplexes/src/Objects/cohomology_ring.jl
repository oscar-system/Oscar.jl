mutable struct SimplicialCohomologyRing{T} <: NCRing
  C::SimplicialCoComplex{T}

  function SimplicialCohomologyRing(K::SimplicialComplex)
    return SimplicialCohomologyRing(SimplicialCoComplex(K))
  end
  function SimplicialCohomologyRing(C::SimplicialCoComplex{T}) where {T}
    return new{T}(C)
  end
end

mutable struct SimplicialCohomologyRingElem{T}
  parent::SimplicialCohomologyRing{T}
  homog_elem::Union{SubQuoModuleElem{T}, Nothing}
  homog_deg::Union{Int, Nothing}
  coeff::Dict{Int, SubquoModuleElem{T}}

  # Constructor for homogeneous elements
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRingElem{T}, 
      p::Int, v::SubquoModuleElem{T}
    ) where {T}
    return new{T}(A, v, p)
  end

  # Constructor for mixed elements
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRingElem{T}, 
      coeff::Dict{Int, SubquoModuleElem{T}}
    ) where {T}
    return new{T}(A, nothing, nothing, coeff)
  end

  # Constructor for the zero
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRingElem{T}
    ) where {T}
    return new{T}(A, nothing, nothing)
  end
end

