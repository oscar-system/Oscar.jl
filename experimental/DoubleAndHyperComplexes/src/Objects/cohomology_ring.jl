mutable struct SimplicialCohomologyRing{T} <: NCRing
  C::SimplicialCoComplex

  function SimplicialCohomologyRing(K::SimplicialComplex)
    return SimplicialCohomologyRing(SimplicialCoComplex(K))
  end
  function SimplicialCohomologyRing(C::SimplicialCoComplex)
    T = elem_type(base_ring(C))
    return new{T}(C)
  end
end

simplicial_co_complex(A::SimplicialCohomologyRing) = A.C

mutable struct SimplicialCohomologyRingElem{T}
  parent::SimplicialCohomologyRing{T}
  homog_elem::Union{SubquoModuleElem{T}, Nothing}
  homog_deg::Union{Int, Nothing}
  coeff::Dict{Int, SubquoModuleElem{T}}

  # Constructor for homogeneous elements
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRing{T}, 
      p::Int, v::SubquoModuleElem{T}
    ) where {T}
    return new{T}(A, v, p)
  end

  # Constructor for mixed elements
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRing{T}, 
      coeff::Dict{Int, SubquoModuleElem{T}}
    ) where {T}
    return new{T}(A, nothing, nothing, coeff)
  end

  # Constructor for the zero
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRing{T}
    ) where {T}
    return new{T}(A, nothing, nothing)
  end
end

parent(a::SimplicialCohomologyRingElem) = a.parent
zero(A::SimplicialCohomologyRing) = SimplicialCohomologyRingElem(A)
zero(a::SimplicialCohomologyRingElem) = zero(parent(a))

function one(A::SimplicialCohomologyRing)
  SimplicialCohomologyRingElem(A, 0, sum(gens(homology(simplicial_co_complex(A), 0)[1])))
end

function deepcopy_internal(d::IdDict, a::SimplicialCohomologyRingElem)
  result = parent(a)()
  if !isnothing(a.homog_elem)
    result.homog_elem = deepcopy_internal(d, a.homog_elem)
    result.homog_deg = copy(a.homog_deg)
    return result
  end
  result.coeff = deepcopy_internal(d, a.coeff)
  return result
end

function +(a::SimplicialCohomologyRingElem, b::SimplicialCohomologyRingElem)
  @assert parent(a) === parent(b) "parent mismatch"
  is_zero(a) && return copy(b)
  is_zero(b) && return copy(a)
  result = parent(a)() # unassigned zero element
  if isnothing(a.homog_elem)
    result.coeff = copy(a.coeff) # the result will not be homogeneously stored 
    !isdefined(a.coeff) && return copy(b) # a was the zero element in the end
    if isnothing(b.homog_elem)
      for (q, v) in b.coeff
        w = get(result.coeff, q, nothing)
        if isnothing(w)
          result.coeff[q] = v
        else
          res = w + v
          if is_zero(res)
            delete!(result.coeff, q)
          else
            result.coeff[q] = res
          end
        end
      end
      return result
    else # b is homogeneous
      q = b.homog_deg
      w = get(result.coeff, q, nothing)
      if isnothing(w)
        result.coeff[q] = b.homog_elem
      else
        res = w + b.homog_elem
        if is_zero(res)
          delete!(result.coeff, q)
        else
          result.coeff[q] = res
        end
      end
    end
  else # a is homogeneous
    q = a.homog_deg
    if isnothing(b.homog_elem) # b is not homogeneous
      result.coeff = copy(b.coeff) # result will not be homogeneously stored
      w = get(result.coeff, q, nothing)
      if isnothing(w)
        result.coeff[q] = a.homog_elem
      else
        res = w + a.homog_elem
        if iszero(res)
          delete!(result.coeff, q)
        else
          result.coeff[q] = res
        end
      end
      return result
    else # b is homogeneous
      p = b.homog_deg
      if p == q
        result.homog_deg = p
        result.homog_elem = a.homog_elem + b.homog_elem
        return result
      else
        result.coeff = Dict{Int, SubquoModuleElem{T}}([q => a.homog_elem, p => b.homog_elem])
        return result
      end
    end
  end
end



