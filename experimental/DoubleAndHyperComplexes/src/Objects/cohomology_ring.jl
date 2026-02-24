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

mutable struct SimplicialCohomologyRingElem{T} <: NCRingElem
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
(A::SimplicialCohomologyRing)() = zero(A)

function one(A::SimplicialCohomologyRing)
  SimplicialCohomologyRingElem(A, 0, sum(gens(homology(simplicial_co_complex(A), 0)[1])))
end

function is_zero(a::SimplicialCohomologyRingElem)
  if isnothing(a.homog_elem)
    !isdefined(a, :coeff) && return true
    isempty(a.coeff) && return true
    return all(iszero(b) for (_, b) in a.coeff)
  end
  return is_zero(a.homog_elem)
end

function deepcopy_internal(a::SimplicialCohomologyRingElem, d::IdDict)
  result = parent(a)()
  if !isnothing(a.homog_elem)
    result.homog_elem = deepcopy_internal(a.homog_elem, d)
    result.homog_deg = deepcopy(a.homog_deg)
    return result
  end
  result.coeff = deepcopy_internal(a.coeff, d)
  return result
end

function +(a::SimplicialCohomologyRingElem{T}, b::SimplicialCohomologyRingElem{T}) where {T}
  @assert parent(a) === parent(b) "parent mismatch"
  is_zero(a) && return deepcopy(b)
  is_zero(b) && return deepcopy(a)
  result = parent(a)() # unassigned zero element
  if isnothing(a.homog_elem)
    result.coeff = deepcopy(a.coeff) # the result will not be homogeneously stored 
    !isdefined(a, :coeff) && return deepcopy(b) # a was the zero element in the end
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
      return result
    end
  else # a is homogeneous
    q = a.homog_deg
    if isnothing(b.homog_elem) # b is not homogeneous
      result.coeff = deepcopy(b.coeff) # result will not be homogeneously stored
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

is_homogeneous_normalized(a::SimplicialCohomologyRingElem) = !isnothing(a.homog_elem)

# extract homogeneous parts
homogeneous_parts(a::SimplicialCohomologyRingElem) = Set(isnothing(a.homog_elem) ? [SimplicialCohomologyRingElem(parent(a), i, m) for (i,m) in pairs(a.coeff)] : [a])

# distribute over homogeneous parts
*(a::SimplicialCohomologyRingElem, b::SimplicialCohomologyRingElem) = (
   isnothing(a.homog_elem) && isnothing(a.homog_elem) ? sum(mul_homog(ah, bh) for ah in homogeneous_parts(a) for bh in homogeneous_parts(b); init = zero(a)) :
   isnothing(a.homog_elem)                            ? sum(mul_homog(ah, b)  for ah in homogeneous_parts(a); init = zero(a)) :
   isnothing(b.homog_elem)                            ? sum(mul_homog(a, bh)  for bh in homogeneous_parts(b); init = zero(a)) :
                                                        mul_homog(a, b)
)

# Multiplication on homogeneous parts
function mul_homog(a, b)
  @req parent(a) === parent(b) "parent mismatch"
  p = a.homog_deg
  q = b.homog_deg
  A = parent(a)
  C = simplicial_co_complex(A)
  H = homology(C, p+q)[1]
  K = simplicial_complex(C)
  cochain = zero(ambient_free_module(H))
  f = Dict(s => i for (i,s) in enumerate(faces(K, p+q)))
  for (ga, ca) in coordinates(repres(a.homog_elem)), (gb, cb) in coordinates(repres(b.homog_elem))
    sa = faces(K, p)[ga]
    sb = faces(K, q)[gb]
    s = union(sa, sb)
    (maximum(sa) == minimum(sb) && s in keys(f)) || continue
    cochain += ca * cb * gens(C[p+q])[f[s]]
  end
  return SimplicialCohomologyRingElem(A, p+q, H(cochain))
end

# Equality for homogeneous elements is straight foward; for inhomogeneous, do it by sets of homogeneous parts
function ==(a::SimplicialCohomologyRingElem, b::SimplicialCohomologyRingElem)
  if is_homogeneous_normalized(a) && is_homogeneous_normalized(b)
    return a.homog_elem == b.homog_elem && a.homog_deg == b.homog_deg
  else
    return homogeneous_parts(a) == homogeneous_parts(b)
  end
end

import Base.hash
# Either hash the homogeneous information, OR the inhomogeneous information (for homogeneous elements, `coeff` is #undef)
hash(a::SimplicialCohomologyRingElem, h::UInt) = is_homogeneous_normalized(a) ? hash(a.homog_elem, hash(a.homog_deg, h)) : hash(a.coeff, h)
