mutable struct SimplicialCohomologyRing{T} <: NCRing
  C::SimplicialCoComplex
  graded_parts::Vector{SubquoModule{T}}

  function SimplicialCohomologyRing(K::SimplicialComplex)
    return SimplicialCohomologyRing(SimplicialCoComplex(K))
  end
  function SimplicialCohomologyRing(C::SimplicialCoComplex)
    T = elem_type(base_ring(C))
    return new{T}(C)
  end
end

simplicial_co_complex(A::SimplicialCohomologyRing) = A.C

function graded_parts(A::SimplicialCohomologyRing)
  if !isdefined(A, :graded_parts)
    K = simplicial_complex(A)
    A.graded_parts = [homology(A.C, i)[1] for i in 0:dim(K)]
  end
  return A.graded_parts
end

function graded_part(A::SimplicialCohomologyRing, i::Int)
  return graded_parts(A)[i+1]
end

mutable struct SimplicialCohomologyRingElem{T} <: NCRingElem
  parent::SimplicialCohomologyRing{T}
  # Elements can be represented in three ways:
  #  1) no fields assigned (the zero element)
  #  2) `homog_elem` and `homog_deg` assigned (homogeneously stored element)
  #  3) a dictionary `p => v` where the degree points to the cohomology 
  #     class for the degree `p` part of this element.
  # The values of `coeff` must not be zero. The zero element in 
  # the ring must be stored as as the empty element. 
  #
  # The distinction whether or not an element is stored as a homogeneous 
  # element can be made by whether or not the field `.homog_elem` is set 
  # to `nothing`. The field `.coeff` is not assigned if it's not needed. 
  homog_elem::Union{SubquoModuleElem{T}, Nothing}
  homog_deg::Union{Int, Nothing}
  coeff::Dict{Int, SubquoModuleElem{T}}

  # Constructor for homogeneous elements
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRing{T}, 
      p::Int, v::SubquoModuleElem{T}
    ) where {T}
    @assert parent(v) === graded_part(A, p)
    return new{T}(A, v, p)
  end

  # Constructor for mixed elements
  function SimplicialCohomologyRingElem(
      A::SimplicialCohomologyRing{T}, 
      coeff::Dict{Int, SubquoModuleElem{T}};
      check::Bool=true
    ) where {T}
    @check all(parent(v) === graded_part(A, p) for (p, v) in coeff) "degree/element incompatibility"
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
      isempty(result.coeff) && return zero(A)
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
      isempty(result.coeff) && return zero(A)
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
      isempty(result.coeff) && return zero(A)
      return result
    else # b is homogeneous
      p = b.homog_deg
      if p == q
        result.homog_deg = p
        result.homog_elem = a.homog_elem + b.homog_elem
        is_zero(result.homog_elem) && return zero(A)
        return result
      else
        result.coeff = Dict{Int, SubquoModuleElem{T}}([q => a.homog_elem, p => b.homog_elem])
        isempty(result.coeff) && return zero(A)
        return result
      end
    end
  end
end

is_homogeneous_normalized(a::SimplicialCohomologyRingElem) = !isnothing(a.homog_elem)

is_homogeneous_denormalized(a::SimplicialCohomologyRingElem) = isnothing(a.homog_elem) && isone(length(a.coeff))

is_homogeneous(a::SimplicialCohomologyRingElem) = iszero(a) || is_homogeneous_normalized(a) || is_homogeneous_denormalized(a)

function degree(a::SimplicialCohomologyRingElem)
  if is_homogeneous_normalized(a)
    return a.homog_deg
  elseif is_homogeneous_denormalized(a)
    return only(keys(a.coeff))
  else
    throw(ArgumentError("not homogeneous"))
  end
end

-(a::SimplicialCohomologyRingElem) = is_homogeneous_normalized(a) ? SimplicialCohomologyRingElem(parent(a), a.homog_deg, -a.homog_elem) : sum(-b for b in homogeneous_parts(a))

-(a::SimplicialCohomologyRingElem, b::SimplicialCohomologyRingElem) = a + (-b)

# extract homogeneous parts. We do a set, because we need it for ==
homogeneous_parts(a::SimplicialCohomologyRingElem) = Set(isnothing(a.homog_elem) ? [SimplicialCohomologyRingElem(parent(a), i, m) for (i,m) in pairs(a.coeff)] : [a])

# distribute over homogeneous parts
*(a::SimplicialCohomologyRingElem, b::SimplicialCohomologyRingElem) = (
   is_homogeneous_normalized(a) && is_homogeneous_normalized(b) ? sum(mul_homog(ah, bh) for ah in homogeneous_parts(a) for bh in homogeneous_parts(b); init = zero(a)) :
   is_homogeneous_normalized(a)                                 ? sum(mul_homog(ah, b)  for ah in homogeneous_parts(a); init = zero(a)) :
   is_homogeneous_normalized(b)                                 ? sum(mul_homog(a, bh)  for bh in homogeneous_parts(b); init = zero(a)) :
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

# Parent and element types
elem_type(::Type{SimplicialCohomologyRing{T}}) where {T} = SimplicialCohomologyRingElem{T}
parent_type(::Type{SimplicialCohomologyRingElem{T}}) where {T} = SimplicialCohomologyRing{T}


# Base ring and base ring type
base_ring(C::SimplicialCohomologyRing) = base_ring(simplicial_co_complex(C))
base_ring_type(::Type{SimplicialCohomologyRingElem{T}}) where {T} = parent_type(T)

# Equality for homogeneous elements is straight foward; for inhomogeneous, do it by sets of homogeneous parts
function ==(a::SimplicialCohomologyRingElem, b::SimplicialCohomologyRingElem)
  if is_homogeneous_normalized(a) && is_homogeneous_normalized(b)
    return a.homog_elem == b.homog_elem && a.homog_deg == b.homog_deg
  else
    return homogeneous_parts(a) == homogeneous_parts(b)
  end
end

function Base.hash(a::SimplicialCohomologyRingElem, u::UInt)
  A = parent(a)
  r = 0x60125ab8e8cd44ca
  if isdefined(a, :homog_elem)
    r = xor(r, hash(a.homog_deg, u))
    r = xor(r, hash(a.homog_elem, u))
    return r
  end

  !isdefined(a, :coeff) && return xor(r, u)
  
  # we need to catch the case where homogeneous elements use the dictionary
  if isone(length(a.coeff))
    p, c = only(a.coeff)
    r = xor(r, hash(p, u))
    r = xor(r, hash(c))
    return r
  end

  return hash(a.coeff, u)
end

