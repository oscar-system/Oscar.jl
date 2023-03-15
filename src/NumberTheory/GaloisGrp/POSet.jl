#TODO: do this properly, with indexed types and such
# Partially Ordered Set ...
mutable struct POSet{T}
  "list of all elements in the POset"
  elem::Vector{T}
  "for two elements that are comparable, return -1, 0, 1"
  cmp::Function #(::T, ::T) -> -1, 0, 1
  "take two elems and test if they are comparable, returns Bool"
  can_cmp::Function #(::T, ::T) -> Bool

  function POSet{S}(elem::Vector{S}, can_cmp::Function, cmp::Function) where {S}
    r = new{S}()
    r.elem = elem
    r.cmp = cmp
    r.can_cmp = can_cmp
    return r
  end
end

function Base.show(io::IO, P::POSet)
  print(io, "POset on $(P.elem)")
end

function Base.push!(P::POSet{T}, a::T; check::Bool = true) where {T}
  if check
    @assert !any(x->P.can_cmp(x, a) && P.cmp(x, a) == 0, P.elem)
  end
  push!(P.elem, a)
end

function Base.in(a::T, P::POSet{T}) where T
  return any(x->P.can_cmp(x, a) && P.cmp(x, a) == 0, P.elem)
end

function Base.findall(a::T, P::POSet{T}) where T
  return [i for i=1:length(P.elem) if P.can_cmp(P.elem[i], a) && P.cmp(P.elem[i], a) == 0]
end

function Base.length(P::POSet)
  return length(P.elem)
end

POSet(elem::Vector{S}, can_cmp::Function, cmp::Function) where {S} = POSet{S}(elem, can_cmp, cmp)

struct POSetElem{T}
  e::Int
  p::POSet{T}
  function POSetElem(L::POSet{S}, i::Int) where {S}
    return new{S}(i, L)
  end
  function POSetElem(L::POSet{S}, e::Any) where {S}
    return new{S}(findfirst(i->L.can_cmp(e, L.elem[i]) && L.cmp(e, L.elem[i]) == 0, 1:length(L.elem)), L)
  end
end

function (P::POSet{T})(a::T) where {T}
  return POSetElem(P, a)
end
Base.getindex(P::POSet, i::Int) = POSetElem(P, i)

Oscar.parent(E::POSetElem) = E.p
Oscar.data(E::POSetElem) = parent(E).elem[E.e]

function Base.show(io::IO, E::POSetElem)
  print(io, "$(E.e)-th element")
end

function Base.iterate(P::POSet)
  length(P) == 0 && return nothing
  return POSetElem(P, 1), 1
end

function Base.iterate(P::POSet, i::Int)
  i >= length(P) && return nothing
  return POSetElem(P, i+1), i+1
end

Base.eltype(P::POSet{T}) where T = POSetElem{T}

function Base.isless(e::POSetElem{T}, f::POSetElem{T}) where T
  P = parent(e)
  @assert P == parent(f)
  @assert P.can_cmp(P.elem[e.e], P.elem[f.e])
  return P.cmp(P.elem[e,e], P.elem[f,e]) < 0
end

function maximal_elements(L::POSet)
  d = Dict{typeof(L.elem[1]), Vector{Int}}()
  for i = 1:length(L.elem)
    e = L.elem[i]
    d[e] = Int[]
    for j=1:length(L.elem)
      if i == j
        continue
      end
      if L.can_cmp(e, L.elem[j]) && L.cmp(e, L.elem[j]) > 0
        push!(d[e], j)
      end
    end
  end
  return [L(x) for (x,v) = d if length(v) == 0]
end

function minimal_elements(L::POSet)
  d = Dict{typeof(L.elem[1]), Vector{Int}}()
  for i = 1:length(L.elem)
    e = L.elem[i]
    d[e] = Int[]
    for j=1:length(L.elem)
      if i == j
        continue
      end
      if L.can_cmp(e, L.elem[j]) && (L.cmp(e, L.elem[j]) < 0)
        push!(d[e], j)
      end
    end
  end
  return [L(x) for (x,v) = d if length(v) == 0]
end


