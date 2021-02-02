#TODO: do this properly, with indexed types and such
# Partially Ordered Set ...
mutable struct POSet{T}
  elem::Array{T, 1}
  "for two elements that are comparable, return -1, 0, 1"
  cmp::Function #(::T, ::T) -> -1, 0, 1
  "take two elems and test if they are comparable, returns Bool"
  can_cmp::Function #(::T, ::T) -> Bool

  function POSet{S}(elem::Array{S, 1}, can_cmp::Function, cmp::Function) where {S}
    r = new{S}()
    r.elem = elem
    r.cmp = cmp
    r.can_cmp = can_cmp
    return r
  end
end

POSet(elem::Array{S, 1}, can_cmp::Function, cmp::Function) where {S} = POSet{S}(elem, can_cmp, cmp)

struct POSetElem{T}
  e::Int
  p::POSet{T}
  function POSetElem(L::POSet{S}, i::Int) where {S}
    return new{S}(i, L)
  end
  function POSetElem(L::POSet{S}, e::Any) where {S}
    return new{S}(findfirst(x->L>can_cmp(e, L.elem[i]) && L.cmp(e, L.elem[i]) == 0, 1:length(L.elem)), L)
  end
end

function maximal_elements(L::POSet)
  d = Dict{typeof(L.elem[1]), Array{Int, 1}}()
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
  return [x for (x,v) = d if length(v) == 0]
end

function minimal_elements(L::POSet)
  d = Dict{typeof(L.elem[1]), Array{Int, 1}}()
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
  return [x for (x,v) = d if length(v) == 0]
end


