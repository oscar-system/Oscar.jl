"""
    IndexedSet{T}

A hybrid data structure combining properties of both a vector and a set.

By implementation you can push things to its end and iterate it sequentially;
but duplicates are simply dropped (like in a set). Therefore, given any element you can
test in O(1) if it is contained in the set, and find out at which position.

This structure is suited well for applications such as orbit algorithms, 
where efficient membership checks and indexed retrieval are beneficial.

```jldoctest
julia> s = IndexedSet(["ball", "flower", "house"])
Indexed set ["ball", "flower", "house"]

julia> "ball" in s
true

julia> s("flower")   # get the position
2

julia> s[2]
"flower"

julia> "stone" in s
false

julia> s("stone")    # the "position" of objects now in there is 0
0

julia> push!(s, "ball", "stone")    # only "new" things are actually added
Indexed set ["ball", "flower", "house", "stone"]
```
"""
struct IndexedSet{T}
  data::Vector{T}
  pos::Dict{T,Int}
  
  IndexedSet{T}() where T = new{T}(Vector{T}(),Dict{T,Int}())
end

function IndexedSet(v::Vector{T}) where T
  s = IndexedSet{T}()
  for x in v
    push!(s, x)
  end
  return s
end

function Base.show(io::IO, s::IndexedSet)
  print(io, "Indexed set ", s.data)
end

function Base.push!(s::IndexedSet, x)
  n = length(s.data)+1
  i = get!(s.pos, x, n)
  n == i && push!(s.data, x)
  return s
end

Base.append!(s::IndexedSet{T}, x::Vector{T}) where {T} = foreach(e -> push!(s, e), x)
Base.union!(s::IndexedSet{T}, iters...) where {T} = foreach(e -> append!(s, e), iters)

function Base.in(x, s::IndexedSet)
  return haskey(s.pos, x)
end

Base.hash(s::IndexedSet, h::UInt) = hash(s.data, h) # ^ 0x0179cc9dd07222ac
Base.isequal(a::IndexedSet, b::IndexedSet) = a.data == b.data

# lookup object via index
Base.getindex(s::IndexedSet, i::Int) = s.data[i]

# lookup index via object
(s::IndexedSet{T})(x::T) where T = get(s, x, 0)
Base.get(s::IndexedSet, x::T, def::Int) where T = get(s.pos, x, def)

Base.iterate(s::IndexedSet) = iterate(s.data)
Base.iterate(s::IndexedSet, state) = iterate(s.data, state)

Base.isempty(s::IndexedSet) = isempty(s.data)
Base.length(s::IndexedSet) = length(s.data)

Base.IteratorSize(::Type{<:IndexedSet}) = Base.HasShape{1}()
Base.size(s::IndexedSet) = size(s.data)
Base.eltype(s::IndexedSet{T}) where T = T

#=
module Orbit
import ..Foo: IndexedSet
import ..Nemo: Group, gens

function orbit_1(generators, pt::T, action = ^; maxlen = nothing) where T
    orb = IndexedSet([pt])
    for b in orb
        for g in generators
            c = action(b, g)::T
            pos = orb(c)
            if pos === 0
                push!(orb, c)
                if maxlen !== nothing
                  if length(orb) > maxlen
                    println("aborting after reaching orbit of length $maxlen")
                    return orb
                  end
                end
            end
        end
    end
    return orb
end

function orbit_2(G::Group, pt::T, action = ^) where T
    gg = gens(G)
    orb = IndexedSet([pt])
    for b in orb
        for g in gg
            c = action(b, g)::T
            pos = orb(c)
            if pos === 0
                push!(orb, c)
            end
        end
    end
    return orb
end

function orbit_3(pt::T, generators, action) where T
    orb = IndexedSet([pt])
    for b in orb
        for g in generators
            c = action(b, g)::T
            if !(c in orb)
                push!(orb, c)
            end
        end
    end
    return orb
end

function orbit_4(pt::T, generators) where T
    orb = IndexedSet([pt])
    schreier = Int[0]
    for b in orb
        for (i,g) in enumerate(generators)
            c = (b^g)::T
            if !(c in orb)
                push!(orb, c)
                push!(schreier, i)
            end
        end
    end
    return orb, schreier
end

end # Orbit
=#
