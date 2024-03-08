# T=type of the group, S=type of the element
@doc raw"""
    GroupCoset{T<: Group, S <: GAPGroupElem}

Type of group cosets.
Two cosets are equal if, and only if, they are both left (resp. right)
and they contain the same elements.
"""
struct GroupCoset{T<: GAPGroup, S <: GAPGroupElem} 
   G::T                    # big group containing the subgroup and the element
   H::T                    # subgroup
   repr::S                 # element
   side::Symbol            # says if the coset is left or right
   X::GapObj               # GapObj(H*repr)
end

Base.hash(x::GroupCoset, h::UInt) = h # FIXME
Base.eltype(::Type{GroupCoset{T,S}}) where {T,S} = S

function _group_coset(G::GAPGroup, H::GAPGroup, repr::GAPGroupElem, side::Symbol, X::GapObj)
  return GroupCoset{typeof(G), typeof(repr)}(G, H, repr, side, X)
end

function ==(x::GroupCoset, y::GroupCoset)
   return x.X == y.X && x.side == y.side
end

function Base.show(io::IO, ::MIME"text/plain", x::GroupCoset)
  side = x.side === :left ? "Left" : "Right"
  io = pretty(io)
  println(io, "$side coset of ", Lowercase(), x.H)
  print(io, Indent())
  println(io, "with representative ", x.repr)
  print(io, "in ", Lowercase(), x.G)
  print(io, Dedent())
end

function Base.show(io::IO, x::GroupCoset)
  side = x.side === :left ? "Left" : "Right"
  if get(io, :supercompact, false)
    print(io, "$side coset of a group")
  else
    print(io, "$side coset of ")
    io = pretty(io)
    print(IOContext(io, :supercompact => true), Lowercase(), x.H, " with representative ", x.repr)
  end
end


"""
    right_coset(H::Group, g::GAPGroupElem)
    *(H::Group, g::GAPGroupElem)

Return the coset `Hg`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Sym(5)

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Sym(3)

julia> right_coset(H, g)
Right coset of Sym(3)
  with representative (1,3)(2,4,5)
  in Sym(5)
```
"""
function right_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   @req GAPWrap.IsSubset(parent(g).X, H.X) "H is not a subgroup of parent(g)"
   return _group_coset(parent(g), H, g, :right, GAP.Globals.RightCoset(H.X,g.X))
end

"""
    left_coset(H::Group, g::GAPGroupElem)
    *(g::GAPGroupElem, H::Group)

Return the coset `gH`.
!!! note
    Since GAP supports right cosets only, the underlying GAP object of
    `left_coset(H,g)` is the right coset `H^(g^-1) * g`.

# Examples
```jldoctest
julia> g = perm([3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Sym(3)

julia> gH = left_coset(H, g)
Left coset of Sym(3)
  with representative (1,3)(2,4,5)
  in Sym(5)
```
"""
function left_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   @req GAPWrap.IsSubset(parent(g).X, H.X) "H is not a subgroup of parent(g)"
   return _group_coset(parent(g), H, g, :left, GAP.Globals.RightCoset(GAP.Globals.ConjugateSubgroup(H.X,GAP.Globals.Inverse(g.X)),g.X))
end


"""
    is_left(c::GroupCoset)

Return whether the coset `c` is a left coset of its acting domain.
"""
is_left(c::GroupCoset) = c.side == :left

"""
    is_right(c::GroupCoset)

Return whether the coset `c` is a right coset of its acting domain.
"""
is_right(c::GroupCoset) = c.side == :right

Base.:*(H::GAPGroup, g::GAPGroupElem) = right_coset(H,g)
Base.:*(g::GAPGroupElem, H::GAPGroup) = left_coset(H,g)

function Base.:*(c::GroupCoset, y::GAPGroupElem)
   @assert y in c.G "element not in the group"
   if c.side == :right
      return right_coset(c.H, c.repr*y)
   else
      return left_coset(c.H^y, c.repr*y)
   end
end

function Base.:*(y::GAPGroupElem, c::GroupCoset)
   @assert y in c.G "element not in the group"
   if c.side == :left
      return left_coset(c.H, y*c.repr)
   else
      return right_coset(c.H^(y^-1), y*c.repr)
   end
end

function Base.:*(c::GroupCoset, d::GroupCoset)
   @req (c.side == :right && d.side == :left) "Wrong input"
   return double_coset(c.H, c.repr*d.repr, d.H)
end

"""
    acting_domain(C::GroupCoset)

If `C` = `Hx` or `xH`, return `H`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Sym(5)

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Sym(3)

julia> gH = left_coset(H,g)
Left coset of Sym(3)
  with representative (1,3)(2,4,5)
  in Sym(5)

julia> acting_domain(gH)
Sym(3)
```
"""
acting_domain(C::GroupCoset) = C.H

"""
    representative(C::GroupCoset)

If `C` = `Hx` or `xH`, return `x`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Sym(5)

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Sym(3)

julia> gH = left_coset(H, g)
Left coset of Sym(3)
  with representative (1,3)(2,4,5)
  in Sym(5)

julia> representative(gH)
(1,3)(2,4,5)
```
"""
representative(C::GroupCoset) = C.repr


"""
    is_bicoset(C::GroupCoset)

Return whether `C` is simultaneously a right coset and a left coset for the same subgroup `H`.  This 
is the case if and only if the coset representative normalizes the acting domain subgroup.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Sym(5)

julia> H = symmetric_group(4)
Sym(4)

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> gH = left_coset(H, g)
Left coset of Sym(4)
  with representative (1,3)(2,4,5)
  in Sym(5)

julia> is_bicoset(gH)
false

julia> f = perm(G,[2,1,4,3,5])
(1,2)(3,4)

julia> fH = left_coset(H, f)
Left coset of Sym(4)
  with representative (1,2)(3,4)
  in Sym(5)

julia> is_bicoset(fH)
true
```
"""
is_bicoset(C::GroupCoset) = GAPWrap.IsBiCoset(C.X)

"""
    right_cosets(G::T, H::T; check::Bool=true) where T<: GAPGroup

Return the G-set that describes the right cosets of `H` in `G`.

If `check == false`, do not check whether `H` is a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Sym(4)

julia> H = symmetric_group(3)
Sym(3)

julia> rc = right_cosets(G, H)
Right cosets of
  Sym(3) in
  Sym(4)

julia> collect(rc)
4-element Vector{GroupCoset{PermGroup, PermGroupElem}}:
 Right coset of H with representative ()
 Right coset of H with representative (1,4)
 Right coset of H with representative (1,4,2)
 Right coset of H with representative (1,4,3)
```
"""
function right_cosets(G::T, H::T; check::Bool=true) where T<: GAPGroup
  return GSetBySubgroupTransversal(G, H, :right, check = check)
end

"""
    left_cosets(G::T, H::T; check::Bool=true) where T<: GAPGroup

Return the G-set that describes the left cosets of `H` in `G`.

If `check == false`, do not check whether `H` is a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Sym(4)

julia> H = symmetric_group(3)
Sym(3)

julia> left_cosets(G, H)
Left cosets of
  Sym(3) in
  Sym(4)
```
"""
function left_cosets(G::T, H::T; check::Bool=true) where T<: GAPGroup
  return GSetBySubgroupTransversal(G, H, :left, check = check)
end


@doc raw"""
    SubgroupTransversal{T<: GAPGroup, S<: GAPGroup, E<: GAPGroupElem}

Type of left/right transversals of subgroups in groups.
The elements are encoded via a right transversal object in GAP.
(Note that GAP does not support left transversals.)

Objects of this type are created by [`right_transversal`](@ref) and
[`left_transversal`](@ref).
"""
struct SubgroupTransversal{T<: GAPGroup, S<: GAPGroup, E<: GAPGroupElem} <: AbstractVector{E}
   G::T                    # big group containing the subgroup
   H::S                    # subgroup
   side::Symbol            # says if the transversal is left or right
   X::GapObj               # underlying *right* transversal in GAP
end

function Base.show(io::IO, ::MIME"text/plain", x::SubgroupTransversal)
  side = x.side === :left ? "Left" : "Right"
  println(io, "$side transversal of length $(length(x)) of")
  io = pretty(io)
  print(io, Indent())
  println(io, Lowercase(), x.H, " in")
  print(io, Lowercase(), x.G)
  print(io, Dedent())
end

function Base.show(io::IO, x::SubgroupTransversal)
  side = x.side === :left ? "Left" : "Right"
  if get(io, :supercompact, false)
    print(io, "$side transversal of groups")
  else
    print(io, "$side transversal of ")
    io = pretty(io)
    print(IOContext(io, :supercompact => true), Lowercase(), x.H, " in ", Lowercase(), x.G)
  end
end

Base.hash(x::SubgroupTransversal, h::UInt) = h # FIXME

Base.length(T::SubgroupTransversal) = index(Int, T.G, T.H)

function Base.getindex(T::SubgroupTransversal, i::Int)
  res = group_element(T.G, T.X[i])
  if T.side === :left
    res = inv(res)
  end
  return res
end

# in order to make `T[end]` work
Base.size(T::SubgroupTransversal) = (index(Int, T.G, T.H),)
Base.lastindex(T::SubgroupTransversal) = length(T)

# in order to make `findfirst` and `findall` work
function Base.keys(T::SubgroupTransversal)
    return keys(1:length(T))
end


"""
    right_transversal(G::T, H::T; check::Bool=true) where T<: GAPGroup

Return a vector containing a complete set of representatives for
the right cosets of `H` in `G`.
This vector is not mutable, and it does not store its entries explicitly,
they are created anew with each access to the transversal.

If `check == false`, do not check whether `H` is a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Sym(4)

julia> H = symmetric_group(3)
Sym(3)

julia> T = right_transversal(G, H)
Right transversal of length 4 of
  Sym(3) in
  Sym(4)

julia> collect(T)
4-element Vector{PermGroupElem}:
 ()
 (1,4)
 (1,4,2)
 (1,4,3)
```
"""
function right_transversal(G::T, H::T; check::Bool=true) where T<: GAPGroup
   @req (!check || GAPWrap.IsSubset(G.X, H.X)) "H is not a subgroup of G"
   return SubgroupTransversal{T, T, eltype(T)}(G, H, :right,
              GAP.Globals.RightTransversal(G.X, H.X))
end

"""
    left_transversal(G::T, H::T; check::Bool=true) where T<: Group

Return a vector containing a complete set of representatives for
the left cosets for `H` in `G`.
This vector is not mutable, and it does not store its entries explicitly,
they are created anew with each access to the transversal.

If `check == false`, do not check whether `H` is a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Sym(4)

julia> H = symmetric_group(3)
Sym(3)

julia> T = left_transversal(G, H)
Left transversal of length 4 of
  Sym(3) in
  Sym(4)

julia> collect(T)
4-element Vector{PermGroupElem}:
 ()
 (1,4)
 (1,2,4)
 (1,3,4)
```
"""
function left_transversal(G::T, H::T; check::Bool=true) where T<: GAPGroup
   @req (!check || GAPWrap.IsSubset(G.X, H.X)) "H is not a subgroup of G"
   return SubgroupTransversal{T, T, eltype(T)}(G, H, :left,
              GAP.Globals.RightTransversal(G.X, H.X))
end

Base.IteratorSize(::Type{<:GroupCoset}) = Base.SizeUnknown()
Base.iterate(G::GroupCoset) = iterate(G, GAPWrap.Iterator(G.X))

function Base.iterate(G::GroupCoset, state)
  GAPWrap.IsDoneIterator(state) && return nothing
  i = GAPWrap.NextIterator(state)::GapObj
  return group_element(G.G, i), state
end



@doc raw"""
    GroupDoubleCoset{T<: Group, S <: GAPGroupElem}

Group double coset.
Two double cosets are equal if, and only if, they contain the same elements.
"""
struct GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
# T=type of the group, S=type of the element
   G::T
   H::T
   K::T
   repr::S
   X::GapObj
end

Base.hash(x::GroupDoubleCoset, h::UInt) = h # FIXME
Base.eltype(::Type{GroupDoubleCoset{T,S}}) where {T,S} = S

function ==(x::GroupDoubleCoset, y::GroupDoubleCoset)
   return x.X == y.X
end

function Base.show(io::IO, ::MIME"text/plain", x::GroupDoubleCoset)
  io = pretty(io)
  println(io, "Double coset of ", Lowercase(), x.H)
  print(io, Indent())
  println(io, "and ", Lowercase(), x.K)
  println(io, "with representative ", x.repr)
  print(io, "in ", Lowercase(), x.G)
  print(io, Dedent())
end

function Base.show(io::IO, x::GroupDoubleCoset)
  if get(io, :supercompact, false)
    print(io, "Double coset of a group")
  else
    print(io, "Double coset of ")
    io = pretty(io)
    print(IOContext(io, :supercompact => true), Lowercase(), x.H,
      " and ", Lowercase(), x.K, " with representative ", x.repr)
  end
end


"""
    double_coset(H::Group, x::GAPGroupElem, K::Group)
    *(H::Group, x::GAPGroupElem, K::Group)

Return the double coset `HxK`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Sym(5)

julia> g = perm(G,[3,4,5,1,2])
(1,3,5,2,4)

julia> H = symmetric_group(3)
Sym(3)

julia> K = symmetric_group(2)
Sym(2)

julia> double_coset(H,g,K)
Double coset of Sym(3)
  and Sym(2)
  with representative (1,3,5,2,4)
  in Sym(5)
```
"""
function double_coset(G::T, g::GAPGroupElem{T}, H::T) where T<: GAPGroup
   @req GAPWrap.IsSubset(parent(g).X,G.X) "G is not a subgroup of parent(g)"
   @req GAPWrap.IsSubset(parent(g).X,H.X) "H is not a subgroup of parent(g)"
   return GroupDoubleCoset(parent(g),G,H,g,GAP.Globals.DoubleCoset(G.X,g.X,H.X))
end

Base.:*(H::GAPGroup, g::GAPGroupElem, K::GAPGroup) = double_coset(H,g,K)

"""
    double_cosets(G::T, H::T, K::T; check::Bool=true) where T<: GAPGroup

Return a vector of all the double cosets `HxK` for `x` in `G`.
If `check == false`, do not check whether `H` and `K` are subgroups of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Sym(4)

julia> H = symmetric_group(3)
Sym(3)

julia> K = symmetric_group(2)
Sym(2)

julia> double_cosets(G,H,K)
3-element Vector{GroupDoubleCoset{PermGroup, PermGroupElem}}:
 Double coset of H and K with representative ()
 Double coset of H and K with representative (1,4)
 Double coset of H and K with representative (1,4,3)
```
"""
function double_cosets(G::T, H::T, K::T; check::Bool=true) where T<: GAPGroup
   if !check
      dcs = GAP.Globals.DoubleCosetsNC(G.X,H.X,K.X)
   else
      @assert is_subset(H, G) "H is not a subgroup of G"
      @assert is_subset(K, G) "K is not a subgroup of G"
      dcs = GAP.Globals.DoubleCosets(G.X,H.X,K.X)
   end
   res = Vector{GroupDoubleCoset{T,elem_type(T)}}(undef, length(dcs))
   for i = 1:length(res)
     dc = dcs[i]
     g = group_element(G, GAPWrap.Representative(dc))
     res[i] = GroupDoubleCoset(G,H,K,g,dc)
   end
   return res
   #return [GroupDoubleCoset(G,H,K,group_element(G.X,GAPWrap.Representative(dc)),dc) for dc in dcs]
end


"""
    order(C::Union{GroupCoset,GroupDoubleCoset})

Return the cardinality of the (double) coset `C`.
"""
order(C::Union{GroupCoset,GroupDoubleCoset}) = GAPWrap.Size(C.X)
Base.length(C::Union{GroupCoset,GroupDoubleCoset}) = GAPWrap.Size(C.X)

"""
    rand(rng::Random.AbstractRNG = Random.GLOBAL_RNG, C::Union{GroupCoset,GroupDoubleCoset})

Return a random element of the (double) coset `C`,
using the random number generator `rng`.
"""
Base.rand(C::Union{GroupCoset,GroupDoubleCoset}) = Base.rand(Random.GLOBAL_RNG, C)

function Base.rand(rng::Random.AbstractRNG, C::Union{GroupCoset,GroupDoubleCoset})
  s = GAP.Globals.Random(GAP.wrap_rng(rng), C.X)
  return group_element(C.G, s)
end

"""
    representative(C::GroupDoubleCoset)

Return a representative `x` of the double coset `C` = `HxK`.
"""
representative(C::GroupDoubleCoset) = C.repr

"""
    left_acting_group(C::GroupDoubleCoset)

Given a double coset `C` = `HxK`, return `H`.
"""
left_acting_group(C::GroupDoubleCoset) = C.H

"""
    right_acting_group(C::GroupDoubleCoset)

Given a double coset `C` = `HxK`, return `K`.
"""
right_acting_group(C::GroupDoubleCoset) = C.K

Base.IteratorSize(::Type{<:GroupDoubleCoset}) = Base.SizeUnknown()

Base.iterate(G::GroupDoubleCoset) = iterate(G, GAPWrap.Iterator(G.X))

function Base.iterate(G::GroupDoubleCoset, state)
  GAPWrap.IsDoneIterator(state) && return nothing
  i = GAPWrap.NextIterator(state)::GapObj
  return group_element(G.G, i), state
end

"""
    intersect(V::AbstractVector{Union{<: GAPGroup, GroupCoset, GroupDoubleCoset}})

Return a vector containing all elements belonging to all groups and cosets
in `V`.
"""
function intersect(V::AbstractVector{Union{<: GAPGroup, GroupCoset, GroupDoubleCoset}})
   if V[1] isa GAPGroup
      G = V[1]
   else
      G = V[1].G
   end
   l = GAP.Obj([v.X for v in V])
   ints = GAP.Globals.Intersection(l)
   L = Vector{typeof(G)}(undef, length(ints))
   for i in 1:length(ints)
      L[i] = group_element(G,ints[i])
   end

   return L
end
