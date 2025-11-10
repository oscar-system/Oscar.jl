# T=type of the group, S=type of the element
@doc raw"""
    GroupCoset{TG <: GAPGroup, TH <: GAPGroup, S <: GAPGroupElem}

Type of right and left cosets of subgroups in groups.

For an element $g$ in a group $G$, and a subgroup $H$ of $G$,
the set $Hg = \{ hg; h \in H \}$ is a right coset of $H$ in $G$,
and the set $gH = \{ gh; h \in H \}$ is a left coset of $H$ in $G$.

- [`group(C::GroupCoset)`](@ref) returns $G$.

- [`acting_group(C::GroupCoset)`](@ref) returns $H$.

- [`representative(C::GroupCoset)`](@ref) returns an element
  (the same element for each call) of `C`.

- [`is_right(C::GroupCoset)`](@ref) and [`is_left(C::GroupCoset)`](@ref)
  return whether `C` is a right or left coset, respectively.

Two cosets are equal if and only if they are both left or right, respectively,
and they contain the same elements.
"""
struct GroupCoset{TG <: GAPGroup, TH <: GAPGroup, S <: GAPGroupElem}
   G::TG                   # big group containing the subgroup and the element
   H::TH                   # subgroup (may have a different type)
   repr::S                 # element
   side::Symbol            # says if the coset is left or right
   X::Ref{GapObj}          # GapObj(H*repr)

   function GroupCoset(G::TG, H::TH, representative::S, side::Symbol) where {TG <: GAPGroup, TH <: GAPGroup, S <:GAPGroupElem}
     return new{TG, TH, S}(G, H, representative, side, Ref{GapObj}())
   end
end

GAP.@install function GapObj(obj::GroupCoset)
  if !isassigned(obj.X)
    g = GapObj(representative(obj))
    if is_right(obj)
      obj.X[] = GAPWrap.RightCoset(GapObj(obj.H), g)
    else
      obj.X[] = GAPWrap.RightCoset(GAPWrap.ConjugateSubgroup(GapObj(obj.H), GAPWrap.Inverse(g)), g)
    end
  end
  return obj.X[]::GapObj
end

Base.hash(x::GroupCoset, h::UInt) = h # FIXME
Base.eltype(::Type{GroupCoset{TG, TH, S}}) where {TG, TH, S} = S

function ==(C1::GroupCoset, C2::GroupCoset)
  C1 === C2 && return true
  H = C1.H
  right = is_right(C1)
  (right == is_right(C2) && C1.G == C2.G && H == C2.H ) || return false
  if right
    # Hx == Hy if x/y in H
    return representative(C1) / representative(C2) in H
  else
    # xH == yH if x\y in H
    return representative(C1) \ representative(C2) in H
  end
end

function Base.show(io::IO, ::MIME"text/plain", x::GroupCoset)
  side = is_left(x) ? "Left" : "Right"
  io = pretty(io)
  println(io, "$side coset of ", Lowercase(), x.H)
  print(io, Indent())
  println(io, "with representative ", representative(x))
  print(io, "in ", Lowercase(), x.G)
  print(io, Dedent())
end

function Base.show(io::IO, x::GroupCoset)
  side = is_left(x) ? "Left" : "Right"
  if is_terse(io)
    print(io, "$side coset of a group")
  else
    print(io, "$side coset of ")
    io = pretty(io)
    print(terse(io), Lowercase(), x.H, " with representative ", representative(x))
  end
end


"""
    right_coset(H::Group, g::GAPGroupElem)
    *(H::Group, g::GAPGroupElem)

Return the coset `Hg`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> right_coset(H, g)
Right coset of symmetric group of degree 3
  with representative (1,3)(2,4,5)
  in symmetric group of degree 5
```
"""
function right_coset(H::GAPGroup, g::GAPGroupElem)
   G = parent(g)
   @req GAPWrap.IsSubset(GapObj(G), GapObj(H)) "H is not a subgroup of parent(g)"
   return GroupCoset(G, H, g, :right)
end

"""
    left_coset(H::Group, g::GAPGroupElem)
    *(g::GAPGroupElem, H::Group)

Return the coset `gH`.
!!! note
    Since GAP supports right cosets only, the underlying GAP object of
    `left_coset(H,g)`, if assigned, is the right coset `H^(g^-1) * g`.

# Examples
```jldoctest
julia> g = perm([3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> gH = left_coset(H, g)
Left coset of symmetric group of degree 3
  with representative (1,3)(2,4,5)
  in symmetric group of degree 5
```
"""
function left_coset(H::GAPGroup, g::GAPGroupElem)
   G = parent(g)
   @req GAPWrap.IsSubset(GapObj(G), GapObj(H)) "H is not a subgroup of parent(g)"
   return GroupCoset(G, H, g, :left)
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
   yy = c.G(y)
   if is_right(c)
      return right_coset(c.H, representative(c)*yy)
   else
      return left_coset(c.H^y, representative(c)*yy)
   end
end

function Base.:*(y::GAPGroupElem, c::GroupCoset)
   yy = c.G(y)
   if is_left(c)
      return left_coset(c.H, yy*representative(c))
   else
      return right_coset(c.H^(y^-1), yy*representative(c))
   end
end

function Base.:*(c::GroupCoset, d::GroupCoset)
   @req (is_right(c) && is_left(d)) "Wrong input"
   return double_coset(c.H, representative(c)*representative(d), d.H)
end


"""
    group(C::GroupCoset)

Return the group `G` that is the parent of all elements in `C`.
That is, `C` is a left or right coset of a subgroup of `G` in `G`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = sylow_subgroup(G, 2)[1]
Permutation group of degree 5 and order 8

julia> C = right_coset(H, gen(G, 1))
Right coset of permutation group of degree 5 and order 8
  with representative (1,2,3,4,5)
  in symmetric group of degree 5

julia> group(C) == G
true
```
"""
group(C::GroupCoset) = C.G


"""
    acting_group(C::GroupCoset)

Return the group `H` such that `C` is `Hx` (if `C` is a right coset)
or `xH` (if `C` is a left coset), for an element `x` in `C`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> C = right_coset(H, gen(G, 1))
Right coset of symmetric group of degree 3
  with representative (1,2,3,4,5)
  in symmetric group of degree 5

julia> acting_group(C) == H
true
```
"""
acting_group(C::GroupCoset) = C.H

"""
    representative(C::GroupCoset)

Return an element `x` in `group(C)` such that
`C` = `Hx` (if `C` is a right coset)
or `xH` (if `C` is a left coset).

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> Hg = right_coset(H, g)
Right coset of symmetric group of degree 3
  with representative (1,3)(2,4,5)
  in symmetric group of degree 5

julia> representative(Hg)
(1,3)(2,4,5)
```
"""
representative(C::GroupCoset) = C.repr


"""
    is_bicoset(C::GroupCoset)

Return whether `C` is simultaneously a right coset and a left coset
for the same subgroup `H`.
This is the case if and only if the coset representative normalizes
`acting_group(C)`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(4)
Symmetric group of degree 4

julia> g = perm(G,[3,4,1,5,2])
(1,3)(2,4,5)

julia> gH = left_coset(H, g)
Left coset of symmetric group of degree 4
  with representative (1,3)(2,4,5)
  in symmetric group of degree 5

julia> is_bicoset(gH)
false

julia> f = perm(G,[2,1,4,3,5])
(1,2)(3,4)

julia> fH = left_coset(H, f)
Left coset of symmetric group of degree 4
  with representative (1,2)(3,4)
  in symmetric group of degree 5

julia> is_bicoset(fH)
true
```
"""
is_bicoset(C::GroupCoset) = GAPWrap.IsBiCoset(GapObj(C))

"""
    right_cosets(G::GAPGroup, H::GAPGroup; check::Bool=true)

Return the G-set that describes the right cosets of `H` in `G`.

If `check == false`, do not check whether `H` is a subgroup of `G`.

Use [`right_transversal`](@ref) to compute the vector of coset representatives.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Symmetric group of degree 4

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> rc = right_cosets(G, H)
Right cosets of
  symmetric group of degree 3 in
  symmetric group of degree 4

julia> collect(rc)
4-element Vector{GroupCoset{PermGroup, PermGroup, PermGroupElem}}:
 Right coset of H with representative ()
 Right coset of H with representative (1,4)
 Right coset of H with representative (1,4,2)
 Right coset of H with representative (1,4,3)
```
"""
function right_cosets(G::GAPGroup, H::GAPGroup; check::Bool=true)
  return GSetBySubgroupTransversal(G, H, :right, check = check)
end

"""
    left_cosets(G::GAPGroup, H::GAPGroup; check::Bool=true)

Return the G-set that describes the left cosets of `H` in `G`.

If `check == false`, do not check whether `H` is a subgroup of `G`.

Use [`left_transversal`](@ref) to compute the vector of coset representatives.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Symmetric group of degree 4

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> left_cosets(G, H)
Left cosets of
  symmetric group of degree 3 in
  symmetric group of degree 4
```
"""
function left_cosets(G::GAPGroup, H::GAPGroup; check::Bool=true)
#T _check_compatible(G, H) ?
  return GSetBySubgroupTransversal(G, H, :left, check = check)
end


@doc raw"""
    SubgroupTransversal{T<: GAPGroup, S<: GAPGroup, E<: GAPGroupElem}

Type of left/right transversals of subgroups in groups.

For a group $G$ and a subgroup $H$ of $G$, $T$ is a right
(resp. left) transversal for $H$ in $G$ if $T$ contains
precisely one element of each right (resp. left) cosets of $H$ in $G$.

Objects of this type are created by [`right_transversal`](@ref) and
[`left_transversal`](@ref).

- [`group(T::SubgroupTransversal)`](@ref) returns $G$.

- [`subgroup(T::SubgroupTransversal)`](@ref) returns $H$.

# Note for developers

The elements are encoded via a right transversal object in GAP.
(Note that GAP does not support left transversals.)
"""
struct SubgroupTransversal{T<: GAPGroup, S<: GAPGroup, E<: GAPGroupElem} <: AbstractVector{E}
   G::T                    # big group containing the subgroup
   H::S                    # subgroup
   side::Symbol            # says if the transversal is left or right
   X::GapObj               # underlying *right* transversal in GAP
end

GAP.@install GapObj(T::SubgroupTransversal) = T.X

function Base.show(io::IO, ::MIME"text/plain", x::SubgroupTransversal)
  side = is_left(x) ? "Left" : "Right"
  println(io, "$side transversal of length $(length(x)) of")
  io = pretty(io)
  print(io, Indent())
  println(io, Lowercase(), x.H, " in")
  print(io, Lowercase(), x.G)
  print(io, Dedent())
end

function Base.show(io::IO, x::SubgroupTransversal)
  side = is_left(x) ? "Left" : "Right"
  if is_terse(io)
    print(io, "$side transversal of groups")
  else
    print(io, "$side transversal of ")
    io = pretty(io)
    print(terse(io), Lowercase(), x.H, " in ", Lowercase(), x.G)
  end
end

is_left(x::SubgroupTransversal) = x.side == :left

is_right(x::SubgroupTransversal) = x.side == :right

Base.hash(x::SubgroupTransversal, h::UInt) = h # FIXME

Base.length(T::SubgroupTransversal) = index(Int, T.G, T.H)

function Base.getindex(T::SubgroupTransversal, i::Int)
  res = group_element(T.G, GapObj(T)[i])
  if is_left(T)
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
    group(T::SubgroupTransversal)

Return the group `G` that contains all of the elements in `T`.
That is, `T` is a left or right transversal of a subgroup of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = sylow_subgroup(G, 2)[1]
Permutation group of degree 5 and order 8

julia> T = right_transversal(G, H)
Right transversal of length 15 of
  permutation group of degree 5 and order 8 in
  symmetric group of degree 5

julia> group(T) == G
true
```
"""
group(T::SubgroupTransversal) = T.G

"""
    subgroup(T::SubgroupTransversal)

Return the group `H` such that `T` is a (left or right)
transversal of `H`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> T = right_transversal(G, H)
Right transversal of length 20 of
  symmetric group of degree 3 in
  symmetric group of degree 5

julia> subgroup(T) == H
true
```
"""
subgroup(T::SubgroupTransversal) = T.H

"""
    right_transversal(G::GAPGroup, H::GAPGroup; check::Bool=true)

Return a vector containing a complete set of representatives for
the right cosets of `H` in `G`.
This vector is not mutable, and it does not store its entries explicitly,
they are created anew with each access to the transversal.

If `check == false`, do not check whether `H` is a subgroup of `G`.

Use [`right_cosets`](@ref) to compute the G-set of right cosets.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Symmetric group of degree 4

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> T = right_transversal(G, H)
Right transversal of length 4 of
  symmetric group of degree 3 in
  symmetric group of degree 4

julia> collect(T)
4-element Vector{PermGroupElem}:
 ()
 (1,4)
 (1,4,2)
 (1,4,3)
```
"""
function right_transversal(G::T1, H::T2; check::Bool=true) where T1 <: GAPGroup where T2 <: GAPGroup
   if check
     @req GAPWrap.IsSubset(GapObj(G), GapObj(H)) "H is not a subgroup of G"
     _check_compatible(G, H)
   end
   return SubgroupTransversal{T1, T2, eltype(T1)}(G, H, :right,
              GAPWrap.RightTransversal(GapObj(G), GapObj(H)))
end

"""
    left_transversal(G::GAPGroup, H::GAPGroup; check::Bool=true)

Return a vector containing a complete set of representatives for
the left cosets for `H` in `G`.
This vector is not mutable, and it does not store its entries explicitly,
they are created anew with each access to the transversal.

If `check == false`, do not check whether `H` is a subgroup of `G`.

Use [`left_cosets`](@ref) to compute the G-set of left cosets.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Symmetric group of degree 4

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> T = left_transversal(G, H)
Left transversal of length 4 of
  symmetric group of degree 3 in
  symmetric group of degree 4

julia> collect(T)
4-element Vector{PermGroupElem}:
 ()
 (1,4)
 (1,2,4)
 (1,3,4)
```
"""
function left_transversal(G::T1, H::T2; check::Bool=true) where T1 <: GAPGroup where T2 <: GAPGroup
   if check
     @req GAPWrap.IsSubset(GapObj(G), GapObj(H)) "H is not a subgroup of G"
     _check_compatible(G, H)
   end
   return SubgroupTransversal{T1, T2, eltype(T1)}(G, H, :left,
              GAPWrap.RightTransversal(GapObj(G), GapObj(H)))
end


@doc raw"""
    GroupDoubleCoset{T<: Group, S <: GAPGroupElem}

Type of double cosets of subgroups in groups.

For an element $g$ in a group $G$, and two subgroups $H$, $K$ of $G$,
the set $HgK = \{ hgk; h \in H, k \in K \}$ is a $H-K$-double coset in $G$.

- [`group(C::GroupDoubleCoset)`](@ref) returns $G$.

- [`left_acting_group(C::GroupDoubleCoset)`](@ref) returns $H$.

- [`right_acting_group(C::GroupDoubleCoset)`](@ref) returns $K$.

- [`representative(C::GroupDoubleCoset)`](@ref) returns an element
  (the same element for each call) of `C`.

Two double cosets are equal if and only if they contain the same elements.
"""
struct GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
# T=type of the group, S=type of the element
   G::T
   H::GAPGroup
   K::GAPGroup
   repr::S
   X::Ref{GapObj}
   size::Ref{ZZRingElem}
   right_coset_reps::Ref{Dict{GAPGroupElem, Tuple{GAPGroupElem, GAPGroupElem}}}

   function GroupDoubleCoset(G::T, H::GAPGroup, K::GAPGroup, representative::S) where {T<: GAPGroup, S<:GAPGroupElem}
     return new{T, S}(G, H, K, representative, Ref{GapObj}(), Ref{ZZRingElem}(),
                      Ref{Dict{GAPGroupElem, Tuple{GAPGroupElem, GAPGroupElem}}}())
   end
end

GAP.@install function GapObj(C::GroupDoubleCoset)
  if !isassigned(C.X)
    C.X[] = GAPWrap.DoubleCoset(GapObj(C.H), GapObj(representative(C)), GapObj(C.K))
  end
  return C.X[]::GapObj
end

Base.hash(x::GroupDoubleCoset, h::UInt) = h # FIXME
Base.eltype(::Type{GroupDoubleCoset{T,S}}) where {T,S} = S

function ==(x::GroupDoubleCoset, y::GroupDoubleCoset)
   # Avoid creating a GAP object if the result is "obvious"
   x === y && return true
   isassigned(x.size) && isassigned(y.size) && order(x) != order(y) && return false
   return GapObj(x) == GapObj(y)
end

function Base.show(io::IO, ::MIME"text/plain", x::GroupDoubleCoset)
  io = pretty(io)
  println(io, "Double coset of ", Lowercase(), x.H)
  print(io, Indent())
  println(io, "and ", Lowercase(), x.K)
  println(io, "with representative ", representative(x))
  print(io, "in ", Lowercase(), x.G)
  print(io, Dedent())
end

function Base.show(io::IO, x::GroupDoubleCoset)
  if is_terse(io)
    print(io, "Double coset of a group")
  else
    print(io, "Double coset of ")
    io = pretty(io)
    print(terse(io), Lowercase(), x.H,
      " and ", Lowercase(), x.K, " with representative ", representative(x))
  end
end


"""
    double_coset(H::Group, x::GAPGroupElem, K::Group)
    *(H::Group, x::GAPGroupElem, K::Group)

Return the double coset `HxK`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> g = perm(G,[3,4,5,1,2])
(1,3,5,2,4)

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> K = symmetric_group(2)
Symmetric group of degree 2

julia> double_coset(H,g,K)
Double coset of symmetric group of degree 3
  and symmetric group of degree 2
  with representative (1,3,5,2,4)
  in symmetric group of degree 5
```
"""
function double_coset(G::GAPGroup, g::GAPGroupElem, H::GAPGroup)
#T what if g is in some subgroup of a group of which G, H are also a subgroup?
   @req GAPWrap.IsSubset(GapObj(parent(g)), GapObj(G)) "G is not a subgroup of parent(g)"
   @req GAPWrap.IsSubset(GapObj(parent(g)), GapObj(H)) "H is not a subgroup of parent(g)"
   return GroupDoubleCoset(parent(g), G, H, g)
end

Base.:*(H::GAPGroup, g::GAPGroupElem, K::GAPGroup) = double_coset(H,g,K)

"""
    double_cosets(G::GAPGroup, H::GAPGroup, K::GAPGroup; check::Bool=true)

Return a vector of all the double cosets `HxK` for `x` in `G`.
If `check == false`, do not check whether `H` and `K` are subgroups of `G`.

# Examples
```jldoctest
julia> G = symmetric_group(4)
Symmetric group of degree 4

julia> H = symmetric_group(3)
Symmetric group of degree 3

julia> K = symmetric_group(2)
Symmetric group of degree 2

julia> double_cosets(G,H,K)
3-element Vector{GroupDoubleCoset{PermGroup, PermGroupElem}}:
 Double coset of H and K with representative ()
 Double coset of H and K with representative (1,4)
 Double coset of H and K with representative (1,4,3)
```
"""
function double_cosets(G::T, H::GAPGroup, K::GAPGroup; check::Bool=true) where T <: GAPGroup
   if check
      @assert is_subset(H, G) "H is not a subgroup of G"
      @assert is_subset(K, G) "K is not a subgroup of G"
   end
   dcs = GAPWrap.DoubleCosetRepsAndSizes(GapObj(G), GapObj(H), GapObj(K))
   res = Vector{GroupDoubleCoset{T, elem_type(T)}}(undef, length(dcs))
   for i in 1:length(res)
     g = group_element(G, dcs[i][1])
     C = GroupDoubleCoset(G, H, K, g)
     n = dcs[i][2]
     C.size[] = ZZRingElem(n)
     res[i] = C
   end
   return res
end

"""
    order(::Type{T} = ZZRingElem, C::Union{GroupCoset,GroupDoubleCoset})

Return the cardinality of the (double) coset `C`,
as an instance of the type `T`.
"""
order(C::Union{GroupCoset,GroupDoubleCoset}) = order(ZZRingElem, C)

function order(::Type{T}, C::GroupCoset) where T <: IntegerUnion
  return T(GAPWrap.Size(GapObj(C)))
end

function order(::Type{T}, C::GroupDoubleCoset) where T <: IntegerUnion
  if !isassigned(C.size)
    C.size[] = ZZRingElem(GAPWrap.Size(GapObj(C)))
  end
  return T(C.size[])::T
end


"""
    rand(rng::Random.AbstractRNG = Random.GLOBAL_RNG, C::Union{GroupCoset,GroupDoubleCoset})

Return a random element of the (double) coset `C`,
using the random number generator `rng`.
"""
Base.rand(C::Union{GroupCoset,GroupDoubleCoset}) = Base.rand(Random.GLOBAL_RNG, C)

function Base.rand(rng::Random.AbstractRNG, C::Union{GroupCoset,GroupDoubleCoset})
  s = GAPWrap.Random(GAP.wrap_rng(rng), GapObj(C))
  return group_element(C.G, s)
end


"""
    group(C::GroupDoubleCoset)

Return the group `G` that is the parent of all elements in `C`.
That is, `C` is a double coset of two subgroups of `G` in `G`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(3); K = symmetric_group(2);

julia> HgK = double_coset(H, gen(G, 1), K)
Double coset of symmetric group of degree 3
  and symmetric group of degree 2
  with representative (1,2,3,4,5)
  in symmetric group of degree 5

julia> group(HgK) == G
true
```
"""
group(C::GroupDoubleCoset) = C.G

"""
    representative(C::GroupDoubleCoset)

Return an element `x` of the double coset `C` = `HxK`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(3); K = symmetric_group(2);

julia> HgK = double_coset(H, gen(G, 1), K)
Double coset of symmetric group of degree 3
  and symmetric group of degree 2
  with representative (1,2,3,4,5)
  in symmetric group of degree 5

julia> representative(HgK)
(1,2,3,4,5)
```
"""
representative(C::GroupDoubleCoset) = C.repr

"""
    left_acting_group(C::GroupDoubleCoset)

Return `H` if `C` = `HxK`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(3); K = symmetric_group(2);

julia> HgK = double_coset(H, gen(G, 1), K)
Double coset of symmetric group of degree 3
  and symmetric group of degree 2
  with representative (1,2,3,4,5)
  in symmetric group of degree 5

julia> left_acting_group(HgK) == H
true
```
"""
left_acting_group(C::GroupDoubleCoset) = C.H

"""
    right_acting_group(C::GroupDoubleCoset)

Return `K` if `C` = `HxK`.

# Examples
```jldoctest
julia> G = symmetric_group(5)
Symmetric group of degree 5

julia> H = symmetric_group(3); K = symmetric_group(2);

julia> HgK = double_coset(H, gen(G, 1), K)
Double coset of symmetric group of degree 3
  and symmetric group of degree 2
  with representative (1,2,3,4,5)
  in symmetric group of degree 5

julia> right_acting_group(HgK) == K
true
```
"""
right_acting_group(C::GroupDoubleCoset) = C.K


############################################################################
#
# iteration over cosets
#
function Base.in(g::GAPGroupElem, C::GroupCoset)
  if is_right(C)
    return g / representative(C) in acting_group(C)
  else
    return g \ representative(C) in acting_group(C)
  end
end

function Base.in(g::GAPGroupElem, C::GroupDoubleCoset)
  if !isassigned(C.right_coset_reps)
    C.right_coset_reps[] = _right_coset_reps(C)
  end
  canon = GAP.Globals.CanonicalRightCosetElement(GapObj(left_acting_group(C)), GapObj(g))
  return haskey(C.right_coset_reps[], group_element(group(C), canon))
end

"""
   _decompose(C::GroupDoubleCoset, x::GAPGroupElem)

Return `flag, u, v` such that `flag` is `true` if `x` is an element
of `C`, and `false` otherwise.

If `flag = true` then `x = u*g*v` holds where `g` is `representative(C)`,
`u` is an element of `left_acting_group(C)`,
and `v` is an element of `right_acting_group(C)`.

# Examples
```jldoctest
julia> G = symmetric_group(5);

julia> H = sylow_subgroup(G, 2)[1]; K = sylow_subgroup(G, 3)[1];

julia> x = gen(G, 1); C = double_coset(H, x, K);

julia> d = Oscar._decompose(C, x^3)
(true, (1,3)(2,4), (1,3,2))

julia> d[2] * x * d[3] == x^3
true

julia> Oscar._decompose(C, x^2)
(false, (), ())
```
"""
function _decompose(C::GroupDoubleCoset, x::GAPGroupElem)
  G = group(C)
  U = left_acting_group(C)
  if !isassigned(C.right_coset_reps)
    C.right_coset_reps[] = _right_coset_reps(C)
  end
  GAP_x = GapObj(x)
  GAP_y = GAP.Globals.CanonicalRightCosetElement(GapObj(U), GAP_x)
  y = group_element(G, GAP_y)
  haskey(C.right_coset_reps[], y) || return false, one(U), one(right_acting_group(C))
  u, v = C.right_coset_reps[][y]
  return true, group_element(U, GAP_x/GAP_y)*u, v
end

# Compute the data for the (constructive) membership test.
function _right_coset_reps(C::GroupDoubleCoset)
  # `C = UxV`
  G = group(C)
  U = left_acting_group(C)
  V = right_acting_group(C)
  Gx = representative(C)
  x = GapObj(Gx)
  GAP_U = GapObj(U)

  # `C` is a disjoint union of right cosets `Ur`
  # where each representative `r` is canonical and has the form `u_r*x*v_r`,
  # with `u_r` in `U` and `v_r` in `V`.
  # `data` stores `(u_r, v_r)` at the key `r`.
  y = GAP.Globals.CanonicalRightCosetElement(GAP_U, x)
  Gy = group_element(G, y)
  data = Dict{GAPGroupElem, Tuple{GAPGroupElem, GAPGroupElem}}(Gy => (group_element(U, y/x), one(V)))
  orb = IndexedSet([Gy])
  for Gr in orb
    r = GapObj(Gr)
    for v in gens(V)
      rv = r*GapObj(v)
      k = GAP.Globals.CanonicalRightCosetElement(GAP_U, rv)
      Gk = group_element(G, k)
      if !(Gk in orb)
        # `k = u*r*v` for some `u` in `U`
        # `r = u_r*x*v_r` means `u_k = u*u_r`, `v_k = v_r*v`.
        u_r, v_r = data[Gr]
        data[Gk] = (group_element(U, k/rv)*u_r, v_r*v)
        push!(orb, Gk)
      end
    end
  end
  return data
end

Base.IteratorSize(::Type{<:GroupCoset{TG, TH, S}}) where {TG, TH, S} = Base.IteratorSize(TH)

# need this function just for the iterator
Base.length(C::Union{GroupCoset,GroupDoubleCoset}) = order(Int, C)

function Base.iterate(C::GroupCoset)
  return iterate(C, iterate(acting_group(C)))
end

function Base.iterate(C::GroupCoset, state)
  state === nothing && return nothing
  G = group(C)
  if is_right(C)
    res = G(state[1]) * representative(C)
  else
    res = representative(C) * G(state[1])
  end
  return res, iterate(acting_group(C), state[2])
end

Base.IteratorSize(::Type{<:GroupDoubleCoset}) = Base.SizeUnknown()
Base.IteratorSize(::Type{GroupDoubleCoset{PermGroup, PermGroupElem}}) = Base.HasLength()

Base.iterate(G::GroupDoubleCoset) = iterate(G, GAPWrap.Iterator(GapObj(G)))

function Base.iterate(G::GroupDoubleCoset, state)
  GAPWrap.IsDoneIterator(state) && return nothing
  i = GAPWrap.NextIterator(state)::GapObj
  return group_element(G.G, i), state
end
