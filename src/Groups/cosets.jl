export
    acting_domain,
    double_coset,
    double_cosets,
    GroupCoset,
    GroupDoubleCoset,
    is_bicoset,
    is_left,
    is_right,
    left_acting_group,
    left_coset,
    left_cosets,
    left_transversal,
    representative,
    right_acting_group,
    right_coset,
    right_cosets,
    right_transversal

# T=type of the group, S=type of the element
"""
    GroupCoset{T<: Group, S <: GAPGroupElem}

Type of group cosets. It is displayed as `H * x` (right cosets) or `x * H`
(left cosets), where `H` is a subgroup of a group `G` and `x` is an element of
`G`. Two cosets are equal if, and only if, they are both left (resp. right)
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


"""
    right_coset(H::Group, g::GAPGroupElem)
    *(H::Group, g::GAPGroupElem)

Return the coset `Hg`.
"""
function right_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAPWrap.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :right, GAP.Globals.RightCoset(H.X,g.X))
end

"""
    left_coset(H::Group, g::GAPGroupElem)
    *(g::GAPGroupElem, H::Group)

Return the coset `gH`.
!!! note
    Since GAP supports right cosets only, the underlying GAP object of
    `left_coset(H,g)` is the right coset `H^(g^-1) * g`.
"""
function left_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAPWrap.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :left, GAP.Globals.RightCoset(GAP.Globals.ConjugateSubgroup(H.X,GAP.Globals.Inverse(g.X)),g.X))
end

function show(io::IO, x::GroupCoset)
   a = String(GAP.Globals.StringViewObj(x.H.X))
   b = String(GAP.Globals.StringViewObj(x.repr.X))
   if x.side == :right
      print(io, "Right coset   $a * $b")
   else
      print(io, "Left coset   $b * $a")
   end
   return nothing
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
   if c.side != :right || d.side != :left
      throw(ArgumentError("Wrong input"))
   end
   return double_coset(c.H, c.repr*d.repr, d.H)
end

"""
    acting_domain(C::GroupCoset)

If `C` = `Hx` or `xH`, return `H`.
"""
acting_domain(C::GroupCoset) = C.H

"""
    representative(C::GroupCoset)

If `C` = `Hx` or `xH`, return `x`.
"""
representative(C::GroupCoset) = C.repr

@deprecate elements(C::GroupCoset) collect(C)

"""
    is_bicoset(C::GroupCoset)

Return whether `C` is simultaneously a right coset and a left coset for the same subgroup `H`.
"""
is_bicoset(C::GroupCoset) = GAPWrap.IsBiCoset(C.X)

"""
    right_cosets(G::Group, H::Group)

Return the vector of the right cosets of `H` in `G`.
"""
function right_cosets(G::GAPGroup, H::GAPGroup)
  L = GAP.Globals.RightCosets(G.X, H.X)
  l = Vector{GroupCoset{typeof(G), elem_type(G)}}(undef, length(L))
  for i = 1:length(l)
    l[i] = _group_coset(G, H, group_element(G, GAP.Globals.Representative(L[i])), :right, L[i])
  end
  return l
end

"""
    left_cosets(G::Group, H::Group)

Return the vector of the left cosets of `H` in `G`.
"""
function left_cosets(G::GAPGroup, H::GAPGroup)
  #L1 = GAP.Globals.RightCosets(G.X, H.X)
  #L = [GAP.Globals.RightCoset(GAP.Globals.ConjugateSubgroup(H.X,GAP.Globals.Representative(L1[i])), GAP.Globals.Representative(L1[i])) for i in 1:length(L1)]
  T = left_transversal(G,H)
  L = [left_coset(H,t) for t in T]
  l = Vector{GroupCoset{typeof(G), elem_type(G)}}(undef, length(L))
  for i = 1:length(l)
    l[i] = _group_coset(G, H, group_element(G, T[i].X), :left, L[i].X)
  end
  return l
end

"""
    right_transversal(G::T, H::T) where T<: Group

Return a vector containing a complete set of representatives for right cosets for `H`.
"""
function right_transversal(G::T, H::T) where T<: GAPGroup
   L = GAP.Globals.RightTransversal(G.X,H.X)
   l = Vector{elem_type(G)}(undef, length(L))
   for i in 1:length(l)
      l[i] = group_element(G,L[i])
   end
   return l
end

"""
    left_transversal(G::T, H::T) where T<: Group

Return a vector containing a complete set of representatives for left cosets for `H`.
"""
function left_transversal(G::T, H::T) where T<: GAPGroup
   return [x^-1 for x in right_transversal(G,H)]
end

Base.IteratorSize(::Type{<:GroupCoset}) = Base.SizeUnknown()
Base.iterate(G::GroupCoset) = iterate(G, GAP.Globals.Iterator(G.X))

function Base.iterate(G::GroupCoset, state)
  if GAPWrap.IsDoneIterator(state)
    return nothing
  end
  i = GAPWrap.NextIterator(state)
  return group_element(G.G, i), state
end



"""
    GroupDoubleCoset{T<: Group, S <: GAPGroupElem}

Group double coset. It is displayed as `H * x * K`, where `H` and `K` are
subgroups of a group `G` and `x` is an element of `G`. Two double cosets are
equal if, and only if, they contain the same elements.
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

function Base.show(io::IO, x::GroupDoubleCoset)
  print(io, String(GAP.Globals.StringViewObj(x.H.X)),
            " * ",
            String(GAP.Globals.StringViewObj(x.repr.X)),
            " * ",
            String(GAP.Globals.StringViewObj(x.K.X)))
end


"""
    double_coset(H::Group, x::GAPGroupElem, K::Group)
    *(H::Group, x::GAPGroupElem, K::Group)

returns the double coset `HxK`.
"""
function double_coset(G::T, g::GAPGroupElem{T}, H::T) where T<: GAPGroup
   # TODO: enforce that G, H have same type
   # TODO: enforce that G, H have common overgroup
   if !GAPWrap.IsSubset(parent(g).X,G.X)
      throw(ArgumentError("G is not a subgroup of parent(g)"))
   end
   if !GAPWrap.IsSubset(parent(g).X,H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return GroupDoubleCoset(parent(g),G,H,g,GAP.Globals.DoubleCoset(G.X,g.X,H.X))
end

Base.:*(H::GAPGroup, g::GAPGroupElem, K::GAPGroup) = double_coset(H,g,K)

"""
    double_cosets(G::T, H::T, K::T; check::Bool=true) where T<: GAPGroup

Return the vector of all the double cosets `HxK` for `x` in `G`.
If `check == false`, do not check whether `H` and `K` are subgroups of `G`.
"""
function double_cosets(G::T, H::T, K::T; check::Bool=true) where T<: GAPGroup
   if !check
      dcs = GAP.Globals.DoubleCosetsNC(G.X,H.X,K.X)
   else
      @assert is_subgroup(G,H)[1] "H is not a subgroup of G"
      @assert is_subgroup(G,K)[1] "K is not a subgroup of G"
      dcs = GAP.Globals.DoubleCosets(G.X,H.X,K.X)
   end
   res = Vector{GroupDoubleCoset{T,elem_type(T)}}(undef, length(dcs))
   for i = 1:length(res)
     dc = dcs[i]
     g = group_element(G, GAP.Globals.Representative(dc))
     res[i] = GroupDoubleCoset(G,H,K,g,dc)
   end
   return res
   #return [GroupDoubleCoset(G,H,K,group_element(G.X,GAP.Globals.Representative(dc)),dc) for dc in dcs]
end

@deprecate elements(C::GroupDoubleCoset) collect(C)

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

if `C` = `HxK`, returns `x`.
"""
representative(C::GroupDoubleCoset) = C.repr

"""
    left_acting_group(C::GroupDoubleCoset)

if `C` = `HxK`, returns `H`
"""
left_acting_group(C::GroupDoubleCoset) = C.H

"""
    right_acting_group(C::GroupDoubleCoset)

if `C` = `HxK`, returns `K`
"""
right_acting_group(C::GroupDoubleCoset) = C.K

Base.IteratorSize(::Type{<:GroupDoubleCoset}) = Base.SizeUnknown()

Base.iterate(G::GroupDoubleCoset) = iterate(G, GAP.Globals.Iterator(G.X))

function Base.iterate(G::GroupDoubleCoset, state)
  if GAPWrap.IsDoneIterator(state)
    return nothing
  end
  i = GAPWrap.NextIterator(state)
  return group_element(G.G, i), state
end

"""
    intersect(V::AbstractVector{Union{T, GroupCoset, GroupDoubleCoset}}) where T <: GAPGroup

Return a vector containing all elements belonging to all groups and cosets
in `V`.
"""
function intersect(V::AbstractVector{Union{T, GroupCoset, GroupDoubleCoset}}) where T <: GAPGroup
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
