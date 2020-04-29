export right_coset, left_coset, representative, elements, isbicoset, right_cosets, left_cosets, right_transversal, left_transversal, double_coset, acting_domain, left_acting_group, right_acting_group

# T=type of the group, S=type of the element
"""
    GroupCoset{T<: Group, S <: GAPGroupElem}
group coset. It is displayed as `H * x` (right cosets) or `x * H` (left cosets), where `H` is a subgroup of a group `G` and `x` is an element of `G`. Two cosets are equal if, and only if, they are both left (resp. right) and they contain the same elements.
"""
mutable struct GroupCoset{T<: GAPGroup, S <: GAPGroupElem} 
   X::T                    # big group containing the subgroup and the element
   H::T                    # subgroup
   repr::S                 # element
   side::Symbol            # says if the coset is left or right
   coset::GapObj           # GapObj(H*repr)
end

function _group_coset(X::GAPGroup, H::GAPGroup, repr::GAPGroupElem, side::Symbol, coset::GapObj)
  return GroupCoset{typeof(X), typeof(repr)}(X, H, repr, side, coset)
end

function ==(x::GroupCoset, y::GroupCoset)
   return x.coset == y.coset && x.side == y.side
end

"""
    right_coset(H::Group, g::GAPGroupElem)
return the coset `Hg`.
"""
function right_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAP.Globals.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :right, GAP.Globals.RightCoset(H.X,g.X))
end

"""
    right_coset(H::Group, g::GAPGroupElem)
return the coset `gH`.
"""
function left_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAP.Globals.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :left, GAP.Globals.RightCoset(GAP.Globals.ConjugateSubgroup(H.X,GAP.Globals.Inverse(g.X)),g.X))
end

function show(io::IO, x::GroupCoset)
   if x.side == :right
      print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.H.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.repr.X)))
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.repr.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.H.X)))
   end
   return nothing
end

"""
    acting_domain(C::GroupCoset)
if `C` = `Hx` or `xH`, returns `H`.
"""
acting_domain(C::GroupCoset) = C.H

"""
    representative(C::GroupCoset)
if `C` = `Hx` or `xH`, returns `x`.
"""
representative(C::GroupCoset) = C.repr

function elements(C::GroupCoset)
  L = GAP.Globals.AsList(C.coset)
  l = Vector{elem_type(C.X)}(undef, length(L))
  for i = 1:length(l)
    l[i] = group_element(C.X, L[i])
  end
  return l
end

"""
    isbicoset(C::GroupCoset)
returns whether `C` is simultaneously a right coset and a left coset for the same subgroup `H`.
"""
isbicoset(C::GroupCoset) = GAP.Globals.IsBiCoset(C.coset)

"""
    right_cosets(G::Group, H::Group)
returns the array of the right cosets of `H` in `G`.
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
returns the array of the left cosets of `H` in `G`.
"""
function left_cosets(G::GAPGroup, H::GAPGroup)
  #L1 = GAP.Globals.RightCosets(G.X, H.X)
  #L = [GAP.Globals.RightCoset(GAP.Globals.ConjugateSubgroup(H.X,GAP.Globals.Representative(L1[i])), GAP.Globals.Representative(L1[i])) for i in 1:length(L1)]
  T = left_transversal(G,H)
  L = [left_coset(H,t) for t in T]
  l = Vector{GroupCoset{typeof(G), elem_type(G)}}(undef, length(L))
  for i = 1:length(l)
    l[i] = _group_coset(G, H, group_element(G, T[i].X), :left, L[i].coset)
  end
  return l
end

"""
    right_transversal(G::T, H::T) where T<: Group
returns an array containing a complete set of representatives for right cosets for `H`.
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
returns an array containing a complete set of representatives for left cosets for `H`.
"""
function left_transversal(G::T, H::T) where T<: GAPGroup
   return [x^-1 for x in right_transversal(G,H)]
end

function Base.iterate(G::GroupCoset)
  L=GAP.Globals.Iterator(G.coset)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return group_element(G.X, i), L
end

function Base.iterate(G::GroupCoset, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return group_element(G.X, i), state
end



"""
    GroupDoubleCoset{T<: Group, S <: GAPGroupElem}
group double coset. It is displayed as `H * x * K`, where `H` and `K` are subgroups of a group `G` and `x` is an element of `G`. Two double cosets are equal if, and only if, they contain the same elements.
"""
# T=type of the group, S=type of the element
mutable struct GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
   X::T
   G::T
   H::T
   repr::S
   coset::GapObj
end

function ==(x::GroupDoubleCoset, y::GroupDoubleCoset)
   return x.coset == y.coset
end

Base.show(io::IO, x::GroupDoubleCoset) =  print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.G.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.repr.X))*" * "*GAP.gap_to_julia(GAP.Globals.StringView(x.H.X)))

"""
    double_coset(H::Group, x::GAPGroupElem, K::Group)
returns the double coset `HxK`.
"""
function double_coset(G::GAPGroup, g::GAPGroupElem, H::GAPGroup)
   if !GAP.Globals.IsSubset(parent(g).X,G.X)
      throw(ArgumentError("G is not a subgroup of parent(g)"))
   end
   if !GAP.Globals.IsSubset(parent(g).X,H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return GroupDoubleCoset(parent(g),G,H,g,GAP.Globals.DoubleCoset(G.X,g.X,H.X))
end

function elements(C::GroupDoubleCoset)
   L=GAP.gap_to_julia(GAP.Globals.AsList(C.coset))
   return elem_type(C.X)[group_element(C.X,x) for x in L]
end

order(C::Union{GroupCoset,GroupDoubleCoset}) = GAP.Globals.Size(C.coset)
Base.length(C::Union{GroupCoset,GroupDoubleCoset}) = GAP.Globals.Size(C.coset)

Base.rand(C::Union{GroupCoset,GroupDoubleCoset}) = group_element(C.X, GAP.Globals.Random(C.coset))

"""
    representative(C::GroupDoubleCoset)
if `C` = `HxK`, returns `x`.
"""
representative(C::GroupDoubleCoset) = C.repr

"""
    left_acting_group(C::GroupDoubleCoset)
if `C` = `HxK`, returns `H`
"""
left_acting_group(C::GroupDoubleCoset) = C.G

"""
    left_acting_group(C::GroupDoubleCoset)
if `C` = `HxK`, returns `K`
"""
right_acting_group(C::GroupDoubleCoset) = C.H

function Base.iterate(G::GroupDoubleCoset)
  L=GAP.Globals.Iterator(G.coset)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return group_element(G.X, i), L
end

function Base.iterate(G::GroupDoubleCoset, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return group_element(G.X, i), state
end
