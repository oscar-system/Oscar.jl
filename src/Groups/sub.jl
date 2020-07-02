import Base.intersect

export
    centralizer,
    centre,
    characteristic_subgroups,
    derived_series,
    derived_subgroup,
    embedding,
    index,
    ischaracteristic_subgroup,
    isnilpotent,
    issolvable,
    issupersolvable,
    maximal_normal_subgroups,
    maximal_subgroups,
    minimal_normal_subgroups,
    normal_subgroups,
    order,
    quo,
    sub

################################################################################
#
#  Subgroup function
#
################################################################################

function _as_subgroup_bare(G::T, H::GapObj) where T
  if T==PermGroup
    H1 = T(H, G.deg)
  else
    H1 = T(H)
  end
  return H1
end

function _as_subgroup(G::T, H::GapObj, ::Type{S}) where { T, S }
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::S -> group_element(G, x.X))
end

function _as_subgroup(G::T, H::GapObj) where T <: GAPGroup
  return _as_subgroup(G, H, elem_type(G))
end

function sub(G::T, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::GAPGroupElem...)
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

"""
    issubgroup(G::T, H::T) where T <: GAPGroup
Return whether `H` is a subgroup of `G`, together with the embedding morphism.
"""
function issubgroup(G::T, H::T) where T <: GAPGroup
   if false in [h in G for h in gens(H)]
      return (false, Nothing)
   else
      return (true, _as_subgroup(G, H.X)[2])
   end
end

"""
    embedding(G::T, H::T) where T <: GAPGroup
Return the embedding morphism of `H` into `G`. It throws ERROR if `H` is not a subgroup of `G`.
"""
function embedding(G::T, H::T) where T <: GAPGroup
   a,f = issubgroup(G,H)
   if !a
      throw(ArgumentError("Not a subgroup"))
   else
      return f
   end
end

###############################################################################
#
#  Index
#
###############################################################################

function index(G::T, H::T) where T <: GAPGroup
  i = GAP.Globals.Index(G.X, H.X)
  return GAP.gap_to_julia(i)
end

###############################################################################
#
#  subgroups computation
#
###############################################################################

# convert a GAP list of subgroups into a vector of Julia groups objects
function _as_subgroups(G::T, subs::GapObj) where T <: GAPGroup
  res = Vector{T}(undef, length(subs))
  for i = 1:length(res)
    res[i] = _as_subgroup_bare(G, subs[i])
  end
  return res
end


"""
   normal_subgroups(G::Group)
Return the list of normal subgroups of `G`, together with their embeddings into `G`.
"""
function normal_subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.NormalSubgroups(G.X))
end

"""
   subgroups(G::Group)
Return the list of subgroups of `G`, together with their embeddings into `G`.
"""
function subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.AllSubgroups(G.X))
end

"""
   maximal_subgroups(G::Group)
Return the list of maximal subgroups of `G`, together with their embeddings into `G`.
"""
function maximal_subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.MaximalSubgroups(G.X))
end

"""
   maximal_normal_subgroups(G::Group)
Return the list of maximal normal subgroups of `G`, together with their embeddings into `G`.
"""
function maximal_normal_subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.MaximalNormalSubgroups(G.X))
end

"""
   minimal_normal_subgroups(G::Group)
Return the list of minimal normal subgroups of `G`, together with their embeddings into `G`.
"""
function minimal_normal_subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.MinimalNormalSubgroups(G.X))
end

"""
   characteristic_subgroups(G::Group)
Return the list of characteristic subgroups of `G`, i.e. the subgroups that are invariant under all automorphisms of `G`, together with their embeddings into `G`.
"""
function characteristic_subgroups(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.CharacteristicSubgroups(G.X))
end

"""
    ischaracteristic_subgroup(G::T, H::T) where T <: Group
Return whether `H` is a characteristic subgroup of `G`, i.e. `H` is invariant under all automorphisms of `G`.
"""
function ischaracteristic_subgroup(G::T, H::T) where T <: GAPGroup
   return GAP.Globals.IsCharacteristicSubgroup(G.X,H.X)
end

"""
    centre(G::Group)
Return the centre of `G`, i.e. the subgroup of all `x` in `G` such that `xy`=`yx` for every `y` in `G`, together with its embedding morphism into `G`.
"""
function centre(G::GAPGroup)
  Z = GAP.Globals.Center(G.X)
  return _as_subgroup(G, Z)
end

"""
   centralizer(G::Group, H::Group)
Return the centralizer of `H` in `G`, i.e. the subgroup of all `g` in `G` such that `gh=hg` for every `h` in `H`, together with its embedding morphism into `G`.
"""
function centralizer(G::T, H::T) where T <: GAPGroup
  C = GAP.Globals.Centralizer(G.X, H.X)
  return _as_subgroup(G, C)
end

"""
   centralizer(G::Group, x::GroupElem) 
Return the centralizer of `x` in `G`, i.e. the subgroup of all `g` in `G` such that `gx=xg`, together with its embedding morphism into `G`.
"""
function centralizer(G::GAPGroup, x::GAPGroupElem)
  C = GAP.Globals.Centralizer(G.X, x.X)
  return _as_subgroup(G, C)
end

centraliser = centralizer

################################################################################
#
#  IsNormal, IsCharacteristic, IsSolvable, IsNilpotent
#
################################################################################

function isnormal(G::T, H::T) where T <: GAPGroup
  return GAP.Globals.IsNormal(G.X, H.X)
end

function ischaracteristic(G::T, H::T) where T <: GAPGroup
  return GAP.Globals.IsCharacteristicSubgroup(G.X, H.X)
end

function issolvable(G::GAPGroup)
  return GAP.Globals.IsSolvable(G.X)
end

function isnilpotent(G::GAPGroup)
  return GAP.Globals.IsNilpotent(G.X)
end

function issupersolvable(G::GAPGroup)
   return GAP.Globals.IsSupersolvableGroup(G.X)
end

################################################################################
#
#  Quotient function
#
################################################################################

function quo(G::FPGroup, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  Q=FPGroup((G.X)/elems_in_gap)
  function proj(x::FPGroupElem)
     return group_element(Q,GAP.Globals.MappedWord(x.X,GAP.Globals.GeneratorsOfGroup(G.X), GAP.Globals.GeneratorsOfGroup(Q.X)))
  end
  return Q, hom(G,Q,proj)
end

"""
    quo(G::T, elements::Vector{S})
return the quotient group `G/H` of type ``FPGroup`` (if `T`=``FPGroup``) or ``PcGroup`` (otherwise), where `H` is the normal closure of `elements` in `G`.
"""
function quo(G::T, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.NormalClosure(G.X,GAP.Globals.Group(elems_in_gap))
  @assert GAP.Globals.IsNormal(G.X, H)
  H1 = T(H)
  return quo(G, H1)
end

"""
    quo(G::T, H::T) where T <: Group
return the quotient group `G/H` of type ``PcGroup``, together with the projection `G` -> `G/H`.
"""
function quo(G::T, H::T) where T <: GAPGroup
  mp = GAP.Globals.NaturalHomomorphismByNormalSubgroup(G.X, H.X)
  cod = GAP.Globals.ImagesSource(mp)
  S = elem_type(G)
  S1 = _get_type(cod)
  codom = S1(cod)
  mp_julia = __create_fun(mp, codom, S)
  return codom, hom(G, codom, mp_julia)
end

function __create_fun(mp, codom, ::Type{S}) where S
  function mp_julia(x::S)
    el = GAP.Globals.Image(mp, x.X)
    return group_element(codom, el)
  end
  return mp_julia
end

################################################################################
#
#  Derived subgroup and derived series
#  
################################################################################

function derived_subgroup(G::GAPGroup)
  H = GAP.Globals.DerivedSubgroup(G.X)
  return _as_subgroup(G, H)
end

function derived_series(G::GAPGroup)
  return _as_subgroups(G, GAP.Globals.DerivedSeries(G.X))
end


################################################################################
#
#  Intersection
#
################################################################################

"""
    intersect(V::T...) where T <: Group
    intersect(V::AbstractVector{T}) where T <: Group
If `V` = [`G_1`, ... , `G_n`], return the group intersection `K` of the groups `G_1`, ..., `G_n`, together with the embeddings `K` -> `G_i`.
"""
function intersect(V::T...) where T<:GAPGroup
   L = GAP.julia_to_gap([G.X for G in V])
   K = GAP.Globals.Intersection(L)
   Embds = [_as_subgroup(G, K)[2] for G in V]
   K = _as_subgroup(V[1], K)[1]
   Arr = Tuple(vcat([K],Embds))
   return Arr
end

function intersect(V::AbstractVector{T}) where T<:GAPGroup
   L = GAP.julia_to_gap([G.X for G in V])
   K = GAP.Globals.Intersection(L)
   Embds = [_as_subgroup(G, K)[2] for G in V]
   K = _as_subgroup(V[1], K)[1]
   Arr = Tuple(vcat([K],Embds))
   return Arr
end


################################################################################
#
#  Conversions between types
#
################################################################################

_get_iso_function(::Type{PermGroup}) = GAP.Globals.IsomorphismPermGroup
_get_iso_function(::Type{FPGroup}) = GAP.Globals.IsomorphismFpGroup
_get_iso_function(::Type{PcGroup}) = GAP.Globals.IsomorphismPcGroup

function isomorphic_group(::Type{T}, G::GAPGroup) where T <: GAPGroup
  f = _get_iso_function(T)
  mp = f(G.X)
  G1 = T(GAP.Globals.ImagesSource(mp))
  fmap = _hom_from_gap_map(G, G1, mp)
  return G1, fmap
end
