export id_hom, trivial_morphism, hom, domain, codomain, image, issurjective, isinjective, 
       isinvertible, isbijective, automorphism_group, sub, quo, kernel, cokernel, haspreimage, isisomorphic,
       center, index, centralizer, order


function Base.show(io::IO, x::GAPGroupHomomorphism)
  print(io, "Group homomorphism from \n")
  println(io, domain(x))
  print(io, "to\n")
  println(io, codomain(x))
end

function ==(f::GAPGroupHomomorphism{S,T}, g::GAPGroupHomomorphism{S,T}) where S where T
   return f.map == g.map
end

Base.:*(f::GAPGroupHomomorphism{S, T}, g::GAPGroupHomomorphism{T, U}) where S where T where U = compose(g, f)

function Base.:inv(f::GAPGroupHomomorphism{S,T}) where S where T
   return GAPGroupHomomorphism{T,S}(codomain(f), domain(f), GAP.Globals.InverseGeneralMapping(f.map))
end

order(f::GAPGroupHomomorphism) = GAP.Globals.Order(f.map)

function Base.:^(f::GAPGroupHomomorphism{S,T}, n::Int64) where S where T
   if n==1  return f  end
   if n<0
      if !isinvertible(f)
         throw(ArgumentError("function not invertible"))
      else
         g = GAP.Globals.Inverse(f.map)
         n = -n
      end
   else
      g = f.map
   end
   if  n==1                             # only in this case it is allowed domain(f) != codomain(f)
      return GAPGroupHomomorphism{T,S}(domain(f),codomain(f),g)
   end
   @assert domain(f) == codomain(f)
   return GAPGroupHomomorphism{S,S}(domain(f),codomain(f),g^n)
end

function compose(g::GAPGroupHomomorphism{T, U}, f::GAPGroupHomomorphism{S, T}) where S where T where U
  dom = domain(f)
  cod = codomain(g)
  @assert codomain(f) == domain(g)
  mp = GAP.Globals.CompositionMapping(g.map, f.map)
  return GAPGroupHomomorphism{S, U}(dom, cod, mp)
end

function id_hom(G::Group)
  return hom(G, G, x -> x)
end

function trivial_morphism(G::Group, H::Group)
  return hom(G, H, x -> one(H))
end

function _hom_from_gap_map(G::Group, H::Group, mp::GapObj)
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

function hom(G::Group, H::Group, img::Function)
  
  #I create the gap function from the julia function
  #The julia function is supposed to be defined on GAPGroupElem
  #We need a function defined on the underlying GapObj
  function gap_fun(x::GapObj)
    el = group_element(G, x)
    img_el = img(el)
    return img_el.X  
  end
  mp = GAP.Globals.GroupHomomorphismByFunction(G.X, H.X, GAP.julia_to_gap(gap_fun))
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

function hom(G::Group, H::Group, gensG::Vector, imgs::Vector)
  vgens = GAP.julia_to_gap(GapObj[x.X for x in gensG])
  vimgs = GAP.julia_to_gap(GapObj[x.X for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

function domain(f::GAPGroupHomomorphism)
  return f.domain
end

function codomain(f::GAPGroupHomomorphism)
  return f.codomain
end

(f::GAPGroupHomomorphism)(x::GAPGroupElem) = image(f, x)

function image(f::GAPGroupHomomorphism, x::GAPGroupElem)
  return group_element(codomain(f), GAP.Globals.Image(f.map,x.X))
end

function issurjective(f::GAPGroupHomomorphism)
  return GAP.Globals.IsSurjective(f.map)
end

function isinjective(f::GAPGroupHomomorphism)
  return GAP.Globals.IsInjective(f.map)
end

function isinvertible(f::GAPGroupHomomorphism)
  return GAP.Globals.IsBijective(f.map)
end

function isbijective(f::GAPGroupHomomorphism)
  return GAP.Globals.IsBijective(f.map)
end

################################################################################
#
#  Image, Kernel, Cokernel
#
################################################################################

function kernel(f::GAPGroupHomomorphism)
  K = GAP.Globals.Kernel(f.map)
  return _as_subgroup(K, domain(f))
end

function image(f::GAPGroupHomomorphism)
  K = GAP.Globals.Image(f.map)
  return _as_subgroup(K, codomain(f))
end

function image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: Group where T <: Group
  H1 = GAP.Globals.Image(H.X)
  return _as_subgroup(H1, codomain(f))
end

function cokernel(f::GAPGroupHomomorphism)
  K, mK = image(f)
  return quo(codomain(f), K)
end

function haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem)
  r = GAP.Globals.PreImagesRepresentative(f.map, x.X)
  if r == GAP.Globals.fail
    return false, one(domain(f))
  else
    return true, group_element(domain(f), r)
  end
end

"""
TODO: document this
"""
function preimage(f::GAPGroupHomomorphism{S, T}, H::T) where S <: Group where T <: Group
  H1 = GAP.Globals.PreImage(f.map, H.X)
  return _as_subgroup(H1, domain(f))
end

################################################################################
#
#  Subgroup function
#
################################################################################

function __as_subgroup(H::GapObj, G::T, ::Type{S}) where { T, S }
  function img(x::S)
    return group_element(G, x.X)
  end
  H1 = T(H)
  return H1, hom(H1, G, img)
end


function _as_subgroup(H::GapObj, G::T) where T <: Group
  S = elem_type(G)
  return __as_subgroup(H, G, S)
end

function sub(G::T, elements::Vector{S}) where T <: Group where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(H, G)
end

function sub(L::GAPGroupElem...)
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

###############################################################################
#
#  Index
#
###############################################################################

function index(G::T, H::T) where T <: Group
  i = GAP.Globals.Index(G.X, H.X)
  return GAP.gap_to_julia(i)
end

###############################################################################
#
#  subgroups computation
#
###############################################################################

function normal_subgroups(G::Group)
  nsubs = GAP.gap_to_julia(GAP.Globals.NormalSubgroups(G.X))
  res = Vector{Tuple{typeof(G), GAPGroupHomomorphism{typeof(G), typeof(G)}}}(undef, length(nsubs))
  for i = 1:length(res)
    N = nsubs[i]
    res[i] = _as_subgroup(N, G)
  end
  return res
end

function subgroups(G::Group)
  subs = GAP.gap_to_julia(GAP.Globals.AllSubgroups(G.X))
  res = Vector{Tuple{typeof(G), GAPGroupHomomorphism{typeof(G), typeof(G)}}}(undef, length(subs))
  for i = 1:length(res)
    N = subs[i]
    res[i] = _as_subgroup(N, G)
  end
  return res
end

function center(G::Group)
  Z = GAP.Globals.Center(G.X)
  return _as_subgroup(Z, G)
end

function centralizer(G::Group, H::Group)
  C = GAP.Globals.Centralizer(G.X, H.X)
  return _as_subgroup(C, G)
end

function centralizer(G::Group, x::GAPGroupElem)
  C = GAP.Globals.Centralizer(G.X, x.X)
  return _as_subgroup(C, G)
end


################################################################################
#
#  IsNormal, IsCharacteristic, IsSolvable, IsNilpotent
#
################################################################################

function isnormal(G::T, H::T) where T <: Group
  return GAP.Globals.IsNormal(G.X, H.X)
end

function ischaracteristic(G::T, H::T) where T <: Group
  return GAP.Globals.IsCharacteristicSubgroup(G.X, H.X)
end

function issolvable(G::Group)
  return GAP.Globals.IsSolvable(G.X)
end

function isnilpotent(G::Group)
  return GAP.Globals.IsNilpotent(G.X)
end

################################################################################
#
#  Quotient function
#
################################################################################

function quo(G::FPGroup, elements::Vector{S}) where T <: Group where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  Q=FPGroup((G.X)/elems_in_gap)
  function proj(x::FPGroupElem)
     return group_element(Q,GAP.Globals.MappedWord(x.X,GAP.Globals.GeneratorsOfGroup(G.X), GAP.Globals.GeneratorsOfGroup(Q.X)))
  end
  return Q, hom(G,Q,proj)
end

function quo(G::T, elements::Vector{S}) where T <: Group where S <: GAPGroupElem
  @assert elem_type(G) == S
  elems_in_gap = GAP.julia_to_gap(GapObj[x.X for x in elements])
  H = GAP.Globals.NormalClosure(GAP.Globals.Group(elems_in_gap))
  @assert GAP.Globals.IsNormal(G.X, H)
  H1 = T(H)
  return quo(G, H1)
end

function quo(G::T, H::T) where T <: Group
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

function derived_subgroup(G::Group)
  H = GAP.Globals.DerivedSubgroup(G.X)
  return _as_subgroup(H, G)
end

function derived_series(G::Group)
  L = GAP.Globals.DerivedSeries(G.X)
  res = Vector{Tuple{typeof(G), GAPGroupHomomorphism}}(undef, length(L))
  for i = 1:length(res)
    res[i] = _as_subgroup(L[i], G)    
  end
  return res
end

################################################################################
#
#  IsIsomorphic
#
################################################################################

function isisomorphic(G::Group, H::Group)
  mp = GAP.Globals.IsomorphismGroups(G.X, H.X)
  if mp == GAP.Globals.fail
    return false, trivial_morphism(G, H)
  else
    return true, _hom_from_gap_map(G, H, mp)
  end
end


################################################################################
#
#  Direct Product
# 
################################################################################

function direct_product(G::Group, H::Group, task::Symbol = :sum)
  @assert task in [:prod, :sum, :both, :none]

  GH_GAP = GAP.Globals.DirectProduct(G.X, H.X)
  T = _get_type(GH_GAP)
  GH = T(GH_GAP)
  if task == :sum
    mp1 = GAP.Globals.Embedding(GH_GAP, 1)
    mp2 = GAP.Globals.Embedding(GH_GAP, 2)
    map1 = _hom_from_gap_map(G, GH, mp1)    
    map2 = _hom_from_gap_map(H, GH, mp2) 
    return GH, map1, map2   
  elseif task == :prod
    mp11 = GAP.Globals.Projection(GH_GAP, 1)
    mp21 = GAP.Globals.Projection(GH_GAP, 2)
    map11 = _hom_from_gap_map(GH, G, mp11)    
    map21 = _hom_from_gap_map(GH, H, mp21) 
    return GH, map11, map21   
  elseif task == :both
    mp1 = GAP.Globals.Embedding(GH_GAP, 1)
    mp2 = GAP.Globals.Embedding(GH_GAP, 2)
    map1 = _hom_from_gap_map(G, GH, mp1)    
    map2 = _hom_from_gap_map(H, GH, mp2) 
    mp11 = GAP.Globals.Projection(GH_GAP, 1)
    mp21 = GAP.Globals.Projection(GH_GAP, 2)
    map11 = _hom_from_gap_map(GH, G, mp11)    
    map21 = _hom_from_gap_map(GH, H, mp21) 
    return GH, map1, map2, map11, map21
  else
    return GH
  end
end

################################################################################
#
#  Wreath Product
#
################################################################################

#To be done properly
function wreath_product(G::Group, H::PermGroup)
  wGH = GAP.Globals.WreathProduct(G, H)
  T = _get_type(wGH)
  return T(wGH)
end


################################################################################
#
#  Automorphism Group
#
################################################################################

function automorphism_group(G::Group)
  AutGAP = GAP.Globals.AutomorphismGroup(G.X)
  return AutomorphismGroup{typeof(G)}(AutGAP, G)
end

function hom(x::AutomorphismGroupElem)
  A = parent(x)
  G = A.G
  return _hom_from_gap_map(G, G, x)
end

(f::AutomorphismGroupElem)(x::GAPGroupElem) = apply_automorphism(f, x)

function apply_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, x::GAPGroupElem{T}; check = true) where T <: Group
  A = parent(f)
  G = parent(x)
  if check
    @assert A.G == G || GAP.Globals.IN(x.X, A.G.X) "Not in the domain of f!"
  end
  return typeof(x)(G, f.X(x.X))
end

function inner_automorphism(g::GAPGroupElem)
  return _hom_from_gap_map(parent(g), parent(g), GAP.Globals.ConjugatorAutomorphism(parent(g).X, g.X))
end


function isinner_automorphism(f::GAPGroupHomomorphism)
  @assert domain(f) == codomain(f) "Not an automorphism!"
  return GAP.Globals.IsInnerAutomorphism(f.map)
end

function isinvariant(f::GAPGroupHomomorphism, H::Group)
  @assert domain(f) == codomain(f) "Not an automorphism!"
  @assert GAP.Globals.IsSubset(codomain(f).X, H.X) "Not a subgroup of the domain"
  return GAP.Globals.Image(f.map, H.X) == H.X
end

function induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupHomomorphism)
  @assert isinvariant(f, kernel(mH)[1]) "The kernel is not invariant under f!"
  map = GAP.Globals.InducedAutomorphism(mH.map, f.map)
  return _hom_from_gap_map(codomain(mH), codomain(mH), map)
end

function restrict_automorphism(f::GAPGroupHomomorphism, H::Group)
  return _hom_from_gap_map(H, H, f.map)
end

################################################################################
#
#  Conversions between types
#
################################################################################

_get_iso_function(::Type{PermGroup}) = GAP.Globals.IsomorphismPermGroup
_get_iso_function(::Type{FPGroup}) = GAP.Globals.IsomorphismFpGroup
_get_iso_function(::Type{PcGroup}) = GAP.Globals.IsomorphismPcGroup

function isomorphic_group(::Type{T}, G::Group) where T <: Group
  f = _get_iso_function(T)
  mp = f(G.X)
  G1 = T(GAP.Globals.ImagesSource(mp))
  fmap = _hom_from_gap_map(G, G1, mp)
  return G1, fmap
end
