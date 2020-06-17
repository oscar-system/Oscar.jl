import Base.intersect

export
    automorphism_group,
    centralizer,
    centre,
    characteristic_subgroups,
    codomain,
    cokernel,
    derived_series,
    derived_subgroup,
    domain,
    embedding,
    haspreimage,
    hom,
    id_hom,
    image,
    index,
    inner_automorphism,
    inner_automorphisms_group,
    isbijective,
    ischaracteristic_subgroup,
    isinjective,
    isinner_automorphism,
    isinvariant,
    isinvertible,
    isisomorphic,
    issurjective,
    kernel,
    maximal_normal_subgroups,
    maximal_subgroups,
    minimal_normal_subgroups,
    normal_subgroups,
    order,
    quo,
    restrict_automorphism,
    restrict_homomorphism,
    sub,
    trivial_morphism,
    wreath_product

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

function Base.inv(f::GAPGroupHomomorphism{S,T}) where S where T
   return GAPGroupHomomorphism{T,S}(codomain(f), domain(f), GAP.Globals.InverseGeneralMapping(f.map))
end

order(f::GAPGroupHomomorphism) = GAP.Globals.Order(f.map)

function Base.:^(f::GAPGroupHomomorphism{S,T}, n::Int64) where S where T
   if n==1
     return f
   end
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

function id_hom(G::GAPGroup)
  return hom(G, G, x -> x)
end

function trivial_morphism(G::GAPGroup, H::GAPGroup)
  return hom(G, H, x -> one(H))
end

function _hom_from_gap_map(G::GAPGroup, H::GAPGroup, mp::GapObj)
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

function hom(G::GAPGroup, H::GAPGroup, img::Function)
  
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

function hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector)
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
Base.:^(x::GAPGroupElem,f::GAPGroupHomomorphism) = image(f,x)

"""
    image(f::GAPGroupHomomorphism, x::GAPGroupElem)
    (f::GAPGroupHomomorphism)(x::GAPGroupElem)
Return `f`(`x`).
"""
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


"""
    isinvariant(f::GAPGroupHomomorphism, H::Group)
    isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)
return whether `f`(`H`) == `H`.
"""
function isinvariant(f::GAPGroupHomomorphism, H::GAPGroup)
  @assert domain(f) == codomain(f) "Not an endomorphism!"
  @assert GAP.Globals.IsSubset(codomain(f).X, H.X) "Not a subgroup of the domain"
  return GAP.Globals.Image(f.map, H.X) == H.X
end

"""
    restrict_homomorphism(f::GAPGroupHomomorphism, H::Group)
    restrict_homomorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: Group
if `H` is invariant under `f`, returns the restriction of `f` to `H`; otherwise it return ERROR.
"""
function restrict_homomorphism(f::GAPGroupHomomorphism, H::GAPGroup)
  @assert isinvariant(f,H) "H is not invariant under f!"
  return _hom_from_gap_map(H, H, GAP.Globals.RestrictedMapping(f.map,H.X))
end

################################################################################
#
#  Image, Kernel, Cokernel
#
################################################################################

"""
    kernel(f::GAPGroupHomomorphism)
Return the kernel of `f`, together with its embedding into `domain`(`f`).
"""
function kernel(f::GAPGroupHomomorphism)
  K = GAP.Globals.Kernel(f.map)
  return _as_subgroup(domain(f), K)
end

"""
    image(f::GAPGroupHomomorphism)
Return the image of `f` as subgroup of `codomain`(`f`), together with the embedding homomorphism.
"""
function image(f::GAPGroupHomomorphism)
  K = GAP.Globals.Image(f.map)
  return _as_subgroup(codomain(f), K)
end

"""
    image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: GAPGroup where T <: GAPGroup
Return `f`(`H`), together with the embedding homomorphism into `codomain`(`f`).
"""
function image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: GAPGroup where T <: GAPGroup
  H1 = GAP.Globals.Image(f.map, H.X)
  return _as_subgroup(codomain(f), H1)
end

"""
    cokernel(f::GAPGroupHomomorphism)
Return the cokernel of `f`.
"""
function cokernel(f::GAPGroupHomomorphism)
  K, mK = image(f)
  return quo(codomain(f), K)
end

"""
TODO: document this
"""
function haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem)
  r = GAP.Globals.PreImagesRepresentative(f.map, x.X)
  if r == GAP.Globals.fail
    return false, one(domain(f))
  else
    return true, group_element(domain(f), r)
  end
end

"""
    preimage(f::GAPGroupHomomorphism{S, T}, H::T) where S <: GAPGroup where T <: GAPGroup
If `H` is a subgroup of the codomain of `f`, return the subgroup `f^-1(H)`, together with its embedding homomorphism into the domain of `f`.
"""
function preimage(f::GAPGroupHomomorphism{S, T}, H::T) where S <: GAPGroup where T <: GAPGroup
  H1 = GAP.Globals.PreImage(f.map, H.X)
  return _as_subgroup(domain(f), H1)
end

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
#  IsIsomorphic
#
################################################################################

"""
    isisomorphic(G::Group, H::Group)
Return (`true`,`f`) if `G` and `H` are isomorphic groups, where `f` is a group isomorphism. Otherwise, return (`false`,`f`), where `f` is the trivial homomorphism.
"""
function isisomorphic(G::GAPGroup, H::GAPGroup)
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

function direct_product(G::GAPGroup, H::GAPGroup, task::Symbol = :sum)
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
#  Wreath Product
#
################################################################################

#To be done properly
function wreath_product(G::GAPGroup, H::PermGroup)
  wGH = GAP.Globals.WreathProduct(G.X, H.X)
  T = _get_type(wGH)
  return T(wGH)
end


################################################################################
#
#  Automorphism Group
#
################################################################################

"""
    automorphism_group(G::Group) -> A::AutomorphismGroup{T}
Return the full automorphism group of `G`. If `f` is an object of type ``GAPGroupHomomorphism`` and it is bijective from `G` to itself, then `A(f)` return the embedding of `f` in `A`. 

Elements of `A` can be multiplied with other elements of `A` or by elements of type ``GAPGroupHomomorphism``; in this last case, the result has type ``GAPGroupHomomorphism``.
"""
function automorphism_group(G::GAPGroup)
  AutGAP = GAP.Globals.AutomorphismGroup(G.X)
  return AutomorphismGroup{typeof(G)}(AutGAP, G)
end

function Base.show(io::IO, A::AutomorphismGroup{T}) where T <: GAPGroup
  print(io, "Aut( "* GAP.gap_to_julia(GAP.Globals.StringView(A.G.X)) *" )")
end

# TODO: why (x::AutomorphismGroupElem) does not work?

"""
    hom(f::GAPGroupElem{AutomorphismGroup{T}}) where T
Return the element f of type ``GAPGroupHomomorphism{T,T}``.
"""
function hom(x::GAPGroupElem{AutomorphismGroup{T}}) where T
  A = parent(x)
  G = A.G
  return _hom_from_gap_map(G, G, x.X)
end

(f::GAPGroupElem{AutomorphismGroup{T}})(x::GAPGroupElem{T}) where T <: GAPGroup = apply_automorphism(f, x)
Base.:^(x::GAPGroupElem{T},f::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = apply_automorphism(f, x)

function (A::AutomorphismGroup{T})(f::GAPGroupHomomorphism{T,T}) where T <: GAPGroup
   @assert domain(f)==A.G && codomain(f)==A.G "f not in A"
   @assert isbijective(f) "f not in A"
   return group_element(A, f.map)
end

function apply_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, x::GAPGroupElem{T}, check=true) where T <: GAPGroup
  A = parent(f)
  G = parent(x)
  if check
    @assert A.G == G || GAP.Globals.IN(x.X, A.G.X) "Not in the domain of f!"
  end
  return typeof(x)(G, GAP.Globals.Image(f.X,x.X))
end

Base.:*(f::GAPGroupElem{AutomorphismGroup{T}}, g::GAPGroupHomomorphism) where T = hom(f)*g
Base.:*(f::GAPGroupHomomorphism, g::GAPGroupElem{AutomorphismGroup{T}}) where T = f*hom(g)

"""
    inner_automorphism(g::GAPGroupElem)
Return the inner automorphism in `automorphism_group(parent(g))` defined by `x` -> `x^g`.
"""
function inner_automorphism(g::GAPGroupElem)
  return _hom_from_gap_map(parent(g), parent(g), GAP.Globals.ConjugatorAutomorphism(parent(g).X, g.X))
end

"""
    isinner_automorphism(f::GAPGroupHomomorphism)
    isinner_automorphism(g::AutomorphismGroupElem)
Return whether `f` is an inner automorphism.
"""
function isinner_automorphism(f::GAPGroupHomomorphism)
  @assert domain(f) == codomain(f) "Not an automorphism!"
  return GAP.Globals.IsInnerAutomorphism(f.map)
end

function isinner_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
  return GAP.Globals.IsInnerAutomorphism(f.X)
end

"""
    inner_automorphisms_group(A::AutomorphismGroup{T})
Return the subgroup of `A` of the inner automorphisms.
"""
function inner_automorphisms_group(A::AutomorphismGroup{T}) where T <: GAPGroup
   AutGAP = GAP.Globals.InnerAutomorphismsAutomorphismGroup(A.X)
   return _as_subgroup(A, AutGAP)
end

function isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T<:GAPGroup
  @assert GAP.Globals.IsSubset(parent(f).G.X, H.X) "Not a subgroup of the domain"
  return GAP.Globals.Image(f.X, H.X) == H.X
end

function induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupHomomorphism)
  @assert isinvariant(f, kernel(mH)[1]) "The kernel is not invariant under f!"
  map = GAP.Globals.InducedAutomorphism(mH.map, f.map)
  return _hom_from_gap_map(codomain(mH), codomain(mH), map)
end

"""
    restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: Group
If `H` is invariant under `f`, returns the restriction of `f` to `H` as automorphism of `H`; otherwise it returns ERROR.
"""
function restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T, A=automorphism_group(H)) where T <: GAPGroup
  @assert isinvariant(f,H) "H is not invariant under f!"
  fh = hom(H,H,gens(H), [f(x) for x in gens(H)])
  return A(fh)
end

function restrict_homomorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: GAPGroup
  return restrict_homomorphism(hom(f),H)
end

## the next function needs a redefinition if G is an AutomorphismGroup
function _as_subgroup(G::AutomorphismGroup{T}, H::GapObj, ::Type{S}) where { T, S }
  function img(x::S)
    return group_element(G, x.X)
  end
  H1 = AutomorphismGroup{T}(H,G.G)
  return H1, hom(H1, G, img)
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
