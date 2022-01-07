export
    defines_automorphism


function _isomorphic_gap_group(A::GrpAbFinGen; T=PcGroup)
  # find independent generators
  if isdiagonal(rels(A))
    exponents = diagonal(rels(A))
    A2 = A
    A2_to_A = identity_map(A)
    A_to_A2 = identity_map(A)
  else
    exponents = elementary_divisors(A)
    A2, A2_to_A = snf(A)
    A_to_A2 = inv(A2_to_A)
  end
  # the isomorphic gap group
  Agap = abelian_group(T, exponents)
  G = GAP.Globals
  # the gap IndependentGenerators may differ from
  # the generators even if the generators are independent
  gensindep = G.IndependentGeneratorsOfAbelianGroup(Agap.X)
  Aindep = abelian_group([G.Order(g) for g in gensindep])

  imgs = [Vector{fmpz}(G.IndependentGeneratorExponents(Agap.X, a.X)) for a in gens(Agap)]
  A2_to_Aindep = hom(A2, Aindep, [Aindep(e) for e in imgs])
  Aindep_to_A2 = inv(A2_to_Aindep)
  Aindep_to_A = compose(Aindep_to_A2, A2_to_A)

  function Agap_to_A(a)
      return Aindep_to_A(_gap_to_oscar(a, Aindep))
  end

  function A_to_Agap(a)
      return _oscar_to_gap(A_to_A2(a), Agap)
  end

  to_gap = Hecke.map_from_func(A_to_Agap, A, Agap)
  to_oscar = Hecke.map_from_func(Agap_to_A, Agap, A)

  return Agap, to_gap, to_oscar
end

function _oscar_to_gap(a::GrpAbFinGenElem, B::GAPGroup)
  gensB = gens(B)
  n = length(a.coeff)
  @assert length(gensB) == n
  img = one(B)
  for i in 1:n
    img *= gensB[i]^a[i]
  end
  return img
end

function _gap_to_oscar(a::Oscar.BasicGAPGroupElem, B::GrpAbFinGen)
  A = parent(a)
  exp = Vector{fmpz}(GAP.Globals.IndependentGeneratorExponents(A.X, a.X))
  return B(exp)
end

"""
    automorphism_group(G::GrpAbFinGen) -> AutomorphismGroup{GrpAbFinGen} 

Return the automorphism group of `G`.
"""
function automorphism_group(G::GrpAbFinGen)
  Ggap, to_gap, to_oscar = _isomorphic_gap_group(G)
  AutGAP = GAP.Globals.AutomorphismGroup(Ggap.X)
  aut = AutomorphismGroup{typeof(G)}(AutGAP, G)
  set_attribute!(aut,:to_gap => to_gap)
  set_attribute!(aut,:to_oscar => to_oscar)
  return aut
end

"""
    group(aut::AutomorphismGroup{GrpAbFinGen}) -> GrpAbFinGen

If `aut == automorphism_group(G)`, return `G`.  
""" 
group(aut::AutomorphismGroup{GrpAbFinGen}) = aut.G

function apply_automorphism(f::AutomorphismGroupElem{GrpAbFinGen}, x::GrpAbFinGenElem, check=true)
  aut = parent(f)
  if check
    @assert parent(x) == aut.G "Not in the domain of f!"
  end
  to_gap = get_attribute(aut,:to_gap)
  to_oscar = get_attribute(aut,:to_oscar)
  xgap = to_gap(x)
  A = parent(f)
  domGap = parent(xgap)
  imgap = typeof(xgap)(domGap, GAP.Globals.Image(f.X,xgap.X))
  return to_oscar(imgap)
end
 
(f::AutomorphismGroupElem{GrpAbFinGen})(x::GrpAbFinGenElem)  = apply_automorphism(f, x, true)
Base.:^(x::GrpAbFinGenElem,f::AutomorphismGroupElem{GrpAbFinGen}) = apply_automorphism(f, x, true)

# the _as_subgroup function needs a redefinition
# to pass on the to_gap and to_oscar attributes to the subgroup
function _as_subgroup(aut::AutomorphismGroup{T}, subgrp::GapObj, ::Type{S}) where { T<:GrpAbFinGen, S }
  function img(x::S)
    return group_element(aut, x.X)
  end
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  subgrp1 = AutomorphismGroup{T}(subgrp, aut.G)
  set_attribute!(subgrp1, :to_gap => to_gap)
  set_attribute!(subgrp1, :to_oscar => to_oscar)
  return subgrp1, hom(subgrp1, aut, img)
end

"""
    hom(f::AutomorphismGroupElem{GrpAbFinGen}) -> GrpAbFinGenMap 

Return the element `f` of type `GrpAbFinGenMap`.
"""
function hom(f::AutomorphismGroupElem{GrpAbFinGen}) 
  A = group(parent(f))
  imgs = [f(a) for a in gens(A)]
  return hom(A, A, imgs)
end


function (aut::AutomorphismGroup{GrpAbFinGen})(f::GrpAbFinGenMap)
  (domain(f) === codomain(f) && domain(f) === group(aut) && isbijective(f)) || error("Map does not define an automorphism of the abelian group.")
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  Agap = domain(to_oscar)
  AA = Agap.X
  function img_gap(x)
    a = to_oscar(group_element(Agap,x))
    b = to_gap(f(a))
    return b.X 
  end
  gene = GAP.Globals.GeneratorsOfGroup(AA)
  img = GAP.julia_to_gap([img_gap(a) for a in gene])
  fgap = GAP.Globals.GroupHomomorphismByImagesNC(AA,AA,img)
  return aut(fgap)
end


function (aut::AutomorphismGroup{GrpAbFinGen})(M::fmpz_mat) 
  defines_automorphism(group(aut),M) || error("Matrix does not define an automorphism of the abelian group.")
  return aut(hom(group(aut),group(aut),M))
end

"""
    matrix(f::AutomorphismGroupElem{GrpAbFinGen}) -> fmpz_mat

Return the underlying matrix of `f` as a module homomorphism.
"""
matrix(f::AutomorphismGroupElem{GrpAbFinGen}) = hom(f).map


"""
    defines_automorphism(G::GrpAbFinGen, M::fmpz_mat) -> Bool

If `M` defines an endomorphism of `G`, return `true` if `M` defines an automorphism of `G`, else `false`.
""" 
defines_automorphism(G::GrpAbFinGen, M::fmpz_mat) = isbijective(hom(G,G,M))

