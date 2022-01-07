export defines_automorphism

@attributes TorQuadMod   # TODO: remove as soon as Hecke is patched

AutGrpAbTor = Union{AutomorphismGroup{GrpAbFinGen},AutomorphismGroup{TorQuadMod}}
AutGrpAbTorElem = Union{AutomorphismGroupElem{GrpAbFinGen},AutomorphismGroupElem{TorQuadMod}}
AbTorElem = Union{GrpAbFinGenElem,TorQuadModElem}

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
  Aindep = abelian_group(fmpz[G.Order(g) for g in gensindep])

  imgs = [Vector{fmpz}(G.IndependentGeneratorExponents(Agap.X, a.X)) for a in gens(Agap)]
  A2_to_Aindep = hom(A2, Aindep, elem_type(Aindep)[Aindep(e) for e in imgs])
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


function apply_automorphism(f::AutGrpAbTorElem, x::AbTorElem, check=true)
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
 
(f::AutGrpAbTorElem)(x::AbTorElem)  = apply_automorphism(f, x, true)
Base.:^(x::AbTorElem,f::AutGrpAbTorElem) = apply_automorphism(f, x, true)

# the _as_subgroup function needs a redefinition
# to pass on the to_gap and to_oscar attributes to the subgroup
function _as_subgroup(aut::AutGrpAbTor, subgrp::GapObj, ::Type{S}) where S
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
function hom(f::AutGrpAbTorElem)
  A = domain(f)
  imgs = [f(a) for a in gens(A)]
  return hom(A, A, imgs)
end


function (aut::AutGrpAbTor)(f::Union{GrpAbFinGenMap,TorQuadModMor};check=false)
  !check || (domain(f) === codomain(f) === domain(aut) && isbijective(f)) || error("Map does not define an automorphism of the abelian group.")
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


function (aut::AutGrpAbTor)(M::fmpz_mat; check=false)
  !check || defines_automorphism(domain(aut),M) || error("Matrix does not define an automorphism of the abelian group.")
  return aut(hom(domain(aut),domain(aut),M); check=check)
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




################################################################################
#
#   Special functions for orthogonal groups of torsion quadratic modules
#
################################################################################


"""
    _orthogonal_group(T::TorQuadMod, gensOT::Vector{fmpz_mat}) -> AutomorphismGroup{TorQuadMod}

Return the subgroup of the orthogonal group of `G` generated by `gensOT`.
"""
function _orthogonal_group(T::TorQuadMod, gensOT::Vector{fmpz_mat}; check=true)
  Ggap, to_gap, to_oscar = _isomorphic_gap_group(abelian_group(T))
  A = abelian_group(T)
  function toA(x)
    return A(x)
  end
  function toT(x)
    return T(x)
  end
  T_to_A = Hecke.map_from_func(toA, T, A )
  A_to_T = Hecke.map_from_func(toT, A, T)
  to_oscar = compose(to_oscar, A_to_T)
  to_gap = compose(T_to_A, to_gap)
  AutGAP = GAP.Globals.AutomorphismGroup(Ggap.X)
  ambient = AutomorphismGroup{typeof(T)}(AutGAP, T)
  set_attribute!(ambient,:to_gap => to_gap)
  set_attribute!(ambient,:to_oscar => to_oscar)
  gens_aut = GapObj([ambient(g, check=check).X for g in gensOT])  # performs the checks
  if check
    # expensive for large groups
    subgrp_gap =GAP.Globals.Subgroup(ambient.X, gens_aut)
  else
    subgrp_gap =GAP.Globals.SubgroupNC(ambient.X, gens_aut)
  end
  aut = AutomorphismGroup{typeof(T)}(subgrp_gap, T)
  set_attribute!(aut,:to_gap => to_gap)
  set_attribute!(aut,:to_oscar => to_oscar)
  return aut
end

function Base.show(io::IO, aut::AutomorphismGroup{TorQuadMod})
  T = domain(aut)
  print(IOContext(io, :compact => true), "Group of isometries of ", T , " generated by ", length(gens(aut)), " elements")
end


"""
    matrix(f::AutomorphismGroupElem{TorQuadMod}) -> fmpz_mat

Return a matrix inducing `f`.
"""
matrix(f::AutomorphismGroupElem{TorQuadMod}) = hom(f).map_ab.map

"""
    defines_automorphism(G::TorQuadMod, M::fmpz_mat) -> Bool

If `M` defines an endomorphism of `G`, return `true` if `M` defines an automorphism of `G`, else `false`.
"""
function defines_automorphism(G::TorQuadMod, M::fmpz_mat)
  g = hom(G, G, M)
  if !isbijective(g)
    return false
  end
  # check that the form is preserved
  B = gens(G)
  n = length(B)
  for i in 1:n
    if Hecke.quadratic_product(B[i]) != Hecke.quadratic_product(g(B[i]))
      return false
    end
    for j in 1:i-1
      if B[i]*B[j] != g(B[i])*g(B[j])
        return false
      end
    end
  end
  return true
end

function Base.show(io::IO, ::MIME"text/plain", f::AutomorphismGroupElem{T}) where T<:TorQuadMod
  D = domain(parent(f))
  print(IOContext(io, :compact => true), "Isometry of ", D, "\n")
  print(io, "defined by")
  print(io, matrix(f))
end

function Base.show(io::IO, f::AutomorphismGroupElem{T}) where T<:TorQuadMod
  print(io, matrix(f))
end


"""
    orthogonal_group(T::TorQuadMod)  -> AutomorphismGroup{TorQuadMod}

Return the full orthogonal group of this torsion quadratic module.
"""
function orthogonal_group(T::TorQuadMod)
  !isdegenerate(T) || error("T must be non-degenerate to compute the full orthogonal group")
  return get_attribute!(T, :orthogonal_group) do
    N,i = normal_form(T)
    j = inv(i)
    gensOT = _compute_gens(N)
    gensOT = [hom(N, N, g) for g in gensOT]
    gensOT = [compose(compose(j,g),i).map_ab.map for g in gensOT]
    return _orthogonal_group(T, gensOT, check=false)
  end
end
