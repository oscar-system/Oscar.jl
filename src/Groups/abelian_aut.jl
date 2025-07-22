const AutGrpAbTor = Union{AutomorphismGroup{FinGenAbGroup},AutomorphismGroup{TorQuadModule}}
const AutGrpAbTorElem = Union{AutomorphismGroupElem{FinGenAbGroup},AutomorphismGroupElem{TorQuadModule}}
const AbTorElem = Union{FinGenAbGroupElem,TorQuadModuleElem}

# function _isomorphic_gap_group(A::FinGenAbGroup; T=PcGroup)
function _isomorphic_gap_group(A::FinGenAbGroup; T=SubPcGroup)
  iso = isomorphism(T, A)
  iso2 = inv(iso)
  return codomain(iso), iso, iso2
end

@doc raw"""
    automorphism_group(G::FinGenAbGroup) -> AutomorphismGroup{FinGenAbGroup}

Return the automorphism group of `G`.

# Examples
```jldoctest
julia> A = abelian_group([2, 3, 4, 4, 4]);

julia> automorphism_group(A)
Automorphism group of
  finitely generated abelian group with 5 generators and 5 relations

julia> S, _ = snf(A);

julia> automorphism_group(S)
Automorphism group of
  Z/2 x (Z/4)^2 x Z/12
```
"""
function automorphism_group(G::FinGenAbGroup)
  Ggap, to_gap, to_oscar = _isomorphic_gap_group(G)
  AutGAP = GAPWrap.AutomorphismGroup(Ggap.X)
  aut = AutomorphismGroup(AutGAP, G)
  set_attribute!(aut, :to_gap => to_gap, :to_oscar => to_oscar)
  return aut
end

function apply_automorphism(f::AutGrpAbTorElem, x::AbTorElem, check::Bool=true)
  aut = parent(f)
  if check
    @assert parent(x) == aut.G "Not in the domain of f!"
  end
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  xgap = to_gap(x)
  A = parent(f)
  domGap = parent(xgap)
  imgap = typeof(xgap)(domGap, GAPWrap.Image(f.X,xgap.X))
  return to_oscar(imgap)::typeof(x)
end

(f::AutGrpAbTorElem)(x::AbTorElem)  = apply_automorphism(f, x, true)
Base.:^(x::AbTorElem,f::AutGrpAbTorElem) = apply_automorphism(f, x, true)

# the _as_subgroup function needs a redefinition
# to pass on the to_gap and to_oscar attributes to the subgroup
function _as_subgroup(aut::AutomorphismGroup{S}, subgrp::GapObj) where S <: Union{TorQuadModule,FinGenAbGroup}
  function img(x::AutomorphismGroupElem{S})
    return group_element(aut, x.X)
  end
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  subgrp1 = AutomorphismGroup{S}(subgrp, aut.G)
  set_attribute!(subgrp1, :to_gap => to_gap, :to_oscar => to_oscar)
  return subgrp1, hom(subgrp1, aut, img)
end

@doc raw"""
    hom(f::AutomorphismGroupElem{FinGenAbGroup}) -> FinGenAbGroupHom

Return the element `f` of type `FinGenAbGroupHom`.

# Examples
```jldoctest
julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [2//3 0; 0 2//9]));

julia> OT = orthogonal_group(T)
Orthogonal group of
  finite quadratic module: Z/3 x Z/9 -> Q/2Z
with 3 generators

julia> f = first(gens(OT))
Isometry of
  finite quadratic module: Z/3 x Z/9 -> Q/2Z
with matrix representation
  [2   0]
  [0   1]

julia> hom(f)
Map
  from finite quadratic module: Z/3 x Z/9 -> Q/2Z
  to finite quadratic module: Z/3 x Z/9 -> Q/2Z
```
"""
function hom(f::AutGrpAbTorElem)
  A = domain(f)
  imgs = elem_type(A)[f(a) for a in gens(A)]
  return hom(A, A, imgs)
end

function (aut::AutGrpAbTor)(f::Union{FinGenAbGroupHom,TorQuadModuleMap}; check::Bool=true)
  !check || (domain(f) === codomain(f) === domain(aut) && is_bijective(f)) || error("Map does not define an automorphism of the abelian group.")
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  Agap = domain(to_oscar)
  AA = Agap.X
  function img_gap(x)
    a = to_oscar(group_element(Agap,x))
    b = to_gap(f(a))
    return b.X
  end
  gene = GAPWrap.GeneratorsOfGroup(AA)
  img = GAP.Obj([img_gap(a) for a in gene])
  fgap = GAP.Globals.GroupHomomorphismByImagesNC(AA,AA,img)
  !check || fgap in aut.X || error("Map does not define an element of the group")
  return aut(fgap)
end

function (aut::AutGrpAbTor)(M::ZZMatrix; check::Bool=true)
  !check || defines_automorphism(domain(aut),M) || error("Matrix does not define an automorphism of the abelian group.")
  return aut(hom(domain(aut),domain(aut),M); check=check)
end

function (aut::AutGrpAbTor)(g::MatrixGroupElem{QQFieldElem, QQMatrix}; check::Bool=true)
  L = relations(domain(aut))
  if check
    B = basis_matrix(L)
    @assert can_solve(B, B*matrix(g),side=:left)
  end
  T = domain(aut)
  g = hom(T, T, elem_type(T)[T(lift(t)*matrix(g)) for t in gens(T)])
  return aut(g, check = false)
end

@doc raw"""
    matrix(f::AutomorphismGroupElem{FinGenAbGroup}) -> ZZMatrix

Return the underlying matrix of `f` as a module homomorphism.

# Examples
```jldoctest
julia> A = abelian_group([3, 9, 12]);

julia> G = automorphism_group(A);

julia> f = first(gens(G))
Automorphism of
  finitely generated abelian group with 3 generators and 3 relations
with matrix representation
  [1   0   0]
  [0   1   0]
  [0   0   7]

julia> matrix(f)
[1   0   0]
[0   1   0]
[0   0   7]
```
"""
matrix(f::AutomorphismGroupElem{FinGenAbGroup}) = matrix(hom(f))


@doc raw"""
    defines_automorphism(G::FinGenAbGroup, M::ZZMatrix) -> Bool

If `M` defines an endomorphism of `G`, return `true` if `M` defines an
automorphism of `G`, else `false`.
"""
defines_automorphism(G::FinGenAbGroup, M::ZZMatrix) = is_bijective(hom(G, G, M))

################################################################################
#
#   Special functions for orthogonal groups of torsion quadratic modules
#
################################################################################

"""
    _orthogonal_group(T::TorQuadModule, gensOT::Vector{ZZMatrix}) -> AutomorphismGroup{TorQuadModule}

Return the subgroup of the orthogonal group of `G` generated by `gensOT`.
"""
function _orthogonal_group(T::TorQuadModule, gensOT::Vector{ZZMatrix}; check::Bool=true)
  A = abelian_group(T)
  As, AstoA = snf(A)
  Ggap, to_gap, to_oscar = _isomorphic_gap_group(As)
  function toAs(x)
    return AstoA\A(x)
  end
  function toT(x)
    return T(AstoA(x))
  end
  T_to_As = MapFromFunc(T, As, toAs, toT)
  to_oscar = compose(to_oscar, inv(T_to_As))
  to_gap = compose(T_to_As, to_gap)
  AutGAP = GAPWrap.AutomorphismGroup(Ggap.X)
  ambient = AutomorphismGroup(AutGAP, T)
  set_attribute!(ambient, :to_gap => to_gap, :to_oscar => to_oscar)
  gens_aut = GapObj([ambient(g, check=check).X for g in gensOT])  # performs the checks
  if check
    # expensive for large groups
    subgrp_gap =GAP.Globals.Subgroup(ambient.X, gens_aut)
  else
    subgrp_gap =GAP.Globals.SubgroupNC(ambient.X, gens_aut)
  end
  aut = AutomorphismGroup(subgrp_gap, T)
  set_attribute!(aut, :to_gap => to_gap, :to_oscar => to_oscar)
  return aut
end

function Base.show(io::IO, ::MIME"text/plain", OT::AutomorphismGroup{TorQuadModule})
  io = pretty(io)
  println(io, "Orthogonal group of", Indent())
  println(io, Lowercase(), OT.G)
  print(io, Dedent(), "with ", ItemQuantity(ngens(OT), "generator"))
end
                                                                                                                                                                                        
function Base.show(io::IO, OT::AutomorphismGroup{TorQuadModule})
  if is_terse(io)
    print(io, "Orthogonal group")
  else
    io = pretty(io)
    print(io, "Orthogonal group of ", Lowercase(), OT.G)
  end
end 

@doc raw"""
    matrix(f::AutomorphismGroupElem{TorQuadModule}) -> ZZMatrix

Return a matrix inducing `f`.

# Examples
```jldoctest
julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [1//12 0; 0 2//9]));

julia> OT = orthogonal_group(T)
Orthogonal group of
  finite quadratic module: Z/3 x Z/36 -> Q/2Z
with 4 generators

julia> f = first(gens(OT))
Isometry of
  finite quadratic module: Z/3 x Z/36 -> Q/2Z
with matrix representation
  [1    0]
  [0   19]

julia> matrix(f)
[1    0]
[0   19]
```
"""
matrix(f::AutomorphismGroupElem{TorQuadModule}) = matrix(hom(f))

@doc raw"""
    defines_automorphism(G::TorQuadModule, M::ZZMatrix) -> Bool

If `M` defines an endomorphism of `G`, return `true` if `M` defines an automorphism of `G`, else `false`.
"""
function defines_automorphism(G::TorQuadModule, M::ZZMatrix)
  g = hom(G, G, M)
  if !is_bijective(g)
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

function Base.show(io::IO, ::MIME"text/plain", f::AutomorphismGroupElem{TorQuadModule})
  D = domain(parent(f))
  io = pretty(io)
  println(io, "Isometry of")
  println(IOContext(io, :compact => true), Indent(), Lowercase(), D, Dedent())
  println(io, "with matrix representation")
  print(io, Indent())
  show(io, MIME"text/plain"(), matrix(f))
  print(io, Dedent())
end

function Base.show(io::IO, f::AutomorphismGroupElem{TorQuadModule})
  print(io, matrix(f))
end


@doc raw"""
    orthogonal_group(T::TorQuadModule)  -> AutomorphismGroup{TorQuadModule}

Return the full orthogonal group of this torsion quadratic module.

# Examples
```jldoctest
julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [1//4 1//2; 1//2 3//4]))
Finite quadratic module
  over integer ring
Abelian group: (Z/4)^2
Bilinear value module: Q/Z
Quadratic value module: Q/2Z
Gram matrix quadratic form:
[1//4   1//2]
[1//2   3//4]

julia> OT = orthogonal_group(T)
Orthogonal group of
  finite quadratic module: (Z/4)^2 -> Q/2Z
with 2 generators

julia> order(OT)
4
```
"""
@attr AutomorphismGroup{TorQuadModule} function orthogonal_group(T::TorQuadModule)
  if is_trivial(abelian_group(T))
    return _orthogonal_group(T, ZZMatrix[identity_matrix(ZZ, ngens(T))], check = false)
  elseif is_semi_regular(T)
    # if T is semi-regular, it is isometric to its normal form for which
    # we know how to compute the isometries.
    N, i = normal_form(T)
    j = inv(i)
    gensOT = _compute_gens(N)
    gensOT = TorQuadModuleMap[hom(N, N, g) for g in gensOT]
    gensOT = ZZMatrix[compose(compose(i,g),j).map_ab.map for g in gensOT]
    unique!(gensOT)
    length(gensOT) > 1 ? filter!(m -> !isone(m), gensOT) : nothing
  elseif iszero(gram_matrix_quadratic(T))
    # in that case, we don't have any conditions regarding the
    # quadratic form, so we have all automorphisms coming
    # from the underlying abelian group
    gensOT = [matrix(g) for g in gens(automorphism_group(abelian_group(T)))]
  else
    # if T is not semi-regular, we distinghuish the cases whether or not
    # it splits its radical quadratic
    i = radical_quadratic(T)[2]
    gensOT = has_complement(i)[1] ? _compute_gens_split_degenerate(T) : _compute_gens_non_split_degenerate(T)
  end
  return _orthogonal_group(T, gensOT, check=false)
end

@doc raw"""
    embedding_orthogonal_group(i::TorQuadModuleMap) -> GAPGroupHomomorphism

Given an embedding $i\colon A \to D$ between two torsion quadratic modules,
such that `A` admits a complement `B` in $D \cong A \oplus B$ to which it is
orthogonal, return the embedding $O(A) \to O(D)$ obtained by extending the
isometries of `A` by the identity on `B`.
"""
function embedding_orthogonal_group(i::TorQuadModuleMap)
  @req is_injective(i) "i must be injective"
  ok, j = has_complement(i)
  @req ok "The domain of i must have a complement in the codomain"
  @req all(v -> i(v[1])*j(v[2]) == 0, Hecke.cartesian_product_iterator([gens(domain(i)), gens(domain(j))], inplace=true)) "The domain of i and its complement must be in orthogonal direct sum"
  A = domain(i)
  B = domain(j)
  D = codomain(i)
  OD = orthogonal_group(D)
  OA = orthogonal_group(A)

  # D = A+B
  gene = data.(union(i.(gens(A)), j.(gens(B))))
  geneOAinOD = elem_type(OD)[]
  for f in gens(OA)
    imgf = data.(union(i.(f.(gens(A))), j.(gens(B))))
    fab = hom(abelian_group(D), abelian_group(D), gene, imgf)
    fD = OD(hom(D, D, fab.map))
    push!(geneOAinOD, fD)
  end
  OAtoOD = hom(OA, OD, geneOAinOD, check = false)
  return OAtoOD::GAPGroupHomomorphism{AutomorphismGroup{TorQuadModule}, AutomorphismGroup{TorQuadModule}}
end

###############################################################################
#
#  Action on injections
#
###############################################################################

@doc raw"""
    is_invariant(f::TorQuadModuleMap, i::TorQuadModuleMap) -> Bool

Given an abelian group morphism $i\colon S \to T$ form a torsion quadratic module
`S` to a torsion quadratic module `T`, and an abelian group endomorphism `f`
of `T`, return whether `f` preserves the image of `i` in `T`, i.e. whether
$f(i(s)) \in i(S)$ for all $s \in S$.
"""
function is_invariant(f::TorQuadModuleMap, i::TorQuadModuleMap)
  @req domain(f) === codomain(f) === codomain(i) "f must be an endomorphism of the target of i"
  U = domain(i)
  for a in gens(U)
    b = f(i(a))
    has_preimage_with_preimage(i, b)[1] || return false
  end
  return true
end

@doc raw"""
    is_invariant(f::AutomorphismGroupElem{TorQuadModule}, i::TorQuadModuleMap)
                                                                        -> Bool

Given an abelian group morphism $i\colon S \to T$ from a torsion quadratic module
`S` to a torsion quadratic module `T`, and an automorphism `f` of `T`, return
whether `f` preserves the image of `i` in `T`, i.e. whether
$f(i(s)) \in i(S)$ for all $s \in S$.
"""
function is_invariant(f::AutomorphismGroupElem{TorQuadModule}, i::TorQuadModuleMap)
  @req domain(parent(f)) === codomain(i) "f must be an automorphism of the target of i"
  return is_invariant(hom(f), i)
end

@doc raw"""
    is_invariant(G::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap)
                                                                        -> Bool

Given an abelian group morphism $i\colon S \to T$ from a torsion quadratic module
`S` to a torsion quadratic module `T`, and a group `G` of automorphisms of `T`,
return whether the image of `i` in `T` is preserved by every element in
`G`, i.e. whether $f(i(s)) \in i(S)$ for all $s \in S$ and all $f \in G$
"""
function is_invariant(G::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap)
  @req domain(G) === codomain(i) "G must consist of automorphisms of the target of i"
  return all(f -> is_invariant(f, i), gens(G))
end

@doc raw"""
    restrict_endomorphism(f::TorQuadModule, i::TorQuadModuleMap)
                                                            -> TorQuadModuleMap

Given an abelian group embedding $i\colon S \to T$ of a torsion quadratic
module `S` in a torsion quadratic module `T`, and an abelian group endomorphism
of `T`, return the restriction of `f` to `S`.

If `S` is not invariant under the action of `f`, then an error is thrown.
"""
function restrict_endomorphism(f::TorQuadModuleMap, i::TorQuadModuleMap; check::Bool = true)
  @req !check || is_injective(i) "i must be an injection"
  @req domain(f) === codomain(f) === codomain(i) "f must be an endomorphism of the target of i"
  imgs = TorQuadModuleElem[]
  U = domain(i)
  for a in gens(U)
    b = f(i(a))
    ok, c = has_preimage_with_preimage(i, b)
    @req ok "The domain of i is not invariant under the action of f"
    push!(imgs, c)
  end
  return hom(U, U, imgs)
end

@doc raw"""
    restrict_automorphism(f::AutomorphismGroupElem{TorQuadModule},
                          i::TorQuadModuleMap)   -> TorQuadModuleMap

Given an abelian group embedding $i\colon S \to T$ of a torsion quadratic
module `S` in a torsion quadratic module `T`, and an automorphism `f` of `T`,
return the restriction of `f` to `S`.

If `S` is not invariant under the action of `f`, then an error is thrown.
"""
function restrict_automorphism(f::AutomorphismGroupElem{TorQuadModule}, i::TorQuadModuleMap; check::Bool = true)
  @req !check || is_injective(i) "i must be an injection"
  @req domain(parent(f)) === codomain(i) "f must be an automorphism of the target of i"
  return restrict_endomorphism(hom(f), i, check = false)
end

@doc raw"""
    restrict_automorphism_group(G::AutomorphismGroup{TorQuadModule},
                                i::TorQuadModuleMap; check::Bool = true)
                                            -> AutomorphismGroup{TorQuadModule},
                                               GAPGroupHomomorphism

Given an embedding $i\colon S \to T$ of a torsion quadratic module `S` in a
torsion quadratic module `T`, and a group `G` of automorphisms of `T`, return
the group of automorphisms `H` of `S` generated by the restrictions of the
elements in `G` to `S`, together with the restriction map $G \to H$.

If `S` is not invariant under the action of `G`, then an error is thrown.

By default, the function checks whether `i` is injective and whether `i`
is a torsion quadratic module morphism. One can disable these checks
by setting `check = false`.
"""
function restrict_automorphism_group(G::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap; check::Bool = true)

  if check
    @req is_injective(i) "i must be an injection"
    @req modulus_bilinear_form(domain(i)) == modulus_bilinear_form(codomain(i)) "The bilinear forms of the domain and the codomain of i must take values in the same torsion module"
    @req modulus_quadratic_form(domain(i)) == modulus_quadratic_form(codomain(i)) "The quadratic forms of the domain and the codomain of i must take values in the same torsion module"
    @req all(a -> all(b -> a*b == i(a)*i(b), gens(domain(i))), gens(domain(i))) "i must preserve bilinear products"
    @req all(a -> Hecke.quadratic_product(a) == Hecke.quadratic_product(i(a)), gens(domain(i))) "i must preserve quadratic products"
  end

  @req domain(G) === codomain(i) "G must consist of automorphisms of the target of i"
  restr = ZZMatrix[]
  for f in gens(G)
    g = try restrict_automorphism(f, i, check = false)
        catch e throw(ArgumentError("The domain of i is not invariant under the action of G"))
        end
    push!(restr, g.map_ab.map)
  end
  H = _orthogonal_group(domain(i), unique(restr), check = check)
  res = hom(G, H, gens(G), H.(restr), check = check)
  return H, res
end

###############################################################################
#
#  Action on submodules
#
###############################################################################

@doc raw"""
    is_conjugate_with_data(O::AutomorphismGroup{TorQuadModule},
                           i::TorQuadModuleMap,
                           j::TorQuadModuleMap)
                                  -> Bool, AutomorphismGroupElem{TorQuadModule}

Return whether the images of `i` and `j` in the domain of `O` lies in the same
orbit of `O` under its action on subgroups.
If yes, return `(true, g)` where `g` is an element of `O` mapping the image of `i`
to the image of `j`.

Otherwise return `(false, nothing)`.

# Examples
```jldoctest
julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [2//3 0; 0 2//5]))
Finite quadratic module
  over integer ring
Abelian group: Z/15
Bilinear value module: Q/Z
Quadratic value module: Q/2Z
Gram matrix quadratic form:
[4//15]

julia> OT = orthogonal_group(T)
Orthogonal group of
  finite quadratic module: Z/15 -> Q/2Z
with 2 generators

julia> T3inT = primary_part(T, 3)[2]
Map
  from finite quadratic module: Z/3 -> Q/2Z
  to finite quadratic module: Z/15 -> Q/2Z

julia> T5inT = primary_part(T, 5)[2]
Map
  from finite quadratic module: Z/5 -> Q/2Z
  to finite quadratic module: Z/15 -> Q/2Z

julia> is_conjugate_with_data(OT, T3inT, T5inT)
(false, [1])
```
"""
function is_conjugate_with_data(O::AutomorphismGroup{TorQuadModule},
                                i::TorQuadModuleMap,
                                j::TorQuadModuleMap)
  @req domain(O) === codomain(i) === codomain(j) "Wrong parents"

  to_gap = get_attribute(O, :to_gap)
  Agap = codomain(to_gap)
  Hgap, _ = sub(Agap, elem_type(Agap)[to_gap(i(a)) for a in gens(domain(i))])
  G = GSetByElements(O, on_subgroups, typeof(Hgap)[Hgap])
  Kgap, _ = sub(Agap, elem_type(Agap)[to_gap(j(a)) for a in gens(domain(j))])
  return is_conjugate_with_data(G, Hgap, Kgap)
end

@doc raw"""
    stabilizer(O::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap)
                                           -> AutomorphismGroup{TorQuadModule},
                                              GAPGroupHomomorphism

Return the stabilizer of the image of `i` under the action of `O` on the subgroups
of the codomain of `i`.

# Examples
```jldoctest
julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [2//3 0; 0 2//5]));

julia> OT = orthogonal_group(T)
Orthogonal group of
  finite quadratic module: Z/15 -> Q/2Z
with 2 generators

julia> T3inT = primary_part(T, 3)[2]
Map
  from finite quadratic module: Z/3 -> Q/2Z
  to finite quadratic module: Z/15 -> Q/2Z

julia> S, _ = stabilizer(OT, T3inT)
(Orthogonal group of finite quadratic module: Z/15 -> Q/2Z, Hom: orthogonal group -> orthogonal group)

julia> order(S)
4
```
"""
function stabilizer(O::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap)
  to_gap = get_attribute(O, :to_gap)
  Agap = codomain(to_gap)
  Hgap, _ = sub(Agap, elem_type(Agap)[to_gap(i(a)) for a in gens(domain(i))])
  return stabilizer(O, Hgap.X, on_subgroups)
end
