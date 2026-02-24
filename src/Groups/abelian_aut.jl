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
  aut = AutomorphismGroup(AutGAP, G, true)
  set_attribute!(aut, :to_gap => to_gap, :to_oscar => to_oscar)
  return aut
end

function apply_automorphism(f::AutGrpAbTorElem, x::AbTorElem, check::Bool=true)
  aut = parent(f)
  @req !check || parent(x) == aut.G "Not in the domain of f!"
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  xgap = to_gap(x)
  A = parent(f)
  domGap = parent(xgap)
  imgap = typeof(xgap)(domGap, GAPWrap.Image(f.X, xgap.X))
  return to_oscar(imgap)::typeof(x)
end

(f::AutGrpAbTorElem)(x::SubPcGroupElem) = group_element(parent(x),GAPWrap.Image(GapObj(f), GapObj(x)))
Base.:^(x::SubPcGroupElem,f::AutGrpAbTorElem) = f(x)

(f::AutGrpAbTorElem)(x::AbTorElem) = apply_automorphism(f, x, true)
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
    hom(f::AutomorphismGroupElem{TorQuadModule}) -> TorQuadModuleMap

Return the underlying homomorphism of ``f``.

# Examples
```jldoctest
julia> A = abelian_group([2, 3, 4]);

julia> G = automorphism_group(A);

julia> f = first(gens(G))
Automorphism of
  finitely generated abelian group with 3 generators and 3 relations
with matrix representation
  [1   0   0]
  [0   1   0]
  [0   0   1]

julia> hom(f)
Map
  from finitely generated abelian group with 3 generators and 3 relations
  to finitely generated abelian group with 3 generators and 3 relations

julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [2//3 0; 0 2//9]));

julia> OT = orthogonal_group(T)
Group of isometries of
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

function (aut::AutGrpAbTor)(f::Union{FinGenAbGroupHom, TorQuadModuleMap}; check::Bool=true)
  @req !check || (domain(f) === codomain(f) === domain(aut) && is_bijective(f)) "Map does not define an automorphism of the abelian group"
  to_gap = get_attribute(aut, :to_gap)
  to_oscar = get_attribute(aut, :to_oscar)
  Agap = domain(to_oscar)
  AA = Agap.X
  function img_gap(x)
    a = to_oscar(group_element(Agap, x))
    b = to_gap(f(a))
    return b.X
  end
  gene = GAPWrap.GeneratorsOfGroup(AA)
  img = GAP.Obj([img_gap(a) for a in gene])
  fgap = GAPWrap.GroupHomomorphismByImagesNC(AA, AA, img)

  @req !check || fgap in aut.X "Map does not define an element of the group"
  return aut(fgap)
end

function (aut::AutomorphismGroup{TorQuadModule})(f::FinGenAbGroupHom; check::Bool=true)
  T = domain(aut)
  @req !check || (domain(f) === codomain(f) === abelian_group(T)) "Map does not define a morphism of the abelian group"
  fT = TorQuadModuleMap(T, T, f)
  return aut(fT; check)
end

function (aut::AutGrpAbTor)(M::ZZMatrix; check::Bool=true)
  @req !check || defines_automorphism(domain(aut), M) "Matrix does not define an automorphism of the abelian group"
  return aut(hom(domain(aut), domain(aut), M); check)
end

function (aut::AutGrpAbTor)(g::MatGroupElem{QQFieldElem, QQMatrix}; check::Bool=true)
  L = relations(domain(aut))
  if check
    B = basis_matrix(L)
    @assert can_solve(B, B*matrix(g); side=:left)
  end
  T = domain(aut)
  g = hom(T, T, elem_type(T)[T(lift(t)*matrix(g)) for t in gens(T)])
  return aut(g; check = false)
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
    _orthogonal_group(
      T::TorQuadModule,
      gensOT::Vector{U},
    ) where U <: Union{FinGenAbGroupHom, TorQuadModuleMap, ZZMatrix} -> AutomorphismGroup{TorQuadModule}

Return the subgroup of the orthogonal group of `G` generated by `gensOT`.
"""
function _orthogonal_group(
  T::TorQuadModule,
  gensOT::Vector{U};
  check::Bool=true,
) where U <: Union{FinGenAbGroupHom, TorQuadModuleMap, ZZMatrix}
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
  gens_aut = GapObj([ambient(g; check).X for g in gensOT])  # performs the checks
  if check
    # expensive for large groups
    subgrp_gap = GAP.Globals.Subgroup(ambient.X, gens_aut)
  else
    subgrp_gap = GAP.Globals.SubgroupNC(ambient.X, gens_aut)
  end
  aut = AutomorphismGroup(subgrp_gap, T)
  set_attribute!(aut, :to_gap => to_gap, :to_oscar => to_oscar)
  return aut
end

function Base.show(io::IO, ::MIME"text/plain", OT::AutomorphismGroup{TorQuadModule})
  io = pretty(io)
  println(io, "Group of isometries of", Indent())
  println(io, Lowercase(), OT.G)
  print(io, Dedent(), "with ", ItemQuantity(ngens(OT), "generator"))
end

function Base.show(io::IO, OT::AutomorphismGroup{TorQuadModule})
  if is_terse(io)
    print(io, "Group of isometries")
  else
    io = pretty(io)
    print(io, "Group of isometries of ", Lowercase(), OT.G)
  end
end 

@doc raw"""
    matrix(f::AutomorphismGroupElem{TorQuadModule}) -> ZZMatrix

Return a matrix inducing `f`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> T = torsion_quadratic_module(matrix(QQ, 2, 2, [1//12 0; 0 2//9]));

julia> OT = orthogonal_group(T)
Group of isometries of
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
Group of isometries of
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
    gensON_mat = unique(_compute_gens(N))
    gensON = TorQuadModuleMap[hom(N, N, g) for g in gensON_mat]
    gensOT = TorQuadModuleMap[i * g * j for g in gensON]
    length(gensOT) > 1 ? filter!(!isone∘matrix, gensOT) : nothing
    return _orthogonal_group(T, gensOT; check=false)
  elseif iszero(gram_matrix_quadratic(T))
    # in that case, we don't have any conditions regarding the
    # quadratic form, so we have all automorphisms coming
    # from the underlying abelian group
    return _orthogonal_group(T, hom.(gens(automorphism_group(abelian_group(T)))); check=false)
  else
    # if T is not semi-regular, we distinghuish the cases whether or not
    # it splits its radical quadratic
    i = radical_quadratic(T)[2]
    gensOT_mat = has_complement(i)[1] ? _compute_gens_split_degenerate(T) : _compute_gens_non_split_degenerate(T)
    return _orthogonal_group(T, gensOT_mat; check=false)
  end
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

function induce_endomorphism(f::TorQuadModuleMap, i::TorQuadModuleMap; check::Bool = true)
  @req !check || is_surjective(i) "i must be a surjection"
  @req domain(f) === codomain(f) === domain(i) "f must be an endomorphism of the domain of i"
  if check 
    K, j = kernel(i)
    @req all(has_preimage_with_preimage(j,f(j(k)))[1] for k in gens(K)) "f must preserve the kernel of i"
  end
  imgs = TorQuadModuleElem[]
  U = codomain(i)
  for a in gens(U)
    ok, c = has_preimage_with_preimage(i, a)
    @assert ok 
    push!(imgs, i(f(c)))
  end
  return hom(U, U, imgs)
end

function induce_automorphism(f::AutomorphismGroupElem{TorQuadModule}, i::TorQuadModuleMap; check::Bool = true)
  return induce_endomorphism(hom(f), i; check)
end 

function induce_automorphism_group(G::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap; check::Bool = true)
  @req domain(G) === domain(i) "G must consists of automorphisms of the domain of i"
  indu = TorQuadModuleMap[]
  for f in gens(G)
    g =  induce_automorphism(f, i; check)
    push!(indu, g)
  end
  H = _orthogonal_group(codomain(i), unique(indu); check)
  res = hom(G, H, gens(G), H.(indu); check)
  return H, res
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
  restr = TorQuadModuleMap[]
  for f in gens(G)
    g = try restrict_automorphism(f, i, check = false)
        catch e throw(ArgumentError("The domain of i is not invariant under the action of G"))
        end
    push!(restr, g)
  end
  H = _orthogonal_group(domain(i), restr; check)
  res = hom(G, H, gens(G), H.(restr); check)
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
Group of isometries of
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
  
function __cokernel(f::TorQuadModuleMap)
  # assumes same ambient space and therefore the underscore
  A = domain(f)
  B = codomain(f)
  BmodA = torsion_quadratic_module(cover(B),cover(A))
  return BmodA, hom(B, BmodA, [BmodA(lift(i)) for i in gens(B)])
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
Group of isometries of
  finite quadratic module: Z/15 -> Q/2Z
with 2 generators

julia> T3inT = primary_part(T, 3)[2]
Map
  from finite quadratic module: Z/3 -> Q/2Z
  to finite quadratic module: Z/15 -> Q/2Z

julia> S, _ = stabilizer(OT, T3inT)
(Group of isometries of finite quadratic module: Z/15 -> Q/2Z, Hom: group of isometries -> group of isometries)

julia> order(S)
4
```
"""
function stabilizer(O::AutomorphismGroup{TorQuadModule}, i::TorQuadModuleMap)
  @req domain(O)===codomain(i) "Domain of automorphism group must agree with codomain of inclusion." 
  A = domain(i)
  C = codomain(i)
  a = elementary_divisors(domain(i))
  if length(a)==0
    return O, id_hom(O) 
  end
  if order(O)<500 || (!is_prime(a[end]) && order(O)<5000)
    # don't do fancy stuff if the group is small
    to_gap = get_attribute(O, :to_gap)
    Agap = codomain(to_gap)
    Hgap, _ = sub(Agap, elem_type(Agap)[to_gap(i(a)) for a in gens(domain(i))])
    st2 =  stabilizer(O, Hgap.X, on_subgroups)
    return st2
  end
  if !is_snf(domain(O))
    # work with a Smith normal form
    # because otherwise GAP will choke on non-invertible matrices
    D = domain(O)
    Snf,toD = snf(D)
    Os,toOs = restrict_automorphism_group(O, toD)
    j = i*inv(toD)
    stab,inc = stabilizer(Os,j)
    return preimage(toOs,stab)
  end 
  if exponent(C) > exponent(A)
    expA = exponent(A)
    Ck, ck = kernel(hom(C,C, [expA*x for x in gens(C)]))
    Ak, ak = sub(Ck, [ck\i(x) for x in gens(A)])
    Ok,iOk = restrict_automorphism_group(O,ck; check=false)
    Sk, isk = stabilizer(Ok, ak)
    st = preimage(iOk, Sk)
    #@assert order(st[1])==order(st2[1])
    return st 
  end
  n = elementary_divisors(C)[end]
  # base case of the recursion n = p prime
  if is_prime(n)
    p = n
    B = matrix(i.map_ab)
    mats = GapObj([GapObj(matrix(x)) for x in gens(O)])
    st = _stab_via_fin_field(O, mats, B, n)
    #@assert order(st[1])==order(st2[1])
    return st
  end
  # for prime power order work with the F_p vector space 
  # (A + p^k*C) / (p^(k+1)C + pA)
  # where k = 0 ... v
  #
  # Suppose (by induction) that f ( A ) ⊆ A + p^k C and that moreover f stabilizes the image of
  # A in ( A + p^k C ) / ( p A + p^{k + 1} C ) , which is ( A + p^{k + 1} C ) / ( p A + p^{k + 1} C ) .
  # Then f ( a + p A + p^{k + 1} C ) = f ( a ) + p A + p^{k + 1} C ∈ ( A + p^{k + 1} C ) / ( p A + p^{k + 1} C ) , 
  # i.e. f ( a ) ∈ A + p^{k + 1} C .
  #
  # By the preceeding: if k is such that p^{k + 1}C = 0 , then f ( A ) ⊆ A.
  fl, v , p = is_prime_power_with_data(n)
  if fl
    A = domain(i)
    C = codomain(i)
    S = O
    pA = [i(p*x) for x in gens(A)]
    for k in 0:v
      # (A + p^k*C) / (p^(k+1)C + pA)
      piC = [p^k*x for x in gens(C)]
      D,iD = sub(C, append!(piC,i.(gens(A))))
      SD, iSD = restrict_automorphism_group(S, iD; check=false)
      E, iE =sub(D, iD.\append!([p^(k+1)*x for x in gens(C)],pA))
      K, iK = __cokernel(iE)
      SK, toSK = induce_automorphism_group(SD,iK; check=false)
      B,j = sub(K, [iK(iD\(i(x))) for x in gens(A)])
      S,_ = stabilizer(SK,j)
      S,_ = preimage(toSK, S)
      S,iS = preimage(iSD, S)
    end
    st = S,iS
    #@assert order(st[1])==order(st2[1])
    return st
  end
  # For composite order iterate over primary parts.
  S = O
  iS = id_hom(O)
  for p in prime_divisors(n)
    Ap, ip = primary_part(domain(i), p)
    Cp, jp = primary_part(codomain(i), p)
    ApinCp, Ap_to_Cp = sub(Cp, [jp\i(ip(x)) for x in gens(Ap)])
    Op,iOp = restrict_automorphism_group(S,jp;check=false)
    Sp, _ = stabilizer(Op, Ap_to_Cp)
    S,iS1 = preimage(iOp, Sp)
    iS = iS1*iS
  end
  st = S,iS
  #@assert order(st[1])==order(st2[1])
  return st 
end

# We adapt the strategy in stabilizer to compute a transporter.
function _is_conjugate_with_data(O::AutomorphismGroup{TorQuadModule}, i1::TorQuadModuleMap, i2::TorQuadModuleMap)
  @req domain(O)===codomain(i1) "Domain of automorphism group must agree with codomain of inclusion." 
  @req domain(O)===codomain(i2) "Domain of automorphism group must agree with codomain of inclusion."
  #i1: A1 -> C
  #i2: A2 -> C
  A1 = domain(i1)
  A2 = domain(i2)
  C = codomain(i1)
  eldiv_1 = elementary_divisors(A1)
  eldiv_1 == elementary_divisors(A2) || return false, one(O)
  length(eldiv_1) == 0 && return true, one(O)
  expA = exponent(A1)
  
  # don't do fancy stuff if the group is small
  if order(O)<80 || (eldiv_1[1]!=eldiv_1[end] && order(O)<130)
    @vprint :Isometry 4 "transporter via GAP group order $(order(O))"
    to_gap = get_attribute(O, :to_gap)
    Agap = codomain(to_gap)
    H1gap, _ = sub(Agap, elem_type(Agap)[to_gap(i1(a)) for a in gens(domain(i1))])
    H2gap, _ = sub(Agap, elem_type(Agap)[to_gap(i2(a)) for a in gens(domain(i2))])
    X = orbit(O, on_subgroups, H1gap)
    fl, g = is_conjugate_with_data(X, H1gap, H2gap)
    @vprintln :Isometry 4 " done"
    return fl, g
  end
  
  if exponent(C) > expA
    # We can always reduce to an action on Ck < C (see next line)
    Ck, ck = kernel(hom(C,C, [expA*x for x in gens(C)]))
    A1k, a1k = sub(Ck, [ck\i1(x) for x in gens(A1)])
    A2k, a2k = sub(Ck, [ck\i2(x) for x in gens(A2)])
    Ok, iOk = restrict_automorphism_group(O, ck; check=false)
    flag, _transporter = _is_conjugate_with_data(Ok, a1k, a2k)
    !flag && return flag, one(O)
    transporter = iOk\_transporter
    return flag, transporter
  end
  
  BL1 = matrix(abelian_group_homomorphism(i1))
  BL2 = matrix(abelian_group_homomorphism(i2))
  # The vector space case is the base case of the recursion.
  if is_prime(expA)
    p = expA
    @vprint :Isometry 4 "transporter elementary $p "
    # prime
    mats = matrix.(gens(O))
    if length(mats)==0
      # catch a corner case which may make gap choke
      return A1==A2, one(O)
    end
    R = fpField(UInt(p))
    BL1mod = change_base_ring(R, BL1)
    BL2mod = change_base_ring(R, BL2)
    r = rref!(BL1mod)
    rref!(BL2mod)
    if nrows(BL1mod) > r
      BL1mod = BL1mod[1:r,:]
    end 
    if nrows(BL2mod) > r
      BL2mod = BL2mod[1:r,:]
    end 

    G = O
    Gnice = GAP.Globals.NiceObject(GapObj(G))
    # FIXME: direct conversion from fpMatrix to GAP matrix seems to be missing?
    BL1mod_gap = GapObj(lift(BL1mod)) * GAP.Globals.Z(GapObj(p))^0
    BL2mod_gap = GapObj(lift(BL2mod)) * GAP.Globals.Z(GapObj(p))^0
    GAP.Globals.ConvertToMatrixRep(BL1mod_gap)
    GAP.Globals.ConvertToMatrixRep(BL2mod_gap)
    mats_gap = GapObj([GapObj(x) for x in mats])
    img = GAP.Globals.RepresentativeAction(Gnice, BL1mod_gap, BL2mod_gap, GAP.Globals.GeneratorsOfGroup(Gnice), mats_gap, GAP.Globals.OnSubspacesByCanonicalBasis)
    img == GAP.Globals.fail && error("") && return false, one(G)
    mono = GAP.Globals.NiceMonomorphism(GapObj(G))
    _g = G(GAP.Globals.PreImage(mono, img))
    @vprintln :Isometry 4 " done"
    if get_assertion_level(:Isometry) > 0
      to_gap = get_attribute(G, :to_gap)
      Agap = codomain(to_gap)
      H1gap, _ = sub(Agap, elem_type(Agap)[to_gap(i1(a)) for a in gens(domain(i1))])
      H2gap, _ = sub(Agap, elem_type(Agap)[to_gap(i2(a)) for a in gens(domain(i2))])
      @assert on_subgroups(H1gap,_g)==H2gap
    end
    return true, _g
  end
  n = expA
  R, iR = residue_ring(ZZ, Int(n); cached=false)
  fl, v, p = is_prime_power_with_data(expA)
  # we can work with the Howell form
  # if D is a free R module
  eldivC = elementary_divisors(C)
  if eldivC[1] == eldivC[end]
    @vprint :Isometry 4 "transporter via Howell"
    BL1mod = change_base_ring(R, BL1)
    BL2mod = change_base_ring(R, BL2)
    k1 = ncols(BL1) - nrows(BL1)
    if k1 > 0
      # howell form requires a square matrix
      BL1mod = vcat(BL1mod, zero_matrix(R, k1, ncols(BL1)))
    end
    k2 = ncols(BL2) - nrows(BL2)
    if k2 > 0
      # howell form requires a square matrix
      BL2mod = vcat(BL2mod, zero_matrix(R, k2, ncols(BL2)))
    end
    howell_form!(BL1mod)
    howell_form!(BL2mod)
    X = orbit(O, on_howell_form, BL1mod)
    flag, _g = is_conjugate_with_data(X, BL1mod, BL2mod)
    @vprintln :Isometry 4 " done"
  elseif fl
    @vprintln :Isometry 4 "transporter via prime power"
    S = O
    iS = id_hom(S)
    transporter = one(O)
    pA1 = [i1(p*x) for x in gens(A1)]
    _i2 = i2
    # We move A1 to A2 step by step going through the quotients 
    # (A1 + p^k C) / (p^(k+1)C + pA1)
    # and compute their stabilizers
    # Note that at each step _A2 := _i2(A2)
    # is contained in A1 + p^k C
    # The basic idea is the same as in 
    # stabilizer(::AutomorphismGroup,::TorQuadModMor) above
    for k in 0:v
      piC = [p^k*x for x in gens(C)]
      D, iD = sub(C, append!(piC, i1.(gens(A1))))
      SD, iSD = restrict_automorphism_group(S, iD; check=false)
      E, iE =sub(D, iD.\append!([p^(k+1)*x for x in gens(C)],pA1))
      K, iK = __cokernel(iE)
      SK, toSK = induce_automorphism_group(SD, iK; check=false)
      B1,j1 = sub(K, [iK(iD\(i1(x))) for x in gens(A1)])
      B2,j2 = sub(K, [iK(iD\(_i2(x))) for x in gens(A2)])
      flag, _transporter = _is_conjugate_with_data(SK, j1, j2)
      !flag && return false, one(O)
      _transporter = iS(iSD\(toSK\_transporter))
      @assert parent(_transporter) == O
      transporter = _transporter * transporter
      _i2 = _i2*hom(inv(_transporter))
      toperm = isomorphism(PermGroup, domain(toSK))
      if k < v
        # not needed in the last iteration
        S,_ = stabilizer(SK, j1)
        @vtime preimage(toSK*inv(toperm), S)
        @vtime :Isometry 4 S,_ = preimage(toSK, S)
        @vtime :Isometry 4 S,iS1 = preimage(iSD, S)
        iS = iS1*iS
        @assert codomain(iS)==O
      end
    end
    _g = transporter
  else 
    @vprintln :Isometry 4 "transporter via composite"
    # for composite order map A1 to A2 prime by prime
    # compute stabilizers to preserve the progress
    _O = O
    _g = one(O)
    _phi = id_hom(O)
    _i2 = i2
    for p in prime_divisors(expA)
      R, iR = residue_ring(ZZ, Int(n); cached=false)
      A1p, i1p = primary_part(domain(i1), p)
      A2p, i2p = primary_part(domain(_i2), p)
      Cp, jp = primary_part(C, p)
      @assert is_snf(Cp)
      A1pinCp, A1p_to_Cp = sub(Cp, [jp\i1(i1p(x)) for x in gens(A1p)])
      A2pinCp, A2p_to_Cp = sub(Cp, [jp\_i2(i2p(x)) for x in gens(A2p)])
      
      Op, iOp = restrict_automorphism_group(_O, jp; check=false)
      @vtime :Isometry 4 flag, transporter = _is_conjugate_with_data(Op, A1p_to_Cp, A2p_to_Cp)
      !flag && return false, one(O)
      @vtime :Isometry 4 _transporter = _phi(preimage(iOp, transporter))
      _g = _transporter*_g
      stab_p, istab_p = stabilizer(Op, A1p_to_Cp)
      # permutation group seems to be a bit faster 
      # ... but still slowish!
      to_Operm = isomorphism(PermGroup, _O)
      to_Op_perm = isomorphism(PermGroup, Op)
      @vtime :Isometry 4 _tmp = preimage(inv(to_Operm)*iOp*to_Op_perm, to_Op_perm(stab_p)[1])[1]
      @vtime :Isometry 4 _O,i_O = preimage(to_Operm,_tmp)
      _phi = i_O*_phi
      _i2 = _i2*hom(inv(_transporter))
      @assert codomain(_phi) === O
      @vprintln :Isometry 4 "done $p"
    end
  end
  if get_assertion_level(:Isometry) > 0
    # confirm a positive result
    to_gap = get_attribute(O, :to_gap)
    Agap = codomain(to_gap)
    H1gap, _ = sub(Agap, elem_type(Agap)[to_gap(i1(a)) for a in gens(domain(i1))])
    H2gap, _ = sub(Agap, elem_type(Agap)[to_gap(i2(a)) for a in gens(domain(i2))])
    if flag
      @assert on_subgroups(H1gap,_g)==H2gap    
    end
  end
  return flag, _g
end
