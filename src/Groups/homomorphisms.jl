export
    automorphism_group,
    codomain,
    cokernel,
    domain,
    haspreimage,
    hom,
    id_hom,
    image,
    induced_automorphism,
    inner_automorphism,
    inner_automorphisms_group,
    isbijective,
    isinjective,
    isinner_automorphism,
    isinvariant,
    isinvertible,
    isisomorphic,
    isisomorphic_with_map,
    isomorphic_fp_group,
    isomorphic_pc_group,
    isomorphic_perm_group,
    isomorphism,
    issurjective,
    kernel,
    order,
    restrict_automorphism,
    restrict_homomorphism,
    simplified_fp_group,
    trivial_morphism

function Base.show(io::IO, x::GAPGroupHomomorphism)
  print(io, "Group homomorphism from \n")
  show(IOContext(io, :compact => true), domain(x))
  print(io, "\nto\n")
  show(IOContext(io, :compact => true), codomain(x))
end

function ==(f::GAPGroupHomomorphism{S,T}, g::GAPGroupHomomorphism{S,T}) where S where T
   return f.map == g.map
end

Base.:*(f::GAPGroupHomomorphism{S, T}, g::GAPGroupHomomorphism{T, U}) where S where T where U = compose(f, g)

function Base.inv(f::GAPGroupHomomorphism{S,T}) where S where T
   @assert GAPWrap.IsBijective(f.map) "f is not bijective"
   return GAPGroupHomomorphism(codomain(f), domain(f), GAP.Globals.InverseGeneralMapping(f.map))
end

order(f::GAPGroupHomomorphism) = GAPWrap.Order(f.map)

function Base.:^(f::GAPGroupHomomorphism{S,T}, n::Int64) where S where T
   if n==1
     return f
   else
      @assert domain(f) == codomain(f) "Domain and codomain do not coincide"
      return GAPGroupHomomorphism(domain(f),codomain(f),(f.map)^n)
   end
end

function compose(f::GAPGroupHomomorphism{S, T}, g::GAPGroupHomomorphism{T, U}) where {S, T, U}
  dom = domain(f)
  cod = codomain(g)
  @assert codomain(f) == domain(g)
  mp = GAP.Globals.CompositionMapping(g.map, f.map)
  return GAPGroupHomomorphism(dom, cod, mp)
end

#(g::GAPGroupHomomorphism{T, U})(f::GAPGroupHomomorphism{S, T}) where {S,T,U} = compose(f,g)

"""
    id_hom(G::GAPGroup)

Return the identity homomorphism on the group `G`.
"""
function id_hom(G::GAPGroup)
  return GAPGroupHomomorphism(G, G, GAP.Globals.IdentityMapping(G.X))
end

"""
    trivial_morphism(G::GAPGroup, H::GAPGroup = G)

Return the homomorphism from `G` to `H` sending every element of `G` into the
identity of `H`.
"""
function trivial_morphism(G::GAPGroup, H::GAPGroup = G)
  return hom(G, H, x -> one(H))
end

"""
    hom(G::GAPGroup, H::GAPGroup, f::Function)

Return the group homomorphism defined by the function `f`.
"""
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
  return GAPGroupHomomorphism(G, H, mp)
end

# This method is used for those embeddings that are
# the identity on the GAP side.
function hom(G::GAPGroup, H::GAPGroup, img::Function, preimg::Function; is_known_to_be_bijective::Bool = false)
  function gap_fun(x::GapObj)
    el = group_element(G, x)
    img_el = img(el)
    return img_el.X
  end

  function gap_pre_fun(x::GapObj)
    el = group_element(H, x)
    preimg_el = preimg(el)
    return preimg_el.X
  end

  if is_known_to_be_bijective
    mp = GAP.Globals.GroupHomomorphismByFunction(G.X, H.X, GAP.julia_to_gap(gap_fun), GAP.julia_to_gap(gap_pre_fun))
  else
    mp = GAP.Globals.GroupHomomorphismByFunction(G.X, H.X, GAP.julia_to_gap(gap_fun), false, GAP.julia_to_gap(gap_pre_fun))
  end

  return GAPGroupHomomorphism(G, H, mp)
end

"""
    hom(G::GAPGroup, H::GAPGroup, gensG::Vector = gens(G), imgs::Vector; check::Bool = true)

Return the group homomorphism defined by `gensG`[`i`] -> `imgs`[`i`] for every
`i`. In order to work, the elements of `gensG` must generate `G`.

If `check` is set to `false` then it is not checked whether the mapping
defines a group homomorphism.
"""
function hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector; check::Bool = true)
  vgens = GapObj([x.X for x in gensG])
  vimgs = GapObj([x.X for x in imgs])
  if check
    mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  else
    mp = GAP.Globals.GroupHomomorphismByImagesNC(G.X, H.X, vgens, vimgs)
  end
  if mp == GAP.Globals.fail throw(ArgumentError("Invalid input")) end
  return GAPGroupHomomorphism(G, H, mp)
end

function hom(G::GAPGroup, H::GAPGroup, imgs::Vector; check::Bool = true)
  return hom(G, H, gens(G), imgs; check)
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
  return group_element(codomain(f), GAPWrap.Image(f.map,x.X))
end

"""
    preimage(f::GAPGroupHomomorphism, x::GAPGroupElem)

Return an element `y` in the domain of `f` with the property `f(y) == x`.
See [`haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem; check::Bool = true)`](@ref)
for a check whether `x` has such a preimage.
"""
function preimage(f::GAPGroupHomomorphism, x::GAPGroupElem)
  fl, p = haspreimage(f, x)
  @assert fl
  return p
end

"""
    issurjective(f::GAPGroupHomomorphism)

Return whether `f` is surjective.
"""
function issurjective(f::GAPGroupHomomorphism)
  return GAPWrap.IsSurjective(f.map)
end

"""
    isinjective(f::GAPGroupHomomorphism)

Return whether `f` is injective.
"""
function isinjective(f::GAPGroupHomomorphism)
  return GAPWrap.IsInjective(f.map)
end

"""
    isinvertible(f::GAPGroupHomomorphism)

Return whether `f` is invertible.
"""
function isinvertible(f::GAPGroupHomomorphism)
  return GAPWrap.IsBijective(f.map)
end

"""
    isbijective(f::GAPGroupHomomorphism)

Return whether `f` is bijective.
"""
function isbijective(f::GAPGroupHomomorphism)
  return GAPWrap.IsBijective(f.map)
end


"""
    isinvariant(f::GAPGroupHomomorphism, H::Group)
    isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return whether `f(H) == H` holds.
An exception is thrown if `domain(f)` and `codomain(f)` are not equal
or if `H` is not contained in `domain(f)`.
"""
function isinvariant(f::GAPGroupHomomorphism, H::GAPGroup)
  @assert domain(f) == codomain(f) "Not an endomorphism!"
  @assert GAPWrap.IsSubset(domain(f).X, H.X) "Not a subgroup of the domain"
  return GAPWrap.Image(f.map, H.X) == H.X
end

"""
    restrict_homomorphism(f::GAPGroupHomomorphism, H::Group)
    restrict_homomorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: Group

Return the restriction of `f` to `H`.
An exception is thrown if `H` is not a subgroup of `domain(f)`.
"""
function restrict_homomorphism(f::GAPGroupHomomorphism, H::GAPGroup)
  # We have to check whether `H` is really a subgroup of `f.domain`,
  # since `GAP.Globals.RestrictedMapping` does not check this.
  # (The GAP documentation does not claim anything about the result
  # in the case that `H` is not a subgroup of `f.domain`,
  # and in fact just the given map may be returned.)
  @assert issubgroup(domain(f), H)[1] "Not a subgroup!"
  return GAPGroupHomomorphism(H, f.codomain, GAP.Globals.RestrictedMapping(f.map,H.X)::GapObj)
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
  K = GAP.Globals.Kernel(f.map)::GapObj
  return _as_subgroup(domain(f), K)
end

"""
    image(f::GAPGroupHomomorphism)

Return the image of `f` as subgroup of `codomain`(`f`), together with the embedding homomorphism.
"""
function image(f::GAPGroupHomomorphism)
  K = GAPWrap.Image(f.map)
  return _as_subgroup(codomain(f), K)
end

"""
    image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: GAPGroup where T <: GAPGroup
    (f::GAPGroupHomomorphism{S, T})(H::S)

Return `f`(`H`), together with the embedding homomorphism into `codomain`(`f`).
"""
function image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: GAPGroup where T <: GAPGroup
  H1 = GAPWrap.Image(f.map, H.X)
  return _as_subgroup(codomain(f), H1)
end

(f::GAPGroupHomomorphism{S, T})(H::S) where S <: GAPGroup where T <: GAPGroup = image(f,H)

"""
    cokernel(f::GAPGroupHomomorphism)

Return the cokernel of `f`, that is, the quotient of the codomain of `f`
by the normal closure of the image.
"""
function cokernel(f::GAPGroupHomomorphism)
  K, mK = image(f)
  C = codomain(f)
  return quo(C, normal_closure(C, K)[1])
end

"""
    haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem; check::Bool = true)

Return (`true`, `y`) if there exists `y` in `domain(f)`
such that `f`(`y`) = `x` holds;
otherwise, return (`false`, `o`) where `o` is the identity of `domain(f)`.

If `check` is set to `false` then the test whether `x` is an element of
`image(f)` is omitted.
"""
function haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem; check::Bool = true)
  # `GAP.Globals.PreImagesRepresentative` does not promise anything
  # if the given element is not in the codomain of the map.
# check && ! (x in codomain(f)) && return false, one(domain(f))
#TODO:
# Apparently the documentation of `GAP.Globals.PreImagesRepresentative`
# is wrong in the situation that `x` is not in the *image* of `f`,
# the function can then run into an error or return some group element,
# see https://github.com/gap-system/gap/issues/4088.
# Until this problem gets fixed on the GAP side, we perform a membership test
# before calling `GAP.Globals.PreImagesRepresentative`.
  check && ! (x in image(f)[1]) && return false, one(domain(f))
  r = GAP.Globals.PreImagesRepresentative(f.map, x.X)::GapObj
  if r == GAP.Globals.fail
    return false, one(domain(f))
  else
    return true, group_element(domain(f), r)
  end
end

"""
    preimage(f::GAPGroupHomomorphism{S, T}, H::T) where S <: GAPGroup where T <: GAPGroup

If `H` is a subgroup of the codomain of `f`, return the subgroup `f^-1(H)`,
together with its embedding homomorphism into the domain of `f`.
"""
function preimage(f::GAPGroupHomomorphism{S, T}, H::T) where S <: GAPGroup where T <: GAPGroup
  H1 = GAP.Globals.PreImage(f.map, H.X)::GapObj
  return _as_subgroup(domain(f), H1)
end


################################################################################
#
#  IsIsomorphic
#
################################################################################

"""
    isisomorphic_with_map(G::Group, H::Group)

Return (`true`,`f`) if `G` and `H` are isomorphic groups, where `f` is a group
isomorphism. Otherwise, return (`false`,`f`), where `f` is the trivial
homomorphism.

# Examples
```jldoctest
julia> isisomorphic_with_map(symmetric_group(3), dihedral_group(6))
(true, Group homomorphism from
Sym( [ 1 .. 3 ] )
to
<pc group of size 6 with 2 generators>)
```
"""
function isisomorphic_with_map(G::GAPGroup, H::GAPGroup)
  mp = GAP.Globals.IsomorphismGroups(G.X, H.X)::GapObj
  if mp === GAP.Globals.fail
    return false, trivial_morphism(G, H)
  else
    return true, GAPGroupHomomorphism(G, H, mp)
  end
end

"""
    isisomorphic(G::Group, H::Group)

Return `true` if `G` and `H` are isomorphic groups, and `false` otherwise.

# Examples
```jldoctest
julia> isisomorphic(symmetric_group(3), dihedral_group(6))
true
```
"""
function isisomorphic(G::GAPGroup, H::GAPGroup)
  mp = GAP.Globals.IsomorphismGroups(G.X, H.X)::GapObj
  return mp !== GAP.Globals.fail
end

"""
    isomorphism(G::Group, H::Group)

Return a group isomorphism between `G` and `H` if they are isomorphic groups.
Otherwise throw an exception.

# Examples
```jldoctest
julia> isomorphism(symmetric_group(3), dihedral_group(6))
Group homomorphism from
Sym( [ 1 .. 3 ] )
to
<pc group of size 6 with 2 generators>
```
"""
function isomorphism(G::GAPGroup, H::GAPGroup)
  mp = GAP.Globals.IsomorphismGroups(G.X, H.X)::GapObj
  mp === GAP.Globals.fail && throw(ArgumentError("the groups are not isomorphic"))
  return GAPGroupHomomorphism(G, H, mp)
end


################################################################################
#
#  Conversions between types
#
################################################################################

_get_iso_function(::Type{PermGroup}) = GAP.Globals.IsomorphismPermGroup
_get_iso_function(::Type{FPGroup}) = GAP.Globals.IsomorphismFpGroup
_get_iso_function(::Type{PcGroup}) = GAP.Globals.IsomorphismPcGroup


"""
    isomorphism(::Type{T}, G::GAPGroup) where T <: Union{FPGroup, PcGroup, PermGroup}

Return an isomorphism from `G` to a group of type `T`.
An exception is thrown if no such isomorphism exists.

Isomorphisms are cached in `G`, subsequent calls of `isomorphism` with the
same `T` yield identical results.

If only the image of such an isomorphism is needed, use `T(G)`.

# Examples
```jldoctest
julia> G = dihedral_group(6)
<pc group of size 6 with 2 generators>

julia> iso = isomorphism(PermGroup, G)
Group homomorphism from
<pc group of size 6 with 2 generators>
to
Group([ (1,2)(3,6)(4,5), (1,3,5)(2,4,6) ])

julia> PermGroup(G)
Group([ (1,2)(3,6)(4,5), (1,3,5)(2,4,6) ])

julia> codomain(iso) === ans
true
```
"""
function isomorphism(::Type{T}, G::GAPGroup) where T <: Union{FPGroup, PcGroup, PermGroup}
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, G, :isomorphisms)
   return get!(isos, T) do
     fun = _get_iso_function(T)
     f = fun(G.X)::GapObj
     f == GAP.Globals.fail && throw(ArgumentError("Could not convert group into a group of type $T"))
     H = T(GAP.Globals.ImagesSource(f)::GapObj)
     return GAPGroupHomomorphism(G, H, f)
   end::GAPGroupHomomorphism{typeof(G), T}
end


"""
    isomorphism(::Type{GrpAbFinGen}, G::GAPGroup)

Return a map from `G` to an isomorphic (additive) group of type `GrpAbFinGen`.
An exception is thrown if `G` is not abelian or not finite.
"""
function isomorphism(::Type{GrpAbFinGen}, G::GAPGroup)
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, G, :isomorphisms)
   return get!(isos, GrpAbFinGen) do
     isabelian(G) || throw(ArgumentError("the group is not abelian"))
     isfinite(G) || throw(ArgumentError("the group is not finite"))
#T this restriction is not nice

     indep = GAP.Globals.IndependentGeneratorsOfAbelianGroup(G.X)::GapObj
     orders = [GAPWrap.Order(x) for x in indep]
     n = length(indep)
     A = abelian_group(GrpAbFinGen, orders)

     f(g) = A(Vector{fmpz}(GAPWrap.IndependentGeneratorExponents(G.X, g.X)))

     finv = function(g::elem_type(GrpAbFinGen))
       res = GAPWrap.One(G.X)
       for i in 1:n
         res = res * indep[i]^GAP.Obj(g.coeff[i])
       end
       return group_element(G, res)
     end

     return MapFromFunc(f, finv, G, A)
   end::MapFromFunc{typeof(G), GrpAbFinGen}
end

"""
    isomorphism(::Type{T}, A::GrpAbFinGen) where T <: Union{GAPGroup, GrpAbFinGen}

Return an isomorphism from `A` to a group of type `T`.
An exception is thrown if no such isomorphism exists or if `A` is not finite.
"""
function isomorphism(::Type{T}, A::GrpAbFinGen) where T <: GAPGroup
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, A, :isomorphisms)
   return get!(isos, T) do
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
     G = abelian_group(T, exponents)
     # `GAP.Globals.GeneratorsOfGroup(G.X)` consists of independent elements
     # of the orders in `exponents`.
     # `GAP.Globals.IndependentGenerators(G.X)` chooses generators
     # that may differ from these generators,
     # and that belong to the exponent vectors returned by
     # `GAPWrap.IndependentGeneratorExponents(G.X, g)`.
     # `GAP.Globals.GeneratorsOfGroup(G.X)` corresponds to `gens(A2)`,
     # we let `hom` compute elements in `A2` that correspond to
     # `GAP.Globals.IndependentGenerators(G.X)`.
     Ggens = Vector{GapObj}(GAP.Globals.GeneratorsOfGroup(G.X)::GapObj)
     gensindep = GAP.Globals.IndependentGeneratorsOfAbelianGroup(G.X)::GapObj
     Aindep = abelian_group(fmpz[GAP.Globals.Order(g) for g in gensindep])

     imgs = [Vector{fmpz}(GAPWrap.IndependentGeneratorExponents(G.X, a)) for a in Ggens]
     A2_to_Aindep = hom(A2, Aindep, elem_type(Aindep)[Aindep(e) for e in imgs])
     Aindep_to_A = compose(inv(A2_to_Aindep), A2_to_A)
     n = length(exponents)

     f = function(a::elem_type(GrpAbFinGen))
       exp = A_to_A2(a)
       img = GAP.Globals.One(G.X)
       for i in 1:n
         img = img * Ggens[i]^GAP.Obj(exp[i])
       end
       return group_element(G, img)
     end

     finv = function(g)
       exp = Vector{fmpz}(GAPWrap.IndependentGeneratorExponents(G.X, g.X))
       return Aindep_to_A(Aindep(exp))
     end

     return MapFromFunc(f, finv, A, G)
   end::MapFromFunc{GrpAbFinGen, T}
end

function isomorphism(::Type{GrpAbFinGen}, A::GrpAbFinGen)
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, A, :isomorphisms)
   return get!(isos, GrpAbFinGen) do
     return identity_map(A)
   end::AbstractAlgebra.Generic.IdentityMap{GrpAbFinGen}
end

"""
    FPGroup(G::T) where T <: Union{GAPGroup, GrpAbFinGen}
    GrpAbFinGen(G::T) where T <: GAPGroup
    PcGroup(G::T) where T <: Union{GAPGroup, GrpAbFinGen}
    PermGroup(G::T) where T <: Union{GAPGroup, GrpAbFinGen}

Return a group of type `T` that is isomorphic with `G`.
If one needs the isomorphism then
[isomorphism(::Type{T}, G::GAPGroup) where T <: Union{FPGroup, PcGroup, PermGroup}](@ref)
can be used instead.
"""
function (::Type{S})(G::T) where {S <: Union{GrpAbFinGen, GAPGroup}, T <: GAPGroup}
   return codomain(isomorphism(S, G))
end

function (::Type{T})(G::GrpAbFinGen) where T <: GAPGroup
   return codomain(isomorphism(T, G))
end


################################################################################
#
# provide and deprecate the old syntax
#
function _isomorphic_perm_group(G::GAPGroup)
   f = isomorphism(PermGroup, G)
   return codomain(f), f
end

@deprecate isomorphic_perm_group(G::GAPGroup) _isomorphic_perm_group(G)

function _isomorphic_pc_group(G::GAPGroup)
   f = isomorphism(PcGroup, G)
   return codomain(f), f
end

@deprecate isomorphic_pc_group(G::GAPGroup) _isomorphic_pc_group(G)

function _isomorphic_fp_group(G::GAPGroup)
   f = isomorphism(FPGroup, G)
   return codomain(f), f
end

@deprecate isomorphic_fp_group(G::GAPGroup) _isomorphic_fp_group(G)

function _isomorphic_group(::Type{T}, G::GAPGroup) where T <: GAPGroup
  fmap = isomorphism(T, G)
  return codomain(fmap), fmap
end

@deprecate isomorphic_group(::Type{T}, G::GAPGroup) where T <: GAPGroup _isomorphic_group(T, G)


"""
    simplified_fp_group(G::FPGroup)

Return a group `H` of type `FPGroup` and an isomorphism `f` from `G` to `H`, where
the presentation of `H` was obtained from the presentation of `G` by applying Tietze
transformations in order to reduce it with respect to the number of
generators, the number of relators, and the relator lengths.

# Examples
```jldoctest
julia> F = free_group(3)
<free group on the generators [ f1, f2, f3 ]>

julia> G = quo(F, [gen(F,1)])[1]
<fp group of size infinity on the generators [ f1, f2, f3 ]>

julia> simplified_fp_group(G)[1]
<fp group of size infinity on the generators [ f2, f3 ]>
```
"""
function simplified_fp_group(G::FPGroup)
   f = GAP.Globals.IsomorphismSimplifiedFpGroup(G.X)
   H = FPGroup(GAPWrap.Image(f))
   # TODO: remove the next line once https://github.com/gap-system/gap/pull/4810
   # is deployed to Oscar
   GAP.Globals.UseIsomorphismRelation(G.X, H.X)
   return H, GAPGroupHomomorphism(G,H,f)
end


################################################################################
#
#  Automorphism Group
#
################################################################################

"""
    automorphism_group(G::Group) -> A::AutomorphismGroup{T}

Return the full automorphism group of `G`. If `f` is an object of type
`GAPGroupHomomorphism` and it is bijective from `G` to itself, then `A(f)`
return the embedding of `f` in `A`.

Elements of `A` can be multiplied with other elements of `A` or by elements
of type `GAPGroupHomomorphism`; in this last case, the result has type
`GAPGroupHomomorphism`.
"""
function automorphism_group(G::GAPGroup)
  AutGAP = GAP.Globals.AutomorphismGroup(G.X)::GapObj
  return AutomorphismGroup(AutGAP, G)
end

function Base.show(io::IO, A::AutomorphismGroup{T}) where T <: GAPGroup
  print(io, "Aut( "* String(GAP.Globals.StringView(A.G.X)) * " )")
end

"""
  domain(A::AutomorphismGroup) -> Group

Return the domain of this group of automorphisms.
"""
domain(A::AutomorphismGroup) = A.G

"""
    domain(f::AutomorphismGroupElem) -> Group

Return the domain of this automorphism.
"""
domain(f::AutomorphismGroupElem) = domain(parent(f))

"""
    hom(f::GAPGroupElem{AutomorphismGroup{T}}) where T

Return the element f of type `GAPGroupHomomorphism{T,T}`.
"""
function hom(x::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
  A = parent(x)
  G = A.G
  return GAPGroupHomomorphism(G, G, x.X)
end

(f::GAPGroupElem{AutomorphismGroup{T}})(x::GAPGroupElem) where T <: GAPGroup = apply_automorphism(f, x, true)
Base.:^(x::GAPGroupElem{T},f::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = apply_automorphism(f, x, true)
#Base.:^(f::GAPGroupElem{AutomorphismGroup{T}},g::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = g^-1*f*g

function (A::AutomorphismGroup{T})(f::GAPGroupHomomorphism{T,T}) where T <: GAPGroup
   @assert domain(f)==A.G && codomain(f)==A.G "f not in A"
   @assert isbijective(f) "f not in A"
   return group_element(A, f.map)
end

function apply_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, x::GAPGroupElem, check::Bool=true) where T <: GAPGroup
  A = parent(f)
  G = parent(x)
  if check
    @assert A.G == G || x.X in A.G.X "Not in the domain of f!"      #TODO Do we really need the IN check?
  end
  return typeof(x)(G, GAPWrap.Image(f.X,x.X))
end

Base.:*(f::GAPGroupElem{AutomorphismGroup{T}}, g::GAPGroupHomomorphism) where T = hom(f)*g
Base.:*(f::GAPGroupHomomorphism, g::GAPGroupElem{AutomorphismGroup{T}}) where T = f*hom(g)

"""
    inner_automorphism(g::GAPGroupElem)

Return the inner automorphism in `automorphism_group(parent(g))` defined by `x` -> `x^g`.
"""
function inner_automorphism(g::GAPGroupElem)
  return GAPGroupHomomorphism(parent(g), parent(g), GAP.Globals.ConjugatorAutomorphism(parent(g).X, g.X))
end

"""
    isinner_automorphism(f::GAPGroupHomomorphism)
    isinner_automorphism(f::GAPGroupElem{AutomorphismGroup{T}})

Return whether `f` is an inner automorphism.
"""
function isinner_automorphism(f::GAPGroupHomomorphism)
  @assert domain(f) == codomain(f) "Not an automorphism!"
  return GAPWrap.IsInnerAutomorphism(f.map)
end

function isinner_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
  return GAPWrap.IsInnerAutomorphism(f.X)
end

"""
    inner_automorphisms_group(A::AutomorphismGroup{T})

Return the subgroup of `A` of the inner automorphisms.
"""
function inner_automorphisms_group(A::AutomorphismGroup{T}) where T <: GAPGroup
   AutGAP = GAP.Globals.InnerAutomorphismsAutomorphismGroup(A.X)
   return _as_subgroup(A, AutGAP)
end

"""
    isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return whether `f`(`H`) == `H`.
"""
function isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T<:GAPGroup
  @assert GAPWrap.IsSubset(parent(f).G.X, H.X) "Not a subgroup of the domain"
  return GAPWrap.Image(f.X, H.X) == H.X
end

"""
    induced_automorphism(f::GAPGroupHomomorphism, g::GAPGroupHomomorphism)
    induced_automorphism(f::GAPGroupHomomorphism, g::GAPGroupElem{AutomorphismGroup{T}})

Return the automorphism `h` of the image of `f` such that `h`(`f`) ==
`f`(`g`), where `g` is an automorphism of a group `G` and `f` is a group
homomorphism defined over `G` such that the kernel of `f` is invariant under
`g`
"""
function induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupHomomorphism)
  @assert isinvariant(mH, kernel(f)[1]) "The kernel is not invariant under g!"
  map = GAP.Globals.InducedAutomorphism(f.map, mH.map)
  A = automorphism_group(image(f)[1])
  return A(GAPGroupHomomorphism(codomain(f), codomain(f), map))
end

induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = induced_automorphism(f,hom(mH))

"""
    restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return the restriction of `f` to `H` as an automorphism of `H`.
An exception is thrown if `H` is not invariant under `f`.
"""
function restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T, A=automorphism_group(H)) where T <: GAPGroup
  @assert isinvariant(f,H) "H is not invariant under f!"
  fh = hom(H, H, gens(H), [f(x) for x in gens(H)], check = false)
  return A(fh)
end

function restrict_homomorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: GAPGroup
  return restrict_homomorphism(hom(f),H)
end

function _as_subgroup_bare(G::AutomorphismGroup{T}, H::GapObj) where T
  return AutomorphismGroup{T}(H, G.G)
end
