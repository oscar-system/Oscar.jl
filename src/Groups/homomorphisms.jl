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
  mp = GAPWrap.GroupHomomorphismByFunction(G.X, H.X, GAP.Obj(gap_fun))
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
    mp = GAPWrap.GroupHomomorphismByFunction(G.X, H.X, GAP.Obj(gap_fun), GAP.Obj(gap_pre_fun))
  else
    mp = GAPWrap.GroupHomomorphismByFunction(G.X, H.X, GAP.Obj(gap_fun), false, GAP.Obj(gap_pre_fun))
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
  @req mp !== GAP.Globals.fail "Invalid input"
  return GAPGroupHomomorphism(G, H, mp)
end

function hom(G::GAPGroup, H::GAPGroup, imgs::Vector; check::Bool = true)
  return hom(G, H, gens(G), imgs; check)
end

# Map `G::GAPGroup` to `A::GrpAbFinGen` by prescribing images.
# Return a composition of homomorphisms `G -> G/G' -> B -> A`,
# not a `GAPGroupHomomorphism`.
function hom(G::GAPGroup, A::GrpAbFinGen, gensG::Vector, imgs::Vector{GrpAbFinGenElem}; check::Bool = true)
  # map G to G/G'
  (q, map1) = quo(G, derived_subgroup(G)[1])

  # map G/G' to an isomorphic additive group B
  iso = isomorphism(GrpAbFinGen, q)
  B = codomain(iso)

  # map B to A as prescribed
  if length(gensG) == 0
    map2 = hom([zero(B)], [zero(A)], check = check)
  else
    map2 = hom([iso(map1(x)) for x in gensG], imgs, check = check)
  end

  # create the composition
  return compose(map1, compose(iso, map2))
end

function hom(G::GAPGroup, A::GrpAbFinGen, imgs::Vector{GrpAbFinGenElem}; check::Bool = true)
  return hom(G, A, gens(G), imgs; check)
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
    is_surjective(f::GAPGroupHomomorphism)

Return whether `f` is surjective.
"""
function is_surjective(f::GAPGroupHomomorphism)
  return GAPWrap.IsSurjective(f.map)
end

"""
    is_injective(f::GAPGroupHomomorphism)

Return whether `f` is injective.
"""
function is_injective(f::GAPGroupHomomorphism)
  return GAPWrap.IsInjective(f.map)
end

"""
    is_invertible(f::GAPGroupHomomorphism)

Return whether `f` is invertible.
"""
function is_invertible(f::GAPGroupHomomorphism)
  return GAPWrap.IsBijective(f.map)
end

"""
    is_bijective(f::GAPGroupHomomorphism)

Return whether `f` is bijective.
"""
function is_bijective(f::GAPGroupHomomorphism)
  return GAPWrap.IsBijective(f.map)
end


"""
    is_invariant(f::GAPGroupHomomorphism, H::Group)
    is_invariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return whether `f(H) == H` holds.
An exception is thrown if `domain(f)` and `codomain(f)` are not equal
or if `H` is not contained in `domain(f)`.
"""
function is_invariant(f::GAPGroupHomomorphism, H::GAPGroup)
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
  @assert is_subset(H, domain(f)) "Not a subgroup!"
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
    is_isomorphic_with_map(G::Group, H::Group)

Return (`true`,`f`) if `G` and `H` are isomorphic groups, where `f` is a group
isomorphism. Otherwise, return (`false`,`f`), where `f` is the trivial
homomorphism.

# Examples
```jldoctest
julia> is_isomorphic_with_map(symmetric_group(3), dihedral_group(6))
(true, Group homomorphism from
Sym( [ 1 .. 3 ] )
to
<pc group of size 6 with 2 generators>)
```
"""
function is_isomorphic_with_map(G::GAPGroup, H::GAPGroup)
  mp = GAP.Globals.IsomorphismGroups(G.X, H.X)::GapObj
  if mp === GAP.Globals.fail
    return false, trivial_morphism(G, H)
  else
    return true, GAPGroupHomomorphism(G, H, mp)
  end
end

function is_isomorphic_with_map(G::GAPGroup, H::GrpGen)
  HtoP = isomorphism(PermGroup, H)
  P = codomain(HtoP)
  fl, GtoP = is_isomorphic_with_map(G, P)
  if !fl
    return false, nothing
  else
    return true, GtoP * inv(HtoP)
  end
end

function is_isomorphic_with_map(G::GrpGen, H::GAPGroup)
  GtoP = isomorphism(PermGroup, G)
  P = codomain(GtoP)
  fl, PtoH = is_isomorphic_with_map(P, H)
  if !fl
    return false, nothing
  else
    return true, GtoP * PtoH
  end
end

"""
    is_isomorphic(G::Group, H::Group)

Return `true` if `G` and `H` are isomorphic groups, and `false` otherwise.

# Examples
```jldoctest
julia> is_isomorphic(symmetric_group(3), dihedral_group(6))
true
```
"""
function is_isomorphic(G::GAPGroup, H::GAPGroup)
  mp = GAP.Globals.IsomorphismGroups(G.X, H.X)::GapObj
  return mp !== GAP.Globals.fail
end

function isisomorphic(G::GAPGroup, H::GrpGen)
  P = PermGroup(H)
  return isisomorphic(G, P)
end

function isisomorphic(G::GrpGen, H::GAPGroup)
  P = PermGroup(G)
  return isisomorphic(P, H)
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
  @req mp !== GAP.Globals.fail "the groups are not isomorphic"
  return GAPGroupHomomorphism(G, H, mp)
end

function isomorphism(G::GAPGroup, H::GrpGen)
  HtoP = isomorphism(PermGroup, H)
  P = codomain(HtoP)
  fl, GtoP = is_isomorphic_with_map(G, P)
  @req fl "the groups are not isomorphic"
  return GtoP * inv(HtoP)
end

function isomorphism(G::GrpGen, H::GAPGroup)
  GtoP = isomorphism(PermGroup, G)
  P = codomain(GtoP)
  fl, PtoH = is_isomorphic_with_map(P, H)
  @req fl "the groups are not isomorphic"
  return GtoP * PtoH
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
   isos = get_attribute!(Dict{Type, Any}, G, :isomorphisms)::Dict{Type, Any}
   return get!(isos, T) do
     fun = _get_iso_function(T)
     f = fun(G.X)::GapObj
     @req f !== GAP.Globals.fail "Could not convert group into a group of type $T"
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
   isos = get_attribute!(Dict{Type, Any}, G, :isomorphisms)::Dict{Type, Any}
   return get!(isos, GrpAbFinGen) do
     @req is_abelian(G) "the group is not abelian"
     @req is_finite(G) "the group is not finite"
#T this restriction is not nice

     indep = GAP.Globals.IndependentGeneratorsOfAbelianGroup(G.X)::GapObj
     orders = ZZRingElem[GAPWrap.Order(x) for x in indep]
     n = length(indep)
     A = abelian_group(GrpAbFinGen, orders)

     f(g) = A(Vector{ZZRingElem}(GAPWrap.IndependentGeneratorExponents(G.X, g.X)))

     finv = function(g::elem_type(GrpAbFinGen))
       res = GAPWrap.One(G.X)
       for i in 1:n
         res = res * indep[i]^GAP.Obj(g.coeff[i])
       end
       return group_element(G, res)
     end

     return GroupIsomorphismFromFunc(f, finv, G, A)
   end::GroupIsomorphismFromFunc{typeof(G), GrpAbFinGen}
end

"""
    isomorphism(::Type{T}, A::GrpAbFinGen) where T <: Union{GAPGroup, GrpAbFinGen}

Return an isomorphism from `A` to a group of type `T`.
An exception is thrown if no such isomorphism exists or if `A` is not finite.
"""
function isomorphism(::Type{T}, A::GrpAbFinGen) where T <: GAPGroup
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, A, :isomorphisms)::Dict{Type, Any}
   return get!(isos, T) do
     # find independent generators
     if is_diagonal(rels(A))
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
     # `GAPWrap.GeneratorsOfGroup(G.X)` consists of independent elements
     # of the orders in `exponents`.
     # `GAP.Globals.IndependentGenerators(G.X)` chooses generators
     # that may differ from these generators,
     # and that belong to the exponent vectors returned by
     # `GAPWrap.IndependentGeneratorExponents(G.X, g)`.
     # `GAPWrap.GeneratorsOfGroup(G.X)` corresponds to `gens(A2)`,
     # we let `hom` compute elements in `A2` that correspond to
     # `GAP.Globals.IndependentGenerators(G.X)`.
     Ggens = Vector{GapObj}(GAPWrap.GeneratorsOfGroup(G.X)::GapObj)
     if length(Ggens) < length(exponents)
       # It may happen that GAP omits the generators of order 1. Insert them.
       @assert length(Ggens) + length(filter(x -> x == 1, exponents)) ==
               length(exponents)
       o = one(G).X
       newGgens = []
       pos = 1
       for i in 1:length(exponents)
         if exponents[i] == 1
           push!(newGgens, o)
         else
           push!(newGgens, Ggens[pos])
           pos = pos+1
         end
       end
       Ggens = newGgens
     end
     gensindep = GAP.Globals.IndependentGeneratorsOfAbelianGroup(G.X)::GapObj
     Aindep = abelian_group(ZZRingElem[GAPWrap.Order(g) for g in gensindep])

     imgs = [Vector{ZZRingElem}(GAPWrap.IndependentGeneratorExponents(G.X, a)) for a in Ggens]
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
       exp = Vector{ZZRingElem}(GAPWrap.IndependentGeneratorExponents(G.X, g.X))
       return Aindep_to_A(Aindep(exp))
     end

     return GroupIsomorphismFromFunc(f, finv, A, G)
   end::GroupIsomorphismFromFunc{GrpAbFinGen, T}
end

####
mutable struct GroupIsomorphismFromFunc{R, T} <: Map{R, T, Hecke.HeckeMap, MapFromFunc}
    map::MapFromFunc{R, T}
end

function GroupIsomorphismFromFunc{R, T}(f, g, D::R, C::T) where {R, T}
  return GroupIsomorphismFromFunc{R, T}(MapFromFunc(f, g, D, C))
end

function GroupIsomorphismFromFunc(f, g, D, C)
  return GroupIsomorphismFromFunc{typeof(D), typeof(C)}(f, g, D, C)
end

# install the same methods as for `MapFromFunc`,
# see `Hecke.jl/src/Map/MapType.jl`

domain(f::GroupIsomorphismFromFunc) = domain(f.map)

codomain(f::GroupIsomorphismFromFunc) = codomain(f.map)

image_function(f::GroupIsomorphismFromFunc) = image_function(f.map)

preimage_function(f::GroupIsomorphismFromFunc) = preimage_function(f.map)

image(f::GroupIsomorphismFromFunc, x) = image(f.map, x)

preimage(f::GroupIsomorphismFromFunc, y) = preimage(f.map, y)

function Base.show(io::IO, M::GroupIsomorphismFromFunc)
  Base.show(io, M.map)
end

Base.inv(M::GroupIsomorphismFromFunc) = GroupIsomorphismFromFunc(inv(M.map))

# additional methods

is_bijective(f::GroupIsomorphismFromFunc) = true

kernel(f::GroupIsomorphismFromFunc) = trivial_subgroup(domain(f))

function images(f::GroupIsomorphismFromFunc{R, T}, G::R) where R where T
  D = domain(f)
  C = codomain(f)
  imgs = eltype(typeof(C))[]
  for x in gens(G)
    if parent(x) === D
      push!(imgs, f(x))
    else
      push!(imgs, f(D(x)))
    end
  end
  return sub(codomain(f), imgs)
end

####

# compute the kernel of a composition of maps, with domain a `GAPGroup`,
# where the kernels of the composed maps can be computed

function kernel(comp::AbstractAlgebra.Generic.CompositeMap{T, GrpAbFinGen}) where T <: GAPGroup
  map1 = comp.map1
  map2 = comp.map2

  if map2 isa GrpAbFinGenMap
    ker2 = kernel(map2, false)::Tuple{GrpAbFinGen, GrpAbFinGenMap}
  else
    ker2 = kernel(map2)
  end
  ker2gens = elem_type(domain(map2))[ker2[2](x) for x in gens(ker2[1])]
  preimages = elem_type(domain(map1))[preimage(map1, x) for x in ker2gens]
  ker1 = kernel(map1)

  # Compute generators of the kernel of `map2`,
  # take their preimages under `map1`,
  # form the closure with the kernel of `map1`
  G = domain(comp)
  K = sub(G, vcat(elem_type(domain(map1))[ker1[2](x) for x in gens(ker1[1])], preimages))
end

####

function isomorphism(::Type{GrpAbFinGen}, A::GrpAbFinGen)
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, A, :isomorphisms)::Dict{Type, Any}
   return get!(isos, GrpAbFinGen) do
     return identity_map(A)
   end::AbstractAlgebra.Generic.IdentityMap{GrpAbFinGen}
end

# We need not find independent generators in order to create
# a presentation of a fin. gen. abelian group.
function isomorphism(::Type{FPGroup}, A::GrpAbFinGen)
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, A, :isomorphisms)::Dict{Type, Any}
   return get!(isos, FPGroup) do
      G = free_group(ngens(A); eltype = :syllable)
      R = rels(A)
      s = vcat(elem_type(G)[i*j*inv(i)*inv(j) for i = gens(G) for j = gens(G) if i != j],
           elem_type(G)[prod([gen(G, i)^R[j,i] for i=1:ngens(A) if !iszero(R[j,i])], init = one(G)) for j=1:nrows(R)])
      F, mF = quo(G, s)
      @assert is_finite(A) == is_finite(F)
      is_finite(A) && @assert order(A) == order(F)
      return MapFromFunc(
        y->F([i => y[i] for i=1:ngens(A)]),
        x->sum([w.second*gen(A, w.first) for w = syllables(x)], init = zero(A)),
        A, F)
   end::MapFromFunc{GrpAbFinGen, FPGroup}
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

# Now for GrpGen

# Remove once Hecke lower bound is >= 0.14.1
@attributes GrpGen

function isomorphism(::Type{T}, A::GrpGen) where T <: GAPGroup
   # Known isomorphisms are cached in the attribute `:isomorphisms`.
   isos = get_attribute!(Dict{Type, Any}, A, :isomorphisms)::Dict{Type, Any}
   return get!(isos, T) do
     S = symmetric_group(order(A))
     gensA = gens(A)
     newgens = elem_type(S)[]
     for g in gensA
       j = g.i
       p = S(A.mult_table[j, :])
       push!(newgens, p)
     end

     GP, m = sub(S, newgens)

     GPgens = map(x -> preimage(m, x), newgens)

     tospin = [(gensA[i], GPgens[i]) for i in 1:length(GPgens)]
     cl = Hecke.closure(tospin, (x, y) -> (x[1] * y[1], x[2] * y[2]),
                        eq = (x, y) -> x[1] == y[1])

     fwd = Dict{elem_type(A), elem_type(GP)}(c[1] => c[2] for c in cl)
     bwd = Dict{elem_type(GP), elem_type(A)}(c[2] => c[1] for c in cl)

     if T === PermGroup
       f = function(a::elem_type(GrpGen))
         return fwd[a]
       end

       finv = function(g)
         return bwd[g]
       end
       return MapFromFunc(f, finv, A, GP)
     else
       m = isomorphism(T, GP)

       f = function(a::elem_type(GrpGen))
         return m(fwd[a])
       end

       finv = function(g)
         return bwd[preimage(m, g)]
       end

       return MapFromFunc(f, finv, A, codomain(m))
     end
   end::MapFromFunc{GrpGen, T}
end

function (::Type{T})(G::GrpGen) where T <: GAPGroup
   return codomain(isomorphism(T, G))
end


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

Groups of automorphisms over a group `G` have parametric type `AutomorphismGroup{T}`, where `T` is the type of `G`. 
# Examples
```jldoctest
julia> S = symmetric_group(3)
Sym( [ 1 .. 3 ] )

julia> typeof(S)
PermGroup

julia> A = automorphism_group(S)
Aut( Sym( [ 1 .. 3 ] ) )

julia> typeof(A)
AutomorphismGroup{PermGroup}
```

The evaluation of the automorphism `f` in the element `x` is analogous to the homomorphism evaluation: 
it can be obtained by typing either `f(x)` or `x^f`.

```jldoctest
julia> S = symmetric_group(4)
Sym( [ 1 .. 4 ] )

julia> A = automorphism_group(S)
Aut( Sym( [ 1 .. 4 ] ) )

julia> x = perm(S,[2,1,4,3])
(1,2)(3,4)

julia> f = A[2]
Pcgs([ (3,4), (2,4,3), (1,4)(2,3), (1,3)(2,4) ]) -> [ (2,3), (2,4,3), (1,3)(2,4), (1,2)(3,4) ]

julia> f(x)
(1,4)(2,3)

julia> x^f
(1,4)(2,3)
```

It is possible to turn an automorphism `f` into a homomorphism by typing `hom(f)`.

```jldoctest
julia> S = symmetric_group(4)
Sym( [ 1 .. 4 ] )

julia> A = automorphism_group(S)
Aut( Sym( [ 1 .. 4 ] ) )

julia> f = A[2]
Pcgs([ (3,4), (2,4,3), (1,4)(2,3), (1,3)(2,4) ]) -> [ (2,3), (2,4,3), (1,3)(2,4), (1,2)(3,4) ]

julia> typeof(f)
AutomorphismGroupElem{PermGroup} (alias for Oscar.BasicGAPGroupElem{AutomorphismGroup{PermGroup}})

julia> typeof(hom(f))
GAPGroupHomomorphism{PermGroup, PermGroup}
```
The converse is also possible: if `g` is a bijective homomorphism from the group `G` to itself and `A` 
is the automorphism group of `G`, then the instruction `A(g)` returns `g` as automorphism of `G`. 
This is the standard way to explicitly build an automorphism (another way, available for inner 
automorphisms, is shown in Section [Inner_automorphisms](@ref inner_automorphisms)).

# Examples
```jldoctest
julia> S = symmetric_group(4)
Sym( [ 1 .. 4 ] )

julia> a = perm(S,[2,1,4,3])
(1,2)(3,4)

julia> f = hom(S,S,x ->x^a)
Group homomorphism from 
Sym( [ 1 .. 4 ] )
to
Sym( [ 1 .. 4 ] )

julia> A = automorphism_group(S)
Aut( Sym( [ 1 .. 4 ] ) )

julia> A(f)
MappingByFunction( Sym( [ 1 .. 4 ] ), Sym( [ 1 .. 4 ] ), <Julia: gap_fun> )
```

Elements of `A` can be multiplied with other elements of `A` or by elements
of type `GAPGroupHomomorphism`; in this last case, the result has type
`GAPGroupHomomorphism`.

# Examples
```jldoctest
julia> S = symmetric_group(4);

julia> A = automorphism_group(S);

julia> g = hom(S,S,x->x^S[1]);

julia> g in A
false

julia> au = A(g);

julia> au in A
true

julia> g == hom(au)
true

julia> x = cperm(S,[1,2,3]);

julia> au(x)
(2,3,4)

julia> g(x) == au(x)
true
```

In Oscar it is possible to multiply homomorphisms and automorphisms (whenever it makes sense); in such cases, the output is always a variable of type `GAPGroupHomomorphism{S,T}`.
```jldoctest
julia> S = symmetric_group(4)
Sym( [ 1 .. 4 ] )

julia> A = automorphism_group(S)
Aut( Sym( [ 1 .. 4 ] ) )

julia> g = hom(S,S,x->x^S[1])
Group homomorphism from 
Sym( [ 1 .. 4 ] )
to
Sym( [ 1 .. 4 ] )

julia> f = A(g)
MappingByFunction( Sym( [ 1 .. 4 ] ), Sym( [ 1 .. 4 ] ), <Julia: gap_fun> )

julia> typeof(g*f)
GAPGroupHomomorphism{PermGroup, PermGroup}
```
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
   @assert is_bijective(f) "f not in A"
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
    is_inner_automorphism(f::GAPGroupHomomorphism)
    is_inner_automorphism(f::GAPGroupElem{AutomorphismGroup{T}})

Return whether `f` is an inner automorphism.
"""
function is_inner_automorphism(f::GAPGroupHomomorphism)
  @assert domain(f) == codomain(f) "Not an automorphism!"
  return GAPWrap.IsInnerAutomorphism(f.map)
end

function is_inner_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
  return GAPWrap.IsInnerAutomorphism(f.X)
end

"""
    inner_automorphism_group(A::AutomorphismGroup{T})

Return the subgroup of `A` of the inner automorphisms.
"""
function inner_automorphism_group(A::AutomorphismGroup{T}) where T <: GAPGroup
   AutGAP = GAP.Globals.InnerAutomorphismsAutomorphismGroup(A.X)
   return _as_subgroup(A, AutGAP)
end

"""
    is_invariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return whether `f`(`H`) == `H`.
"""
function is_invariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T<:GAPGroup
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
  @assert is_invariant(mH, kernel(f)[1]) "The kernel is not invariant under g!"
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
  @assert is_invariant(f,H) "H is not invariant under f!"
  fh = hom(H, H, gens(H), [f(x) for x in gens(H)], check = false)
  return A(fh)
end

function restrict_homomorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: GAPGroup
  return restrict_homomorphism(hom(f),H)
end

function _as_subgroup_bare(G::AutomorphismGroup{T}, H::GapObj) where T
  return AutomorphismGroup{T}(H, G.G)
end
