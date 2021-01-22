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
    isomorphic_fp_group,
    isomorphic_pc_group,
    isomorphic_perm_group,
    issurjective,
    kernel,
    order,
    restrict_automorphism,
    restrict_homomorphism,
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

Base.:*(f::GAPGroupHomomorphism{S, T}, g::GAPGroupHomomorphism{T, U}) where S where T where U = compose(g, f)

function Base.inv(f::GAPGroupHomomorphism{S,T}) where S where T
   @assert GAP.Globals.IsBijective(f.map) "f is not bijective"
   return GAPGroupHomomorphism{T,S}(codomain(f), domain(f), GAP.Globals.InverseGeneralMapping(f.map))
end

order(f::GAPGroupHomomorphism) = GAP.Globals.Order(f.map)

function Base.:^(f::GAPGroupHomomorphism{S,T}, n::Int64) where S where T
   if n==1
     return f
   else
      @assert domain(f) == codomain(f) "Domain and codomain do not coincide"
      return GAPGroupHomomorphism{S,S}(domain(f),codomain(f),(f.map)^n)
   end
end

function compose(g::GAPGroupHomomorphism{T, U}, f::GAPGroupHomomorphism{S, T}) where S where T where U
  dom = domain(f)
  cod = codomain(g)
  @assert codomain(f) == domain(g)
  mp = GAP.Globals.CompositionMapping(g.map, f.map)
  return GAPGroupHomomorphism{S, U}(dom, cod, mp)
end

(g::GAPGroupHomomorphism{T, U})(f::GAPGroupHomomorphism{S, T}) where {S,T,U} = compose(g,f)

"""
    id_hom(G::GAPGroup)

Return the identity homomorphism on the group `G`.
"""
function id_hom(G::GAPGroup)
  return hom(G, G, x -> x)
end

"""
    trivial_morphism(G::GAPGroup, H::GAPGroup)

Return the homomorphism from `G` to `H` sending every element of `G` into the identity of `H`. If `H` is not specified, it is taken equal to `G`.
"""
function trivial_morphism(G::GAPGroup, H::GAPGroup)
  return hom(G, H, x -> one(H))
end

function trivial_morphism(G::GAPGroup)
  return hom(G, G, x -> one(G))
end

function _hom_from_gap_map(G::GAPGroup, H::GAPGroup, mp::GapObj)
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

"""
    hom(G::GAPGroup, H::GAPGroup, f::Function)

Return the group homomorphism defined by the function ``f``.
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
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

"""
    hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector)

Return the group homomorphism defined by `gensG`[`i`] -> `imgs`[`i`] for every `i`. In order to work, the elements of `gensG` must generate `G`.
"""
function hom(G::GAPGroup, H::GAPGroup, gensG::Vector, imgs::Vector)
  vgens = GAP.julia_to_gap(GapObj[x.X for x in gensG])
  vimgs = GAP.julia_to_gap(GapObj[x.X for x in imgs])
  mp = GAP.Globals.GroupHomomorphismByImages(G.X, H.X, vgens, vimgs)
  if mp == GAP.Globals.fail throw(ArgumentError("Invalid input")) end
  return GAPGroupHomomorphism{typeof(G), typeof(H)}(G, H, mp)
end

"""
    domain(f::GAPGroupHomomorphism)

Return the domain of `f`.
"""
function domain(f::GAPGroupHomomorphism)
  return f.domain
end

"""
    codomain(f::GAPGroupHomomorphism)

Return the codomain of `f`.
"""
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

"""
    issurjective(f::GAPGroupHomomorphism)

Return whether `f` is surjective.
"""
function issurjective(f::GAPGroupHomomorphism)
  return GAP.Globals.IsSurjective(f.map)
end

"""
    isinjective(f::GAPGroupHomomorphism)

Return whether `f` is injective.
"""
function isinjective(f::GAPGroupHomomorphism)
  return GAP.Globals.IsInjective(f.map)
end

"""
    isinvertible(f::GAPGroupHomomorphism)

Return whether `f` is invertible.
"""
function isinvertible(f::GAPGroupHomomorphism)
  return GAP.Globals.IsBijective(f.map)
end

"""
    isbijective(f::GAPGroupHomomorphism)

Return whether `f` is bijective.
"""
function isbijective(f::GAPGroupHomomorphism)
  return GAP.Globals.IsBijective(f.map)
end


"""
    isinvariant(f::GAPGroupHomomorphism, H::Group)
    isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return whether `f`(`H`) == `H`.
"""
function isinvariant(f::GAPGroupHomomorphism, H::GAPGroup)
  @assert domain(f) == codomain(f) "Not an endomorphism!"
  @assert GAP.Globals.IsSubset(codomain(f).X, H.X) "Not a subgroup of the domain"
  return GAP.Globals.Image(f.map, H.X) == H.X
end

"""
    restrict_homomorphism(f::GAPGroupHomomorphism, H::Group)
    restrict_homomorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T <: Group

Return the restriction of `f` to `H`; otherwise it return ERROR.
"""
function restrict_homomorphism(f::GAPGroupHomomorphism, H::GAPGroup)
  return _hom_from_gap_map(H, f.codomain, GAP.Globals.RestrictedMapping(f.map,H.X))
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
    (f::GAPGroupHomomorphism{S, T})(H::S)

Return `f`(`H`), together with the embedding homomorphism into `codomain`(`f`).
"""
function image(f::GAPGroupHomomorphism{S, T}, H::S) where S <: GAPGroup where T <: GAPGroup
  H1 = GAP.Globals.Image(f.map, H.X)
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
    haspreimage(f::GAPGroupHomomorphism, x::GAPGroupElem)

Return (`true`,`y`) if there exists `y` such that `f`(`y`) = `x`; otherwise, return (`false`,`1`).
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

"""
    isomorphic_perm_group(G::GAPGroup)

Return a permutation group `H` and an isomorphism `f` from `G` to `H`.

If `G` is infinite, then no such isomorphism exists and an exception is thrown.
"""
function isomorphic_perm_group(G::GAPGroup)
   f = GAP.Globals.IsomorphismPermGroup(G.X)
   f!=GAP.Globals.fail || throw(ArgumentError("Could not convert group into a permutation group"))
   H = GAP.Globals.Image(f)
   n = GAP.Globals.NrMovedPoints(H)
   H = PermGroup(H,n)
   return H, _hom_from_gap_map(G,H,f)
end

"""
    isomorphic_pc_group(G::GAPGroup)

Return a group `H` of type `PcGroup` and an isomorphism `f` from `G` to `H`.

If `G` is infinite or not solvable, then no such isomorphism exists and an exception is thrown.
"""
function isomorphic_pc_group(G::GAPGroup)
   f = GAP.Globals.IsomorphismPcGroup(G.X)
   f!=GAP.Globals.fail || throw(ArgumentError("Could not convert group into a group of type PcGroup"))
   H = PcGroup(GAP.Globals.Image(f))
   return H, _hom_from_gap_map(G,H,f)
end

"""
    isomorphic_fp_group(G::GAPGroup)

Return a group `H` of type `FPGroup` and an isomorphism `f` from `G` to `H`.
"""
function isomorphic_fp_group(G::GAPGroup)
   f = GAP.Globals.IsomorphismFpGroup(G.X)
   f!=GAP.Globals.fail || throw(ArgumentError("Could not convert group into a group of type FPGroup"))
   H = FPGroup(GAP.Globals.Image(f))
   return H, _hom_from_gap_map(G,H,f)
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

"""
    hom(f::GAPGroupElem{AutomorphismGroup{T}}) where T

Return the element f of type ``GAPGroupHomomorphism{T,T}``.
"""
function hom(x::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup
  A = parent(x)
  G = A.G
  return _hom_from_gap_map(G, G, x.X)
end

(f::GAPGroupElem{AutomorphismGroup{T}})(x::GAPGroupElem) where T <: GAPGroup = apply_automorphism(f, x)
Base.:^(x::GAPGroupElem,f::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = apply_automorphism(f, x)
Base.:^(f::GAPGroupElem{AutomorphismGroup{T}},g::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = g^-1*f*g

function (A::AutomorphismGroup{T})(f::GAPGroupHomomorphism{T,T}) where T <: GAPGroup
   @assert domain(f)==A.G && codomain(f)==A.G "f not in A"
   @assert isbijective(f) "f not in A"
   return group_element(A, f.map)
end

function apply_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, x::GAPGroupElem, check=true) where T <: GAPGroup
  A = parent(f)
  G = parent(x)
  if check
    @assert A.G == G || GAP.Globals.IN(x.X, A.G.X) "Not in the domain of f!"      #TODO Do we really need the IN check?
  end
  if x isa MatrixGroupElem
     return MatrixGroupElem(G,GAP.Globals.Image(f.X,x.X))
  else
     return typeof(x)(G, GAP.Globals.Image(f.X,x.X))
  end
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
    isinner_automorphism(f::GAPGroupElem{AutomorphismGroup{T}})

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

"""
    isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

Return whether `f`(`H`) == `H`.
"""
function isinvariant(f::GAPGroupElem{AutomorphismGroup{T}}, H::T) where T<:GAPGroup
  @assert GAP.Globals.IsSubset(parent(f).G.X, H.X) "Not a subgroup of the domain"
  return GAP.Globals.Image(f.X, H.X) == H.X
end

"""
    induced_automorphism(f::GAPGroupHomomorphism, g::GAPGroupHomomorphism)
    induced_automorphism(f::GAPGroupHomomorphism, g::GAPGroupElem{AutomorphismGroup{T}})

Return the automorphism `h` of the image of `f` such that `h`(`f`) == `f`(`g`), where `g` is an automorphism of a group `G` and `f` is a group homomorphism defined over `G` such that the kernel of `f` is invariant under `g`
"""
function induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupHomomorphism)
  @assert isinvariant(mH, kernel(f)[1]) "The kernel is not invariant under g!"
  map = GAP.Globals.InducedAutomorphism(f.map, mH.map)
  A = automorphism_group(image(f)[1])
  return A(_hom_from_gap_map(codomain(f), codomain(f), map))
end

induced_automorphism(f::GAPGroupHomomorphism, mH::GAPGroupElem{AutomorphismGroup{T}}) where T <: GAPGroup = induced_automorphism(f,hom(mH))

"""
    restrict_automorphism(f::GAPGroupElem{AutomorphismGroup{T}}, H::T)

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


