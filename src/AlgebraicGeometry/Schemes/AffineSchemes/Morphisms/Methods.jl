###########################################################
# (1) The fiber product of two morphisms of affine schemes
###########################################################

@doc raw"""
    fiber_product(f::AbsAffineSchemeMor, g::AbsAffineSchemeMor)

For morphisms ``f : X → Z`` and ``g : Y → Z`` return the fiber
product ``X×Y`` over ``Z`` together with its two canonical projections.

Whenever you have another set of maps `a: W → X` and `b : W → Y` forming
a commutative square with `f` and `g`, you can use
`induced_map_to_fiber_product` to  create the resulting map `W → X×Y`.
"""
function fiber_product(
    f::AbsAffineSchemeMor,
    g::AbsAffineSchemeMor
  )
  Y = domain(f)
  X = codomain(f)
  X == codomain(g) || error("maps need to have the same codomain")
  Z = domain(g)
  YxZ, pY, pZ = product(Y, Z)
  RX = ambient_coordinate_ring(X)
  W = subscheme(YxZ, [pullback(pY)(pullback(f)(x)) - pullback(pZ)(pullback(g)(x)) for x in gens(RX)])
  return W, restrict(pY, W, Y, check=false), restrict(pZ, W, Z, check=false)
end

# Whenever one of the maps, say f, in a fiber product is a `PrincipalOpenEmbedding`
# then the fiber product is only the restriction of g to g^{-1}(image(f)).
# This can be computed much easier and, in particular, without introducing
# extra variables: One just pulls back the `complement_equations` for `f` to
# the domain of `g`.
#
# We need this more simple procedure for refinements of coverings. In particular,
# it is important that the resulting fiber product is a PrincipalOpenSubset of
# the domain of `g` (or the domain of `f` when it's the other way around), so that
# the ancestry-tree for patches is preserved.
function fiber_product(f::PrincipalOpenEmbedding, g::AbsAffineSchemeMor)
  @assert codomain(f) === codomain(g) "codomains are not the same"
  A = domain(f)
  B = domain(g)
  C = codomain(f)
  h = complement_equations(f)
  pbh = pullback(g).(h)
  result = PrincipalOpenSubset(B, pbh)
  ff = PrincipalOpenEmbedding(morphism(result, B, gens(OO(result)), check=false), pbh, check=false)
  f_res_inv = inverse_on_image(f)
  gg = compose(restrict(g, result, image(f), check=false), f_res_inv)
  return result, gg, ff
end

function fiber_product(f::AbsAffineSchemeMor, g::PrincipalOpenEmbedding)
  result, ff, gg = fiber_product(g, f)
  return result, gg, ff
end

# additional method to remove method ambiguity
function fiber_product(f::PrincipalOpenEmbedding, g::PrincipalOpenEmbedding)
  @assert codomain(f) === codomain(g) "codomains are not the same"
  A = domain(f)
  B = domain(g)
  C = codomain(f)
  h = complement_equations(f)
  pbh = pullback(g).(h)
  result = PrincipalOpenSubset(B, pbh)
  ff = PrincipalOpenEmbedding(morphism(result, B, gens(OO(result)), check=false), pbh, check=false)
  f_res_inv = inverse_on_image(f)
  gg = compose(restrict(g, result, image(f), check=false), f_res_inv)
  hg = complement_equations(g)
  pbhg = pullback(f).(hg)
  return result, PrincipalOpenEmbedding(gg, pbhg, check=false), ff
end

@doc raw"""
    induced_map_to_fiber_product(
        a::AbsAffineSchemeMor, b::AbsAffineSchemeMor,
        f::AbsAffineSchemeMor, g::AbsAffineSchemeMor;
        fiber_product::Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:AbsAffineSchemeMor}=fiber_product(f, g)
      )

In a commutative diagram
```
          b
   W ------------.
   |             |
   |             V
  a|    X x Y -->Y
   |      |      | g
   |      V      V
   `----->X----> Z
             f
```
this computes the canonical map `W -> X x Y`.
"""
function induced_map_to_fiber_product(
    a::AbsAffineSchemeMor, b::AbsAffineSchemeMor,
    f::AbsAffineSchemeMor, g::AbsAffineSchemeMor;
    fiber_product::Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:AbsAffineSchemeMor}=fiber_product(f, g),
    check::Bool=true
  )
  # All checks are done here. The actual computations are carried out
  # in an internal method.
  X = domain(f)
  Y = domain(g)
  Z = codomain(f)
  @assert codomain(g) === Z
  XxY = fiber_product[1]
  gg = fiber_product[2]
  ff = fiber_product[3]
  @assert codomain(ff) === Y
  @assert codomain(gg) === X

  W = domain(a)
  @assert W === domain(b)
  @assert codomain(a) === X
  @assert codomain(b) === Y
  @check compose(a, f) == compose(b, g) "maps do not commute"
  @check compose(ff, g) == compose(gg, f) "maps do not commute"
  return _induced_map_to_fiber_product(a, b, f, g, fiber_product=fiber_product, check=check)
end

function _induced_map_to_fiber_product(
    a::AbsAffineSchemeMor, b::AbsAffineSchemeMor,
    f::AbsAffineSchemeMor, g::AbsAffineSchemeMor;
    fiber_product::Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:AbsAffineSchemeMor}=fiber_product(f, g),
    check::Bool=true
  )
  # The ambient scheme of XxY is the actual product of X and Y
  # over Spec(k), the coefficient ring. If it is not, then
  # this is due to special dispatch which has to also be caught
  # with a special method for this function here.
  XxY = fiber_product[1]
  gg = fiber_product[2]
  ff = fiber_product[3]
  X = domain(f)
  Y = domain(g)
  W = domain(a)
  @check gens(OO(XxY)) == vcat(pullback(gg).(gens(OO(X))), pullback(ff).(gens(OO(Y)))) "variables must be pullbacks of variables on the factors"

  img_gens = vcat(_images_of_generators(pullback(a)), _images_of_generators(pullback(b)))
  return morphism(W, XxY, img_gens, check=check)
end

# When the fiber product was created from at least one `PrincipalOpenEmbedding`,
# then the construction did not proceed via the `product` of `X` and `Y`.
# In this case, the induced map must be created differently.
function _induced_map_to_fiber_product(
    a::AbsAffineSchemeMor, b::AbsAffineSchemeMor,
    f::PrincipalOpenEmbedding, g::AbsAffineSchemeMor;
    fiber_product::Tuple{<:AbsAffineScheme, <:AbsAffineSchemeMor, <:PrincipalOpenEmbedding}=fiber_product(f, g),
    check::Bool=true
  )
  # XxY is a principal open subset of Y.
  XxY = fiber_product[1]
  W = domain(a)
  Y = codomain(b)
  return morphism(W, XxY, _images_of_generators(pullback(b)), check=check)
end

function _induced_map_to_fiber_product(
    a::AbsAffineSchemeMor, b::AbsAffineSchemeMor,
    f::AbsAffineSchemeMor, g::PrincipalOpenEmbedding;
    fiber_product::Tuple{<:AbsAffineScheme, <:PrincipalOpenEmbedding, <:AbsAffineSchemeMor}=fiber_product(f, g),
    check::Bool=true
  )
  # XxY is a principal open subset of X.
  XxY = fiber_product[1]
  X = domain(f)
  W = domain(a)
  return morphism(W, XxY, _images_of_generators(pullback(a)), check=check)
end

# additional method to remove ambiguity
function _induced_map_to_fiber_product(
    a::AbsAffineSchemeMor, b::AbsAffineSchemeMor,
    f::PrincipalOpenEmbedding, g::PrincipalOpenEmbedding;
    fiber_product::Tuple{<:AbsAffineScheme, <:PrincipalOpenEmbedding, <:PrincipalOpenEmbedding}=fiber_product(f, g),
    check::Bool=true
  )
  W = domain(a)
  # XxY is a principal open subset of Y.
  XxY = fiber_product[1]
  Y = domain(g)
  return morphism(W, XxY, _images_of_generators(pullback(b)), check=check)
end

### Some helper functions

_images_of_generators(f::Map) = f.(gens(domain(f)))
_images_of_generators(f::MPolyAnyMap) = _images(f)
_images_of_generators(f::MPolyLocalizedRingHom) = _images(restricted_map(f))
_images_of_generators(f::MPolyQuoLocalizedRingHom) = _images(restricted_map(f))

function _restrict_domain(f::AbsAffineSchemeMor, D::PrincipalOpenSubset; check::Bool=true)
  D === domain(f) && return f
  !_has_coefficient_map(pullback(f)) && ambient_scheme(D) === domain(f) && return morphism(D, codomain(f), [OO(D)(x; check) for x in _images_of_generators(pullback(f))]; check=check)

  @check is_subscheme(D, domain(f)) "domain incompatible"
  !_has_coefficient_map(pullback(f)) && return morphism(D, codomain(f), [OO(D)(x; check) for x in _images_of_generators(pullback(f))], check=check)
  return morphism(D, codomain(f), coefficient_map(pullback(f)), [OO(D)(x; check) for x in _images_of_generators(pullback(f))], check=check)
end

function _restrict_domain(f::AbsAffineSchemeMor, D::AbsAffineScheme; check::Bool=true)
  D === domain(f) && return f
  inc = inclusion_morphism(D, domain(f); check)
  return compose(inc, f)
end

function _restrict_codomain(f::AbsAffineSchemeMor, D::PrincipalOpenSubset; check::Bool=true)
  D === codomain(f) && return f
  if ambient_scheme(D) === codomain(f)
    @check is_unit(pullback(f)(complement_equation(D))) "complement equation does not pull back to a unit"
    !_has_coefficient_map(pullback(f)) && return morphism(domain(f), D, [OO(domain(f))(x; check) for x in _images_of_generators(pullback(f))]; check=check)
    return morphism(domain(f), D, coefficient_map(pullback(f)), [OO(domain(f))(x; check) for x in _images_of_generators(pullback(f))], check=check)
  end
  @check is_subscheme(D, codomain(f)) "codomain incompatible"
  @check is_subscheme(domain(f), preimage(f, D))
  !_has_coefficient_map(pullback(f)) && return morphism(domain(f), D, [OO(domain(f))(x; check) for x in _images_of_generators(pullback(f))]; check=check)
  return morphism(domain(f), D, coefficient_map(pullback(f)), [OO(domain(f))(x; check) for x in pullback(f).(gens(OO(codomain(f))))], check=check)
end

# Some missing constructors
morphism(A::AbsAffineScheme, B::AbsAffineScheme, coeff_map::Any, img_gens::Vector; check::Bool=true) = morphism(A, B, hom(OO(B), OO(A), coeff_map, img_gens; check); check)

function hom(L::MPolyLocRing, P::NCRing, coeff_map::Any, img_gens::Vector; check::Bool=true)
  R = base_ring(L)
  phi = hom(R, P, coeff_map, img_gens; check)
  return MPolyLocalizedRingHom(L, P, phi; check)
end

function hom(L::MPolyQuoLocRing, P::NCRing, coeff_map::Any, img_gens::Vector; check::Bool=true)
  R = base_ring(L)
  phi = hom(R, P, coeff_map, img_gens; check)
  return MPolyQuoLocalizedRingHom(L, P, phi; check)
end

function _restrict_codomain(f::AbsAffineSchemeMor, D::AbsAffineScheme; check::Bool=true)
  @check is_subscheme(D, codomain(f)) "codomain incompatible"
  @check is_subscheme(domain(f), preimage(f, D)) "new domain is not contained in preimage of codomain"
  !_has_coefficient_map(pullback(f)) && return morphism(domain(f), D, [OO(domain(f))(x; check) for x in _images_of_generators(pullback(f))]; check)
  return morphism(domain(f), D, coefficient_map(pullback(f)), [OO(domain(f))(x; check) for x in pullback(f).(gens(OO(codomain(f))))], check=check)
end

@doc raw"""
    restrict(f::AbsAffineSchemeMor, D::AbsAffineScheme, Z::AbsAffineScheme; check::Bool=true)

This method restricts the domain of the morphism ``f``
to ``D`` and its codomain to ``Z``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> restrict(identity_map(X), Y, Y) == identity_map(Y)
true
```
"""
function restrict(f::AbsAffineSchemeMor, D::AbsAffineScheme, Z::AbsAffineScheme; check::Bool=true)
  interm = _restrict_domain(f, D; check)
  return _restrict_codomain(interm, Z; check)
end

function Base.:(==)(f::AbsAffineSchemeMor, g::AbsAffineSchemeMor)
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return pullback(f) == pullback(g)
end

###########################################################
# (2) The direct product of two affine schemes
###########################################################

# First the product of the ambient spaces. Documented below.
function product(X::AbsAffineScheme{BRT, RT}, Y::AbsAffineScheme{BRT, RT};
    change_var_names_to::Vector{String}=["", ""]
  ) where {BRT, RT<:MPolyRing}
  K = OO(X)
  L = OO(Y)
  # V = localized_ring(K)
  # W = localized_ring(L)
  k = base_ring(K)
  k == base_ring(L) || error("varieties are not defined over the same base ring")

  m = ngens(K)
  n = ngens(L)
  new_symb = Symbol[]
  if length(change_var_names_to[1]) == 0
    new_symb = symbols(K)
  else
    new_symb = Symbol.([change_var_names_to[1]*"$i" for i in 1:ngens(L)])
  end
  if length(change_var_names_to[2]) == 0
    new_symb = vcat(new_symb, symbols(L))
  else
    new_symb = vcat(new_symb, Symbol.([change_var_names_to[2]*"$i" for i in 1:ngens(L)]))
  end
  KL, z = polynomial_ring(k, new_symb)
  XxY = spec(KL)
  pr1 = morphism(XxY, X, gens(KL)[1:m], check=false)
  pr2 = morphism(XxY, Y, gens(KL)[m+1:m+n], check=false)
  return XxY, pr1, pr2
end

@doc raw"""
    product(X::AbsAffineScheme, Y::AbsAffineScheme)

Return a triple ``(X×Y, p₁, p₂)`` consisting of the product ``X×Y`` over
the common base ring ``𝕜`` and the two projections ``p₁ : X×Y → X`` and
``p₂ : X×Y → Y``.
"""
function product(X::AbsAffineScheme, Y::AbsAffineScheme;
    change_var_names_to::Vector{String}=["", ""]
  )
  # take the product of the ambient spaces and restrict
  base_ring(X) == base_ring(Y) || error("schemes are not defined over the same base ring")
  A = ambient_space(X)
  B = ambient_space(Y)
  AxB,prA, prB  = product(A, B, change_var_names_to=change_var_names_to)
  XxY = intersect(preimage(prA, X, check=false), preimage(prB, Y,check=false))
  prX = restrict(prA, XxY, X, check=false)
  prY = restrict(prB, XxY, Y, check=false)
  return XxY, prX, prY
end



#=
function product(X::StdAffineScheme, Y::StdAffineScheme;
    change_var_names_to::Vector{String}=["", ""]
  )
  K = OO(X)
  L = OO(Y)
  V = localized_ring(K)
  W = localized_ring(L)
  R = base_ring(K)
  S = base_ring(L)
  k = base_ring(R)
  k == base_ring(S) || error("varieties are not defined over the same field")

  m = ngens(R)
  n = ngens(S)
  new_symb = Symbol[]
  if length(change_var_names_to[1]) == 0
    new_symb = symbols(R)
  else
    new_symb = Symbol.([change_var_names_to[1]*"$i" for i in 1:ngens(R)])
  end
  if length(change_var_names_to[2]) == 0
    new_symb = vcat(new_symb, symbols(S))
  else
    new_symb = vcat(new_symb, Symbol.([change_var_names_to[2]*"$i" for i in 1:ngens(S)]))
  end
  RS, z = polynomial_ring(k, new_symb)
  inc1 = hom(R, RS, gens(RS)[1:m], check=false)
  inc2 = hom(S, RS, gens(RS)[m+1:m+n], check=false)
  IX = ideal(RS, inc1.(gens(modulus(underlying_quotient(OO(X))))))
  IY = ideal(RS, inc2.(gens(modulus(underlying_quotient(OO(Y))))))
  UX = MPolyPowersOfElement(RS, inc1.(denominators(inverted_set(OO(X)))))
  UY = MPolyPowersOfElement(RS, inc2.(denominators(inverted_set(OO(Y)))))
  XxY = spec(RS, IX + IY, UX*UY)
  pr1 = morphism(XxY, X, gens(RS)[1:m], check=false)
  pr2 = morphism(XxY, Y, gens(RS)[m+1:m+n], check=false)
  return XxY, pr1, pr2
end
=#




########################################
# (4) Equality
########################################

function ==(f::AffineSchemeMorType, g::AffineSchemeMorType) where {AffineSchemeMorType<:AbsAffineSchemeMor}
  X = domain(f)
  X == domain(g) || return false
  codomain(f) == codomain(g) || return false
  return all(OO(X)(x) == OO(X)(y) for (x, y) in zip(_images_of_generators(pullback(f)), _images_of_generators(pullback(g))))
end



########################################
# (5) Display
########################################

# Since the morphism is given in terms of pullback on the local coordinates,
# we need to adapt the printing to have everything aligned.
function Base.show(io::IO, ::MIME"text/plain", f::AbsAffineSchemeMor)
  io = pretty(io)
  X = domain(f)
  cX = coordinates(X)
  Y = codomain(f)
  cY = coordinates(Y)
  co_str = String[]
  str = "["*join(cX, ", ")*"]"
  kX = length(str)                                  # Length coordinates domain
  push!(co_str, str)
  str = "["*join(cY, ", ")*"]"
  kY = length(str)                                  # Length coordinates codomain
  push!(co_str, str)
  k = max(length.(co_str)...)                       # Maximum to estimate offsets
  println(io, "Affine scheme morphism")
  print(io, Indent(), "from ")
  print(io, co_str[1]*" "^(k-kX+2))                 # Consider offset for alignment
  println(IOContext(io, :show_coordinates => false), Lowercase(), X)
  print(io, "to   ")
  print(io, co_str[2]*" "^(k-kY+2))                 # Consider offset for alignment
  print(IOContext(io, :show_coordinates => false), Lowercase(), Y)
  x = coordinates(codomain(f))
  # If there are no coordinates, we do not print anything (since the target is
  # empty then)
  if length(x) > 0
    println(io)
    print(io, Dedent(), "given by the pullback function")
    pf = pullback(f)
    print(io, Indent())
    for i in 1:length(x)
      println(io)
      print(io, "$(x[i]) -> $(pf(x[i]))")
    end
  end
  print(io, Dedent())
end

function Base.show(io::IO, f::AbsAffineSchemeMor)
  if is_terse(io)
    print(io, "Affine scheme morphism")
  else
    io = pretty(io)
    print(io, "Hom: ")
    print(io, Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
  end
end

########################################################################
# (6) Base change
########################################################################

@doc raw"""
    base_change(phi::Any, f::AbsAffineSchemeMor)
        domain_map::AbsAffineSchemeMor=base_change(phi, domain(f))[2],
        codomain_map::AbsAffineSchemeMor=base_change(phi, codomain(f))[2]
      )

For a morphism ``f : X → Y`` between two schemes over a `base_ring` ``𝕜``
and a ring homomorphism ``φ : 𝕜 → 𝕂`` this returns a triple
`(b₁, F, b₂)` consisting of the maps in the commutative diagram
```
             f
    X        →    Y
    ↑ b₁          ↑ b₂
  X×ₖSpec(𝕂) → Y×ₖSpec(𝕂)
             F
```
The optional arguments `domain_map` and `codomain_map` can be used
to specify the morphisms `b₁` and `b₂`, respectively.
"""
function base_change(phi::Any, f::AbsAffineSchemeMor;
    domain_map::AbsAffineSchemeMor=base_change(phi, domain(f))[2],
    codomain_map::AbsAffineSchemeMor=base_change(phi, codomain(f))[2]
  )
  X = domain(f)
  Y = codomain(f)
  XX = domain(domain_map)
  YY = domain(codomain_map)
  pbf = pullback(f)
  pb1 = pullback(domain_map)
  pb2 = pullback(codomain_map)

  R = OO(Y)
  S = OO(X)
  RR = OO(YY)
  SS = OO(XX)

  img_gens = [pb1(pbf(x)) for x in gens(R)]
  # For the pullback of F no explicit coeff_map is necessary anymore
  # since both rings in domain and codomain have the same (extended/reduced)
  # coefficient ring by now.
  pbF = hom(RR, SS, img_gens, check=false)
  return domain_map, morphism(XX, YY, pbF, check=false), codomain_map
end

function base_change(phi::Any, f::ClosedEmbedding;
    domain_map::AbsAffineSchemeMor=base_change(phi, domain(f))[2],
    codomain_map::AbsAffineSchemeMor=base_change(phi, codomain(f))[2]
  )
  @assert codomain(codomain_map) === codomain(f)
  @assert codomain(domain_map) === domain(f)
  g = underlying_morphism(f)
  _, bc_g, _ = base_change(phi, g; domain_map, codomain_map)
  I = image_ideal(f)
  @assert base_ring(I) === OO(codomain(f))
  @assert base_ring(I) === OO(codomain(f))
  #bc_I = ideal(OO(codomain(bc_g)), pullback(codomain_map).(gens(I)))
  bc_I = pullback(codomain_map)(I)
  @assert domain(bc_g) === domain(domain_map)
  @assert codomain(bc_g) === domain(codomain_map)
  return domain_map, ClosedEmbedding(bc_g, bc_I; check=false), codomain_map
end

function _register_birationality!(f::AbsAffineSchemeMor,
    g::AbsAffineSchemeMor, ginv::AbsAffineSchemeMor)
  set_attribute!(g, :inverse, ginv)
  set_attribute!(ginv, :inverse, g)
  return _register_birationality(f, g)
end

function _register_birationality!(f::AbsAffineSchemeMor,
    g::AbsAffineSchemeMor
  )
  set_attribute!(f, :is_birational, true)
  set_attribute!(f, :iso_on_open_subset, g)
end

