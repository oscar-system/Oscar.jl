export pushforward
# This struct associates to a morphism f of type `MorphismType` 
# the mathematical symbol f^* which, in one way or the other, 
# denotes some pullback functor. 
#
# Which functor exactly this should be, is determined only from 
# the context in which this symbol appears. Hence, it is no 
# rigorous mathematical object, but really only a symbol. 
#
# It is introduced to facilitate writing 
#
#   pbf = pullback(f) # f some morphism f : X → Y
#   N = pbf(M) # M some object on Y
# 
# rather than 
# 
#   N = pullback(f, M)
#
# In particular, this makes the pullback a univariate 
# function which might come in handy when passing it as an argument 
# to other functions. For instance, for a collection of objects 
# A = [M₁, M₂, … ] on Y you can now write 
#
#   B = pbf.(A)
#
# rather than
#
#   B = (a->pullback(f, a).(A)
#
# or
#
#   B = [pullback(f, a) for a in A]
struct UniversalPullbackSymbol{MorphismType}
  f::MorphismType
end

# The following refers a call of f_star(M) to a bivariate pullback-method.
# The latter is what really needs to be implemented by the programmer. 
# Also, all parent checks, etc. are expected to happen there.
function (f_star::UniversalPullbackSymbol)(M::Any)
  return pullback(f_star.f, M)
end

# We can do the same thing with pushforward and other symbols.
struct UniversalPushforwardSymbol{MorphismType}
  f::MorphismType
end

# The following refers a call of f_star(M) to a bivariate pushforward-method.
# The latter is what really needs to be implemented by the programmer. 
# Also, all parent checks, etc. are expected to happen there.
function (f_lower_star::UniversalPushforwardSymbol)(M::Any)
  return pushforward(f_lower_star.f, M)
end

@doc raw"""
    pullback(f::CoveredSchemeMorphism)

Return a function `phi` which takes any reasonable 
argument `M` associated to the `codomain` of `f` and 
produces ``f*(M)``.

**Note:** Internally, this simply calls `pullback(f, M)`. 
Hence, that method needs to be implemented.
"""
function pullback(f::AbsCoveredSchemeMorphism)
  return UniversalPullbackSymbol{typeof(f)}(f)
end

@doc raw"""
    pushforward(f::CoveredSchemeMorphism)

Return a function `phi` which takes any reasonable 
argument `M` associated to the `domain` of `f` and 
produces ``f_*(M)``.

**Note:** Internally, this simply calls `pushforward(f, M)`. 
Hence, that method needs to be implemented.
"""
function pushforward(f::AbsCoveredSchemeMorphism)
  return UniversalPushforwardSymbol{typeof(f)}(f)
end

function pullback(f::AbsCoveredSchemeMorphism, II::AbsIdealSheaf)
  X = domain(f)
  Y = codomain(f)
  scheme(II) === Y || error("ideal sheaf is not defined on the codomain of the function")
  return PullbackIdealSheaf(f, II)
end

@doc raw"""
    total_transform(f::AbsSimpleBlowupMorphism, II::IdealSheaf)

Compute the total transform of an ideal sheaf along a blowup.

In particular, this applies in the toric setting. However, note that
currently (October 2023), ideal sheaves are only supported on smooth
toric varieties.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> bl = blow_up(P2, [1, 1])
Toric blowup morphism

julia> S = cox_ring(P2);

julia> x, y, z = gens(S);

julia> I = ideal_sheaf(P2, ideal([x*y]))
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_1_1*x_2_1)
  2: Ideal (x_2_2)
  3: Ideal (x_1_3)

julia> total_transform(bl, I)
Sheaf of ideals
  on normal toric variety
with restrictions
  1: Ideal (x_1_1*x_2_1^2)
  2: Ideal (x_1_2^2*x_2_2)
  3: Ideal (x_2_3)
  4: Ideal (x_1_4)
```
"""
function total_transform(f::AbsSimpleBlowupMorphism, II::AbsIdealSheaf)
  return pullback(f, II)
end

function total_transform(f::AbsBlowupMorphism, II::AbsIdealSheaf)
  return pullback(f, II)
end

function pullback(f::CompositeCoveredSchemeMorphism, C::EffectiveCartierDivisor)
  result = C
  for g in reverse(maps(f))
    result = pullback(g, result)
  end
  return result
end

function pullback(f::AbsCoveredSchemeMorphism, C::EffectiveCartierDivisor)
  X = domain(f)
  Y = codomain(f)
  Y === ambient_scheme(C) || error("divisor must be defined on the codomain of the map")
  # The challenge is that phi has two coverings cov1 → cov2 on which it is defined. 
  # The covering cov3 on which C is principalized might be different from cov2. 
  # Thus, we need to first pass to a common refinement cov' of cov2 and cov3, 
  # restrict f to that, obtain a new underlying covering morphism psi cov1' → cov', 
  # and pull back along psi.
  phi = covering_morphism(f)
  triv_dict = IdDict{AbsAffineScheme, RingElem}()
  E, a, b = common_refinement(codomain(phi), trivializing_covering(C))
  psi = restrict(f, E)
  psi = compose(psi, b)
  triv_cov = domain(psi)
  for U in patches(triv_cov)
    V = codomain(psi[U])
    pbgens = pullback(psi[U]).(C(V))

    # Do the sanity checks
    length(pbgens) == 1 || error("cartier divisor is not principal on this patch")
    g = first(pbgens)
    # See the Stacks project on cartier divisors
    is_zero_divisor(g) && error("pullback of the local equation is a zero divisor; pullback of the cartier divisor is not defined")
    triv_dict[U] = g
  end

  return EffectiveCartierDivisor(X, triv_dict, trivializing_covering=triv_cov, check=false)
end

function pullback(f::AbsCoveredSchemeMorphism, C::CartierDivisor)
  R = coefficient_ring(C)
  result = CartierDivisor(domain(f), R)
  for (D, c) in coefficient_dict(C)
    result += c*pullback(f, D)
  end
  return result
end

function pullback(f::AbsCoveredSchemeMorphism, CC::Covering)
  psi = restrict(f, CC)
  return domain(psi)
end

function pullback(f::AbsCoveredSchemeMorphism, M::AbsCoherentSheaf)
  return PullbackSheaf(f, M)
end

@doc raw"""
    restrict(f::CoveredSchemeMorphism, DD::Covering)

Given a morphism of `CoveredScheme`s ``f : X → Y``, described via ``φ : C → D`` 
on `Covering`s `C` of ``X`` and `D` of ``Y``, and a refinement `DD ≤ D`, 
compute the preimages ``Uᵢ ∩ f⁻¹(Vⱼ)`` for ``Uᵢ`` the `patches` in `C` and 
``Vⱼ`` those of `DD` and restrictions of ``f`` to these new patches. 
Return the resulting `CoveringMorphism` ``ψ : C ∩ f⁻¹(DD) → DD``.
"""
function restrict(f::AbsCoveredSchemeMorphism, DD::Covering)
  X = domain(f)
  Y = codomain(f)
  phi = covering_morphism(f)
  C = domain(phi)
  D = codomain(phi)
  
  # Check the cache
  res_cache = restriction_cache(f)
  haskey(res_cache, DD) && return res_cache[DD]

  # First check the trivial case that we can just compose.
  if haskey(refinements(Y), (D, DD))
    return compose(phi, refinements(Y)[(D, DD)])
  end

  # Then check whether we can nevertheless build the composition manually
  success, psi = is_refinement(D, DD)
  success && return compose(phi, psi)

  # DD needs to be a refinement of D; otherwise quit.
  all(x->has_ancestor(y->any(z->(z===y), patches(D)), x), patches(DD)) || error("second argument needs to be a refinement of the codomain covering on which the first argument is defined")

  OOX = OO(X)
  OOY = OO(Y)
  # We need to do the following:
  # The covering DD has patches Vⱼ in the codomain Y of f.
  # Their preimages must be taken in every patch Uᵢ of X in 
  # the domain's covering for phi. 
  res_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(domain(phi))
    V = codomain(phi[U])
    for W in patches(DD)
      # In case V and W do not happen to lay in the same affine 
      # chart, the preimage of W in U will be empty. So we can 
      # skip the computation in that case.
      success, par = _have_common_ancestor(V, W)
      if success
        # Reconstruct W as a PrincipalOpenSubset of the affine chart directly
        # (Note that we might have gone through some calls to `simplify(...)`, 
        # so this really is a non-trivial step)
        iso_W_flat = _flatten_open_subscheme(W, par)
        W_flat = codomain(iso_W_flat)
        # To cheaply construct the preimage of W in U, just pull back the 
        # complement equation.
        h = complement_equation(codomain(iso_W_flat))
        UW = PrincipalOpenSubset(U, pullback(phi[U])(OOY(par, V)(h)))
        # Manually assemble the restriction of phi to this new patch
        ff = morphism(UW, W_flat, 
                     hom(OO(W_flat), OO(UW), OO(UW).(pullback(phi[U]).(OOY(par, V).(gens(OO(par))))), check=false),
                     check=false
                    )
        f_res = compose(ff, inverse(iso_W_flat))
        res_dict[UW] = f_res
      end
    end
  end
  new_domain = Covering(collect(keys(res_dict)), IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}())
  inherit_gluings!(new_domain, domain(phi))
  _register!(new_domain, X)

  psi = CoveringMorphism(new_domain, DD, res_dict, check=false)
  restriction_cache(f)[DD] = psi
  return psi
end

function _register!(C::Covering, X::AbsCoveredScheme)
  push!(coverings(X), C)
  refinements(X)[(C, default_covering(X))] = _canonical_map(C, default_covering(X))
  return C
end

function _canonical_map(C::Covering, D::Covering)
  map_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(C)
    f, _ = _find_chart(U, D)
    map_dict[U] = f
  end
  phi = CoveringMorphism(C, D, map_dict, check=false)
  return phi
end

# Several objects on the codomain of f : X → Y might use the same covering 
# in their internals. Now when pulling them back, we do not want to recreate 
# the necessary refinement on X again and again. In particular, because 
# we do not want the pulled back modules/ideals/functions... to live on different 
# affine patches in X while their originals were defined on the same patches in Y.
@attr IdDict{<:Covering, <:CoveringMorphism} function restriction_cache(f::AbsCoveredSchemeMorphism)
  res_dict = IdDict{Covering, CoveringMorphism}()
  phi = covering_morphism(f)
  res_dict[codomain(phi)] = phi
  return res_dict
end

struct InheritGluingData
  orig::Covering
  X::AbsAffineScheme
  Y::AbsAffineScheme
end

function _compute_gluing(gd::InheritGluingData)
  X = gd.X
  Y = gd.Y
  C = gd.orig

  success, Z = _have_common_ancestor(X, Y)
  if success
    # This is the easy case: Gluing within one chart. 
    # Keep in mind, however, that we might have gone through some simplify(...) calls, 
    # so we can not assume everything to be happening in the same ambient_ring.
    iso_X = _flatten_open_subscheme(X, Z)
    iso_Y = _flatten_open_subscheme(Y, Z)
    h_X = complement_equation(codomain(iso_X))
    h_Y = complement_equation(codomain(iso_Y))
    XY = PrincipalOpenSubset(X, pullback(iso_X)(OO(codomain(iso_X))(h_Y)))
    YX = PrincipalOpenSubset(Y, pullback(iso_Y)(OO(codomain(iso_Y))(h_X)))
    if iszero(h_X*h_Y)
      # Gluing along the empty set. This is trivial.
      g = morphism(YX, XY, hom(OO(XY), OO(YX), [zero(OO(YX)) for i in 1:ngens(OO(XY))], check=false), check=false)
      f = morphism(XY, YX, hom(OO(YX), OO(XY), [zero(OO(XY)) for i in 1:ngens(OO(YX))], check=false), check=false)
      
      return SimpleGluing(X, Y, f, g, check=false)
    end

    XYZ = PrincipalOpenSubset(Z, h_X*h_Y)

    x_img = gens(OO(X))
    x_img = pullback(inverse(iso_X)).(x_img)
    x_img = OO(XYZ).(x_img)
    phi = restrict(iso_Y, YX, XYZ, check=false)
    x_img = pullback(phi).(x_img)
    g = morphism(YX, XY, hom(OO(XY), OO(YX), x_img, check=false), check=false)
    
    y_img = gens(OO(Y))
    y_img = pullback(inverse(iso_Y)).(y_img)
    y_img = OO(XYZ).(y_img)
    psi = restrict(iso_X, XY, XYZ, check=false)
    y_img = pullback(psi).(y_img)
    f = morphism(XY, YX, hom(OO(YX), OO(XY), y_img, check=false), check=false)
    return SimpleGluing(X, Y, f, g, check=false)
  end

  # As the easy case would have been caught before, we are now facing an inherited 
  # gluing across charts: X ↪ A ⊃ U ≅ V ⊂ B ↩ Y. 
  # We need to compute the intersection of X with Y along the identifications of U and V
  # and cook up the SimpleGluing from that.
  iso_X = _flatten_open_subscheme(X, C)
  iso_Y = _flatten_open_subscheme(Y, C)
  A = ambient_scheme(codomain(iso_X))
  B = ambient_scheme(codomain(iso_Y))
  G = C[A, B] # The original gluing needed
  U, V = gluing_domains(G)
  f, g = gluing_morphisms(G)
  U isa PrincipalOpenSubset && ambient_scheme(U) === A || error("incorrect intermediate result")
  V isa PrincipalOpenSubset && ambient_scheme(V) === B || error("incorrect intermediate result")

  UX = intersect(U, codomain(iso_X))
  UX isa PrincipalOpenSubset && ambient_scheme(UX) === A || error("incorrect intermediate result")
  h_Y = pullback(f)(complement_equation(codomain(iso_Y)), check=false)
  UXY = PrincipalOpenSubset(codomain(iso_X), 
                            OO(codomain(iso_X))(lifted_numerator(h_Y)*complement_equation(U), check=false))
  UXY isa PrincipalOpenSubset && ambient_scheme(UXY) === codomain(iso_X) || error("incorrect intermediate output")
  UY = PrincipalOpenSubset(U, h_Y)
  XY = PrincipalOpenSubset(X, pullback(iso_X)(complement_equation(UXY)))

  VY = intersect(V, codomain(iso_Y))
  VY isa PrincipalOpenSubset && ambient_scheme(VY) === B || error("incorrect intermediate result")
  h_X = pullback(g)(complement_equation(codomain(iso_X)), check=false)
  VX = PrincipalOpenSubset(V, h_X)
  VYX = PrincipalOpenSubset(codomain(iso_Y), 
                            OO(codomain(iso_Y))(lifted_numerator(h_X)*complement_equation(V), check=false))
  VYX isa PrincipalOpenSubset && ambient_scheme(VYX) === codomain(iso_Y) || error("incorrect intermediate output")
  YX = PrincipalOpenSubset(Y, pullback(iso_Y)(complement_equation(VYX)))

  fres = restrict(f, UXY, VYX, check=false)
  gres = restrict(g, VYX, UXY, check=false)

  x_img = gens(OO(X))
  x_img = pullback(inverse(iso_X)).(x_img)
  x_img = [OO(UXY)(x, check=false) for x in x_img]
  x_img = pullback(gres).(x_img)
  phi = restrict(iso_Y, YX, VYX, check=false)
  x_img = pullback(phi).(x_img)
  gg = morphism(YX, XY, hom(OO(XY), OO(YX), x_img, check=false), check=false)

  y_img = gens(OO(Y))
  y_img = pullback(inverse(iso_Y)).(y_img)
  y_img = [OO(VYX)(y, check=false) for y in y_img]
  y_img = pullback(fres).(y_img)
  psi = restrict(iso_X, XY, UXY, check=false)
  y_img = pullback(psi).(y_img)
  ff = morphism(XY, YX, hom(OO(YX), OO(XY), y_img, check=false), check=false)

  return SimpleGluing(X, Y, ff, gg, check=false)
end

number_of_generators(Q::MPolyQuoLocRing) = number_of_generators(base_ring(Q))

@doc raw"""
    inherit_gluings!(ref::Covering, orig::Covering)

For a refinement `ref` of a `Covering` `orig` add all missing gluings 
in `ref` as inherited gluings from those in `orig`.
"""
function inherit_gluings!(ref::Covering, orig::Covering)
  for U in patches(ref)
    for V in patches(ref)
      if !haskey(gluings(ref), (U, V))
        gluings(ref)[(U, V)] = LazyGluing(U, V, 
                                            InheritGluingData(orig, U, V)
                                           )
      end
    end
  end
  return ref
end

### Generic pullback and pushforward for composite maps
pushforward(f::Generic.CompositeMap, a::Any) = pushforward(map2(f), pushforward(map1(f), a))
pullback(f::Generic.CompositeMap, a::Any) = pullback(map1(f), pullback(map2(f), a))


### Strands of graded modules
function _coordinates_in_monomial_basis(v::T, b::Vector{T}) where {B <: MPolyRingElem, T <: FreeModElem{B}}
  F = parent(v)
  R = base_ring(F)
  kk = coefficient_ring(R)
  result = SRow(kk)
  iszero(v) && return result
  pos = [findfirst(==(m), b) for m in monomials(v)]
  vals = collect(coefficients(v))
  return SRow(kk, pos, vals)
end

### Take a complex of graded modules and twist all 
# modules so that the (co-)boundary maps become homogeneous
# of degree zero. 
function _make_homogeneous(C::ComplexOfMorphisms{T}) where {T<:SparseFPModule}
  R = base_ring(C[first(range(C))])
  is_standard_graded(R) || error("ring must be standard graded")
  all(k->base_ring(C[k])===R, range(C)) || error("terms in complex must have the same base ring")
  all(k->is_graded(C[k]), range(C)) || error("complex must be graded")
  all(k->C[k] isa FreeMod, range(C)) || error("terms in complex must be free")

  new_maps = Map[]
  new_chains = [C[first(range(C))]]
  offset = zero(grading_group(R))
  for i in range(C)
    i == last(range(C)) && break
    phi = map(C, i)
    dom = domain(phi)
    cod = codomain(phi)
    x = gens(dom)
    if !isempty(x)
      delta = degree(phi(first(x))) - degree(first(x))
      @assert all(u->degree(phi(u)) - degree(u) == delta, x[2:end]) "map is not homogeneous"
      offset = offset + delta
    end
    new_chain = twist(cod, offset)
    imgs = phi.(gens(domain(phi)))
    imgs = [new_chain(coordinates(y)) for y in imgs]
    new_map = hom(last(new_chains), new_chain, imgs)
    push!(new_maps, new_map)
    push!(new_chains, new_chain)
  end
  
  return ComplexOfMorphisms(T, new_maps, typ = :chain, seed=last(range(C)))
end

### Take a complex of free ℤ-graded modules with (co-)boundary 
# maps of degree zero and return the complex given by all the 
# degree d parts.
function strand(C::ComplexOfMorphisms{T}, d::Int) where {T<:SparseFPModule}
  R = base_ring(C[first(range(C))])
  is_standard_graded(R) || error("ring must be standard graded")
  all(k->base_ring(C[k])===R, range(C)) || error("terms in complex must have the same base ring")
  all(k->is_graded(C[k]), range(C)) || error("complex must be graded")
  all(k->C[k] isa FreeMod, range(C)) || error("terms in complex must be free")

  kk = coefficient_ring(R)
  res_maps = Map[]

  i = first(range(C))
  mons = collect(all_monomials(C[i], d))
  chains = [FreeMod(kk, length(mons))]
  
  for i in range(C)
    i == last(range(C)) && break # How else to skip that???
    new_mons = collect(all_monomials(C[i-1], d))
    push!(chains, FreeMod(kk, length(new_mons)))
    imgs = elem_type(last(chains))[]
    phi = map(C, i)
    M = last(chains)
    for x in mons
      c = _coordinates_in_monomial_basis(phi(x), new_mons)
      v = sum(x*M[i] for (i, x) in c; init=zero(M))
      push!(imgs, v)
    end
    push!(res_maps, hom(chains[end-1], chains[end], imgs))
    mons = new_mons
  end
  return ComplexOfMorphisms(typeof(domain(first(res_maps))), res_maps, typ = :chain, seed=last(range(C)), check=false)
end

### Get a subcomplex in a specific range. 
function getindex(c::ComplexOfMorphisms{T}, r::UnitRange) where {T}
  first(r) in range(c) || error("range out of bounds")
  last(r) in range(c) || error("range out of bounds")
  if is_chain_complex(c)
    maps = [map(c, i) for i in last(r):-1:first(r)+1]
    return ComplexOfMorphisms(T, maps, seed=first(r), typ=:chain, check=false)
  else
    maps = [map(c, i) for i in first(r):last(r)-1]
    return ComplexOfMorphisms(T, maps, seed=first(r), typ=:cochain, check=false)
  end
end


