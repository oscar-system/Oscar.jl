export BlowupMorphism
export center
export exceptional_divisor
export projection

########################################################################
# An abstract type for blowdown morphisms.
#
# This should also comprise sequences of simple blowups leading
# to a partial or full resolution of singularities. The interface
# is specified below.
########################################################################
abstract type AbsBlowdownMorphism{DomainType<:AbsCoveredScheme,
                                  CodomainType<:AbsCoveredScheme,
                                  BaseMorphismType<:Nothing,
                                  BlowdownMorphismType
   } <: AbsCoveredSchemeMorphism{
                                 DomainType,
                                 CodomainType,
                                 BaseMorphismType,
                                 BlowdownMorphismType
                                }
end

# The interface inherits all functionality from AbsCoveredSchemeMorphism.
# This is extended by the following:
@doc raw"""
    exceptional_divisor(f::AbsBlowdownMorphism)

Return a `CartierDivisor` on the `domain` of `f` which is the
preimage of the vanishing locus of the `center` of `f`.

!!! note If `f` is not a classical blowup, but, for instance, a small contraction, then the exceptional locus need not be a divisor. See `exceptional_locus` for such cases.
"""
function exceptional_divisor(f::AbsBlowdownMorphism)
  error("not implemented")
end

@doc raw"""
    exceptional_locus(f::AbsBlowdownMorphism)

Return an `AbsAlgebraicCycle` on the `domain` of `f` which is the preimage
of the `center` of `f` in a reasonable sense. In particular, `f` is an isomorphism
onto its image outside the support of its `exceptional_locus`. See also `exceptional_divisor`.
"""
function exceptional_locus(f::AbsBlowdownMorphism)
  error("not implemented")
end

@doc raw"""
    center(f::AbsBlowdownMorphism)

Return an `IdealSheaf` on the `codomain` of `f` such that on the complement
of the vanishing locus of that ideal sheaf `f` is an isomorphism.
The support of the `exceptional_locus` of `f` coincides with the vanishing
locus of the pullback of the `center`.

!!! note For classical blowup constructions this will be the center which has been blown up. For more exotic blowups, however, this might be more special.
"""
function center(f::AbsBlowdownMorphism)
  error("not implemented")
end

@doc raw"""
    strict_transform(f::AbsBlowdownMorphism, a::Any)

Take the closure of the `pullback` of `a` along `f` on the complement
of the `exceptional_locus` in a mathematically reasonable sense.
Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function strict_transform(f::AbsBlowdownMorphism, a::Any)
  error("not implemented")
end

@doc raw"""
    total_transform(f::AbsBlowdownMorphism, a::Any)

Take the `pullback` or preimage of `a` along `f` in a mathematically
reasonable sense.
Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function total_transform(f::AbsBlowdownMorphism, a::Any)
  error("not implemented")
end

########################################################################
# An abstract type for classical blowups of ideal sheaves.
#
# This can either be a BlowupMorphism as below, but also a special
# toric morphism induced by fan subdivisions.
########################################################################
abstract type AbsSimpleBlowdownMorphism{DomainType<:AbsCoveredScheme,
                                     CodomainType<:AbsCoveredScheme,
                                     BaseMorphismType<:Nothing,
                                     BlowdownMorphismType
    } <: AbsBlowdownMorphism{
                             DomainType,
                             CodomainType,
                             BaseMorphismType,
                             BlowdownMorphismType
                            }
end

@doc raw"""
    exceptional_divisor(f::AbsSimpleBlowdownMorphism)

Return a `CartierDivisor` on the `domain` of `f` which coincides
with the pullback of the `center` of `f`.
"""
function exceptional_divisor(f::AbsSimpleBlowdownMorphism)
  error("not implemented")
end

@doc raw"""
    exceptional_locus(f::AbsBlowdownMorphism)

Return the `WeilDivisor` of the `exceptional_divisor`.
"""
function exceptional_locus(f::AbsSimpleBlowdownMorphism)
  error("not implemented")
end

@doc raw"""
    center(f::AbsSimpleBlowdownMorphism)

Return an `IdealSheaf` on the `codomain` of `f` the blowup of which
leads to `f`.
"""
function center(f::AbsSimpleBlowdownMorphism)
  error("not implemented")
end

@doc raw"""
    strict_transform(f::AbsSimpleBlowdownMorphism, a::Any)

Take the closure of the `pullback` of `a` along `f` on the complement
of the `exceptional_divisor` in a mathematically reasonable sense.

Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function strict_transform(f::AbsSimpleBlowdownMorphism, a::Any)
  error("not implemented")
end

@doc raw"""
    total_transform(f::AbsSimpleBlowdownMorphism, a::Any)

Take the `pullback` or preimage of `a` along `f` in a mathematically
reasonable sense.

Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function total_transform(f::AbsSimpleBlowdownMorphism, a::Any)
  error("not implemented")
end

@doc raw"""
    controlled_transform(f::AbsSimpleBlowdownMorphism, a::Any, k::Int)

Take the `pullback` or preimage of `a` along `f` saturated by `k` times
the `exceptional_divisor` of `f` in a mathematically reasonable sense.

Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function controlled_transform(f::AbsSimpleBlowdownMorphism, a::Any, k::Int)
  error("not implemented")
end

########################################################################
# BlowupMorphism
#
# A datastructure to maintain all information necessary to effectively
# handle blowups. This is work in progress and will one day serve as
# a building block for sequences of blowups
########################################################################

@doc raw"""
    BlowupMorphism

A datastructure to encode blowups of covered schemes in some sheaves of ideals.

It is described as a morphism from the new scheme to the blown-up scheme, with
information about its center (i.e. the ideal sheaves blown-up in the bottom
scheme) and its exceptional locus (i.e. the preimage of the center under the
blowup).

# Examples
```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> A3 = spec(R)
Spectrum
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field

julia> I = ideal(R, [x,y,z])
Ideal generated by
  x
  y
  z

julia> bl = blow_up(A3, I)
Blowup
  of scheme over QQ covered with 1 patch
    1b: [x, y, z]   affine 3-space
  in sheaf of ideals with restriction
    1b: Ideal (x, y, z)
with domain
  scheme over QQ covered with 3 patches
    1a: [(s1//s0), (s2//s0), x]   scheme(0, 0, 0)
    2a: [(s0//s1), (s2//s1), y]   scheme(0, 0, 0)
    3a: [(s0//s2), (s1//s2), z]   scheme(0, 0, 0)
and exceptional divisor
  effective cartier divisor defined by
    sheaf of ideals with restrictions
      1a: ideal(x)
      2a: ideal(y)
      3a: ideal(z)

julia> E = exceptional_divisor(bl)
Effective cartier divisor
  on scheme over QQ covered with 3 patches
    1: [(s1//s0), (s2//s0), x]   scheme(0, 0, 0)
    2: [(s0//s1), (s2//s1), y]   scheme(0, 0, 0)
    3: [(s0//s2), (s1//s2), z]   scheme(0, 0, 0)
defined by
  sheaf of ideals with restrictions
    1: ideal(x)
    2: ideal(y)
    3: ideal(z)

julia> Z = center(bl)
Sheaf of ideals
  on scheme over QQ covered with 1 patch
    1: [x, y, z]   affine 3-space
with restriction
  1: Ideal (x, y, z)
```
"""
@attributes mutable struct BlowupMorphism{
     DomainType<:AbsCoveredScheme, # Not a concrete type in general because this is lazy
     CodomainType<:AbsCoveredScheme,
     BaseMorphismType # Nothing in case of no base change
   } <: AbsSimpleBlowdownMorphism{
                                  DomainType,
                                  CodomainType,
                                  BaseMorphismType,
                                  BlowupMorphism
                                 }
  projective_bundle::CoveredProjectiveScheme
  codomain::CodomainType   # in general a CoveredScheme
  center::IdealSheaf      # on codomain
  projection::AbsCoveredSchemeMorphism
  domain::AbsCoveredScheme # in general a CoveredScheme
  exceptional_divisor::EffectiveCartierDivisor

  function BlowupMorphism(
      IP::CoveredProjectiveScheme,
      I::IdealSheaf
    )
    X = base_scheme(IP)
    X === scheme(I) || error("ideal sheaf not compatible with blown up variety")
    return new{AbsCoveredScheme, typeof(X), Nothing}(IP, X, I)
  end
end

### Forward the essential functionality
underlying_morphism(phi::BlowupMorphism) = projection(phi)

function domain(p::BlowupMorphism)
  if !isdefined(p, :domain)
    p.domain = covered_scheme(p.projective_bundle)
    simplify!(p.domain)        # if simplify hangs, no other computation would go through anyway
  end
  return p.domain
end

codomain(p::BlowupMorphism) = p.codomain
center(p::BlowupMorphism) = p.center

function projection(p::BlowupMorphism)
  if !isdefined(p, :projection)
    p.projection = covered_projection_to_base(p.projective_bundle)
  end
  return p.projection
end

# TODO: Find better name!
covered_projective_scheme(p::BlowupMorphism) = p.projective_bundle

@doc raw"""
    exceptional_divisor(p::BlowupMorphism)

For a `BlowupMorphism` ``p : Y → X`` coming from the blowup of an
`IdealSheaf` ``ℐ`` on X, return the `EffectiveCartierDivisor` ``E``
on ``Y`` associated to the (relative) tautological bundle ``𝒪(1)``.

On a pair of charts ``V → U`` of the `covered_scheme` of the
`projection` of ``p`` this returns the pullback of the `i`-th
generator of ``ℐ(U)`` when ``V`` is the `i-1`-st canonical chart
of the local blowup over ``U``.
"""
function exceptional_divisor(p::BlowupMorphism)
  if !isdefined(p, :exceptional_divisor)
    error("exceptional divisor needs to be cached during construction")
    # The exceptional divisor must be created
    # and set during the construction of the BlowupMorphism.
  end
  return p.exceptional_divisor
end

@doc raw"""
    strict_transform(p::BlowupMorphism, inc::CoveredClosedEmbedding)

For a `BlowupMorphism` ``p : Y → X`` and a `CoveredClosedEmbedding`
``ι : Z ↪ X``, compute the strict transform ``Z'`` of ``Z`` along ``p`` and
return a triple ``(Z', j, π)`` containing the `CoveredClosedEmbedding`
``j : Z' ↪ Y`` and the induced projection ``π : Z' → Z``.
"""
function strict_transform(p::AbsSimpleBlowdownMorphism, inc::CoveredClosedEmbedding)
  Y = domain(p)
  X = codomain(p)
  Z = domain(inc)
  codomain(inc) === X || error("maps must have the same codomain")
  I_trans = strict_transform(p, image_ideal(inc))
  inc_Z_trans = CoveredClosedEmbedding(Y, I_trans,
                                       covering=simplified_covering(Y), # Has been set by the previous call
                                       check=false)
  inc_dom_cov = covering_morphism(inc_Z_trans)
  inc_cod_cov = covering_morphism(inc)

  Z_trans = domain(inc_Z_trans)
  pr_res = restrict(projection(p), inc_Z_trans, inc)

  if has_attribute(p, :isomorphism_on_open_subset) # will only happen when p is a BlowupMorphism
    OOX = OO(X)
    OOY = OO(Y)
    p_iso = isomorphism_on_open_subset(p)
    U = domain(p_iso)::PrincipalOpenSubset
    V = codomain(p_iso)::PrincipalOpenSubset
    V_amb = ambient_scheme(V)
    U_amb = ambient_scheme(U)

    # We have the following diagram:
    #             inc_dom
    #     Z_trans    ↪     Y   ⊃ U_amb ⊃ U
    #
    #  pr_res↓             ↓ p           ↓ p_iso
    #        Z       ↪     X   ⊃ V_amb ⊃ V
    #             inc_cod
    #
    # Given all the refinements that potentially had to be done, we have
    # to do the following.
    #  1. Find a `patch` `U_sub` of the `domain` of `inc_dom_cov` for which
    #    inc_dom : U_sub -> W has a codomain `W` which has `U_amb` as an
    #    ancestor. Since `U_amb` is one of the `affine_charts` of `Y`, this will work.

    k = findfirst(x->has_ancestor(y->y===U_amb, codomain(inc_dom_cov[x])), patches(domain(inc_dom_cov)))
    U_sub = patches(domain(inc_dom_cov))[k]

    #  2. pr_res : U_amb -> V_amb has some codomain such that there exists
    #    an ancestor `VV` in the `domain` of `inc_dom_cov` such that
    #    inc_dom : VV -> W' has a codomain with `V_amb` as an ancestor.
    #  3. outside the `complement_equation`s of `U` and `V` the projection
    #    was an isomorphism. Then the restriction of `pr_res` to those
    #    complements in `U_sub` and `V_sub` also is.
    U_sub_res = PrincipalOpenSubset(U_sub,
                   pullback(inc_dom_cov[U_sub])(
                       OOX(U_amb, codomain(inc_dom_cov[U_sub]))(
                           complement_equation(U)
                         )
                     )
                 )

    pr_res_cov = covering_morphism(pr_res)
    pr_sub = pr_res_cov[U_sub]

    V_sub = codomain(pr_sub)
    @assert has_ancestor(x->any(y->y===x, patches(domain(inc_cod_cov))), V_sub)
    V_sub_inc, comp_eqns = _find_chart(V_sub, domain(inc_cod_cov))
    OOZ = OO(Z)
    dummy_dom = codomain(V_sub_inc)
    dummy_map = inc_cod_cov[dummy_dom]
    dummy_cod = codomain(dummy_map)
    V_sub_res = PrincipalOpenSubset(V_sub,
                                    OOZ(dummy_dom, V_sub)(pullback(dummy_map)(OOY(V_amb, dummy_cod)(complement_equation(V)))))
    result = restrict(pr_sub, U_sub_res, V_sub_res)
    # TODO: Obtain the inverse another way?
    @assert is_isomorphism(result)
    set_attribute!(pr_res, :isomorphism_on_open_subset, result)
  end
  return Z_trans, inc_Z_trans, pr_res
end

@doc """
    strict_transform(p::BlowupMorphism, I::IdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
strict transform of ``I`` on ``Y``.
"""
function strict_transform(p::AbsSimpleBlowdownMorphism, I::IdealSheaf)
  Istrict,_ =_do_transform(p, I, -1)
  return Istrict
end

@doc """
    weak_transform(p::BlowupMorphism, I::IdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
weak transform ``J`` of ``I`` on ``Y``, i.e. an `IdealSheaf` satisfying ``E^m J = p^*I`` with ``m``
maximal and ``E`` the 'IdealSheaf' of the exceptional divisor of ``p``.
"""
function weak_transform(p::AbsSimpleBlowdownMorphism, I::IdealSheaf)
  Iweak,_ =_do_transform(p,I,0)
  return Iweak
end

@doc """
    weak_transform_with_multiplicity((p::BlowupMorphism, I::IdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
weak transform ``J`` of ``I`` on ``Y`` and the multiplicity ``m`` of the exceptional divisor, i.e.
the maximal ``m`` such that ``E^m J = p^*I``, where ``E`` denotes the `IdealSheaf` of the exceptional
divisor of ``p``.
"""
function weak_transform_with_multiplicity(p::AbsSimpleBlowdownMorphism, I::IdealSheaf)
  Iweak, multi = _do_transform(p,I,0)
  return Iweak,multi
end

@doc """
    controlled_transform(p::BlowupMorphism, I::IdealSheaf, b::Int)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
controlled transform of ``I`` on ``Y`` with control ``b``,i.e. an `IdealSheaf` ``J`` such that
``E^b J = p^*I`` where ``E``denotes the `IdealSheaf` of the exceptional divisor.
"""
function controlled_transform(p::AbsSimpleBlowdownMorphism, I::IdealSheaf, b::Int)
  Icontrol,_ = _do_transform(p,I,b)
  return Icontrol
end

##########################################################################################################
## central internal method for strict, weak and controlled transforms of IdealSheafs and subschemes
##########################################################################################################
function _do_transform(p::AbsSimpleBlowdownMorphism, I::IdealSheaf, method::Int=-1)
## method: -1  strict transform
##          0  weak transform
##         b>0  controlled transform with control b>0
##         < -1 error
  method > -2  || error("unknown method of transform", method)

  ## initializations and sanity checks for p
  X = scheme(I)
  Y = domain(p)
  X === codomain(p) || error("ideal sheaf is not defined on the codomain of the morphism")
  IE = ideal_sheaf(exceptional_divisor(p))
  ID = IdDict{AbsAffineScheme, Ideal}()

  ## get our hands on the coverings -- using simplified covering for CY
  p_cov_temp = covering_morphism(p)
  CX = codomain(p_cov_temp)
  CY = domain(p_cov_temp)
  CY_simp = (CY === default_covering(Y) ? simplified_covering(Y) : CY)
  phi = (CY === default_covering(Y) ? Y[CY_simp,CY] : identity_map(CY_simp))
  p_cov = compose(phi,p_cov_temp)    # blow up using simplified covering

  ## do the transform on the charts
  b = -2                               # safe initialization of multiplicity return value
  bmin = -2                            # safe initialization of minimal multiplicity for weak transform in different charts
  for U in patches(CY_simp)
    V = codomain(p_cov[U])             # affine patch on X
    Iorig_chart = I(V)                 # I on this patch
    Itotal_chart = ideal(OO(U), pullback(p_cov[U]).(gens(Iorig_chart)))
                                       # total transform on Chart

    ## don't try saturating with respect to an empty set
    ## not expensive for Cartier divisor, GB is cached (?) after first computation
    if is_one(IE(U))
      ID[U] = Itotal_chart
      continue
    end

    IE_chart = IE(U)

    ## do different methods according to integer argument method
    if method == -1
      Itrans_chart,btemp = saturation_with_index(Itotal_chart, IE_chart)                      # strict
      b = max(b,btemp)
    elseif method == 0
      Itrans_chart,btemp = iterated_quotients(Itotal_chart,IE_chart, method)             # weak
      if b == -2
         b = btemp
      end
      bmin = min(b,btemp)
    else
      Itrans_chart,b = iterated_quotients(Itotal_chart,IE_chart, method)                 # controlled
    end
    ID[U] = Itrans_chart
  end

  bmin == -2 || bmin==b || error("different multiplicities in different charts, use controlled transform with control ",bmin)
  b > -2 || error("no patches in CY_simp!!!")
  I_trans = IdealSheaf(Y,ID,check=false)
  return I_trans,b
end

##########################################################################################################
## Handle Cartier divisors separately, as ideal quotients are quotients of ring elements in this case
##########################################################################################################
@doc """
  strict_transform(p::BlowupMorphism, C::EffectiveCartierDivisor)

For a `BlowupMorphism`  ``p : Y → X`` and an `EffectiveCartierDivisor` ``C`` on ``X`` return the
strict transform of ``C`` on ``Y``.
"""
function strict_transform(p::AbsSimpleBlowdownMorphism, C::EffectiveCartierDivisor)
  X = scheme(C)
  Y = domain(p)
  X === codomain(p) || error("cartier divisor is not defined on the codomain of the morphism")
  E = exceptional_divisor(p)
  ID = IdDict{AbsAffineScheme, RingElem}()

  ## get our hands on the coverings -- trivializing covering for C leading the way
  CX = trivializing_covering(C)
  pr_refined = restrict(p, CX)::CoveringMorphism
  CY = domain(pr_refined)

  ## do the transform on the charts
  multEInC = -1
  for U in patches(CY)
    V = codomain(pr_refined[U])        # affine patch on X

    ## determine single generator of Cartier divisor C on V
    length(C(V)) == 1 || error("ideal of divisor is not principal")
                       # sanity check -- we are on a trivializing covering after all!
    h_orig = C(V)[1]
    h_total = pullback(pr_refined[U]).(h_orig)
    if is_unit(h_total)
      ID[U] = one(OO(U))
      continue
    end

    ## determine single generator of Cartier divisor E on U
    length(E(U)) == 1 || error("exceptional divisor is not principal")
                       # sanity check -- default covering of Y is already trivializing for E!
    e = E(U)[1]
    if is_unit(e)
      ID[U] = h_total
      continue
    end

    ## find correct multiplicity for sanity check on result of iterated division
    ## iterated division only reliable, if multiplicity has expected value
    ## philosophy: if sanity check holds, iterated division ensures principality
    ##             if not, we throw an error with a good explanation (for now)
    if multEInC == -1
      _,multEInC = saturation_with_index(ideal(OO(U),[h_total]),ideal(OO(U),[e]))
    end

    ## now it is just division of ring elements
    epower = e^multEInC
    good, h_strict = divides(h_total,epower)
    bad,_ = divides(h_total, e*epower)
    (good && !bad) ||error("setting not suitable for iterated division -- use strict transform on IdealSheaf instead")

    ## fill in data of C_strict
    ID[U] = h_strict
  end

  ## we are good to go now
  C_strict = EffectiveCartierDivisor(Y, ID, check=false)
  return C_strict
end

function strict_transform(p::AbsSimpleBlowdownMorphism, C::CartierDivisor)
  X = codomain(p)
  Y = domain(p)
  X === scheme(C) || error("cartier divisor not defined on the codomain of the map")
  kk = coefficient_ring(C)
  result = CartierDivisor(Y, kk)
  for c in components(C)
    result = result + C[c] * strict_transform(p, c)
  end
  return result
end

@doc raw"""
    restrict(f::AbsCoveredSchemeMorphism,
        inc_dom::CoveredClosedEmbedding,
        inc_cod::CoveredClosedEmbedding;
        check::Bool=true
      )

For a diagram

  Z' ↪ Y
       ↓ f
  Z ↪  X
with `inc_dom` and `inc_cod` the respective horizontal maps
we assume ``f(Z') ⊂ Z``, compute and return the restriction ``f : Z' → Z``.
"""
function restrict(f::AbsCoveredSchemeMorphism,
    inc_dom::CoveredClosedEmbedding,
    inc_cod::CoveredClosedEmbedding;
    check::Bool=true
  )
  f_cov = covering_morphism(f)
  inc_dom_cov = covering_morphism(inc_dom)
  inc_cod_cov = covering_morphism(inc_cod)

  # We need to do the following.
  # - Pass to a common refinement ref_cod in X that both
  #   f and inc_cod can restrict to.
  # - Pass to a common refinement in Y
  ref_cod, a, b = _register!(common_refinement(codomain(f_cov), codomain(inc_cod_cov)), codomain(f))
  inc_cod_ref = restrict(inc_cod, ref_cod)
  f_res = restrict(f, ref_cod)
  ref_dom, aa, bb = _register!(common_refinement(domain(f_res), codomain(inc_dom_cov)), domain(f))
  inc_dom_ref = restrict(inc_dom, ref_dom)
  inc_dom_ref = compose(inc_dom_ref, aa)
  # Collecting the maps for the restricted projection here
  map_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(domain(inc_dom_ref))
    q_res = compose(inc_dom_ref[U], f_res[codomain(inc_dom_ref[U])])
    V = codomain(q_res)
    g = maps_with_given_codomain(inc_cod_ref, V)
    if !isone(length(g))
      error()
    end
    pre_V = domain(first(g))
    map_dict[U] = restrict(q_res, domain(q_res), pre_V, check=false)
  end
  psi = CoveringMorphism(domain(inc_dom_ref), domain(inc_cod_ref), map_dict, check=false)
  return CoveredSchemeMorphism(domain(inc_dom), domain(inc_cod), psi)
end

function _register!(data::Tuple{<:Covering, <:CoveringMorphism, <:CoveringMorphism},
    X::AbsCoveredScheme
  )
  push!(coverings(X), data[1])
  refinements(X)[(domain(data[2]), codomain(data[2]))] = data[2]
  refinements(X)[(domain(data[3]), codomain(data[3]))] = data[3]
  return data
end

function maps_with_given_codomain(phi::CoveringMorphism, V::AbsAffineScheme)
  result = Vector{AbsAffineSchemeMor}()
  for U in keys(morphisms(phi))
    floc = morphisms(phi)[U]
    codomain(floc) === V || continue
    push!(result, floc)
  end
  return result
end

##############################################################################
# show functions for Blowup morphisms
##############################################################################
function Base.show(io::IO, Bl::AbsSimpleBlowdownMorphism)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Blowup morphism")
  else
    print(io, "Blow-down: ", Lowercase(), domain(Bl))
    print(io, " -> ", Lowercase(), codomain(Bl))
  end
end

function show(io::IO, ::MIME"text/plain", Bl::AbsSimpleBlowdownMorphism)
  ## data of the original scheme
  X0 = codomain(Bl)
  C0 = get_attribute(X0, :simplified_covering, default_covering(X0))

  ## data of the blown up scheme
  X1 = domain(Bl)
  C1 = get_attribute(X1, :simplified_covering, default_covering(X1))

  ## data of the blowing up itself
  ED = exceptional_divisor(Bl)
  C_X0 = center(Bl)

  ## create the output
  io = pretty(io)
  println(io, "Blowup")
  print(io, Indent(), "of ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => C0, :label => "b"), X0)
  println(io)

  print(io, "in ", Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => C0, :label => "b"), C_X0)
  println(io, Dedent())

  println(io, "with domain")
  print(io, Indent(), Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => C1, :label => "a"), X1)
  println(io, Dedent())

  println(io, "and exceptional divisor")
  print(io, Indent(), Lowercase())
  show(IOContext(io, :show_semi_compact => true, :covering => C1, :label => "a"), ED)
  print(io, Dedent())
end

@attr AbsCoveredSchemeMorphism function isomorphism_on_complement_of_center(f::BlowupMorphism)
  iso_dict = get_attribute(f, :isos_on_complement_of_center)
  p = projection(f)
  X = domain(f)
  Y = codomain(f)
  dom_cov = Covering([U for U in keys(iso_dict)])
  inherit_gluings!(dom_cov, default_covering(X))
  cod_cov = Covering([codomain(p) for p in values(iso_dict)])
  inherit_gluings!(cod_cov, default_covering(Y))
  XU = CoveredScheme(dom_cov)
  YV = CoveredScheme(cod_cov)
  p_res_cov = CoveringMorphism(dom_cov, cod_cov, iso_dict, check=false)
  p_res = CoveredSchemeMorphism(XU, YV, p_res_cov)

  # Assemble the inverse
  iso_inv_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for (U, q) in iso_dict
    V = codomain(q)
    iso_inv_dict[V] = inverse(q)
  end
  p_res_inv_cov = CoveringMorphism(cod_cov, dom_cov, iso_inv_dict, check=false)
  p_res_inv = CoveredSchemeMorphism(YV, XU, p_res_inv_cov)

  set_attribute!(p_res, :inverse, p_res_inv)
  set_attribute!(p_res_inv, :inverse, p_res)
  return p_res
end

### Some functionality used by function fields

# If set this attribute shall return some isomorphism on some open subsets of the domain
# and codomain of phi. If both are irreducible, this automatically implies that both are
# dense, but we do not check for this.
@attr AbsAffineSchemeMor function isomorphism_on_open_subset(f::BlowupMorphism)
  X = domain(f)
  Y = codomain(f)
  iso_dict = get_attribute(f, :isos_on_complement_of_center)
  pr_res = first(values(iso_dict))
  U = domain(pr_res)
  V = codomain(pr_res)
  iso_U = _flatten_open_subscheme(U, default_covering(X))
  U_flat = codomain(iso_U)
  iso_V = _flatten_open_subscheme(V, default_covering(Y))
  V_flat = codomain(iso_V)
  phi = morphism(U_flat, V_flat, pullback(inverse(iso_U)).(pullback(pr_res).(pullback(iso_V).(gens(OO(V_flat))))), check=false)
  phi_inv = morphism(V_flat, U_flat, pullback(inverse(iso_V)).(pullback(inverse(pr_res)).(pullback(iso_U).(gens(OO(U_flat))))), check=false)
  set_attribute!(phi, :inverse, phi_inv)
  set_attribute!(phi_inv, :inverse, phi)
  return phi
end

function pullback(f::BlowupMorphism, g::VarietyFunctionFieldElem)
  X = domain(projection(f))
  Y = codomain(projection(f))
  FX = function_field(X)
  FY = function_field(Y)
  phi = isomorphism_on_open_subset(f)
  U = ambient_scheme(domain(phi))
  V = ambient_scheme(codomain(phi))
  parent(g) === FY || error("element does not belong to the correct field")
  h = g[V]
  pbg = fraction(pullback(phi)(OO(V)(numerator(h))))//fraction(pullback(phi)(OO(V)(denominator(h))))
  return FX.(pbg)
end

function pushforward(f::BlowupMorphism, g::VarietyFunctionFieldElem)
  X = domain(projection(f))
  Y = codomain(projection(f))
  FX = function_field(X)
  FY = function_field(Y)
  phi = inverse(isomorphism_on_open_subset(f))
  U = ambient_scheme(domain(phi))
  V = ambient_scheme(codomain(phi))
  parent(g) === FX || error("element does not belong to the correct field")
  h = g[V]
  pfg = fraction(pullback(phi)(OO(V)(numerator(h))))//fraction(pullback(phi)(OO(V)(denominator(h))))
  return FY.(pfg)
end

@attr AbsAffineSchemeMor function isomorphism_on_open_subset(phi::AbsCoveredSchemeMorphism)
  error("attribute not found; this needs to be set manually in general")
end

function compose(f::AbsSimpleBlowdownMorphism, g::AbsSimpleBlowdownMorphism)
  return composite_map(f, g)
end

function compose(f::AbsSimpleBlowdownMorphism, g::AbsCoveredSchemeMorphism)
  return composite_map(f, g)
end

function compose(f::AbsCoveredSchemeMorphism, g::AbsSimpleBlowdownMorphism)
  return composite_map(f, g)
end

