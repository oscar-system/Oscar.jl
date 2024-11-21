########################################################################
# `AbsDesingMor` and `AbsBlowupMorphism`
#
# An abstract type for blowup morphisms.
#
# This should also comprise sequences of simple blowups leading
# to a partial or full resolution of singularities. The interface
# is specified below.
########################################################################

### See experimental/Schemes/src/BlowupMorphismTypes.jl for the definition of 
# `AbsDesingMor` and `AbsBlowupMorphism`.

# The interface inherits all functionality from AbsCoveredSchemeMorphism.
# This is extended by the following:
@doc raw"""
    exceptional_divisor(f::AbsBlowupMorphism)

Return a `CartierDivisor` on the `domain` of `f` which is the
preimage of the vanishing locus of the `center` of `f`.

!!! note If `f` is not a classical blowup, but, for instance, a small contraction, then the exceptional locus need not be a divisor. See `exceptional_locus` for such cases.
"""
function exceptional_divisor(f::AbsBlowupMorphism)
  error("not implemented")
end

@doc raw"""
    exceptional_locus(f::AbsBlowupMorphism)

Return an `AbsAlgebraicCycle` on the `domain` of `f` which is the preimage
of the `center` of `f` in a reasonable sense. In particular, `f` is an isomorphism
onto its image outside the support of its `exceptional_locus`. See also `exceptional_divisor`.
"""
function exceptional_locus(f::AbsBlowupMorphism)
  error("not implemented")
end

@doc raw"""
    center(f::AbsBlowupMorphism)

Return an `AbsIdealSheaf` on the `codomain` of `f` such that on the complement
of the vanishing locus of that ideal sheaf `f` is an isomorphism.
The support of the `exceptional_locus` of `f` coincides with the vanishing
locus of the pullback of the `center`.

!!! note For classical blowup constructions this will be the center which has been blown up. For more exotic blowups, however, this might be more special.
"""
function center(f::AbsBlowupMorphism)
  error("not implemented")
end

@doc raw"""
    strict_transform(f::AbsBlowupMorphism, a::Any)

Take the closure of the `pullback` of `a` along `f` on the complement
of the `exceptional_locus` in a mathematically reasonable sense.
Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function strict_transform(f::AbsBlowupMorphism, a::Any)
  error("not implemented")
end

@doc raw"""
    total_transform(f::AbsBlowupMorphism, a::Any)

Take the `pullback` or preimage of `a` along `f` in a mathematically
reasonable sense.
Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function total_transform(f::AbsBlowupMorphism, a::Any)
  error("not implemented")
end

########################################################################
# `AbsSimpleBlowupMorphism`
# An abstract type for classical blowups of ideal sheaves.
#
# This can either be a BlowupMorphism as below, but also a special
# toric morphism induced by fan subdivisions.
########################################################################

### See BlowupMorphismTypes.jl for the definition of 
# `AbsSimpleBlowupMorphism`

@doc raw"""
    exceptional_divisor(f::AbsSimpleBlowupMorphism)

Return a `CartierDivisor` on the `domain` of `f` which coincides
with the pullback of the `center` of `f`.
"""
function exceptional_divisor(f::AbsSimpleBlowupMorphism)
  error("not implemented")
end

@doc raw"""
    exceptional_locus(f::AbsBlowupMorphism)

Return the `WeilDivisor` of the `exceptional_divisor`.
"""
function exceptional_locus(f::AbsSimpleBlowupMorphism)
  error("not implemented")
end

@doc raw"""
    center(f::AbsSimpleBlowupMorphism)

Return an `AbsIdealSheaf` on the `codomain` of `f` the blowup of which
leads to `f`.
"""
function center(f::AbsSimpleBlowupMorphism)
  error("not implemented")
end

@doc raw"""
    strict_transform(f::AbsSimpleBlowupMorphism, a::Any)

Take the closure of the `pullback` of `a` along `f` on the complement
of the `exceptional_divisor` in a mathematically reasonable sense.

Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function strict_transform(f::AbsSimpleBlowupMorphism, a::Any)
  error("not implemented")
end

@doc raw"""
    total_transform(f::AbsSimpleBlowupMorphism, a::Any)

Take the `pullback` or preimage of `a` along `f` in a mathematically
reasonable sense.

Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function total_transform(f::AbsSimpleBlowupMorphism, a::Any)
  error("not implemented")
end

@doc raw"""
    controlled_transform(f::AbsSimpleBlowupMorphism, a::Any, k::Int)

Take the `pullback` or preimage of `a` along `f` saturated by `k` times
the `exceptional_divisor` of `f` in a mathematically reasonable sense.

Depending on `a` this either returns the corresponding object on
the `domain` of `f`, or a tuple consisting of the object and eventually
induced maps.
"""
function controlled_transform(f::AbsSimpleBlowupMorphism, a::Any, k::Int)
  error("not implemented")
end

########################################################################
# `BlowupMorphism`
#
# A datastructure to maintain all information necessary to effectively
# handle blowups. This is work in progress and will one day serve as
# a building block for sequences of blowups
########################################################################

### See experimental/Schemes/src/BlowupMorphismTypes.jl for the concrete type 
# `BlowupMorphism`.

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
`AbsIdealSheaf` ``ℐ`` on X, return the `EffectiveCartierDivisor` ``E``
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
function strict_transform(p::AbsSimpleBlowupMorphism, inc::CoveredClosedEmbedding)
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
  pr_res = restrict(projection(p), inc_Z_trans, inc; check=false)

  return Z_trans, inc_Z_trans, pr_res
end

@doc """
    strict_transform(p::BlowupMorphism, I::AbsIdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `AbsIdealSheaf` ``I`` on ``X`` return the
strict transform of ``I`` on ``Y``.
"""
function strict_transform(p::AbsSimpleBlowupMorphism, I::AbsIdealSheaf)
  return StrictTransformIdealSheaf(p, I)
  Istrict,_ =_do_transform(p, I, -1)
  return Istrict
end

@doc """
    weak_transform(p::BlowupMorphism, I::AbsIdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `AbsIdealSheaf` ``I`` on ``X`` return the
weak transform ``J`` of ``I`` on ``Y``, i.e. an `AbsIdealSheaf` satisfying ``E^m J = p^*I`` with ``m``
maximal and ``E`` the 'AbsIdealSheaf' of the exceptional divisor of ``p``.
"""
function weak_transform(p::AbsSimpleBlowupMorphism, I::AbsIdealSheaf)
  Iweak,_ =_do_transform(p,I,0)
  return Iweak
end

@doc """
    weak_transform_with_multiplicity((p::BlowupMorphism, I::AbsIdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `AbsIdealSheaf` ``I`` on ``X`` return the
weak transform ``J`` of ``I`` on ``Y`` and the multiplicity ``m`` of the exceptional divisor, i.e.
the maximal ``m`` such that ``E^m J = p^*I``, where ``E`` denotes the `AbsIdealSheaf` of the exceptional
divisor of ``p``.
"""
function weak_transform_with_multiplicity(p::AbsSimpleBlowupMorphism, I::AbsIdealSheaf)
  Iweak, multi = _do_transform(p,I,0)
  return Iweak,multi
end

@doc """
    controlled_transform(p::BlowupMorphism, I::AbsIdealSheaf, b::Int)

For a `BlowupMorphism`  ``p : Y → X`` and an `AbsIdealSheaf` ``I`` on ``X`` return the
controlled transform of ``I`` on ``Y`` with control ``b``,i.e. an `AbsIdealSheaf` ``J`` such that
``E^b J = p^*I`` where ``E``denotes the `AbsIdealSheaf` of the exceptional divisor.
"""
function controlled_transform(p::AbsSimpleBlowupMorphism, I::AbsIdealSheaf, b::Int)
  Icontrol,_ = _do_transform(p,I,b)
  return Icontrol
end

##########################################################################################################
## central internal method for strict, weak and controlled transforms of AbsIdealSheafs and subschemes
##########################################################################################################
function _do_transform(p::AbsSimpleBlowupMorphism, I::AbsIdealSheaf, method::Int=-1)
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
function strict_transform(p::AbsSimpleBlowupMorphism, C::EffectiveCartierDivisor)
  return strict_transform_with_multiplicity(p,C)[1]
end

function strict_transform_with_multiplicity(p::AbsSimpleBlowupMorphism, C::EffectiveCartierDivisor)
  X = ambient_scheme(C)
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
    (good && !bad) ||error("setting not suitable for iterated division -- use strict transform on AbsIdealSheaf instead")

    ## fill in data of C_strict
    ID[U] = h_strict
  end

  ## we are good to go now
  C_strict = EffectiveCartierDivisor(Y, ID, check=false)
  return C_strict,multEInC
end

function strict_transform(p::AbsSimpleBlowupMorphism, C::CartierDivisor)
  X = codomain(p)
  Y = domain(p)
  X === ambient_scheme(C) || error("cartier divisor not defined on the codomain of the map")
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

  # Build up the common refinement 
  success, dom_ref = is_refinement(codomain(inc_dom_cov), domain(f_cov))
  @assert success "restriction not implemented for this constellation of refinements"
  inc_dom_f = codomain(inc_dom_cov) === domain(f_cov) ? compose(inc_dom_cov, f_cov) : compose(compose(inc_dom_cov, dom_ref), f_cov)

  success, cod_ref_map = is_refinement(codomain(inc_dom_f), codomain(inc_cod_cov))
  inc_dom_f_ref = inc_dom_f # initialize the variable
  if !success
    success, cod_ref_map = is_refinement(codomain(inc_cod_cov), codomain(inc_dom_f))
    @assert success "restriction not implemented for this constellation of refinements"
    dom_ref2, to_inc_dom_f, to_inc_cod_cov = fiber_product(inc_dom_f, cod_ref_map)
    inc_dom_f_ref = to_inc_cod_cov
  else
    inc_dom_f_ref = compose(inc_dom_f, cod_ref_map)
  end

  # Collecting the maps for the restricted projection here
  map_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  for U in patches(domain(inc_dom_f_ref))
    q_res = inc_dom_f_ref[U]
    V = codomain(q_res)
    g = maps_with_given_codomain(inc_cod_cov, V)
    if !isone(length(g))
      error()
    end
    pre_V = domain(first(g))
    map_dict[U] = restrict(q_res, domain(q_res), pre_V; check)
  end
  psi = CoveringMorphism(domain(inc_dom_f_ref), domain(inc_cod_cov), map_dict; check)
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
function Base.show(io::IO, Bl::AbsSimpleBlowupMorphism)
  io = pretty(io)
  if is_terse(io)
    print(io, "Blowup morphism")
  else
    print(io, "Blow-up: ", Lowercase(), domain(Bl))
    print(io, " -> ", Lowercase(), codomain(Bl))
  end
end

function show(io::IO, ::MIME"text/plain", Bl::AbsSimpleBlowupMorphism)
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

function compose(f::AbsSimpleBlowupMorphism, g::AbsSimpleBlowupMorphism)
  return composite_map(f, g)
end

function compose(f::AbsSimpleBlowupMorphism, g::AbsCoveredSchemeMorphism)
  return composite_map(f, g)
end

function compose(f::AbsCoveredSchemeMorphism, g::AbsSimpleBlowupMorphism)
  return composite_map(f, g)
end

### See experimental/Schemes/src/BlowupMorphismTypes.jl for declaration of `BlowUpSequence`.

########################################################################
# Convenience method                                                   #
########################################################################

function blow_up(m::AbsCoveredScheme, I::AbsIdealSheaf; coordinate_name = "e")
  @assert m === scheme(I) "Incompatible scheme and ideal sheaf"
  return blow_up(I)
end


### See experimental/Schemes/src/BlowupMorphismTypes.jl for the declaration of `StrictTransformIdealSheaf`.

morphism(I::StrictTransformIdealSheaf) = I.morphism
original_ideal_sheaf(I::StrictTransformIdealSheaf) = I.orig
underlying_presheaf(I::StrictTransformIdealSheaf) = I.underlying_presheaf

function produce_object_on_affine_chart(I::StrictTransformIdealSheaf, U::AbsAffineScheme)
  f = morphism(I)
  X = domain(f)
  Y = codomain(f)
  J = original_ideal_sheaf(I)
  @assert any(x->x===U, affine_charts(X))
  if f isa ToricBlowupMorphism
    # This is not actually an exceptional divisor of a blowup along an ideal sheaf.
    # This is the prime Weil divisor corresponding to the added/chosen ray.
    # This is the exceptional divisor of a blowup along a certain Rees algebra.
    E = exceptional_prime_divisor(f)
  else
    E = exceptional_divisor(f)
  end
  IE = ideal_sheaf(E)
  # We assume that the covering morphism has the default_covering of X as its domain.
  @assert domain(covering_morphism(f)) === default_covering(X) "not implemented for this covering"
  f_loc = covering_morphism(f)[U]
  V = codomain(f_loc)
  IE_loc = IE(U)
  # It is usually better to pass to the simplified covering to do the computations
  simp_cov = simplified_covering(X)
  U_simp = first([V for V in patches(simp_cov) if original(V) === U])
  a, b = identification_maps(U_simp)
  tot = pullback(f_loc)(J(V))
  if isone(ngens(IE_loc))
    result = _iterative_saturation(pullback(a)(tot), elem_type(OO(U_simp))[pullback(a)(u) for (u, _) in factor(lifted_numerator(first(gens(IE_loc))))])
    return pullback(b)(result)
  else
    result, _ = saturation_with_index(pullback(a)(tot), pullback(a)(IE_loc))
    return result
  end
end

@attr Bool function is_prime(I::StrictTransformIdealSheaf)
  is_subset(center(morphism(I)), radical(original_ideal_sheaf(I))) && return false # It's the unit ideal sheaf in this case
  return is_prime(original_ideal_sheaf(I))
end

function cheap_sub_ideal(I::StrictTransformIdealSheaf, U::AbsAffineScheme)
  II = pullback_ideal_sheaf(I)
  return cheap_sub_ideal(II, U)
end

@attr PullbackIdealSheaf function pullback_ideal_sheaf(I::StrictTransformIdealSheaf)
  f = morphism(I)
  J = original_ideal_sheaf(I)
  return PullbackIdealSheaf(f, J)
end

