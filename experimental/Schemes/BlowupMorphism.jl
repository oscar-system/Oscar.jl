export BlowupMorphism
export center
export exceptional_divisor
export projection


########################################################################
# BlowupMorphism 
#
# A datastructure to maintain all information necessary to effectively 
# handle blowups. This is work in progress and will one day serve as 
# a building blow for sequences of blowups
########################################################################
mutable struct BlowupMorphism{
                              CodomainType<:AbsCoveredScheme
                             } # TODO: Derive this from AbsCoveredSchemeMorphism ? 
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
    return new{typeof(X)}(IP, X, I)
  end
end

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
function strict_transform(p::BlowupMorphism, inc::CoveredClosedEmbedding)
  Y = domain(p)
  X = codomain(p)
  Z = domain(inc)
  codomain(inc) === X || error("maps must have the same codomain")
  I_trans = strict_transform(p, image_ideal(inc))
  inc_Z_trans = CoveredClosedEmbedding(Y, I_trans, 
                                       covering=simplified_covering(Y), # Has been set by the previous call
                                       check=false)
  inc_cov = covering_morphism(inc_Z_trans)

  Z_trans = domain(inc_Z_trans)
  pr_res = restrict(projection(p), inc_Z_trans, inc)
  return Z_trans, inc_Z_trans, pr_res
end

@doc """
    strict_transform(p::BlowupMorphism, I::IdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
strict transform of ``I`` on ``Y``.
"""
function strict_transform(p::BlowupMorphism, I::IdealSheaf)
  Istrict,_ =_do_transform(p,I,-1)
  return Istrict
end

@doc """
    weak_transform(p::BlowupMorphism, I::IdealSheaf)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
weak transform ``J`` of ``I`` on ``Y``, i.e. an `IdealSheaf` satisfying ``E^m J = p^*I`` with ``m``
maximal and ``E`` the 'IdealSheaf' of the exceptional divisor of ``p``.
"""
function weak_transform(p::BlowupMorphism, I::IdealSheaf)
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
function weak_transform_with_multiplicity(p::BlowupMorphism, I::IdealSheaf)
  Iweak, multi = _do_transform(p,I,0)
  return Iweak,multi
end

@doc """
    controlled_transform(p::BlowupMorphism, I::IdealSheaf, b::Int)

For a `BlowupMorphism`  ``p : Y → X`` and an `IdealSheaf` ``I`` on ``X`` return the
controlled transform of ``I`` on ``Y`` with control ``b``,i.e. an `IdealSheaf` ``J`` such that
``E^b J = p^*I`` where ``E``denotes the `IdealSheaf` of the exceptional divisor.
"""
function controlled_transform(p::BlowupMorphism, I::IdealSheaf, b::Int)
  Icontrol,_ = _do_transform(p,I,b)
  return Icontrol
end

##########################################################################################################
## central internal method for strict, weak and controlled transforms of IdealSheafs and subschemes
##########################################################################################################
function _do_transform(p::BlowupMorphism, I::IdealSheaf, method::Int=-1)
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
  ID = IdDict{AbsSpec, Ideal}()

  ## get our hands on the coverings -- using simplified covering for CY
  pr = projection(p)
  p_cov_temp = covering_morphism(pr)
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
function strict_transform(p::BlowupMorphism, C::EffectiveCartierDivisor)
  X = scheme(C)
  Y = domain(p) 
  X === codomain(p) || error("cartier divisor is not defined on the codomain of the morphism")
  E = exceptional_divisor(p)
  ID = IdDict{AbsSpec, RingElem}()

  ## get our hands on the coverings -- trivializing covering for C leading the way
  CX = trivializing_covering(C)
  pr = projection(p)
  pr_refined = restrict(pr,CX)::CoveringMorphism
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
    if isunit(h_total)
      ID[U] = one(OO(U))
      continue
    end

    ## determine single generator of Cartier divisor E on U
    length(E(U)) == 1 || error("exceptional divisor is not principal")
                       # sanity check -- default covering of Y is already trivializing for E!
    e = E(U)[1]
    if isunit(e)
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

function strict_transform(p::BlowupMorphism, C::CartierDivisor)
  X = codomain(p)
  Y = domain(p) 
  X === scheme(C) || error("cartier divisor not defined on the codomain of the map")
  kk = coefficient_ring(C)
  result = CartierDivisor(Y, kk)
  for c in components(C)
    result = result + C[c]*strict_transform(p, c)
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
  map_dict = IdDict{AbsSpec, AbsSpecMor}()
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

function maps_with_given_codomain(phi::CoveringMorphism, V::AbsSpec)
  result = Vector{AbsSpecMor}()
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
function Base.show(io::IO,Bl::BlowupMorphism)

## data of the original scheme
  X0 = codomain(Bl)
  C0 = (has_attribute(X0, :simplified_covering) ? simplified_covering(X0) : default_covering(X0))
  n0 = npatches(C0)

## data of the blown up scheme
  X1 = domain(Bl)
  C1 = (has_attribute(X1, :simplified_covering) ? simplified_covering(X1) : default_covering(X1))
  n1 = npatches(C1)

## create the output
  println(io,"Blow up of a Covered Scheme with ",n0," Charts leading to a Covered Scheme with ",n1," Charts")
end

@doc raw"""
  show_details(Bl::BlowupMorphism)

For a `BlowupMorphism` ``p : Y → X`` display domain, codomain, center and exceptional divisor in a detailed view.
# Examples
```jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"];

julia> A3 = Spec(R);

julia> I = ideal(R,[x,y,z]);

julia> bl = blow_up(A3,I)
Blow up of a Covered Scheme with 1 Charts leading to a Covered Scheme with 3 Charts

julia> show_details(bl)
Blow up of a Covered Scheme with 1 Charts at a 0-dimensional Center leading to a Covered Scheme with 3 Charts.

=====================================
Affine charts of the original scheme:
Spec of Multivariate Polynomial Ring in x, y, z over Rational Field

=====================================
Affine charts of the blown up scheme:
Spec of Localization of Quotient of Multivariate Polynomial Ring in (s1//s0), (s2//s0), x over Rational Field by ideal(0, 0, 0) at the multiplicative set powers of QQMPolyRingElem[1]

Spec of Localization of Quotient of Multivariate Polynomial Ring in (s0//s1), (s2//s1), y over Rational Field by ideal(0, 0, 0) at the multiplicative set powers of QQMPolyRingElem[1]

Spec of Localization of Quotient of Multivariate Polynomial Ring in (s0//s2), (s1//s2), z over Rational Field by ideal(0, 0, 0) at the multiplicative set powers of QQMPolyRingElem[1]

=====================================
Data of center:
Ideal Sheaf on Covered Scheme with 1 Charts:

Chart 1:
   ideal(x, y, z)

=====================================
Exceptional divisor:
Effective Cartier Divisor on Covered Scheme with 3 Charts:

Chart 1:
   ideal in Localization of Quotient of Multivariate Polynomial Ring in (s1//s0), (s2//s0), x, y, z over Rational Field by ideal(-(s1//s0)*z + (s2//s0)*y, (s2//s0)*x - z, (s1//s0)*x - y) at the multiplicative set powers of QQMPolyRingElem[1] generated by [x]

Chart 2:
   ideal in Localization of Quotient of Multivariate Polynomial Ring in (s0//s1), (s2//s1), x, y, z over Rational Field by ideal((s2//s1)*y - z, -(s0//s1)*z + (s2//s1)*x, -(s0//s1)*y + x) at the multiplicative set powers of QQMPolyRingElem[1] generated by [y]

Chart 3:
   ideal in Localization of Quotient of Multivariate Polynomial Ring in (s0//s2), (s1//s2), x, y, z over Rational Field by ideal(-(s1//s2)*z + y, -(s0//s2)*z + x, -(s0//s2)*y + (s1//s2)*x) at the multiplicative set powers of QQMPolyRingElem[1] generated by [z]

```
"""
function show_details(Bl::BlowupMorphism)
   show_details(stdout, Bl)
end

function show_details(io::IO, Bl::BlowupMorphism)
## data of the original scheme
  X0 = codomain(Bl)
  C0 = (has_attribute(X0, :simplified_covering) ? simplified_covering(X0) : default_covering(X0))
  n0 = npatches(C0)

## data of the blown up scheme
  X1 = domain(Bl)
  C1 = (has_attribute(X1, :simplified_covering) ? simplified_covering(X1) : default_covering(X1))
  n1 = npatches(C1)

## data of the blowing up itself
  ED = exceptional_divisor(Bl)
  C_X0 = Bl.center

## create the output
  println(io,"Blow up of a Covered Scheme with ",n0," Charts at a ",dim(C_X0),"-dimensional Center leading to a Covered Scheme with ",n1," Charts.\n")

  println(io,"=====================================")
  println(io,"Affine charts of the original scheme:")
  for U in patches(C0)
    println(io,U,"\n")
  end

  println(io,"=====================================")
  println(io,"Affine charts of the blown up scheme:")
  for V in patches(C1)
    println(io,V,"\n")
  end

  println(io,"=====================================")
  println(io,"Data of center:")
  show_details(io,Bl.center)

  println(io,"=====================================")
  println(io,"Exceptional divisor:")
  show_details(io,exceptional_divisor(Bl))
end

