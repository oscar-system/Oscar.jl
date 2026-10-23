##############################################################################
#
# Conversion to and from Singular: in particular, some rings are
# special as they exist natively in Singular and thus should be used
#
##############################################################################
#
# Singular's polynomial rings are not recursive:
# 1. iso_oscar_singular_poly_ring(R::Ring) tries to create a Singular.PolyRing (with
#    elements of type Singular.spoly) isomorphic to R
# 2. iso_oscar_singular_coeff_ring(R::Ring) tries to create a ring isomorphic to R that is
#    acceptable to Singular.jl as 'coefficients'
#
# a real native Singular polynomial ring with Singular's native QQ as coefficients:
#  codomain(iso_oscar_singular_poly_ring(QQ[t])) => Singular.PolyRing{Singular.n_Q}
#
# Singular's native Fp(5):
#  codomain(iso_oscar_singular_coeff_ring(GF(5))) => Singular.N_ZpField
#
# Singular wrapper of the Oscar type QQPolyRingElem:
#  codomain(iso_oscar_singular_coeff_ring((QQ[t])) => Singular.N_Ring{QQPolyRingElem}
#
# even more wrappings of the immutable Oscar type FpFieldElem:
#  codomain(iso_oscar_singular_coeff_ring(GF(ZZRingElem(5)))) => Singular.N_Field{Singular.FieldElemWrapper{FpField, FpFieldElem}}

"""
    iso_oscar_singular_coeff_ring(R::Ring) -> Map

Return a ring isomorphism of `R` onto a native Singular ring.
"""
iso_oscar_singular_coeff_ring

"""
    iso_oscar_singular_poly_ring(R::Ring, ...; kw...) -> Map

Given a polynomial ring `R[x_1,...x_n]` return a ring isomorphism onto a native
Singular ring `S[x_1,...,x_n]` with `R` isomorphic to `S`.
"""
iso_oscar_singular_poly_ring

#############################################################################
#
# Coefficient rings
#
abstract type OscarSingularCoefficientRingMap{D, C} <: Map{D, C, Any, Any} end

domain(f::OscarSingularCoefficientRingMap) = f.R

codomain(f::OscarSingularCoefficientRingMap) = f.S

# fallback
function image(f::OscarSingularCoefficientRingMap, x)
  parent(x) !== domain(f) && error("Element not in domain")
  return codomain(f)(x)
end

(f::OscarSingularCoefficientRingMap)(x) = image(f, x)

function preimage(f::OscarSingularCoefficientRingMap, x)
  parent(x) !== codomain(f) && error("Element not in codomain")
  return domain(f)(x)
end

#############################################################################
#
# generic catchall 
#
struct OscarSingularCoefficientRingMapGeneric{D,C} <: OscarSingularCoefficientRingMap{D,C} 
  R::D
  S::C
end

function iso_oscar_singular_coeff_ring(F::AbstractAlgebra.Ring)
  return OscarSingularCoefficientRingMapGeneric(F, Singular.CoefficientRing(F))
end

# image(f, a) done by parent call overloading

function preimage(f::OscarSingularCoefficientRingMapGeneric, a::Singular.n_unknown)
  parent(a) !== codomain(f) && error("Element not in codomain")
  b = Singular.libSingular.julia(Singular.libSingular.cast_number_to_void(a.ptr))
  if b isa Singular.FieldElemWrapper || b isa Singular.RingElemWrapper
    # handle immutable ring elements (such as QQAbFieldElem) which are
    # put in a wrapper on the Singular side
    return b.data::elem_type(domain(f))
  end
  return b::elem_type(domain(f))
end

# ZZ
iso_oscar_singular_coeff_ring(R::ZZRing) = OscarSingularCoefficientRingMapGeneric(R, Singular.Integers())

# QQ
iso_oscar_singular_coeff_ring(R::QQField) = OscarSingularCoefficientRingMapGeneric(R, Singular.Rationals())

# prime field
iso_oscar_singular_coeff_ring(R::fpField) = OscarSingularCoefficientRingMapGeneric(R, Singular.Fp(characteristic(R)))

# ZZ/nZZ, n small and big
iso_oscar_singular_coeff_ring(R::Union{zzModRing, ZZModRing}) = OscarSingularCoefficientRingMapGeneric(R, Singular.residue_ring(Singular.Integers(), BigInt(modulus(R)))[1]
)

# QQ(a)
function iso_oscar_singular_coeff_ring(R::AbsSimpleNumField)
  minpoly = defining_polynomial(R)
  Qa = parent(minpoly)
  SQa, (Sa,) = Singular.FunctionField(Singular.QQ, _variables_for_singular(symbols(Qa)))
  Sminpoly = SQa(coeff(minpoly, 0))
  var = one(SQa)
  for i in 1:degree(minpoly)
      var = mul!(var, Sa)
      Sminpoly = addmul!(Sminpoly, SQa(coeff(minpoly, i)), var)
  end
  SK, _ = Singular.AlgebraicExtensionField(SQa, Sminpoly)
  return OscarSingularCoefficientRingMapGeneric(R, SK)
end

# GF(p, n), small p
function iso_oscar_singular_coeff_ring(F::fqPolyRepField) 
  minpoly = modulus(F)
  Fa = parent(minpoly)
  SFa, (Sa,) = Singular.FunctionField(Singular.Fp(characteristic(F)),
                                                    _variables_for_singular(symbols(Fa)))
  Sminpoly = SFa(coeff(minpoly, 0))
  var = one(SFa)
  for i in 1:degree(minpoly)
      var = mul!(var, Sa)
      Sminpoly = addmul!(Sminpoly, SFa(coeff(minpoly, i)), var)
  end
  SF, _ = Singular.AlgebraicExtensionField(SFa, Sminpoly)
  return OscarSingularCoefficientRingMapGeneric(F, SF)
end

# Conversion between an Oscar finite field `K` and a Singular algebraic
# extension `SK` of F_p, both presented by the same minimal polynomial.
# `c(i)` returns the i-th coordinate in the power basis as an integer.
function _singular_algext_elem(SK::Singular.N_AlgExtField, c, d::Int)
  Sa = gen(SK)
  res = SK(c(0))
  var = one(SK)
  for i in 1:d-1
    var = mul!(var, Sa)
    res = addmul!(res, SK(c(i)), var)
  end
  return res
end

function _oscar_finite_field_elem(K::Field, a::Singular.n_algExt)
  SF = parent(Singular.modulus(parent(a)))
  SFa = SF(a)
  numSa = Singular.n_transExt_to_spoly(numerator(SFa))
  denSa = first(AbstractAlgebra.coefficients(Singular.n_transExt_to_spoly(denominator(SFa))))
  @assert isone(denSa)
  res = zero(K)
  Ka = gen(K)
  for (c, e) in zip(AbstractAlgebra.coefficients(numSa), AbstractAlgebra.exponent_vectors(numSa))
    res += K(Int(c))*Ka^e[1]
  end
  return res
end

function image(f::OscarSingularCoefficientRingMapGeneric{fqPolyRepField}, a::fqPolyRepFieldElem)
  parent(a) !== domain(f) && error("Element not in domain")
  return _singular_algext_elem(codomain(f), i -> coeff(a, i), degree(domain(f)))
end

function preimage(f::OscarSingularCoefficientRingMapGeneric{fqPolyRepField}, a::Singular.n_algExt)
  parent(a) !== codomain(f) && error("Element not in codomain")
  return _oscar_finite_field_elem(domain(f), a)
end

#############################################################################
#
# Finite field (FqField)
#
struct OscarSingularCoefficientRingMapFqField{D} <: OscarSingularCoefficientRingMap{FqField, D}
  R::FqField
  S::D
  iso::MapFromFunc{FqField, FqField}

  OscarSingularCoefficientRingMapFqField(R, S) = new{typeof(S)}(R, S)
  OscarSingularCoefficientRingMapFqField(R, S, iso) = new{typeof(S)}(R, S, iso)
end

function _absolute_field(F::FqField)
  if is_absolute(F)
    error("don't use me, you are already absolute")
  end
  Kx = parent(defining_polynomial(F))
  Fabs, QQtoFabs = Nemo._residue_field(defining_polynomial(F); absolute = true, check = false)
  Fabsx, = polynomial_ring(Fabs, :x; cached = false)
  return Fabs, MapFromFunc(Fabs, F, a -> F(preimage(QQtoFabs, a)), b -> begin a = Fabs(); Nemo.set!(a, b); a end)
end

# FqField (aka fq_default from flint)
function iso_oscar_singular_coeff_ring(F::FqField)
  # we are way beyond type stability, so just do what you want
  if !is_absolute(F)
    Fabs, FabstoF = _absolute_field(F)
  else
    Fabs = F
  end

  if absolute_degree(F) == 1
    S = Singular.Fp(characteristic(F))
  elseif nbits(characteristic(F)) <= 29
    minpoly = modulus(F)
    Fa = parent(minpoly)
    SFa, (Sa,) = Singular.FunctionField(Singular.Fp(characteristic(F)),
                                        _variables_for_singular(symbols(Fa)))
    Sminpoly = SFa(lift(ZZ, coeff(minpoly, 0)))
    var = one(SFa)
    for i in 1:degree(minpoly)
        var = mul!(var, Sa)
        Sminpoly = addmul!(Sminpoly, SFa(lift(ZZ, coeff(minpoly, i))), var)
    end
    S, _ = Singular.AlgebraicExtensionField(SFa, Sminpoly)
  else
    S = Singular.CoefficientRing(F)
  end

  if !is_absolute(F)
    return OscarSingularCoefficientRingMapFqField(F, S, FabstoF)
  end
  return OscarSingularCoefficientRingMapFqField(F, S)
end

function image(f::OscarSingularCoefficientRingMapFqField, a::FqFieldElem)
  parent(a) !== domain(f) && error("Element not in domain")

  if codomain(f) isa Singular.N_ZpField
    return codomain(f)(lift(ZZ, a))
  end

  if codomain(f) isa Singular.N_ZnRing
    return codomain(f)(a)
  end

  if codomain(f) isa Singular.N_UnknownSingularCoefficientRing
    absolute_degree(domain(f)) == 1 || error("Cannot convert non-prime finite field element to unknown Singular coefficient ring")
    return codomain(f)(BigInt(lift(ZZ, a)))
  end

  if codomain(f) isa Singular.N_Field
    return codomain(f)(a)
  end

  b = isdefined(f, :iso) ? preimage(f.iso, a) : a
  return _singular_algext_elem(codomain(f), i -> lift(ZZ, coeff(b, i)), degree(parent(b)))
end

function preimage(f::OscarSingularCoefficientRingMapFqField, a::Singular.n_FieldElem)
  parent(a) !== codomain(f) && error("Element not in codomain")
  return Singular.libSingular.julia(Singular.libSingular.cast_number_to_void(a.ptr))::FqFieldElem
end

function preimage(f::OscarSingularCoefficientRingMapFqField, a::Singular.n_Zp)
  parent(a) !== codomain(f) && error("Element not in codomain")
  return domain(f)(Int(a))
end

function preimage(f::OscarSingularCoefficientRingMapFqField, a::Singular.n_Zn)
  parent(a) !== codomain(f) && error("Element not in codomain")
  return domain(f)(a)
end

function preimage(f::OscarSingularCoefficientRingMapFqField, a::Singular.n_unknownsingularcoefficient)
  parent(a) !== codomain(f) && error("Element not in codomain")
  absolute_degree(domain(f)) == 1 || error("Cannot convert unknown Singular coefficient to non-prime finite field element")
  S = parent(a)
  z = GC.@preserve a S Singular.libSingular.n_GetMPZ(a.ptr, S.ptr)
  return domain(f)(ZZ(z))
end

function preimage(f::OscarSingularCoefficientRingMapFqField, a::Singular.n_algExt)
  parent(a) !== codomain(f) && error("Element not in codomain")

  isdefined(f, :iso) || return _oscar_finite_field_elem(domain(f), a)
  return image(f.iso, _oscar_finite_field_elem(domain(f.iso), a))
end

#############################################################################
#
# fraction field of polynomial rings over QQ and Fp
#
struct OscarSingularCoefficientRingMapFractionField{U, V, W, X} <: OscarSingularCoefficientRingMap{U, V}
  R::U
  S::V
  g::W
  Spoly::X #= we create univariate polynomials over S during the conversion =#
end

function _map_oscar_singular_univariate(Rx, g, f::Singular.spoly)
  @assert base_ring(Rx) === domain(g)
  @assert ngens(parent(f)) == 1
  return Rx(elem_type(domain(g))[preimage(g, c) for c in coefficients_of_univariate(f)])
end

function iso_oscar_singular_coeff_ring(F::AbstractAlgebra.Generic.FracField{<:PolyRingElem{T}}) where {T <: Union{FqFieldElem, QQFieldElem}}
  R = base_ring(F)
  g = iso_oscar_singular_coeff_ring(base_ring(R))
  S, = Singular.FunctionField(codomain(g), [var(R)])
  Sx, = polynomial_ring(S, :x; cached = false)
  return OscarSingularCoefficientRingMapFractionField(F, S, g, Sx)
end

function preimage(f::OscarSingularCoefficientRingMapFractionField{<:AbstractAlgebra.Generic.FracField{<:PolyRingElem}}, a::Singular.n_transExt)
  parent(a) !== codomain(f) && error("Element not in codomain")
  F = domain(f)
  R = base_ring(F)
  n, d = Singular.n_transExt_to_spoly.([numerator(a), denominator(a)]; cached = false)
  return F(_map_oscar_singular_univariate(R, f.g, n), _map_oscar_singular_univariate(R, f.g, d))
end

function image(f::OscarSingularCoefficientRingMapFractionField{<:AbstractAlgebra.Generic.FracField{<:PolyRingElem}}, a)
  parent(a) !== domain(f) && error("Element not in domain")
  F = base_ring(base_ring(domain(f))) # the F in domain(f) = F(X)
  K = codomain(f)
  @assert Singular.transcendence_degree(K) == 1 "wrong number of generators"
  t, = Singular.transcendence_basis(K)
  # It must be a transcendental extension of Q or Fp (other things are not supported)
  n = map_coefficients(numerator(a); parent = f.Spoly) do x
    if F isa FinField
      K(lift(ZZ, x))
    else
      K(f.g(x))
    end
  end
  d = map_coefficients(denominator(a); cached = false) do x
    if F isa FinField
      K(lift(ZZ, x))
    else
      K(f.g(x))
    end
  end
  return divexact(n(t), d(t))
end

# maps for fraction field of multivariate polynomial ring
function iso_oscar_singular_coeff_ring(F::AbstractAlgebra.Generic.FracField{<:MPolyRingElem{T}}) where {T <: Union{FqFieldElem, QQFieldElem}}
  R = base_ring(F)
  g = iso_oscar_singular_coeff_ring(base_ring(R))
  S, = Singular.FunctionField(codomain(g),_variables_for_singular(symbols(R)))
  Sx, = polynomial_ring(S, nvars(R), :x; cached = false)
  return OscarSingularCoefficientRingMapFractionField(F, S, g, Sx)
end

function image(f::OscarSingularCoefficientRingMapFractionField{<:AbstractAlgebra.Generic.FracField{<:MPolyRingElem}}, a)
  parent(a) !== domain(f) && error("Element not in domain")
  x = Singular.transcendence_basis(codomain(f))
  F = base_ring(base_ring(domain(f))) # the F in domain(f) = F(X)
  K = codomain(f)
  map_coeff_map = x -> begin
    if F isa FinField
      K(lift(ZZ, x))
    else
      K(f.g(x))
    end
  end
  return divexact(evaluate(map_coefficients(map_coeff_map,
                                            numerator(a); parent = f.Spoly), x
                          ),
                  evaluate(map_coefficients(map_coeff_map,
                                            denominator(a); parent = f.Spoly), x
                          )
                 )
end

function preimage(f::OscarSingularCoefficientRingMapFractionField, a::Singular.n_transExt)
  parent(a) !== codomain(f) && error("Element not in codomain")
  F = domain(f)
  R = base_ring(F)
  n, d = Singular.n_transExt_to_spoly.([numerator(a), denominator(a)]; cached = false)
  return F(map_coefficients(x -> preimage(f.g, x), n; parent = R),
           map_coefficients(x -> preimage(f.g, x), d; parent = R))
end

#############################################################################
#
# rational function field
#
struct OscarSingularCoefficientRingMapRationalFunctionField{D, C, W} <: OscarSingularCoefficientRingMap{D, C}
  R::D
  S::C
  g::W
end

function iso_oscar_singular_coeff_ring(R::Generic.RationalFunctionField)
  g = iso_oscar_singular_coeff_ring(R.fraction_field)
  return OscarSingularCoefficientRingMapRationalFunctionField(R, codomain(g), g)
end

function image(f::OscarSingularCoefficientRingMapRationalFunctionField, a)
  parent(a) !== domain(f) && error("Element not in domain")
  return f.g(a.d)
end

function preimage(f::OscarSingularCoefficientRingMapRationalFunctionField, a)
  parent(a) !== codomain(f) && error("Element not in codomain")
  return domain(f)(preimage(f.g, a))
end

# All Singular rings created for one Oscar ring must share their coefficient
# ring, so isos between them stay composable (see `_iso_with_codomain`).
const _coeff_ring_iso_cache = AbstractAlgebra.WeakKeyIdDict{Any, Any}()

function _cached_iso_oscar_singular_coeff_ring(R::Ring)
  # immutable rings (QQ, ZZ) cannot be weak keys; their Singular rings are cached anyway
  ismutable(R) || return iso_oscar_singular_coeff_ring(R)
  return get!(() -> iso_oscar_singular_coeff_ring(R), _coeff_ring_iso_cache, R)
end

#############################################################################
#
# Singular polynomial ring
#
struct OscarSingularPolyRingMap{D, C, W} <: Map{D, C, Any, Any}
  R::D
  S::C
  f::W
end

domain(f::OscarSingularPolyRingMap) = f.R

codomain(f::OscarSingularPolyRingMap) = f.S

(f::OscarSingularPolyRingMap)(x) = image(f, x)

_coeff_iso(f::OscarSingularPolyRingMap) = f.f

# some helper function shared with singlar_poly_ring

function _create_singular_poly_ring(S, Rx; keep_ordering::Bool = false)
  if keep_ordering
    Sx = Singular.polynomial_ring(S,
              _variables_for_singular(symbols(Rx)),
              ordering = internal_ordering(Rx),
              cached = false)[1]
  else
    Sx = Singular.polynomial_ring(S,
              _variables_for_singular(symbols(Rx)),
              cached = false)[1]
  end
  return Sx
end

function _create_singular_poly_ring(S, Rx, ord::Symbol)
  Sx =  Singular.polynomial_ring(S,
              _variables_for_singular(symbols(Rx)),
              ordering = ord,
              cached = false)[1]
  return Sx
end

function _create_singular_poly_ring(S, Rx, ord::Singular.sordering)
  Sx =  Singular.polynomial_ring(S,
              _variables_for_singular(symbols(Rx)),
              ordering = ord,
              cached = false)[1]
  return Sx
end

function _create_singular_poly_ring(S, Rx, ord::MonomialOrdering)
  Sx =  Singular.polynomial_ring(S,
              _variables_for_singular(symbols(Rx)),
              ordering = singular(ord),
              cached = false)[1]
  return Sx
end

function iso_oscar_singular_poly_ring(Rx::MPolyRing; keep_ordering::Bool = false)
  fcoeff = _cached_iso_oscar_singular_coeff_ring(base_ring(Rx))
  S = codomain(fcoeff)
  Sx = _create_singular_poly_ring(S, Rx; keep_ordering)
  return OscarSingularPolyRingMap(Rx, Sx, fcoeff)
end

function iso_oscar_singular_poly_ring(Rx::MPolyRing, ord::Union{Symbol, Singular.sordering, MonomialOrdering})
  fcoeff = _cached_iso_oscar_singular_coeff_ring(base_ring(Rx))
  S = codomain(fcoeff)
  Sx = _create_singular_poly_ring(S, Rx, ord)
  return OscarSingularPolyRingMap(Rx, Sx, fcoeff)
end

function image(f::OscarSingularPolyRingMap, a)
  parent(a) !== domain(f) && error("Element not in domain")
  g = MPolyBuildCtx(codomain(f))
  for (c, e) = zip(Nemo.coefficients(a), Nemo.exponent_vectors(a))
    push_term!(g, f.f(c), e)
  end
  return finish(g)
end

function preimage(f::OscarSingularPolyRingMap, a; check::Bool = true)
  check && (parent(a) === codomain(f) || error("Element not in codomain"))
  g = MPolyBuildCtx(domain(f))
  for (c, e) = Base.Iterators.zip(AbstractAlgebra.coefficients(a), AbstractAlgebra.exponent_vectors(a))
    push_term!(g, preimage(f.f, c), e)
  end
  return finish(g)
end

# Singular rings for the same Oscar ring differ only in their ordering and share
# the coefficient ring, so an iso onto one yields an iso onto the others.
function _iso_with_codomain(f::OscarSingularPolyRingMap, Sx::Singular.PolyRing)
  codomain(f) === Sx && return f
  @req base_ring(Sx) === base_ring(codomain(f)) "coefficient rings must agree"
  @req symbols(Sx) == symbols(codomain(f)) "variables must agree"
  return OscarSingularPolyRingMap(domain(f), Sx, f.f)
end

# iso from `R` onto `Sx`, which must have been created by `iso_oscar_singular_poly_ring(R, ...)`
function _iso_onto(R::MPolyRing, Sx::Singular.PolyRing)
  fcoeff = _cached_iso_oscar_singular_coeff_ring(coefficient_ring(R))
  @req base_ring(Sx) === codomain(fcoeff) "Singular ring was not created for this ring"
  @req ngens(Sx) == ngens(R) "number of variables must agree"
  return OscarSingularPolyRingMap(R, Sx, fcoeff)
end

#############################################################################
#
# Quotient rings
#
struct OscarSingularPolyRingQuoMap{D, C, W} <: Map{D, C, Any, Any}
  R::D
  S::C
  f::W #= isomorphism of the underlying "base rings" =#
end

domain(f::OscarSingularPolyRingQuoMap) = f.R

codomain(f::OscarSingularPolyRingQuoMap) = f.S

function (f::OscarSingularPolyRingQuoMap)(a)
  return image(f, a)
end

# cached, as it is needed for every conversion of an element of `R`
function _iso_oscar_singular_poly_ring(R::MPolyQuoRing)
  return get_attribute!(R, :iso_oscar_singular_poly_ring) do
    _groebner_basis(R)
    return OscarSingularPolyRingQuoMap(R, R.SQR, _iso_onto(base_ring(R), R.SQR))
  end::OscarSingularPolyRingQuoMap
end

function image(f::OscarSingularPolyRingQuoMap, b::MPolyQuoRingElem)
  @assert parent(b) === domain(f)
  a = b.f
  # we take a lift and map it into the singular quotient ring
  # by applying the coefficient map, which is f.f.f
  return map_coefficients(f.f.f, a; parent = codomain(f))
end

# for some reason this is used
function image(f::OscarSingularPolyRingQuoMap, b::MPolyRingElem)
  @assert parent(b) === base_ring(domain(f))
  return image(f, domain(f)(b))
end

function preimage(f::OscarSingularPolyRingQuoMap, a::Singular.spoly; check::Bool = true)
  check && (parent(a) === codomain(f) || error("Element not in codomain"))
  return domain(f)(preimage(f.f, a; check = false))
end

function _iso_with_codomain(f::OscarSingularPolyRingQuoMap, Sx::Singular.PolyRing)
  @req codomain(f) === Sx "quotient rings support only one Singular ring"
  return f
end

function _iso_onto(Q::MPolyQuoRing, Sx::Singular.PolyRing)
  return _iso_with_codomain(iso_oscar_singular_poly_ring(Q), Sx)
end

iso_oscar_singular_poly_ring(Q::MPolyQuoRing; keep_ordering::Bool = false) = _iso_oscar_singular_poly_ring(Q)
iso_oscar_singular_poly_ring(Q::MPolyQuoRing, ordering::MonomialOrdering) = _iso_oscar_singular_poly_ring(Q)
