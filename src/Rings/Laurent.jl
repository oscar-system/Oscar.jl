import AbstractAlgebra.Generic: LaurentMPolyWrapRing, LaurentMPolyWrap
import AbstractAlgebra: LaurentMPolyRing, LaurentMPolyRingElem

@attributes mutable struct LaurentMPolyAnyMap{D, C} <: Map{D, C, Map, LaurentMPolyAnyMap}
  R::D
  S::C
  image_of_gens
  data # map from _polyringquo(R) -> C
  
  function LaurentMPolyAnyMap(R::D, S::C, image_of_gens) where {D, C}
    @assert all(x -> parent(x) === S, image_of_gens)
    return new{D, C}(R, S, image_of_gens)
  end
end

mutable struct LaurentMPolyIdeal{T} <: Ideal{T}
  R
  gens::Vector{T}
  data # ideal of of _polyringquo(R)

  function LaurentMPolyIdeal(R::LaurentMPolyRing, gens::Vector)
    @assert all(x -> parent(x) === R, gens)
    return new{eltype(gens)}(R, gens)
  end
end

mutable struct _LaurentMPolyBackend{D, C, M}
  R::D
  Q::C
  inv::M
  _gens_cache

  function _LaurentMPolyBackend(R::D, Q::C) where {D, C}
    _inv = hom(Q, R, vcat(gens(R), inv.(gens(R))))
    return new{D, C, typeof(_inv)}(R, Q, _inv)
  end
end

# Ideals and maps for multivariate Laurent polynomial rings
# 
# Tommy Hofmann:
# If R = K[x1^+-,...,xn^+-] is a such a ring, then this is isomorphic as a K-algebra
# to the affine algebra Raff = K[x1,...,xn,y1,...,yn]/(x1y1 - 1,...,xnyn - 1)
# We piggyback on quotient rings for all serious computations. More precisely
#
# g = _polyringquo(R)
#
# constructs a map g : R -> Raff, which supports
# - image(g, x)/g(x) for x an element or ideal of R or a map R -> ?
# - preimage(g, x) for x an element or ideal of Raff
#
# In this way we can move everything to quotient rings

################################################################################
#
#  Conversions to quotient rings
#
################################################################################

function _polyringquo(R::LaurentMPolyWrapRing)
  get_attribute!(R, :polyring) do
    n = nvars(R)
    C = base_ring(R)
    Cx, x, xinv = polynomial_ring(C,"x#" => 1:n, "x#^-1" => 1:n; cached = false)
    I = ideal(Cx, [x[i]*xinv[i] - 1 for i in 1:n])
    Q, = quo(Cx, I)
    return _LaurentMPolyBackend(R, Q)
  end::_LaurentMPolyBackend  # TODO: make the type fully concrete
end

function _evaluate_gens_cache(f::_LaurentMPolyBackend{D, C}) where {D, C}
  if !isdefined(f, :_gens_cache)
    f._gens_cache = gens(codomain(f))[1:nvars(domain(f))]
  end
  return f._gens_cache::Vector{elem_type(C)}
end

domain(f::_LaurentMPolyBackend) = f.R

codomain(f::_LaurentMPolyBackend) = f.Q

(f::_LaurentMPolyBackend)(p) = image(f, p)

# conversion of elements
# computes an exponent vector e (possible negative exponents) and a polynomial q
# such that p = x^e * q
# note that q is in the underlying multivariate polynomial ring
# in particular <p> = <q> in the Laurent polynomial ring
_split(p::LaurentMPolyWrap) = p.mindegs, p.mpoly

function image(f::_LaurentMPolyBackend, p::LaurentMPolyWrap)
  @assert parent(p) === domain(f)
  Q = codomain(f)
  n = nvars(domain(f))
  _gens = _evaluate_gens_cache(f)
  # I try to be a bit clever here
  exps, poly = _split(p)
  r = evaluate(poly, _gens)
  v = zeros(Int, ngens(Q))
  for i in 1:length(exps)
    if exps[i] >= 0
      v[i] = exps[i]
    else
      v[n + i] = -exps[i]
    end
  end
  m = Q(monomial(base_ring(Q), v))
  # TODO:
  # One can remove the evaluate with some more tricks, since this is just
  # K[x1,...xn] -> K[x1,...,xn,y1,...yn] -> K[x,y]/(xy - 1)
  return m * evaluate(poly, _gens)
end

function preimage(f::_LaurentMPolyBackend, q::MPolyQuoRingElem)
  @assert parent(q) === codomain(f)
  return f.inv(q)
end

# conversion of ideals
function image(f::_LaurentMPolyBackend{D, C, M}, I::LaurentMPolyIdeal) where {D, C, M}
  @assert base_ring(I) === domain(f)
  if isdefined(I, :data)
    return I.data::ideal_type(C)
  else
    I.data = ideal(codomain(f), f.(gens(I)))
    return I.data::ideal_type(C)
  end
end

function preimage(f::_LaurentMPolyBackend, J::MPolyQuoIdeal)
  @assert base_ring(J) === codomain(f)
  I = ideal(domain(f), map(x -> preimage(f, x), gens(J)))
  I.data = J
  return I
end

# conversion for maps
function image(f::_LaurentMPolyBackend, g::LaurentMPolyAnyMap)
  Q = codomain(f)
  _images = (g.image_of_gens)
  __images = append!(copy(_images), inv.(_images))
  return hom(Q, codomain(g), __images, check = false)
end

################################################################################
#
#  Maps from Laurent polynomial rings
#
################################################################################

domain(f::LaurentMPolyAnyMap) = f.R

codomain(f::LaurentMPolyAnyMap) = f.S

function hom(R::LaurentMPolyRing, S::Ring, images::Vector; check::Bool = true)
  @req length(images) == nvars(R) "Wrong number of images"
  if check
    @req all(is_unit, images) "Images of generators must be units"
  end
  _images = S.(images)
  return LaurentMPolyAnyMap(R, S, _images)
end

function image(f::LaurentMPolyAnyMap{D, C}, g::LaurentMPolyRingElem) where {D, C}
  @req parent(g) === domain(f) "Element not in domain of map"
  return evaluate(g, f.image_of_gens::Vector{elem_type(C)})
end

function preimage(f::LaurentMPolyAnyMap, g)
  @req parent(g) === codomain(f) "Element not in domain of map"
  q = _polyringquo(domain(f))
  qf = q(f)
  return preimage(q, preimage(qf, g))
end

(f::LaurentMPolyAnyMap)(g::LaurentMPolyRingElem) = image(f, g)

# preimage for ideals
function preimage(f::MPolyAnyMap{X, <: LaurentMPolyRing, Y, Z}, I::LaurentMPolyIdeal) where {X<:Union{MPolyRing, MPolyQuoRing}, Y, Z}
  R = domain(f)
  S = codomain(f)
  if coefficient_ring(R) === base_ring(S) && f.(gens(R)) == gens(S)
    # We try to do something clever if f : K[x1,...xn] -> K[x1^+-,...xn^+-] is
    # the natural inclusion
    polygens = map(x -> _split(x)[2], gens(I))
    # These are polynomials generating the same ideal as I
    II = ideal(polygens) # this is an element of the internal polynomial ring underpinning
                         # the Laurent polynomial ring R
    _R = base_ring(II)
    for g in gens(_R)
      II = Oscar._saturation2(II, g*_R)
    end
    # We need to translate to an ideal of R
    return ideal(R, map(x -> map_coefficients(identity, x, parent = R), gens(II)))
  end
  q = _polyringquo(codomain(f))
  qI = q(I) # ideal in the quotient
  # map from domain(f) to quotient
  qf = hom(domain(f), codomain(q), [q(f(g)) for g in gens(domain(f))])
  return preimage(qf, qI)
end

################################################################################
#
#  Ideals
#
################################################################################

base_ring(I::LaurentMPolyIdeal{T}) where {T} = I.R::parent_type(T)

base_ring_type(::Type{LaurentMPolyIdeal{T}}) where {T} = parent_type(T)

gens(I::LaurentMPolyIdeal) = I.gens

@enable_all_show_via_expressify LaurentMPolyIdeal

function AbstractAlgebra.expressify(a::LaurentMPolyIdeal; context = nothing)
  return Expr(:call, :ideal, [AbstractAlgebra.expressify(g, context = context) for g in gens(a)]...)
end

function ideal(R::LaurentMPolyRing, x::Vector)
  return LaurentMPolyIdeal(R, filter!(!iszero, R.(x)))
end

function in(x::LaurentMPolyRingElem, I::LaurentMPolyIdeal)
  R = parent(x)
  if parent(x) !== base_ring(I)
    return false
  end
  f = _polyringquo(R)
  return f(x) in f(I)
end

function +(I::LaurentMPolyIdeal, J::LaurentMPolyIdeal)
  @req base_ring(I) === base_ring(J) "Rings must be equal"
  R = base_ring(I)
  f = _polyringquo(R)
  IpJ = preimage(f, f(I) + f(J))
  return IpJ
end



function quo(R::Oscar.LaurentMPolyWrapRing, I::Oscar.LaurentMPolyIdeal)
  @req R === base_ring(I) "ring and ideal do not match"
  poly_repr = Oscar._polyringquo(R);
  underlying_poly_ring_quo = codomain(poly_repr);
  II = ideal(underlying_poly_ring_quo, poly_repr.(gens(I)));
  Q,phi = quo(underlying_poly_ring_quo, II);
  map_down(f::AbstractAlgebra.Generic.LaurentMPolyWrap) = phi(poly_repr(f));
  lift_up(f::MPolyQuoRingElem) = poly_repr.inv(preimage(phi,f));
  return Q, MapFromFunc(R,Q, map_down, lift_up);
end
