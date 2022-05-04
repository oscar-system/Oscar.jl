module GaloisGrp

using Oscar, Markdown
import Base: ^, +, -, *, ==
import Oscar: Hecke, AbstractAlgebra, GAP
using Oscar: SLPolyRing, SLPoly, SLPolynomialRing, CycleType

export galois_group, slpoly_ring, elementary_symmetric, galois_quotient,
       power_sum, to_elementary_symmetric, cauchy_ideal, galois_ideal, fixed_field,
       extension_field

import Hecke: orbit, fixed_field


function __init__()
  GAP.Packages.load("ferret"; install=true)

  Hecke.add_verbose_scope(:GaloisGroup)
  Hecke.add_verbose_scope(:GaloisInvariant)
  Hecke.add_assert_scope(:GaloisInvariant)
end

"""
A poor mans version of a (range of) tropical rings.

Actually, not even a ring. Defined in terms of operations
 - `mul` used for multiplication
 - `add` for addition
 - `pow` for powering with Int exponent
 - `map` to create elements from (large) integers/ fmpz
 - `name` is only used for printing.
"""
struct BoundRing{T}  <: AbstractAlgebra.Ring
  mul#::(T,T) -> T
  add#::(T,T) -> T
  pow#::(T, Int) -> T
  map#:: R -> T
  name::String
end

function Base.show(io::IO, b::BoundRing{T}) where {T}
  print(io, "$(b.name) for type $T")
end

struct BoundRingElem{T} <: AbstractAlgebra.RingElem
  val::T
  p::BoundRing # the parent
end

function Base.show(io::IO, b::BoundRingElem)
  print(io, "(x <= $(b.val))")
end

function check_parent(a::BoundRingElem, b::BoundRingElem)
  parent(a) == parent(b) || error("Elements must have same parent")
  return true
end

Base.:(==)(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && a.val == b.val
function +(a::BoundRingElem, b::BoundRingElem) 
  check_parent(a, b) 
  c = BoundRingElem(a.p.add(a.val, b.val), a.p)
#  @show a, "+", b, "=", c
  return c
end
-(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && BoundRingElem(a.p.add(a.val, b.val), a.p)
function *(a::BoundRingElem, b::BoundRingElem) 
  check_parent(a, b) 
  c = BoundRingElem(a.p.mul(a.val, b.val), a.p)
#  @show a, "*", b, ":=", c
  return c
end

function *(a::fmpz, b::BoundRingElem) 
  c = BoundRingElem(b.p.mul(b.p.map(a), b.val), b.p)
#  @show a, ":*", b, ":=", c
  return c
end
function ^(a::BoundRingElem, b::Int) 
  c = BoundRingElem(a.p.pow(a.val, b), a.p)
#  @show a, ":^", b, ":=", c
  return c
end
-(a::BoundRingElem) = a

Oscar.parent(a::BoundRingElem) = a.p
value(a::BoundRingElem) = a.val
Base.isless(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && isless(value(a), value(b))

Oscar.parent_type(::BoundRingElem{T}) where T = BoundRing{T}
Oscar.elem_type(::BoundRing{T}) where T = BoundRingElem{T}
Oscar.parent_type(::Type{BoundRingElem{T}}) where T = BoundRing{T}
Oscar.elem_type(::Type{BoundRing{T}}) where T = BoundRingElem{T}

(R::BoundRing{T})(a::fmpz) where T = BoundRingElem{T}(R.map(abs(a)), R)
(R::BoundRing{T})(a::Integer) where T = BoundRingElem{T}(fmpz(a), R)
(R::BoundRing{T})() where T = BoundRingElem{T}(fmpz(0), R)
(R::BoundRing{T})(a::T) where T = BoundRingElem{T}(a, R)
(R::BoundRing{T})(a::BoundRingElem{T}) where T = a
Oscar.one(R::BoundRing) = R(1)
Oscar.zero(R::BoundRing) = R(0)


"""
Intended to be used to get upper bounds under evaluation in function fields (non 
  Archimedean valuaton, degree)

Operations are: 
 - `a+b := max(a, b)`
 - `ab  := a+b`
"""
function max_ring()
  return BoundRing{fmpz}( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x, x->x, "max-ring")
end

"""
Normal ring
"""
function add_ring(;type::Type=fmpz)
  return BoundRing{type}( (x,y) -> x*y, (x,y) -> x+y, (x,y) -> x^y, x->abs(x), "add-ring")
end

#roots rt are power series sum a_n x^n
#we have |a_n| <= r^-n B (n+1)^k for B = x[1], k = x[2]
#and deg(rt) <= x[3] (infinite valuation)
function qt_ring()
  return BoundRing{Tuple{fmpz, Int, fmpq}}( (x,y) -> (x[1]*y[1], x[2]+y[2]+1, x[3]+y[3]),
                                            (x,y) -> (x[1]+y[1], max(x[2], y[2]), max(x[3], y[3])),
                                            (x,y) -> (x[1]^y, y*x[2]+y-1, y*x[3]),
                                                x -> _coerce_qt(x), "qt-ring")
end

_coerce_qt(x::fmpz) = (abs(x), 0, fmpq(0))
_coerce_qt(x::Integer) = (fmpz(x), 0, fmpq(0))
function _coerce_qt(x::fmpz_poly)
  if iszero(x)
    return (fmpz(0), 0, fmpq(0))
  end
  return (maximum(abs, coefficients(x))*(degree(x)+1), 0, fmpq(0))
end

(R::BoundRing{Tuple{fmpz, Int, fmpq}})(a::Tuple{fmpz, Int, fmpq}) = BoundRingElem(a, R)
(R::BoundRing{Tuple{fmpz, Int, fmpq}})(a::Integer) = BoundRingElem(R.map(a), R)
(R::BoundRing{Tuple{fmpz, Int, fmpq}})(a::fmpz) = BoundRingElem(R.map(a), R)
(R::BoundRing{Tuple{fmpz, Int, fmpq}})(a::fmpz_poly) = BoundRingElem(R.map(a), R)

"""
An slpoly evaluated at `cost_ring` elements `0` will count the number
of multiplications involved. A measure of the cost of evaluation at more
interesting scalars.

Operations:
 - `xy := x+y+1`
 - `x+y := x+y`
 - `x^y := x+2*log_2(y)`
 - all constants are mapped to `0`
"""
function cost_ring()
  return BoundRing{fmpz}( (x,y) -> x+y+1, (x,y) -> x+y, (x,y) -> x+2*nbits(y), x->0, "cost-ring")
end

"""
An slpoly evaluated at `degree_ring` elements `1` will bound the total degree
from above.

Operations:
 - `xy := x+y`
 - `x+y := max(x,y)`
 - `x^y := yx`
 - all constants are mapped to `0`
"""
function degree_ring()
  return BoundRing{fmpz}( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x, x->0, "degree-ring")
end

@doc Markdown.doc"""
    cost(I::SLPoly)

Counts the number of multiplications to evaluate `I`, optionally
a Tschirnhaus transformation (`fmpz_poly`) can be passed in as well.
"""
function cost(I::SLPoly)
  n = ngens(parent(I))
  C = cost_ring()
  return value(evaluate(I, [C(0) for i = 1:n]))
end
function cost(I::SLPoly, ts::fmpz_poly)
  n = ngens(parent(I))
  C = cost_ring()
  return value(evaluate(I, [C(0) for i = 1:n]))+n*degree(ts)
end

@doc Markdown.doc"""
    total_degree(I::SLPoly)

Determines an upper bound for the total degree of `I`.
"""
function total_degree(I::SLPoly)
  n = ngens(parent(I))
  C = degree_ring()
  return value(evaluate(I, [C(1) for i = 1:n]))
end

Oscar.mul!(a::BoundRingElem, b::BoundRingElem, c::BoundRingElem) = b*c
Oscar.addeq!(a::BoundRingElem, b::BoundRingElem) = a+b

#my 1st invariant!!!
@doc Markdown.doc"""
    sqrt_disc(a::Vector)

The product of differences ``a[i] - a[j]`` for all indices ``i<j``.    
"""
function sqrt_disc(a::Vector)
  if length(a) == 1
    return one(parent(a[1]))
  end
  return prod([a[i] - a[j] for i = 1:length(a)-1 for j = i+1:length(a)])
end

@doc Markdown.doc"""
    elementary_symmetric(g::Vector, i::Int) -> 

Evaluates the `i`-th elementary symmetric polynomial at the values in `g`.    
The `i`-th elementary symmetric polynomial is the sum over all
products of `i` distinct variables.
"""
function elementary_symmetric(g::Vector, i::Int)
  return sum(prod(g[i] for i = s) for s = Hecke.subsets(Set(1:length(g)), i))
end

@doc Markdown.doc"""
    power_sums(g::Vector, i::Int) ->

Evaluates the `i`-th power sums at the values in `g`, ie. the sum
of the `i`-th power of the values.
"""
function power_sum(g::Vector, i::Int)
  return sum(a^i for a = g)
end

@doc Markdown.doc"""
    discriminant(g::Vector)

Compute the product of all differences of distinct elements in the array.    
"""
function Oscar.discriminant(g::Vector{<:RingElem})
  return prod(a-b for a = g for b = g if a!=b)
end


function slpoly_ring(R::AbstractAlgebra.Ring, n::Int; cached ::Bool = false)
  return SLPolynomialRing(R, [ Symbol("x_$i") for i=1:n], cached = cached)
end

function slpoly_ring(R::AbstractAlgebra.Ring, p::Pair{Symbol, UnitRange{Int}}...; cached::Bool = false)
  return SLPolynomialRing(R, p..., cached = cached)
end

function (R::SLPolyRing)(a::SLPoly)
  parent(a) == R && return a
  error("wrong parent")
end

@doc Markdown.doc"""
    roots_upper_bound(f::fmpz_poly) -> fmpz

An upper upper_bound for the absolute value of the complex roots of the input.    
Uses the Cauchy bound.
"""
function Nemo.roots_upper_bound(f::fmpz_poly)
  a = coeff(f, degree(f))
  b = ceil(fmpz, abs(coeff(f, degree(f)-1)//a))
  for i=0:degree(f)-2
    b = max(b, root(ceil(fmpz, abs(coeff(f, i)//a)), degree(f)-i)+1)
  end
  return 2*b
  return max(fmpz(1), maximum([ceil(fmpz, abs(coeff(f, i)//a)) for i=0:degree(f)]))
end
function Nemo.roots_upper_bound(f::fmpq_poly)
  a = coeff(f, degree(f))
  return max(fmpz(1), maximum([ceil(fmpz, abs(coeff(f, i)//a)) for i=0:degree(f)]))
end

#roots are sums of m distinct roots of f
#from https://doi.org/10.2307/2153295
#Symmetric Functions, m-Sets, and Galois Groups
#by David Casperson and John McKay
@doc Markdown.doc"""
    msum_poly(f::PolyElem, m::Int)

Compute the polynomial with roots sums of `m` roots of `f` using
resultants.
"""
function msum_poly(f::PolyElem, m::Int)
  f = divexact(f, leading_coefficient(f))
  N = binomial(degree(f), m)
  p = Hecke.polynomial_to_power_sums(f, N)
  p = vcat([degree(f)*one(base_ring(f))], p)
  S, a = PowerSeriesRing(base_ring(f), N+1, "a")
  Hfs = S([p[i]//factorial(fmpz(i-1)) for i=1:length(p)], N+1, N+1, 0)
  H = [S(1), Hfs]
  for i=2:m
    push!(H, 1//i*sum((-1)^(h+1)*Hfs(h*a)*H[i-h+1] for h=1:i))
  end
  p = [coeff(H[end], i)*factorial(fmpz(i)) for i=0:N]
  return Hecke.power_sums_to_polynomial(p[2:end])
end

@doc Markdown.doc"""
A `GaloisCtx`, is the context object used for the computation of Galois
groups of (univariate) polynomials. It contains
 - the polynomial
 - an object that can compute the roots in a fixed order up to a given 
   precision. Currently, this is a q-adic field, that is an unramified
   extension of the p-adic numbers.
 - a upper bound on the _size_ of the roots, currently an upper bound on the
   complex absolute value.
 - at the end, the Galois group  

This is constructed implicitly while computing a Galois group and returned
together with the group.

Not type stable, not sure what to do about it:
Depends on type of
 - `f`
However, `f` is hardly ever used. 
"""
mutable struct GaloisCtx{T}
  f::PolyElem
  C::T # a suitable root context
  B::BoundRingElem # a "bound" on the roots, might be "anything"
  G::PermGroup
  rt_num::Dict{Int, Int}
  chn::Vector{Tuple{PermGroup, SLPoly, fmpz_poly, Vector{PermGroupElem}}}
  start::Tuple{Int, Vector{Vector{Vector{Int}}}} # data for the starting group:
    # if start[1] == 1: start[2] is a list of the block systems used
    #             == 2  start[2] is a list of the orbits of the msum-poly, ie. a list of a list of pairs
    #                            where the pairs are Vector{Int} of length 2
  data::Any #whatever else is needed in special cases
  #= the descent chain, recording
   - the group
   - the invariant
   - the tschirnhaus transformation
   - the cosets used
   should probably also record if the step was proven or not
    the starting group  and the block systems used to get them
  =#
  prime::Any #=can be
   - fmpz/ Int: prime number, used over Q
   - NfOrdIdl : prime ideal , used over NfAbs
   - (Int, Int): evaluation point, prime number used over Q(t)
   =#

  function GaloisCtx(f::fmpz_poly, p::Int)
    r = new{Hecke.qAdicRootCtx}()
    r.prime = p
    r.f = f
    r.C = Hecke.qAdicRootCtx(f, p, splitting_field = true)
    r.B = add_ring()(leading_coefficient(f)*roots_upper_bound(f))
    r.chn = Tuple{PermGroup, SLPoly, fmpz_poly, Vector{PermGroupElem}}[]
    return r
  end

  function GaloisCtx(f::fmpq_poly, p::Int)
    d = mapreduce(denominator, lcm, coefficients(f))
    return GaloisCtx(Hecke.Globals.Zx(d*f), p)
  end
  #=
  Roots in F_q[[t]] for q = p^d
    - needs to be tweaked to do Q_q[[t]
    - need q-lifting as well as t-lifting
    - possibly also mul_ks for Q_q[[t]] case
    - needs merging in Hecke
  =#

  function GaloisCtx(f::fmpz_mpoly, shft::Int, p::Int, d::Int)
    f = evaluate(f, [gen(parent(f), 1), gen(parent(f), 2)+shft])
    #f(x, T+t), the roots are power series in T over qAdic(p, d)
    #so basically for f in Qq<<(T+t)>> 
    @assert ngens(parent(f)) == 2
    Qq, _ = QadicField(p, d, 10)
    F, mF = ResidueField(Qq)
    H = Hecke.MPolyFact.HenselCtxFqRelSeries(f, F)
    SQq, _ = PowerSeriesRing(Qq, 2, "s", cached = false)
    SQqt, _ = PolynomialRing(SQq, "t", cached = false)
    mc(f) = map_coefficients(x->map_coefficients(y->setprecision(preimage(mF, y), 1), x, parent = SQq), f, parent = SQqt)
    HQ = Hecke.MPolyFact.HenselCtxFqRelSeries(H.f, map(mc, H.lf), map(mc, H.cf), H.n)
    r = new{Hecke.MPolyFact.HenselCtxFqRelSeries{AbstractAlgebra.Generic.RelSeries{qadic}}}()
    r.prime = (shft, p)
    Qt, t = RationalFunctionField(QQ, "t", cached = false)
    Qts, s = PolynomialRing(Qt, "s", cached = false)
    r.f = evaluate(f, [s, Qts(t)])
    r.C = HQ
    r.chn = Tuple{PermGroup, SLPoly, fmpz_poly, Vector{PermGroupElem}}[]
    
    vl = roots_upper_bound(f)
    r.B = qt_ring()(vl[1])
    r.data = Any[vl[2], shft, false] #false: not simulating C
    @assert typeof(r.data[3]) == Bool
    return r
  end
  function GaloisCtx(T::Type)
    r = new{T}()
    r.chn = Tuple{PermGroup, SLPoly, fmpz_poly, Vector{PermGroupElem}}[]
    return r
  end
end

Base.round(::Type{Int}, q::fmpq) = Int(round(fmpz, q))

function Oscar.prime(C::GaloisCtx{Hecke.MPolyFact.HenselCtxFqRelSeries{Generic.RelSeries{qadic}}})
  return prime(base_ring(base_ring(C.C.lf[1])))
end

function bound_to_precision(G::GaloisCtx{T}, B::BoundRingElem{Tuple{fmpz, Int, fmpq}}, extra=(0, 0)) where {T}
  if isa(extra, Int)
    extra = (extra, min(2, div(extra, 3)))
  end
  C, k, d = B.val
  r = G.data[1]
  #so power series prec need to be floor(Int, d)
  n = floor(fmpz, d+1)
  #padic: we ne |a_i| for i=0:n and |a_i| <= C (i+1)^k/r^i
  #and then log_p()
  #according to the Qt file, a_i is maximal around k/log(r) -1
  if G.data[3] #we're simulating CC, so we don't care about the p-adics (too much)
    b = max(C, fmpz(10)^10)
  elseif isone(r)
    b = C*(n+1)^k
  else
    c = max(1, floor(Int, k/log(r)-1))
    if n<c
      b = C*(n+1)^k//r^n
    else
      b = max(C*(c+1)^k//r^c, C*(c+2)^k//r^(c+1))
    end
  end
  b = max(b, fmpz(1))
  return (clog(floor(fmpz, b), prime(G))+extra[1], Int(n)+extra[2])
end

function bound_to_precision(G::GaloisCtx{T}, B::BoundRingElem{fmpz}, extra::Int = 5) where {T}
  return clog(B.val, G.C.p)+extra
end

function Nemo.roots_upper_bound(f::fmpz_mpoly, t::Int = 0)
  @assert nvars(parent(f)) == 2
  Qs, s = RationalFunctionField(FlintQQ, "t", cached = false)
  Qsx, x = PolynomialRing(Qs, cached = false)
  F = evaluate(f, [x, Qsx(s)])
  dis = numerator(discriminant(F))
  @assert !iszero(dis(t))
  rt = roots(dis, ComplexField(20))
  r = Hecke.lower_bound(minimum([abs(x-t) for x = rt]), fmpz)
  @assert r > 0
  ff = map_coefficients(abs, f)
  C = roots_upper_bound(Hecke.Globals.Zx(map(x->evaluate(x, fmpz[r, 0]), coefficients(ff, 2))))
  C1 = maximum(map(x->evaluate(x, fmpz[r, 0]), coefficients(ff, 2)))
  #the infinite valuation... need Newton
  vl = valuations_of_roots(F)
  return (C+1, 0, maximum(x[1] for x = vl)), r
end

function Base.show(io::IO, GC::GaloisCtx{Hecke.qAdicRootCtx})
  print(io, "Galois Context for $(GC.f) and prime $(GC.C.p)")
end
function Base.show(io::IO, GC::GaloisCtx{<:Hecke.MPolyFact.HenselCtxFqRelSeries})
  print(io, "Galois Context for $(GC.f)")
end


#TODO: change pr to be a "bound_ring_elem": in the Qt case this has to handle
#      both power series prec as well as q-adic...
@doc Markdown.doc"""
    roots(G::GaloisCtx, pr::Int)

The roots of the polynomial used to define the Galois-context in the fixed order
used in the algorithm. The roots are returned up to a precision of `pr`
p-adic digits, thus they are correct modulo ``p^pr``

For non-monic polynomials they roots are scaled by the leading coefficient.
If `raw` is set to true, the scaling is omitted.
The bound in the `GaloisCtx` is also adjusted.
"""
function Hecke.roots(G::GaloisCtx{Hecke.qAdicRootCtx}, pr::Int=5; raw::Bool = false)
  a = Hecke.roots(G.C, pr)::Vector{qadic}
  b = Hecke.expand(a, all = true, flat = false, degs = Hecke.degrees(G.C.H))::Vector{qadic}
  if isdefined(G, :rt_num)
    b = [b[G.rt_num[i]] for i=1:length(G.rt_num)]
  end
  if raw
    return b
  else
    return leading_coefficient(G.f) .* b
  end
end

function Hecke.setprecision(a::Generic.RelSeries, p::Int)
  b = parent(a)(a.coeffs, min(length(a.coeffs), p), p+valuation(a), valuation(a))
end

#TODO: does not really work: the polynomial is not truncated, only
#      the entry for precision is updated: 1 + ... + s^40 + O(s^2)
function Hecke.setprecision!(G::GaloisCtx{<:Hecke.MPolyFact.HenselCtxFqRelSeries}, pr::Tuple{Int, Int})
  c = coeff(G.C.lf[1], 0)
  @assert precision(c) >= pr[2]
  @assert precision(coeff(c, 0)) >= pr[1]
  for f = G.C.lf
    Hecke.set_precision!(f, pr[2])
  end
  for f = G.C.cf  
    Hecke.set_precision!(f, pr[2])
  end
end

function Hecke.roots(G::GaloisCtx{<:Hecke.MPolyFact.HenselCtxFqRelSeries}, pr::Tuple{Int, Int} = (5, 2); raw::Bool = false)
  C = G.C
  while precision(C)[1] < pr[1]
    Hecke.MPolyFact.lift_q(C)
  end
  #TODO: truncate precision where neccessary, working with an insane
  #      series precision is costly
  while precision(C)[2] < pr[2]
    Hecke.MPolyFact.lift(C)
  end
  rt = [-coeff(x, 0) for x = C.lf[1:C.n]]
  rt = map(y->map_coefficients(x->setprecision(x, pr[1]), setprecision(y, pr[2]), parent = parent(y)), rt)
  if isdefined(G, :rt_num)
    rt = [rt[G.rt_num[i]] for i=1:length(G.rt_num)]
  end
  return rt
end

@doc Markdown.doc"""
    upper_bound(G::GaloisCtx, f...)

Given a `GaloisCtx` and some multivariate function, upper_bound the image of `f`
upon evaluation at the roots implicit in `G`.

`f` can be
 - a multivariate polynomial or straight-line polynomial (strictly: any object
   allowing `evaluate`
 - `elementary_symmetric` or `power_sum`, in which case more arguments are
   needed: the array with the values and the index.
   `upper_bound(G, power_sum, A, i)` is equivalent to `upper_bound(G, power_sum(A, i))`
   but more efficient.

In every case a univariate polynomial (over the integers) can be added, it
will act as a Tschirnhaus-transformation, ie. the roots (bounds) implicit
in `G` will first be transformed.
"""
function upper_bound end

function upper_bound(G::GaloisCtx, f)
#  @show :eval, f, G.B, degree(G.f)
  return Oscar.evaluate(f, [G.B for i=1:degree(G.f)])
end

function upper_bound(G::GaloisCtx, f, ts::fmpz_poly)
  B = ts(G.B)
  return Oscar.evaluate(f, [B for i=1:degree(G.f)])
end

function upper_bound(G::GaloisCtx, ::typeof(elementary_symmetric), A::Vector, i::Int, ts::fmpz_poly = gen(Oscar.Hecke.Globals.Zx))
  if ts != gen(Hecke.Globals.Zx)
    A = [ts(x) for x = A]
  end
  B = [upper_bound(G, x) for x = A]
  n = length(B)
  b = sort(B)
  return parent(B[1])(binomial(n, i))*prod(b[max(1, n-i+1):end])
end

function upper_bound(G::GaloisCtx, ::typeof(power_sum), A::Vector, i::Int, ts::fmpz_poly = gen(Oscar.Hecke.Globals.Zx))
  if ts != gen(Hecke.Globals.Zx)
    A = [ts(x) for x = A]
  end
  B = [upper_bound(G, x)^i for x = A]
  return sum(B)
end

function upper_bound(G::GaloisCtx, ::typeof(elementary_symmetric), i::Int, ts::fmpz_poly = gen(Oscar.Hecke.Globals.Zx))
  if ts != gen(Hecke.Globals.Zx)
    b = ts(G.B)
  else
    b = G.B
  end
  n = degree(G.f)
  return parent(b)(binomial(n, i))*b^i
end

function upper_bound(G::GaloisCtx, ::typeof(power_sum), i::Int, ts::fmpz_poly = gen(Oscar.Hecke.Globals.Zx))
  if ts != gen(Hecke.Globals.Zx)
    b = ts(G.B)
  else
    b = G.B
  end
  return parent(b)(degree(G.f))*b^i
end

function Hecke.orbit(G::Oscar.PermGroup, f::MPolyElem)
  s = Set([f])
  while true
    n = Set(x^g for x = s for g = gens(G))
    sn = length(s)
    union!(s, n)
    if length(s) == sn
      break
    end
  end
  return s
end

function Hecke.evaluate(I::SLPoly, p, a::Vector)
  return evaluate(I, [a[p(i)] for i=1:length(a)])
end

probable_orbit(G::Oscar.PermGroup, f::MPolyElem) = orbit(G, f)

"""
`slprogram`s can be compiled into "normal" Julia functions, but there is
some overhead in the compilation itself. By default, aparently nothing is
compiled, so we allow to force this there.

`isPoly` allows the use of inplace operations, as `SLPoly`s result
in programs where intermediate results are used only once.
"""
function compile!(f::SLPoly)
  if !isdefined(f.slprogram, :f)
    Oscar.StraightLinePrograms.compile!(f.slprogram, isPoly = true)
  end
end

#one cannot compare (==) slpoly, no hash either..
#(cannot be done, thus comparison is indirect via evaluation)
#I assume algo can be improved (TODO: Max?)
function probable_orbit(G::Oscar.PermGroup, f::SLPoly; limit::Int = typemax(Int))
  n = ngens(parent(f))
  F = GF(next_prime(2^50))
  p = [rand(F) for i=1:n]
  s = [f]
  v = Set([evaluate(f, p)])
  while true
    nw = []
    for g = gens(G)
      for h = s
        z = evaluate(h^g, p)
        if !(z in v)
          push!(nw, h^g)
          push!(v, z)
          if length(s) + length(nw) > limit
            append!(s, nw)
            return s
          end
        end
      end
    end
    if length(nw) == 0
      return s
    end
    append!(s, nw)
  end
end

#TODO:
#- Bessere Abstraktion um mehr Grundkoerper/ Ringe zu erlauben
#- Bessere Teilkpoerper: ich brauche "nur" maximale
#- sanity-checks
#- "datenbank" fuer Beispiele

#a gimmick, not used in galois groups
@doc Markdown.doc"""
    to_elementary_symmetric(f)

For a multivariate symmetric polynomial `f`, (i.e. `f` is invariant under
permutation of the variables), compute a new polynomial `g` s.th.
`g` evaluated at the elementary symmetric polynomials recovers `f`.

This is using a rather elementary algorithm.
# Example
We recover the Newton-Girard formulas:
```julia-repl
julia> R, x = PolynomialRing(QQ, 3);
julia> d = power_sum(x, 3)
x1^3 + x2^3 + x3^3
julia> g = to_elementary_symmetric(d)
x1^3 - 3*x1*x2 + 3*x3
julia> evaluate(g, [elementary_symmetric(x, i) for i=1:3])
x1^3 + x2^3 + x3^3
```
"""
function to_elementary_symmetric(f)
  S = parent(f)
  n = ngens(S)
  if n == 1 || isconstant(f)
    return f
  end
  T = PolynomialRing(base_ring(S), n-1)[1]
  g1 = to_elementary_symmetric(evaluate(f, vcat(gens(T), [T(0)])))
  es = [elementary_symmetric(gens(S), i) for i=1:n-1]
  f = f - evaluate(g1, es)
  h = divexact(f, elementary_symmetric(gens(S), n))
  g2 = to_elementary_symmetric(h)
  g1 = evaluate(g1, gens(S)[1:n-1])
  return g1 + gen(S, n)*g2
end

function ^(f::SLPoly, p::Oscar.PermGroupElem)
  #TODO: replace by makeing the permutation of the input an internal
  #      operation.

  g = gens(parent(f))
  h = typeof(f)[]
  for i=1:ngens(parent(f))
    push!(h, g[p(i)])
  end
  e = evaluate(f, h)
  if typeof(e) != typeof(f)
    @show "bad case"
    return f
  end
  return e
end

@doc Markdown.doc"""
    isprobably_invariant(g, p) -> Bool

For a multivariate function, mainly an `SLPoly`, test if this is
likely to be invariant under the permutation `p`. Due to the representation
of `SLPoly`s as trees, it is not possible to test this exactly. Instead
`p` is evaluated at random elements in a large finite field.
"""
function isprobably_invariant(g, p)
  R = parent(g)
  k = GF(next_prime(2^20))
  n = ngens(R)
  lp = [rand(k) for i=1:n]
  return evaluate(g, lp) == evaluate(g^p, lp)
end
#TODO: think about the order of arguments!

function isprobably_invariant(p, G::PermGroup)
  R = parent(p)
  k = GF(next_prime(2^20))
  n = ngens(R)
  lp = [rand(k) for i=1:n]
  gp = evaluate(p, lp)
  return all(x->gp == evaluate(p^x, lp), gens(G))
end

function set_orbit(G::PermGroup, H::PermGroup)
  #from Elsenhans
  #https://math.uni-paderborn.de/fileadmin/mathematik/AG-Computeralgebra/inv_transfer_5_homepage.pdf
  #    http://dblp.uni-trier.de/db/journals/jsc/jsc79.html#Elsenhans17
  # https://doi.org/10.1016/j.jsc.2016.02.005

  l = low_index_subgroups(H, 2*degree(G)^2)
  S, g = slpoly_ring(ZZ, degree(G), cached = false)

  sort!(l, lt = (a,b) -> isless(order(b), order(a)))
  for U = l
    O = orbits(U)
    for o in O
      #TODO: should use orbits of Set(o)...
      f = sum(g[o])
      oH = probable_orbit(H, f)
      oG = probable_orbit(G, f, limit = length(oH)+5)
      if length(oH) < length(oG)
        for i = 1:length(o)
          I = sum(x^i for x = oH)
          if isprobably_invariant(I, H) &&
             !isprobably_invariant(I, G)
            @vprint :GaloisInvariant 2 "SetOrbit won\n"
            return true, I
          end
        end
        f = prod(g[o])
        oH = probable_orbit(H, f)
        I = sum(oH)
        if isprobably_invariant(I, H) &&
           !isprobably_invariant(I, G)
          @vprint :GaloisInvariant 2 "SetOrbit won - final attempt\n"
          return true, I
        end
       end
     end
   end
   return false, g[1]
end

@doc Markdown.doc"""
    invariant(G::PermGroup, H::PermGroup)

For a permutation group `G` and a maximal subgroup `H`, find
a (multivariate (`SLPoly`)) `f` with `G`-stabilizer `H`, i.e. 
`f` is invariant under all permutations in `H`, but not invariant under
any other element of `G`.
"""
function invariant(G::PermGroup, H::PermGroup)
  @vprint :GaloisInvariant 1 "Searching G-relative H-invariant\n"
  @vprint :GaloisInvariant 2 "that is a $G-relative $H-invariant\n"

  S, g = slpoly_ring(ZZ, degree(G), cached = false)

  if istransitive(G) && !istransitive(H)
    @vprint :GaloisInvariant 2 "top group transitive, bottom not\n"
    return sum(probable_orbit(H, g[1]))
  end

  if !istransitive(G) 
    @vprint :GaloisInvariant 2 "both groups are intransitive\n"
    OG = [sort(x) for x = orbits(G)]
    OH = [sort(x) for x = orbits(H)]
    d = setdiff(OH, OG)
    if length(d) > 0
      @vprint :GaloisInvariant 2 "groups have different orbits\n"
      return sum(probable_orbit(H, g[d[1][1]]))
    end
    #OH == OG
    for o = OH
      h = action_homomorphism(G, o)
      hG = image(h)[1]
      hH = image(h, H)[1]
      if order(hG) > order(hH)
        @vprint :GaloisInvariant 2 "differ on action on $o, recursing\n"
        @hassert :GaloisInvariant 0 ismaximal(hG, hH)
        I = invariant(hG, hH)
        return evaluate(I, g[o])
      end
    end
    @vprint :GaloisInvariant 2 "going transitive...\n"
    #creating transitive version.
    gs = Oscar.gset(G, [[x[1] for x = OG]])
    os = Oscar.orbits(gs)
    @assert length(os) == 1
    os = os[1]
    h = Oscar.action_homomorphism(os)
    GG = h(G)[1]
    HH = h(H)[1]
    I = invariant(GG, HH)
    ex = 1
    while true
      J = evaluate(I, [sum(g[o])^ex for o = elements(os)])
      if !isprobably_invariant(J, G)
        I = J
        break
      end
      ex += 1
    end
    @hassert :GaloisInvariant 2 isprobably_invariant(I, H)
    @hassert :GaloisInvariant 2 !isprobably_invariant(I, G)

    return I
  end

  if isprimitive(G) && isprimitive(H)
    if isodd(G) && iseven(H)
      @vprint :GaloisInvariant 3 "using sqrt_disc\n"
      return sqrt_disc(g)
    end
    @vtime :GaloisInvariant 2 fl, I = set_orbit(G, H)
    fl && return I
  end

  bG = all_blocks(G)
  bH = all_blocks(H)

  d = setdiff(bH, bG)
  @vprint :GaloisInvariant 2 "Have block systems, they differ at $d\n"
  if length(d) > 0
    @vprint :GaloisInvariant 3  "using F-invar\n"
    return prod(sum(g[b]) for b = block_system(H, d[1]))
  else
    for BB = bG
      #TODO: Max(?): some of this possibly cannot happen anymore
      #      in starting group the sieving has been improved
      #      classical D_sm or similar may no-longer be possible
      B = block_system(H, BB)
      m = length(B)
      l = length(B[1])
      y = [sum(g[b]) for b = B]
      d = [sqrt_disc(g[b]) for b = B]
      D = sqrt_disc(y)

      I = D
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G) #TODO: this can be decided theoretically
        @vprint :GaloisInvariant 3 "using D-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, 1)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G) #TODO: this can be decided theoretically
        @vprint :GaloisInvariant 3 "using s1-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, m)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G) #TODO: this can decided theroetically
        @vprint :GaloisInvariant 3 "using sm-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, 2)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G) #TODO
        @vprint :GaloisInvariant 3 "using s2-invar for $BB\n"
        return I
      end
      I = D*elementary_symmetric(d, m)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G)
        @vprint :GaloisInvariant 3 "using D sm-invar for $BB\n"
        return I
      end
    end

    for BB = bG
      h = action_on_blocks(G, BB)
      B = block_system(G, BB)
      hG = image(h)[1]
      hH = image(h, H)[1]
      if length(hH) < length(hG)
        y = [sum(g[b]) for b = B]
        J = invariant(hG, hH)
        E = evaluate(J, y)
        @vprint :GaloisInvariant 3 "using E-invar for $BB (4.1.2)\n"
        return E
      end

      sG = stabilizer(G, BB, on_sets)[1]
      sH = stabilizer(H, BB, on_sets)[1]
      hG = action_homomorphism(sG, BB)
      ssG = image(hG, sG)[1]
      ssH = image(hG, sH)[1]
      if length(ssH) < length(ssG)
        J = invariant(ssG, ssH)
        C = left_transversal(H, sH)
        gg = g[BB]
        J = evaluate(J, gg)
        F = sum(J^t for t = C)
        @hassert :GaloisInvariant 1 isprobably_invariant(F, H)
        @hassert :GaloisInvariant 1 !isprobably_invariant(F, G)
        @vprint :GaloisInvariant 3 "using F-invar for $BB (4.1.4)\n"
        return F
      end
    end
  end

  @vprint :GaloisInvariant 3 "no special invar found, resorting to generic\n"

  m = prod(gen(S, i)^i for i=1:degree(G)-1)
  return sum(m^s for s = H)
end

@doc Markdown.doc"""
    resolvent(C::GaloisCtx, G::PermGroup, U::PermGroup)

Find a `G`-relative `H`-invariant `I` and form the correspondig resolvent polynomial
``\prod (x-I^t)``
where the product runs over all coset-representatives of `G/U`.
"""
function resolvent(C::GaloisCtx, G::PermGroup, U::PermGroup, extra::Int = 5)
  I = invariant(G, U)
  t = right_transversal(G, U)
  n = length(t)
  rt = roots(C)
  #make square-free (in residue field)
  k, mk = ResidueField(parent(rt[1]))
  k_rt = map(mk, rt)
  ts = find_transformation(k_rt, I, t)

  B = 2*n*evaluate(I, map(ts, [C.B for i = 1:ngens(parent(I))]))^n
  rt = roots(C, bound_to_precision(C, B, extra))
  rt = map(ts, rt)
  rt = [evaluate(I^s, rt) for s = t]
  pr = copy(rt)
  ps = [isinteger(C, B, sum(pr))[2]]
  while length(ps) < n
    pr .*=  rt
    push!(ps, isinteger(C, B, sum(pr))[2])
  end
  if isa(ps[1], fmpz)
    return Hecke.power_sums_to_polynomial(map(fmpq, ps))
  else
    return Hecke.power_sums_to_polynomial(ps)
  end
end

"""
    minpoly(C::GaloisCtx, I, extra::Int = 5)

Computes the minimal polynomial of `I` evaluated at the roots
stored in `C`.
"""
function Hecke.minpoly(C::GaloisCtx, I, extra::Int = 5)
  O = collect(probable_orbit(C.G, I))
  n = length(O)
  rt = roots(C)
  #make square-free (in residue field)
  k, mk = ResidueField(parent(rt[1]))
  k_rt = map(mk, rt)
  ts = gen(Hecke.Globals.Zx)
  while true
    r = map(ts, k_rt)
    s = Set(map(x->evaluate(x, r), O))
    if length(s) < length(O)
      while true
        ts = rand(Zx, 2:rand(2:max(2, length(r))), -4:4) #TODO: try smaller degrees stronger
        if degree(ts) > 0
          break
        end
      end
    else
      break
    end
  end

  B = 2*n*evaluate(I, map(ts, [C.B for i = 1:ngens(parent(I))]))^n
  rt = roots(C, bound_to_precision(C, B, extra))
  rt = map(ts, rt)
  rt = [evaluate(i, rt) for i = O]
  pr = copy(rt)
  ps = [isinteger(C, B, sum(pr))[2]]
  while length(ps) < n
    pr .*=  rt
    push!(ps, isinteger(C, B, sum(pr))[2])
  end
  if isa(ps[1], fmpz)
    return Hecke.power_sums_to_polynomial(map(fmpq, ps))
  else
    return Hecke.power_sums_to_polynomial(ps)
  end
end

Hecke.minpoly(R::FmpzPolyRing, C::GaloisCtx, I, extra::Int = 5) = Hecke.minpoly(C, I, extra = extra)(gen(R))
Hecke.minpoly(R::FmpqPolyRing, C::GaloisCtx, I, extra::Int = 5) = Hecke.minpoly(C, I, extra = extra)(gen(R))

struct GroupFilter
  f::Vector{Tuple{Function, String}}
  GroupFilter() = new([(x->true, "")])
end

#TODO: Max: rank the filter function by cost, cheapest first...
#      possibly needs a third field in the structure
function (F::GroupFilter)(G::PermGroup)
  do_print = Hecke.get_verbose_level(:GaloisGroup) >= 1
  do_all = Hecke.get_verbose_level(:GaloisGroup) >= 2

  res = true
  for (x, s) = F.f
    fl = x(G)
    fl && continue
    if do_all
      if res 
        println("group rejected by:\n\t- $s")
      elseif length(s) > 0
        println("\t- $s")
      end
      res = false
    elseif do_print
      println("group rejected by: $s")
      return false
    else
      return false
    end
  end
  return res
end

function Base.push!(G::GroupFilter, F::Function, s::String="")
  push!(G.f, (F, s))
end

@doc Markdown.doc"""
# `DescentEnv`

Given a permutation group `G`, known to contain the Galois group of some polynomial,
the key step in Stauduhar's method is to systematically test all maximal subgroups
of `G`. This structure is used to organize this process. On creation, all the
maximal subgroups are computed and possibly filtered (if the polynomial
is irreducible, only transitive groups are needed, if the stem field if primitive,
then so is the group...). The filtering is done via a `GroupFilter`.

Furthermore we may compute (and store)
 - a `G`-relative `H`-invariant for any subgroup (`I`)
 - a list of Tschirnhausen transformations used (`T`)
 - a list (`l`) marking which subgroups we have dealt with (`1` means group has been
 processed completely)
"""
mutable struct DescentEnv
  G::PermGroup
  s::Vector{PermGroup}
  I::Dict{Int, SLPoly}
  T::Dict{Int, Vector{fmpz_poly}}
  l::Vector{Int}
  #the coset reps need to be cached as well
  #a "work limit" on the "invariant" function
  #a more select choice of group....

  function DescentEnv(G::PermGroup, f::GroupFilter = GroupFilter())
    s = maximal_subgroup_reps(G)
    r = new()
    r.G = G
    @vprint :GaloisGroup 1 "starting with $(length(s)) maximal subgroup classes\n"
    r.s = filter(f, s)
    @vprint :GaloisGroup 1 "filtering down to $(length(r.s)) maximal subgroup classes\n"
    r.I = Dict{Int, SLPoly}()
    r.T = Dict{Int, Vector{fmpz_poly}}()
    r.l = zeros(Int, length(r.s))
    return r
  end
end

#TODO: instead of iterating, use the cheapest invariant
#      use invariant with a worklevel to limit costs in finding them
#      use the V_4 (index 2) stuff when necessary
function Base.iterate(D::DescentEnv, i::Int=1)
  #r.l[i] == 1: done,  0:: nothing
  all(isone, D.l) && return nothing

  while i<=length(D.s) && D.l[i] == 1
    i += 1
  end

  local I, TS
  if haskey(D.I, i)
    I = D.I[i]
  else
    I = D.I[i] = invariant(D.G, D.s[i])
  end

  if haskey(D.T, i)
    T = D.T[i]
    ts = rand(Hecke.Globals.Zx, 2:degree(D.G), -4:4)
    while degree(ts) < 2 || ts in T
      ts = rand(Hecke.Globals.Zx, 2:degree(D.G), -4:4)
    end
    push!(D.T[i], ts)
    TS = ts
  else
    D.T[i] = [gen(Hecke.Globals.Zx)]
    TS = D.T[i][1]
  end

  return (D.s[i], I, TS), i
end

function Base.push!(D::DescentEnv, i::Int)
  D.l[i] = 0
end

function Base.pop!(D::DescentEnv, i::Int)
  D.l[i] = 1
end

function sum_orbits(K, Qt_to_G, r)
  @vprint :GaloisGroup 1 "group is primitive (no subfields), trying operations on pairs\n"

  #starting group as the stabilizer of the factorisation of the 2-sum poly,
  #the polynomial with roots r_i + r_j, i<j
  #we need this square-free, so we compute this over the finite field
  #actually, we only check that the roots over the finite field are distinct.
  #if not: transform (using ts) and try again.
  
  k = parent(r[1])
  m = Dict{typeof(r[1]), Vector{Int}}()  #the vector has length 2 always
  local ts = gen(Hecke.Globals.Zx)
  c = r
  while true
    for i=1:length(r)-1
      for j=i+1:length(r)
        m[c[i]+c[j]] = [i,j]
      end
    end
    if length(keys(m)) < binomial(degree(K), 2)
      @vprint :GaloisGroup 2 " 2-sum: found duplicate, transforming...\n"
      while true
        ts = rand(Hecke.Globals.Zx, 2:degree(K), -4:4)
        if degree(ts) > 1
          break
        end
      end
      c = map(ts, r)
      empty!(m)
    else
      break
    end
  end

  @vprint :GaloisGroup 1 "have everything, now getting the 2-sum poly\n"
  if gen(Hecke.Globals.Zx) == ts
    @vtime :GaloisGroup 2 g = msum_poly(defining_polynomial(K), 2) #if f has roots a_i, then g has roots a_i+a_j, if G is 2-transitive, this is pointless.
  else
    @vtime :GaloisGroup 2 g = msum_poly(minpoly(ts(gen(K))), 2)
  end
  @vprint :GaloisGroup 1 "... factoring...\n"
  @vtime :GaloisGroup 2 fg  = factor(g)
  @assert all(isone, values(fg.fac))

  O = []
  for f = keys(fg.fac)
    r = roots(map_coefficients(Qt_to_G, f))
    push!(O, [m[x] for x = r])
  end
  @vprint :GaloisGroup 2 "partitions: $O\n"
  return O
end

"""
Given 2 normal subgroups of index 2, there is a 3 index 2 subgroup
automatically, given that the quotient ov G by the intersection
is a V_4 group.

This computes (in a dumb way) the 3rd group
"""
function index2sum(G::PermGroup, H::PermGroup, V::PermGroup)
  @assert index(G, H) == index(G, V) == 2
  U = intersect(H, V)[1]
  @assert index(G, U) == 4
  t = right_transversal(G, U)
  s = [x for x = t if !(x in H || x in V)]
  @assert length(s) == 1
  return sub(G, vcat(gens(U), s))[1]
end

"""
Maos elements of the base ring into the complation (the splitting field).
FOr computations over QQ this does not matter as QQ just coerces into
any q/p-adic field.

FOr computations over number fields this is not true, here an explciti map
is needed.

TODO: check and add a precision parameter
TODO: implement also for Q(t) (even though not used)
"""
function map_coeff(C::GaloisCtx{Hecke.qAdicRootCtx}, x)
  return (C.C.Q[1])(x)
end

@doc Markdown.doc"""
    starting_group(GC::GaloisCtx, K::SimpleNumberField) 

Finds a _starting group_, that is a group `G` as a subgroup of the
symmetric group acting on the roots in the explicit ordering in the
galois context.

If the field is imprimitive, the group is derived from the subfields, otherwise
the factorisation of the 2-sum polynomial is used.

Returns a triple: 
 - the group
 - a filter for all groups that can occur
 - a permutation representing the operation of the Frobenius automorphism
"""
function starting_group(GC::GaloisCtx, K::T; useSubfields::Bool = true) where T <: SimpleNumField
  c = roots(GC, 5, raw = true)

  @vprint :GaloisGroup 1 "computing starting group (upper bound for Galois group)\n"

  if useSubfields
    @vprint :GaloisGroup 1 "computing principal subfields ...\n"
    @vtime :GaloisGroup 2 S = Hecke.principal_subfields(K)
  else
    S = []
  end

  #compute the block system for all subfields...
  bs = Vector{Vector{Int}}[]
  fld = Dict{Vector{Vector{Int}}, Any}()
  for (s, ms) = S
    if degree(s) == degree(K) || degree(s) == 1
      continue
    end
    pr = 5
    local v
    #TODO: maybe also use cycle_types for lower bounds?
    while true
      c = roots(GC, pr, raw = true)
      g = ms(gen(s))
      gg = map_coefficients(x->map_coeff(GC, x), parent(K.pol)(g))
      d = map(gg, c)
      b = Dict{typeof(c[1]), Vector{Int}}()
      for i=1:length(c)
        if Base.haskey(b, d[i])
          push!(b[d[i]], i)
        else
          b[d[i]] = Int[i]
        end
      end
      v = collect(values(b))
      if any(x->length(x) != div(degree(K), degree(s)), v)
        pr *= 2
        if pr > 100
          error("too bad")
        end
      else
        break
      end
    end
    #each subfield defines a block system, each block system results in
    #a wreath product representing the largest subgroup of Sn with this
    #block system: Sl wr Sk
    #However, this wreath prodct has automatically 3 maximal subgroups
    #of index 2 (Eichenlaub Prop 3, p 28 or Geissler Satz 6.8).
    # - Sl wr Ak: possible iff disc(subfield) is a square, so
    #   we can test this and then reduce the starting group 
    #      or
    #   add this as an exclusion: we don't need to test subgroups of this
    #   any more
    # - one of them should be (Sl wr Sk) meet A(lk))
    #   disc(K) is a square
    # - see below in starting_group
    sort!(v, lt = (x,y) -> minimum(x) < minimum(y))
    push!(bs, v)
    fld[v] = (s, ms)
  end

  F = GroupFilter()
  push!(F, istransitive, "not transitive") #poly is defining number field, hence irreducible


  #selecting maximal block systems only...
  if length(bs) > 1
    @vprint :GaloisGroup 2 "group will have (all) block systems: $([x[1] for x = bs])\n"
    L = POSet(bs, (x,y) -> issubset(x[1], y[1]) || issubset(y[1], x[1]),
                  (x,y) -> issubset(x[1], y[1]) - issubset(y[1], x[1]))
    bs = map(Oscar.data, minimal_elements(L))
    @vprint :GaloisGroup 1 "group will have (maximal) block systems: $([x[1] for x = bs])\n"
  end
  GC.start = (1, bs)

  G = symmetric_group(degree(K))
  if degree(K) == 1
    return G, F, one(G)
  end
  S = G

  #we only need groups that allow the subfield to keep the
  #"correct" galois group. Up to isomorphism, the subfield galois group
  #is the operation on the block-system (permuting of blocks)
  for v = bs
    s, ms = fld[v] #it appears that principal_subfields are maximal..
                   #(all others are via field intersection, hence smaller)
    if degree(s) > 2
      #no point in checking C_2
      local gc, bs
      # TODO: infer the subfields of s from the ones in K
      #       the poset should make this trivial.
      @vprint :GaloisGroup 1 "computing group of subfield of degree $(degree(s))\n"
      local g, gc
      if valuation(discriminant(s), GC.prime) > 0
        throw(Hecke.BadPrime(GC.prime))
      end
      g, gc = try
          Hecke.pushindent()
          galois_group(s, prime = GC.prime)
        catch e  #can "legally" be BadPrime, will be dealt with one level up
          Hecke.popindent()
          rethrow()
        end
      Hecke.popindent()

      i = order(g)

      pr = 5
      r = roots(gc, pr, raw = true)
      R = roots(GC,pr, raw = true)

      gg = map_coefficients(x->map_coeff(GC, x), parent(K.pol)(ms(gen(s))))
      d = map(gg, R)
      f, mf = ResidueField(parent(r[1]))
      _F, mF = ResidueField(parent(R[1]))
      mfF = find_morphism(f, _F)
      #we should have
      # - d == r (in the approriate setting)
      # - order(Set(map(mF, d))) == order(Set(map(mf, r))) == degree(s)
      # - Set(map(mfF, map(mf, r))) == Set(map(mF, d))
      # - d is partitioned
      # and g (the small Galois group) operates on the blocks in d
      # Now we need to create the wreath_product(sym(deg(K)), g)
      # that produces THIS block system
      # The ambiguity of the different possible embeddings (mfF)
      # should be harmless as they are linked by the Frobenius
      # and Frobenius in g, so g^Frobenius == g....
      r = map(mf, r)
      r = map(mfF, r)
      @assert parent(R[1]) == parent(d[1])
      d = map(mF, d)
      R = map(mF, R)
      if length(Set(r)) == length(r) &&
         length(Set(d)) == divexact(length(d), length(r))
         length(Set(R)) == length(R) #no problems with precision, 
                                     #we can proceed.
        @assert Set(r) == Set(d)
        #assuming wreath_product(A, B) has block system [[1..l], [l+1..2l]...]
        # the re-ordering is easy: W^s for s = vcat(bs)
        bs = map(x->findall(isequal(x), d), r)
        @assert all(x->length(x) == length(bs[1]), bs)
        W = isomorphic_perm_group(wreath_product(symmetric_group(length(bs[2])), g))[1]
        #should have the block system as above..
        W = W^symmetric_group(degree(W))(vcat(bs...))
        G = intersect(G, W)[1]
      else
        @vprint :GaloisGroup 2 "precision problem in wreath product, using size only\n"
      end
      let _v = v, _i = i 
        push!(F, x->order(image(action_on_blocks(x, _v[1]))[1]) % _i == 0 , "subfield for $(_v[1]) would have wrong group (too small)")
      end
    end
  end

  d = map(frobenius, c)
  si = [findfirst(y->y==x, c) for x = d]

  @vprint :GaloisGroup 1 "found Frobenius element: $si\n"

  in_br = false
  if issquare(discriminant(K))
    G = intersect(G, alternating_group(degree(K)))[1]
    in_br = true
  else
    push!(F, isodd, "has to be odd")
  end

  for b = bs
    #each block system b corresponds to a subfield and induces an upper
    #limit of i
    #  G =  Sym(k) wr Sym(l) 
    #on the Galois group
    #This wreath product has always (at least) 3 index 2 subgroups:
    #  G meet Alt(k*l) 
    #  Sym(k) wr Alt(l)
    #  index2sum of the above
    #We are a subgroup of the 1st iff disc(K) is a square (in_br = true)
    #                         2nd iff disc(k) is a square (in_ar = true)
    #                         3rd iff disc(k)*disc(K) is a square
    #This is used in both directions
    # - in an intersection to limit the starting group (upper bound)
    #   (true case)
    # - in the filter to avoid unneccessary tests
    #   (false case)
    # disc(k) is up to an unknown Tschirnhaus transformation exactly the
    # D(y_1, .. y_k) invar in Geisser/Eichenlaub (case 2)
    # and the product then is one for the 3rd case.
    s = S(vcat(b...))

    wr = PermGroup(wreath_product(symmetric_group(length(b[1])), symmetric_group(length(b))))
    br = intersect(wr, alternating_group(degree(K)))[1]
    wr = wr^s
    br = br^s

    G = intersect(G, wr)[1]

    can_use_wr = length(b) > 2
    #for length(b) == 2, the bottom field is quadratic and thus always S2
    #the wreath product would not be transitive here
    #TODO: the new stuff above shoudl guarantee that the induced action
    #  on a block system can never be too large:
    #    (sym(l) wr gal(k) is constructed above
    #    does this part still make sense????
    if can_use_wr
      ar = PermGroup(wreath_product(symmetric_group(length(b[1])), alternating_group(length(b))))
      ar = ar^s
      cr = index2sum(wr, ar, br)
    end


    in_ar = false
    if issquare(discriminant(fld[b][1]))
      in_ar = true
    end
    #in_ar == true iff galois <= ar
    #in_br == true iff galois <= br
    in_cr = issquare(discriminant(K)*discriminant(fld[b][1]))

    if can_use_wr
      if in_ar 
        G = intersect(G, ar)[1]
      else
        let ar = ar 
          push!(F, x->!issubgroup(ar, x)[1], "subfield is even/odd")
        end
      end
    end
    if in_br
      #G = intersect(G, br)[1] #already done above
    else
      #push!(F, x->!issubgroup(br, x)[1]) #already done above
    end
    if can_use_wr
      if in_cr
        G = intersect(G, cr)[1]
      else
        let cr = cr
          push!(F, x->!issubgroup(cr, x)[1], "third subgroup of wreath product")
        end
      end
    end
  end

  if length(bs) == 0 #primitive case: no subfields, no blocks, primitive group!
    push!(F, isprimitive, "primitivity")
    pc = parent(c[1])
    k, mk = ResidueField(pc)
    O = sum_orbits(K, x->mk(pc(map_coeff(GC, x))), map(mk, c))
    GC.start = (2, O)
    
    #the factors define a partitioning of pairs, the stabiliser of this
    #partition is the largest possible group...
    #code from Max...

    #TODO: wrap this properly
    if hasproperty(GAP.Globals, :ConStabilize)
      tmp = [GAP.Globals.ConStabilize(GAP.Obj(sort(o), recursive=true), GAP.Globals.OnSetsSets) for o in O]
      H = GAP.Globals.Solve(GAP.Obj(vcat(GAP.Globals.ConInGroup(G.X), tmp)))
      G = Oscar._as_subgroup(G, H)[1]
    else
      #TODO: fallback if ferret wasn't loaded for some reason; this should be removed
      # once we have GAP_pkg_ferret
      G = intersect([stabilizer(G, sort(o), on_sets_sets)[1] for o=O]...)[1]
    end
  end

  si = G(si)

  return G, F, si
end

@doc Markdown.doc"""
    find_prime(f::fmpq_poly, extra::Int = 5; prime::Int = 0, pStart::Int = 2*degree(f)) -> p, c

Tries to find a useful prime for the computation of Galois group. Useful means
  - degree of the splitting field not too small
  - degree of the splitting field not too large
Starts searching at `pStart`, returns the prime as well as the cycle types identified.

If `prime` is given, no search is performed.
"""
function find_prime(f::fmpq_poly, extra::Int = 5; prime::Int = 0, pStart::Int = 2*degree(f), filter_prime = x->true, filter_pattern = x->true)
  if prime != 0
    p = prime
    lf = factor(f, GF(p))
    return p, Set([CycleType(map(degree, collect(keys(lf.fac))))])
  end
  if pStart < 0
    error("should no longer happen")
    p = -pStart
    lf = factor(f, GF(p))
    return p, Set([CycleType(map(degree, collect(keys(lf.fac))))])
  end

  d_max = degree(f)^2
  d_min = min(2, d_max)
  p_best = 1
  cnt = 5
  ct = Set{Vector{Int}}()
  ps = Vector{Tuple{Int, Int}}()

  #find a q-adic splitting field of "good degree":
  # - too small, then the Frobenius automorphisms is not containing lots of
  #   information
  # - too large, all is slow.
  # careful: group could be (C_2)^n hence d_min might be small...
  no_p = 0
  for p = Hecke.PrimesSet(pStart, -1)
    filter_prime(p) || continue
    k = GF(p)
    if k(leading_coefficient(f)) == 0
      continue
    end
    lf = factor(f, GF(p))
    if any(x->x>1, values(lf.fac))
      continue
    end
    filter_pattern(lf) || continue
    no_p +=1 
    push!(ct, sort(map(degree, collect(keys(lf.fac))))) # ct = cycle types as vector of cycle lengths
    d = lcm([degree(x) for x = keys(lf.fac)])
    if d <= d_max 
      if d >= d_min
        push!(ps, (p, d))
      else
        no_p += 1
        if no_p > degree(f)
          d_min = 1
        end
      end
    end
    if length(ps) > 2*degree(f)
      break
    end
  end

  @vprint :GaloisGroup 2 "possible primes and degrees: $ps\n"

  #want: degree > degree K/2 - which might not be possible.
  l = [x for x = ps if 2*x[2] <= degree(f)]
  if length(l) > 0
    pd = l[argmax([x[2] for x = l])]
  else
    pd = ps[argmin([x[2] for x = ps])]
  end

  p = p_best = pd[1]
  d_max = pd[2]

  @vprint :GaloisGroup 1 "using prime $p_best with degree $d_max\n"
  @vprint :GaloisGroup 2 "and cycle types $ct\n"
  return p_best, Set(CycleType(x) for x = ct)
end

#given cycle types, try to find a divisor (large) of the transitive
#group. Used for sieving.
# Elsenhans-Klueners:
# Computing subfields of number fields and
# applications to Galois group computations
function order_from_shape(ct::Set{CycleType}, n)
  o1 = fmpz(1)
  #Lemma 15:
  o  = n * reduce(lcm, fmpz[divexact(
           reduce(lcm, fmpz[l for (l,v) = x]),
           reduce(gcd, fmpz[l for (l,v) = x])) for x = ct if length(x) > 0], 
           init = fmpz(1))
  #Lemma 16:         
  for p = PrimesSet(1, n)
    #find cycletypes that would correpond to elt of order a multiple of p
    #the exact order is lcm(x) for x in ct
    ct_p = [x for x = ct if any(t->t[1] % p == 0, x)]
    #so ct_t contains all cycle types belongign to elements where the
    #order is divisible by p
    length(ct_p) == 0 && continue
    #now x in ct_p is is the cycle type of an element where the order is a multiple
    #of p. I need the type of x^f (which is of maximal order p-power order)
    #e.g. a cycle of length 6 will produce 3 of length 2
    ct_pp = Set{CycleType}()
    for x = ct_p
      y = x^Int(ppio(order(x), fmpz(p))[2])
      push!(ct_pp, y^y.s[1][1])
    end
    # Now everything has a fixed point, so we can bound the
    # size of Stab_G(1) from below (we're transitive, so the fixed point
    # is 1)
    # the plan is to take to 2 different largest ones

    mu = [(y.s[end][1], y.s[end][1]*y.s[end][2]) for y = ct_pp]
    sort!(mu, lt = (a,b) -> a[1] < b[1])
    if length(mu) == 1
      p_part = mu[1][1]
      o1 *= p_part
    else
      p_part = maximum([x[1]*y[1] for x = mu for y = mu if x[2] != y[2]], init = maximum(x[1] for x = mu))
      o1 *= p_part
    end
  end
  return lcm(n*o1, o)
end


#https://mathoverflow.net/questions/286316/recognition-of-symmetric-groups-in-gap
#show primitivity by finding a permutation which has a p cycle for n/2<p<n−1
function primitive_by_shape(ct::Set{CycleType}, n::Int)
  for p = PrimesSet(div(n, 2)+1, n-2)
    if any(y->any(x->x[1] % p == 0, y), ct)
      return true
    end
  end
  return false
end

function an_sn_by_shape(ct::Set{CycleType}, n::Int)
  n <= 3 && return true

  ub = n > 8 ? n-3 : n-2
  lp = reduce(lcm, map(fmpz, collect(PrimesSet(div(n, 2)+1, ub))), init = fmpz(1))

  if any(x->any(y ->gcd(y[1], lp) > 1, x.s), ct)
    n != 6 && return true
    #from Elsenhans (ask Max???) need the 5-cycle and 
    #    [1,1,1,3] or [1,2,3] ie. not [3,3]
    # Elsenhans also tries to seperate Sn/An but there I am not certain
    # I think if there is a type proving odd, then we have Sn, hower
    # not finding one does not prove anyhting
    if any(x->(3=>2) in x, ct)
      return true
    end
  end
  return false
end

function Oscar.cycle_structure(x::GroupConjClass{PermGroup, PermGroupElem})
  return cycle_structure(representative(x))
end

function cycle_structures(G::PermGroup)
  r = conjugacy_classes(G)
  return Set(cycle_structure(x) for x = r)
end

#
#TODO
# - subfield lattice directly
# - subfield to block system outsource (also for lattice)
# - more invars
# - re-arrange to make cheap things first, expensive later
# - proof, live and later
# - tschirnhausen transformation: slowly increasing degree and coeffs.
# - for larger degrees: improve starting group by doing Galois groups of subfields
# - more base rings
# - applications: subfields of splitting field (done), towers, solvability by radicals
@doc Markdown.doc"""
    galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true, pStart::Int = 2*degree(K)) -> PermGroup, GaloisCtx

Computes the Galois group of the splitting field of the defining polynomial of `K`.
Currently the polynomial needs to be monic.

The group is returned as an explicit permutation group permuting the roots as contained
in the contex object (the 2nd return value). The roots live in a suitable unramifed
extension of the p-adics.

EXAMPLE
```jldoctest
julia> K, a = cyclotomic_field(5);
julia> G, C = galois_group(K)
julia> G, C = galois_group(K)
(Group([ (1,4,2,3), (1,2)(3,4) ]), Galois Context for x^4 + x^3 + x^2 + x + 1 and prime 19)

julia> describe(G)
"C4"

julia> roots(C, 2)
4-element Vector{qadic}:
 (15*19^0 + 16*19^1 + O(19^2))*a + 9*19^0 + 7*19^1 + O(19^2)
 (4*19^0 + 2*19^1 + O(19^2))*a + 5*19^0 + 9*19^1 + O(19^2)
 (18*19^0 + 18*19^1 + O(19^2))*a + 12*19^0 + O(19^2)
 (19^0 + O(19^2))*a + 11*19^0 + 19^1 + O(19^2)
```
"""
function galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true, pStart::Int = 2*degree(K), prime::Int = 0)

  if prime != 0
    p = prime
    ct = Set{CycleType}()
  else
    p, ct = find_prime(K.pol, pStart = pStart, prime = prime)
  end

  # TODO: otherwise, try to detect here if we are primitive or even 2-transitive
  #       primitive: done.
  #       2-transitive: ????
  # TODO: sieve (filter) by 
  #       - cycle_structure (done)
  #       - primitive/ transitive/ 2-transitive/ ... (partly done)
  #       - blocks (not sure: we're using the subfield group and wreath prods)

  # TODO: this too needs to be made either generic or replicated in the 
  #       relative case
  #       or the RootCtx needs to learn to deal with bad primes
 
  while true
    GC = GaloisCtx(Hecke.Globals.Zx(K.pol), p)
    r = roots(GC, 5)
    if length(r) < degree(K) 
      @vprint :GaloisGroup 1 "in recursion: bad prime detected\n"
      prime == 0 || throw(Hecke.BadPrime(p))
      p, _ = find_prime(K.pol, pStart = p+1)
      continue
    end

    if an_sn_by_shape(ct, degree(K))
      @vprint :GaloisGroup 1 "An/Sn by cycle type\n"
      if issquare(discriminant(K))
        G = alternating_group(degree(K))
      else
        G = symmetric_group(degree(K))
      end
      GC.G = G
      return G, GC
    end

    G, F, si = try
      starting_group(GC, K, useSubfields = useSubfields && !primitive_by_shape(ct, degree(K)))
    catch e
      if isa(e, Hecke.BadPrime)
        @vprint :GaloisGroup 1 "in recursion: bad prime detected\n"
        prime == 0 || rethrow()
        p, _ = find_prime(K.pol, pStart = p+1)
        continue
      else
        rethrow(e)
      end
    end
    # the filter will filter at least on 
    #   size (have a lower bound from cycles)
    #   even/odd
    #   some subgroups of the wreath products
    if degree(K) == 1
      GC.G = G
      return G, GC
    end

    os = order_from_shape(ct, degree(K))
    @vprint :GaloisGroup 2 "using lower bound $os for order to sieve\n"
    push!(F, x->order(x) % os == 0, "group too small")
    if degree(K) < 20 #arbitrary....
      push!(F, x->issubset(ct, cycle_structures(x)), "by cycle structures")
    end

    # TODO: here we know if we are primitive; can we detect 2-transitive (inside starting_group)?

    return descent(GC, G, F, si, extra = extra)
  end
end

@doc Markdown.doc"""
    descent(GC::GaloisCtx, G::PermGroup, F::GroupFilter, si::PermGroupElem; grp_id = transitive_identification, extra::Int = 5)

Performs a generic Stauduhar descent: starting with the group `G` that needs to be a 
supergroup of the Galois group, operating on the roots in `GC`, the context object.

The groups are filtered by `F` and the result needs to contain the permutation `si`.
For verbose output, the groups are printed through `grp_id`.
"""
function descent(GC::GaloisCtx, G::PermGroup, F::GroupFilter, si::PermGroupElem; grp_id = transitive_identification, extra::Int = 5)
  @vprint :GaloisGroup 2 "Have starting group with id $(grp_id(G))\n"

  n = degree(GC.f)

  nG = order(G)
  while true
    D = DescentEnv(G, F)
    @vprint :GaloisGroup 2 "found $(length(D.s)) maximal subgroups\n"

    while true
      d = iterate(D)
      if d === nothing
        break
      end
      s, I, ts = d[1]

      @vprint :GaloisGroup 2 "testing descent $(grp_id(G)) -> $(grp_id(s)) of index $(index(G, s))\n"

      @hassert :GaloisInvariant 1 all(x->isprobably_invariant(I, x), s)
      @hassert :GaloisInvariant 1 any(x->!isprobably_invariant(I, x), G)

      B = upper_bound(GC, I, ts)
      @vprint :GaloisGroup 2 "invariant uses $(cost(I, ts)) multiplications\n"
      @vprint :GaloisGroup 2 "of maximal degree $(total_degree(I)*degree(ts))\n"
      @vprint :GaloisGroup 2 "root upper_bound: $(value(B))\n"

      c = roots(GC, bound_to_precision(GC, B, extra))
      c = map(ts, c)

      local fd
      cnt = 0

      cs = Set{typeof(c[1])}()
      fd = []

      local lt
      if index(G, s) < 100
        @vtime :GaloisGroup 2 lt = right_transversal(G, s)
      elseif isnormal(G, s)
        lt = [one(G)] # I don't know how to get the identity
      else
        @vtime :GaloisGroup 2 lt = short_right_transversal(G, s, si)
        @hassert :GaloisInvariant 2 all(x->si in s^x, lt)
      end
      #x^7-2 triggers wrong descent, occaisonly, if lt is only 10 elems long
      while length(lt) < 20 && length(lt) < index(G, s)
        r = rand(G)
        if all(x->right_coset(s, r) != right_coset(s, x), lt)
          push!(lt, r)
        end
      end

      compile!(I)
      for t = lt
        e = evaluate(I, t, c)
        if e in cs
          @vprint :GaloisGroup 2 " evaluation found duplicate, transforming...\n"
          push!(D, d[2])
          break
        end

        push!(cs, e)
        fl, l = isinteger(GC, B, e)
        if fl
          @vprint :GaloisGroup 2 "found descent at $t and value $l\n"
          push!(fd, t)
        end
      end
      if length(cs) == length(lt)
        pop!(D, d[2])
        if length(fd) > 0
          @vprint :GaloisGroup 2 "descending via $fd\n"
        else
          @vprint :GaloisGroup 2 "no descent\n"
        end
        if length(fd)>0
          push!(GC.chn, (G, I, ts, fd))
          G = intersect([s^x for x = fd]...)[1]
          @assert length(fd) > 1 || order(G) == order(s)
          @assert order(G) < nG
          @vprint :GaloisGroup 2 "descending to $(grp_id(G))\n"
          
          if order(G) == n
            GC.G = G
            return G, GC
          end
          break
        end
      end
    end
    if nG == order(G)
      GC.G = G
      return G, GC
    end
    nG = order(G)
  end      
  GC.G = G
  return G, GC
end

"""
    isinteger(C::GaloisCtx, B::BoundRingElem, v)

For an element `v` representing an integral polynomial evaluated at the
roots stored in `C`, known to be bounded from above by `B`, either return
`true` and an explicit (algebraic) integer in the base ring of the context or
return `false`.
"""
function isinteger end

function isinteger(GC::GaloisCtx{Hecke.qAdicRootCtx}, B::BoundRingElem{fmpz}, e)
  p = GC.C.p
  if e.length<2
    l = coeff(e, 0)
    lz = lift(l)
    lz = Hecke.mod_sym(lz, fmpz(p)^precision(l))
    if abs(lz) < value(B)
      return true, lz
    else
      return false, lz
    end
  else
    return false, fmpz(0)
  end
end

#TODO: also allow for
# p-adic, q-adic
# all locals
# finite fields
# rel. ext
# ...
function extension_field(f::fmpz_poly, n::String = "_a"; cached::Bool = true, check::Bool = true)
  return NumberField(f, n, cached = cached, check = check)
end
function extension_field(f::fmpq_poly, n::String = "_a"; cached::Bool = true, check::Bool = true)
  return NumberField(f, n, cached = cached, check = check)
end

function extension_field(f::Generic.Poly{Generic.Rat{T}}, n::String = "_a";  cached::Bool = true, check::Bool = true) where {T}
  return FunctionField(f, n, cached = cached)
end

function extension_field(f::Generic.Poly{nf_elem}, n::String = "_a";  cached::Bool = true, check::Bool = true) where {T}
  return NumberField(f, n, cached = cached)
end

Hecke.function_field(f::Generic.Poly{Generic.Rat{T}}, n::String = "_a";  cached::Bool = true, check::Bool = true) where {T} = FunctionField(f, n, cached = cached)


@doc Markdown.doc"""
Finds a Tschirnhausen transformation, ie a polynomial in `Zx` s.th.
`I` evaluated at the (roots in ) `r` does not have repetitions.

  ``|\{ I^s(t(r_1), ..., t(r_n)) | s in T\}| = |T|``
"""
function find_transformation(r, I::SLPoly, T::Vector{PermGroupElem})
  Zx = Hecke.Globals.Zx
  ts = gen(Zx)
  while true
    rt = map(ts, r)
    conj = [evaluate(I^t, rt) for t = T]
    if length(Set(conj)) == length(T)
      return ts
    end
    while true
      ts = rand(Zx, 2:rand(2:max(2, length(r))), -4:4) #TODO: try smaller degrees stronger
      if degree(ts) > 0
        break
      end
    end
  end
end

@doc Markdown.doc"""
    fixed_field(GC::GaloisCtx, U::PermGroup, extra::Int = 5)

Given the `GaloisCtx` as returned by a call to `galois_group` and a subgroup
`U` of the Galois group, compute the field fixed by `U` as a simple
extension.
"""
function fixed_field(GC::GaloisCtx, U::PermGroup, extra::Int = 5)
  G = GC.G
  if index(G, U) == 1 # not type stable
    return QQ
  end
  #XXX: seems to be broken for reducible f, ie. intransitive groups

  c = reverse(maximal_subgroup_chain(G, U))
  @vprint :GaloisGroup 2 "using a subgroup chain with orders $(map(order, c))\n"

  I = [invariant(c[i], c[i+1]) for i=1:length(c)-1]
  #the I[i] is a relative primitive element - it may be absolute or not...
  # right_transversal: U*g, thus G>U>V -> V a b where a runs (U/V), b (G/U):
  # G = cup U b, U = cup V a, so G = V a b
  ts = [gen(Hecke.Globals.Zx) for i = I]
  tv  = [right_transversal(c[i], c[i+1]) for i=1:length(c)-1]
  r = roots(GC, bound_to_precision(GC, GC.B))
  k, mk = ResidueField(parent(r[1]))
  r = map(mk, r)
  mu = ones(Int, length(I))
  #need tschirni per invar
  local conj
  ts[1] = find_transformation(r, I[1], tv[1])
  T = tv[1]
  if ts[1] == gen(Hecke.Globals.Zx)
    a = I[1]
  else
    a = evaluate(I[1], map(ts[1], gens(parent(I[1]))))
  end

  for j=2:length(I)
    local nc, rt
    ts[j] = find_transformation(r, I[j], tv[j])
    lT = length(T)
    for k=2:length(tv[j])
      append!(T, [tv[j][k] * t for t =  T[1:lT]])
    end

    mu[j] = 0

    if ts[j] == gen(Hecke.Globals.Zx)
      b = I[j]
    else
      b = evaluate(I[j], map(ts[j], gens(parent(I[1]))))
    end
    while true
      d = Set([evaluate((mu[j]*a+b)^t, r) for t = T])
      if length(d) == length(T)
        break
      end
      mu[j] += 1
    end
    if iszero(mu[j])
      a = b
    else
      a = mu[j]*a+b
    end
  end
  @vprint :GaloisGroup 2  "have primitive element via $mu \n$a\n"

  B = upper_bound(GC, a)
  m = length(T)
  B = m*B^m

  @vtime :GaloisGroup 2 r = roots(GC, bound_to_precision(GC, B, extra))
  compile!(a)
  @vtime :GaloisGroup 2 conj = [evaluate(a, t, r) for t = T]

  fl, val = isinteger(GC, B, sum(conj))
  @assert fl
  ps = [val]
  d = copy(conj)
  while length(ps) < m
#    @show length(ps)
    @vtime :GaloisGroup 2 d .*= conj
    fl, val = isinteger(GC, B, sum(d))
    @assert fl
    push!(ps, val)
  end
  ps = map(base_ring(GC.f), ps)

  k = extension_field(Hecke.power_sums_to_polynomial(ps), check = false, cached = false)[1]
  @assert all(x->isone(denominator(x)), coefficients(k.pol))
  @assert ismonic(k.pol)
  return k
end

"""
    galois_quotient(C::GaloisCtx, Q::PermGroup)

Finds all(?) subfields of the splitting field s.th. the galois
group will be permutation isomorphic to Q.
"""
function galois_quotient(C::GaloisCtx, Q::PermGroup)
  G = C.G
  if order(G) % degree(Q) != 0
    return []
  end
  s = subgroup_reps(G, order = divexact(order(G), degree(Q)))
  res = []
  for U = s
    phi = right_coset_action(G, U)
    if isisomorphic(Q, image(phi)[1])[1]
      push!(res, fixed_field(C, U))
    end
  end
  return res
end

"""
    galois_quotient(C::GaloisCtx, d::Int)

Finds all(?) subfields (up to isomorphism) of the splitting field of degree d
with galois group isomorphic to the original one.

# EXAMPLE
```jldoctest
julia> Qx, x = QQ["x"];
julia> G, C = galois_group(x^3-2);
julia> galois_quotient(C, 6)
1-element Vector{Any}:
 Number field over Rational Field with defining polynomial x^6 + 324*x^4 - 4*x^3 + 34992*x^2 + 1296*x + 1259716

julia> galois_group(ans[1])
(Group([ (), (1,5)(2,4)(3,6), (1,2,3)(4,5,6) ]), Galois Context for x^6 + 324*x^4 - 4*x^3 + 34992*x^2 + 1296*x + 1259716 and prime 13)

julia> isisomorphic(ans[1], G)
true
```
"""
function galois_quotient(C::GaloisCtx, d::Int)
  G = C.G
  if order(G) % d != 0
    return []
  end
  s = subgroup_reps(G, order = divexact(order(G), d))
  res = []
  for U = s
    phi = right_coset_action(G, U)
    if order(image(phi)[1]) == order(G)
      push!(res, fixed_field(C, U))
    end
  end
  return res
end

"""
    galois_quotient(C::GaloisCtx, d::Int, n::Int)

Finds all subfields of the splitting field with galois group the n-th
transitive group in degree d
"""
galois_quotient(C::GaloisCtx, d::Int, n::Int) = 
            galois_quotient(C, transitive_group(d, n))

"""
    galois_quotient(f::PolyElem, p::Vector{Int})

Equivalent to

    galois_quotient(galois_group(f)[2], p[1], p[2])

Finds all subfields of the splitting field of f with galois group
the p[2]-th transitive group of degree p[1]
"""
function galois_quotient(f::PolyElem, p::Vector{Int})
  G, C = galois_group(f)
  @assert length(p) == 2
  return galois_quotient(C, p[1], p[2])
end

#based on 
#  doi:10.1006/jsco.2000.0376
# Using Galois Ideals for Computing Relative Resolvents
# Aubrey & Valibouze
# JSC (2000) 30, 635-651

function cauchy_ideal(f::fmpz_poly; parent::MPolyRing = PolynomialRing(QQ, degree(f), cached = false)[1])
  return cauchy_ideal(f(gen(Hecke.Globals.Qx)), parent=parent)
end

@doc Markdown.doc"""
    cauchy_ideal(f::PolyElem{<:FieldElem})

The coefficients of `f` are the elementary symmetric functions evaluated
at the roots of `f`. The `cauchy_ideal` is the ideal generated
by the differences between the elementary symmetric functions and the
coefficients.

EXAMPLE
```jldoctest
julia> Qx, x = QQ["x"]
julia> cauchy_ideal(x^4-2)
ideal(x4^4 - 2, x3^3 + x3^2*x4 + x3*x4^2 + x4^3, x2^2 + x2*x3 + x2*x4 + x3^2 + x3*x4 + x4^2, x1 + x2 + x3 + x4)
```
"""
function cauchy_ideal(f::PolyElem{<:FieldElem}; parent::MPolyRing = PolynomialRing(base_ring(f), degree(f), cached = false)[1])
  x = gens(parent)
  n = degree(f)
  f = f(x[n])
  c = [f]
  for i=1:n-1
    f = divexact(evaluate(f, vcat(x[1:n-i], [x[n-i]], x[n-i+2:n])) - f, x[n-i]-x[n-i+1])
    push!(c, f)
  end
  return ideal(c)
end

@doc Markdown.doc"""
    galois_ideal(C::GaloisCtx, extra::Int = 5)

The so called Galois ideal is a descriptio of the splitting field of `f` as
a quotient by some maximal ideal. Algebraically, this ideal is an irreducible
component of the Cauchy ideal, the ideal geberated by the elementary symmetric
functions and the coefficients of the polynomial.

EXAMPLE
```jldoctest
julia> Qx, x = QQ["x"]
julia> i = galois_ideal(galois_group(x^4-2)[1])
ideal(x4^4 - 2, x3^3 + x3^2*x4 + x3*x4^2 + x4^3, x2^2 + x2*x3 + x2*x4 + x3^2 + x3*x4 + x4^2, x1 + x2 + x3 + x4, -x1*x2 - x1*x3 - x2*x4 - x3*x4)

julia> k, _ = number_field(i);
julia> length(roots(x^4-2, k))
4
```
"""
function galois_ideal(C::GaloisCtx, extra::Int = 5)
  f = C.f
  id = gens(cauchy_ideal(f))
  R = parent(id[1])
  x = gens(R)
  n = length(x)
  #we need to go down (by hand) to the starting group.
  #This is either C.G (if there was no descent in Stauduhar)
  #or the 1st group in the descent chain (C.chn)...
  #TODO: the subfields use, implicitly, special invariants, so
  #      we should be able to avoid the chain
  if length(C.chn) == 0
    c = maximal_subgroup_chain(symmetric_group(n), C.G)
  else
    c = maximal_subgroup_chain(symmetric_group(n), C.chn[1][1])
  end

  r = roots(C, bound_to_precision(C, C.B))
  k, mk = ResidueField(parent(r[1]))
  r = map(mk, r)

  for i=1:length(c)-1
    I = invariant(c[i+1], c[i]) 
    T = right_transversal(c[i+1], c[i])
    ts = find_transformation(r, I, T)
    B = upper_bound(C, I, ts)
    r = roots(C, bound_to_precision(C, B))
    r = map(ts, r)
    compile!(I)
    for t = T
      e = evaluate(I, t, r)
      fl, v = isinteger(C, B, e)
      if fl
        push!(id, v-evaluate(I, t, map(ts, x)))
        break
      end
    end
  end
  for (_, I, ts, T) = C.chn
    B = upper_bound(C, I, ts)
    r = roots(C, bound_to_precision(C, B, extra))
    r = map(ts, r)
    compile!(I)
    for t = T
      e = evaluate(I, t, r)
      fl, v = isinteger(C, B, e)
      @assert fl
      push!(id, v-evaluate(I, t, x))
    end
  end
  return ideal(id)
end

function Hecke.absolute_primitive_element(K::Oscar.NfNSGen{fmpq, fmpq_mpoly})
  while true
    a = rand(K, -10:10)
    f = minpoly(a)
    if degree(f) == degree(K)
      return a
    end
  end
end

function Hecke.absolute_minpoly(a::Oscar.NfNSGenElem{fmpq, fmpq_mpoly})
  return minpoly(a)
end

#TODO copied from MPolyFact in Hecke....
function find_morphism(k::FqNmodFiniteField, K::FqNmodFiniteField)
   if degree(k) > 1
    phi = Nemo.find_morphism(k, K) #avoids embed - which stores the info
  else
    phi = MapFromFunc(x->K((coeff(x, 0))), y->k((coeff(y, 0))), k, K)
  end
  return phi
end

function blow_up(G::PermGroup, C::GaloisCtx, lf::Vector, con::PermGroupElem=one(G))
  if all(x->x[2] == 1, lf)
    return G, C
  end
  
  n = degree(G)
  st = 0
  mp = Dict{Int, Int}(i => i for i=1:n)

  icon = inv(con)

  gs = map(Vector, gens(G))
  for (g, k) = lf
    for j=2:k
      for i=1:degree(g)
        mp[n+i] = (st+i)^con
      end
      for h = gs
        for i=1:degree(g)
          push!(h, h[(st+i)^con]^icon-st+n)
        end
      end
      n += degree(g)
    end
    st += degree(g)
  end

  C.rt_num = mp
  S = symmetric_group(n)
  GG, _ = sub(S, map(S, gs))

  h = hom(G, GG, gens(G), gens(GG))
  @assert isinjective(h) && issurjective(h)
  return GG, C
end


#TODO: do not move from fmpz_poly to fmpq_poly to fmpz_poly...
function galois_group(f::fmpz_poly; pStart::Int = 2*degree(f), prime::Int = 0)
  return galois_group(f(gen(Hecke.Globals.Qx)), pStart = pStart, prime = prime)
end

@doc Markdown.doc"""
    galois_group(f::PolyElem{<:FieldElem})

Computes the automorphism group of a splitting field of `f` as an explicit
group of permutations of the roots. Furthermore, the `GaloisCtx` is
returned allowing algorithmic access to the splitting field.
"""
function galois_group(f::PolyElem{<:FieldElem}; prime=0, pStart::Int = 2*degree(f))
  lf = [(k,v) for  (k,v) = factor(f).fac]

  if length(lf) == 1
    @vprint :GaloisGroup 1 "poly has only one factor\n"
    G, C = galois_group(number_field(lf[1][1], cached = false)[1])
    return blow_up(G, C, lf)
  end

  @vprint :GaloisGroup 1 "factoring input...gives\n$lf\n"

  lg = [x[1] for x = lf]
  
  g = prod(lg)
  p, ct = find_prime(g, pStart = pStart, prime = prime)

  C = [galois_group(number_field(x, cached = false)[1], prime = p)[2] for x = lg]
  G, emb, pro = inner_direct_product([x.G for x = C], morphisms = true)

  CC = GaloisCtx(g, p)
  rr = roots(CC, 5, raw = true) #raw is neccessary for non-monic case
                                #the scaling factor is the lead coeff
                                #thus not the same for all factors...
  @assert length(Set(rr)) == length(rr)

  d = map(frobenius, rr)
  si = [findfirst(y->y==x, rr) for x = d]

  @vprint :GaloisGroup 1 "found Frobenius element: $si\n"

  k, mk = ResidueField(parent(rr[1]))
  rr = map(mk, rr)
  po = Int[]
  for GC = C
    r = roots(GC, 5, raw = true)
    K, mK = ResidueField(parent(r[1]))
    r = map(mK, r)
    phi = find_morphism(K, k)
    po = vcat(po, [findfirst(x->x == phi(y), rr) for y = r])
  end

  con = (symmetric_group(length(po))(po))
  G = G^con
  F = GroupFilter()
#=
  function fi(x)
    @show [pro[y](x^inv(con))[1] == C[y].G for y=1:length(C)]
    return all(y->pro[y](x^inv(con))[1] == C[y].G, 1:length(C))
  end
  push!(F, fi)
=#  

  push!(F, x->all(y->pro[y](x^inv(con))[1] == C[y].G, 1:length(C)), "subdirect case: wrong projections")
  G, C =  descent(CC, G, F, G(si), grp_id = x->(:-,:-))
  return blow_up(G, C, lf, con)
end

function Nemo.cyclotomic(n::Int, x::fmpq_poly)
  return Nemo.cyclotomic(n, gen(Hecke.Globals.Zx))(x)
end


################################################################################
#
#  Promote rules
#
################################################################################

AbstractAlgebra.promote_rule(::Type{BoundRingElem}, ::Type{fmpz}) = BoundRingElem

AbstractAlgebra.promote_rule(::Type{BoundRingElem}, ::Type{T}) where {T <: Integer} = BoundRingElem

include("Group.jl")
include("POSet.jl")
include("Subfields.jl")
include("SeriesEval.jl")
include("Qt.jl")
include("RelGalois.jl")

end

using .GaloisGrp
export galois_group, slpoly_ring, elementary_symmetric, galois_quotient,
       power_sum, to_elementary_symmetric, cauchy_ideal, galois_ideal, fixed_field, 
       maximal_subgroup_reps, extension_field
       
#=
       M12: 2-transitive, hence msum is a waste
x^12 - 4*x^11 + 4*x^10 + 12*x^9 - 72*x^8 + 168*x^7 - 132*x^6 - 324*x^5 + 1197*x^4 - 1752*x^3 + 1500*x^2 - 672*x + 207 

several index 2 groups, uses generic invar
x^16 + 209*x^14 + 102*x^13 + 18138*x^12 + 13232*x^11 + 814855*x^10 + 673869*x^9 + 19699456*x^8 + 17373605*x^7 + 258261711*x^6 + 283233089*x^5 + 2368948579*x^4 + 2075376015*x^3 + 13009666150*x^2 + 14279830429*x + 37222748001


primitive, 2-transitive
x^10 - 2*x^9 - 263*x^8 + 1136*x^7 + 18636*x^6 - 120264*x^5 - 81916*x^4 + 1314656*x^3 - 290197*x^2 - 3135542*x + 2052019

needs the group theory test for isinvar in line 805 and related
x^16 + 4*x^15 + 20*x^14 + 44*x^13 + 106*x^12 + 120*x^11 + 180*x^10 + 44*x^9 + 134*x^8 - 120*x^7 + 172*x^6 - 32*x^5 + 598*x^4 + 312*x^3 + 616*x^2 + 147
=#
