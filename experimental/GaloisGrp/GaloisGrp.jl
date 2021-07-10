module GaloisGrp

using Oscar, Markdown
import Base: ^, +, -, *
import Oscar: Hecke, AbstractAlgebra, GAP
using Oscar: SLPolyRing, SLPoly, SLPolynomialRing

export galois_group, transitive_group_identification, slpoly_ring, elementary_symmetric,
       power_sum, to_elementary_symmetric, cauchy_ideal, galois_ideal, fixed_field

import Hecke: orbit, fixed_field


function __init__()
  GAP.Packages.load("ferret", install = true)

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
struct BoundRing  <: AbstractAlgebra.Ring
  mul
  add
  pow
  map
  name::String
end

struct BoundRingElem <: AbstractAlgebra.RingElem
  val::fmpz
  p::BoundRing # the parent
end

function Base.show(io::IO, b::BoundRingElem)
  print(io, "x <= $(b.val)")
end

function check_parent(a::BoundRingElem, b::BoundRingElem)
  parent(a) == parent(b) || error("Elements must have same parent")
  return true
end

+(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && BoundRingElem(a.p.add(a.val, b.val), a.p)
-(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && BoundRingElem(a.p.add(a.val, b.val), a.p)
*(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && BoundRingElem(a.p.mul(a.val, b.val), a.p)
^(a::BoundRingElem, b::Int) = BoundRingElem(a.p.pow(a.val, b), a.p)
-(a::BoundRingElem) = a

Oscar.parent(a::BoundRingElem) = a.p
value(a::BoundRingElem) = a.val
Base.isless(a::BoundRingElem, b::BoundRingElem) = check_parent(a, b) && isless(value(a), value(b))

(R::BoundRing)(a::fmpz) = BoundRingElem(R.map(abs(a)), R)
(R::BoundRing)(a::Integer) = BoundRingElem(fmpz(a), R)
(R::BoundRing)(a::BoundRingElem) = a
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
  return BoundRing( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x, x->x, "max-ring")
end

"""
Normal ring
"""
function add_ring()
  return BoundRing( (x,y) -> x*y, (x,y) -> x+y, (x,y) -> x^y, x->x, "add-ring")
end

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
  return BoundRing( (x,y) -> x+y+1, (x,y) -> x+y, (x,y) -> x+2*nbits(y), x->0, "cost-ring")
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
  return BoundRing( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x, x->0, "degree-ring")
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

Oscar.parent_type(::BoundRingElem) = BoundRing
Oscar.elem_type(::BoundRing) = BoundRingElem

#my 1st invariant!!!
@doc Markdown.doc"""
    sqrt_disc(a::Array{<:Any, 1})

The product of differences ``a[i] - a[j]`` for all indices ``i<j``.    
"""
function sqrt_disc(a::Array{<:Any, 1})
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
function elementary_symmetric(g::Array{<:Any, 1}, i::Int)
  return sum(prod(g[i] for i = s) for s = Hecke.subsets(Set(1:length(g)), i))
end

@doc Markdown.doc"""
    power_sums(g::Vector, i::Int) ->

Evaluates the `i`-th power sums at the values in `g`, ie. the sum
of the `i`-th power of the values.
"""
function power_sum(g::Array{<:Any, 1}, i::Int)
  return sum(a^i for a = g)
end

@doc Markdown.doc"""
    discriminant(g::Vector)

Compute the product of all differences of distinct elements in the array.    
"""
function Oscar.discriminant(g::Array{<:RingElem, 1})
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
  chn::Array{Tuple{PermGroup, SLPoly, fmpz_poly, Array{PermGroupElem, 1}}, 1}
  #= the descent chain, recodring
   - the group
   - the invariant
   - the tschirnhaus transformation
   - the cosets used
   should probably also record if the step was proven or not
   and the starting group  
  =#
  function GaloisCtx(f::fmpz_poly, p::Int)
    r = new{Hecke.qAdicRootCtx}()
    r.f = f
    r.C = Hecke.qAdicRootCtx(f, p, splitting_field = true)
    r.B = add_ring()(roots_upper_bound(f))
    r.chn = Tuple{PermGroup, SLPoly, fmpz_poly, Array{PermGroupElem, 1}}[]
    return r
  end
  #=
  Roots in F_q[[t]] for q = p^d
    - needs to be tweaked to do Q_q[[t]
    - need q-lifting as well as t-lifting
    - possibly also mul_ks for Q_q[[t]] case
    - needs merging in Hecke
  =#

  function GaloisCtx(f::fmpz_mpoly, p::Int, d::Int)
    @assert ngens(parent(f)) == 2
    Qq, _ = QadicField(p, d, 10)
    F, mF = ResidueField(Qq)
    H = Hecke.MPolyFact.HenselCtxFqRelSeries(f, F)
    SQq, _ = PowerSeriesRing(Qq, 2, "s", cached = false)
    SQqt, _ = PolynomialRing(SQq, "t", cached = false)
    mc(f) = map_coefficients(x->map_coefficients(y->setprecision(preimage(mF, y), 1), x, parent = SQq), f, parent = SQqt)
    HQ = Hecke.MPolyFact.HenselCtxFqRelSeries(H.f, map(mc, H.lf), map(mc, H.cf), H.n)
    r = new{Hecke.MPolyFact.HenselCtxFqRelSeries{AbstractAlgebra.Generic.RelSeries{qadic}}}()
    Qt, t = PolynomialRing(QQ, "t", cached = false)
    Qts, s = PolynomialRing(Qt, "s", cached = false)
    r.f = evaluate(f, [Qts(t), s])
    r.C = HQ
    #r.B = more complicated: needs degree (inf. val.) bound as well as coeffs
    r.chn = Tuple{PermGroup, SLPoly, fmpz_poly, Array{PermGroupElem, 1}}[]
    return r
  end
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
"""
function Hecke.roots(G::GaloisCtx{Hecke.qAdicRootCtx}, pr::Int)
  a = Hecke.roots(G.C, pr)::Array{qadic, 1}
  return Hecke.expand(a, all = true, flat = false, degs = Hecke.degrees(G.C.H))::Array{qadic, 1}
end
function Hecke.roots(G::GaloisCtx{<:Hecke.MPolyFact.HenselCtxFqRelSeries}, pr)
  C = G.C
  while precision(C)[1] < pr[1]
    Hecke.MPolyFact.lift_q(C)
  end
  while precision(C)[2] < pr[2]
    Hecke.MPolyFact.lift(C)
  end
  return [-coeff(x, 0) for x = C.lf[1:C.n]]
end
#=
too simplistic. RootCtx computes roots in F_q[[t]], but currently not
in Q_q[[t]]

function Hecke.roots(G::GaloisCtx{Main.MPolyFact.RootCtx}, pr::Int)
  a = Main.MPolyFact.root(G.C, 1, 1)
  return a
end
=#


@doc Markdown.doc"""
    upper_bound(G::GaloisCtx, f...)

Given a `GaloisCtx` and some multivariate function, upper_bound the image of `f`
upon evaluation at the roots implicit in `G`.

`f` can be
 - a multivariate polynomial or straight-line polynomial (strictly: any object
   allowing `evaluate`
 - `elementary_symmetric` or `power_sum`, in which case more arguments are
   needed: the array with the values and the index.
   `upper_bound(G, power_sum, A, i)` is equivalent to `bound(G, power_sum(A, i))`
   but more efficient.

In every case a univariate polynomial (over the integers) can be added, it
will act as a Tschirnhaus-transformation, ie. the roots (bounds) implicit
in `G` will first be transformed.
"""
function upper_bound end

function upper_bound(G::GaloisCtx, f)
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
some overhead inthe compilation itself. By default, aparently nothing is
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
#- ansehen der ZykelTypen um Sn/An zu erkennen (as well as to avoid 2-sum)
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
      @show ex += 1
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
      B = block_system(H, BB)
      m = length(B)
      l = length(B[1])
      y = [sum(g[b]) for b = B]
      d = [sqrt_disc(g[b]) for b = B]
      D = sqrt_disc(y)
      I = D
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G)
        @vprint :GaloisInvariant 3 "using D-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, 1)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G)
        @vprint :GaloisInvariant 3 "using s1-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, m)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G)
        @vprint :GaloisInvariant 3 "using sm-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, 2)
      if all(p->isprobably_invariant(I, p), H) &&
         any(p->!isprobably_invariant(I, p), G)
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

      sG = set_stabilizer(G, BB)[1]
      sH = set_stabilizer(H, BB)[1]
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
  rt = roots(C, 5)
  #make square-free (in residue field)
  k, mk = ResidueField(parent(rt[1]))
  k_rt = map(mk, rt)
  ts = find_transformation(k_rt, I, t)

  B = 2*n*evaluate(I, map(ts, [C.B for i = 1:ngens(parent(I))]))^n
  rt = roots(C, clog(value(B), C.C.p)+extra)
  rt = map(ts, rt)
  rt = [evaluate(I^s, rt) for s = t]
  pr = copy(rt)
  ps = fmpq[isinteger(C, B, sum(pr))[2]]
  while length(ps) < n
    pr .*=  rt
    push!(ps, isinteger(C, B, sum(pr))[2])
  end
  return Hecke.power_sums_to_polynomial(ps)
end

struct GroupFilter
  f::Array{Function, 1}
  GroupFilter() = new([x->true])
end

function (F::GroupFilter)(G::PermGroup)
  return all(x->x(G), F.f)
end

function Base.push!(G::GroupFilter, F::Function)
  push!(G.f, F)
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
  s::Array{PermGroup, 1}
  I::Dict{Int, SLPoly}
  T::Dict{Int, Array{fmpz_poly, 1}}
  l::Array{Int, 1}
  #the coset reps need to be cached as well
  #a "work limit" on the "invariant" function
  #a more select choice of group....

  function DescentEnv(G::PermGroup, f::GroupFilter = GroupFilter())
    s = maximal_subgroup_reps(G)
    r = new()
    r.G = G
    r.s = filter(f, s)
    r.I = Dict{Int, SLPoly}()
    r.T = Dict{Int, Array{fmpz_poly, 1}}()
    r.l = zeros(Int, length(r.s))
    return r
  end
end

#TODO: instead of iterating, use the cheapest invariant
#      use invariant with a worklevel to limit costs in finding them
#      use the V_4 (index 2) stuff when neccessary
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

@doc Markdown.doc"""
    starting_group(GC::GaloisCtx, K::AnticNumberField) 

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
function starting_group(GC::GaloisCtx, K::AnticNumberField; useSubfields::Bool = true)
  c = roots(GC, 5)

  @vprint :GaloisGroup 1 "computing starting group (upper bound for Galois group)\n"

  if useSubfields
    @vprint :GaloisGroup 1 "computing principal subfields ...\n"
    @vtime :GaloisGroup 2 S = Hecke.principal_subfields(K)
  else
    S = []
  end

  #compute the block system for all subfields...
  bs = Array{Array{Int, 1}, 1}[]
  for (s, ms) = S
    if degree(s) == degree(K) || degree(s) == 1
      continue
    end
    pr = 5
    local v
    while true
      c = roots(GC, pr)

      g = ms(gen(s))
      d = map(parent(K.pol)(g), c)
      b = Dict{typeof(c[1]), Array{Int, 1}}()
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
    sort!(v, lt = (x,y) -> minimum(x) < minimum(y))
    push!(bs, v)
  end

  F = GroupFilter()
  push!(F, istransitive) #poly is defining number field, hence irreducible


  #selecting maximal block systems only...
  if length(bs) > 1
    @vprint :GaloisGroup 2 "group will have (all) block systems: $([x[1] for x = bs])\n"
    L = POSet(bs, (x,y) -> issubset(x[1], y[1]) || issubset(y[1], x[1]),
                  (x,y) -> issubset(x[1], y[1]) - issubset(y[1], x[1]))
    bs = minimal_elements(L)
    @vprint :GaloisGroup 1 "group will have (maximal) block systems: $([x[1] for x = bs])\n"
  end

  d = map(frobenius, c)
  si = [findfirst(y->y==x, c) for x = d]

  @vprint :GaloisGroup 1 "found Frobenius element: $si\n"

  G = symmetric_group(degree(K))
  if degree(K) == 1
    return G, F, one(G)
  end

  S = G
  for b = bs
    W = isomorphic_perm_group(wreath_product(symmetric_group(length(b[1])), symmetric_group(length(b))))[1]
    s = S(vcat(b...))
    # W^s is largest group having "b" as a block system
    # courtesy of Max...
    G = intersect(G, W^s)[1]
  end

  if issquare(discriminant(K))
    G = intersect(G, alternating_group(degree(K)))[1]
    push!(F, iseven)
  else
    push!(F, isodd)
  end


  if length(bs) == 0 #primitive case: no subfields, no blocks, primitive group!
    push!(F, isprimitive)
    @vprint :GaloisGroup 1 "group is primitive (no subfields), trying operations on pairs\n"

    #starting group as the stabilizer of the factorisation of the 2-sum poly,
    #the polynomial with roots r_i + r_j, i<j
    #we need this square-free, so we compute this over the finite field
    #actually, we only check that the roots over the finite field are distinct.
    #if not: transform (using ts) and try again.
    
    k, mk = ResidueField(parent(c[1]))
    m = Dict{fq_nmod, Tuple{Int, Int}}()
    local ts = gen(Hecke.Globals.Zx)
    while true
      cc = map(mk, c)
      for i=1:length(c)-1
        for j=i+1:length(c)
          m[cc[i]+cc[j]] = (i,j)
        end
      end
      if length(keys(m)) < binomial(degree(K), 2)
        @vprint :GaloisGroup 2 " m-sum: found duplicate, transforming...\n"
        while true
          ts = rand(Hecke.Globals.Zx, 2:degree(K), -4:4)
          if degree(ts) > 1
            break
          end
        end
        c = map(ts, roots(GC, 5))
        empty!(m)
      else
        break
      end
    end

    @vprint :GaloisGroup 1 "have everything, now getting the 2-sum poly\n"
    if gen(Hecke.Globals.Zx) == ts
      @vtime :GaloisGroup 2 g = msum_poly(K.pol, 2)
    else
      @vtime :GaloisGroup 2 g = msum_poly(minpoly(ts(gen(K))), 2)
    end
    @vprint :GaloisGroup 1 "... factoring...\n"
    @vtime :GaloisGroup 2 fg  = factor(g)
    @assert all(isone, values(fg.fac))

    O = []
    for f = keys(fg.fac)
      r = roots(f, k)
      push!(O, [m[x] for x = r])
    end
    #the factors define a partitioning of pairs, the stabiliser of this
    #partition is the largest possible group...
    #code from Max...

    @vprint :GaloisGroup 2 "partitions: $O\n"
    #TODO: wrap this properly
    G = Oscar._as_subgroup(G, GAP.Globals.Solve(GAP.julia_to_gap(vcat(GAP.Globals.ConInGroup(G.X), [GAP.Globals.ConStabilize(GAP.julia_to_gap(sort(o), Val(true)), GAP.Globals.OnSetsSets) for o in O]))))[1]
    #@show G = intersect([stabilizer(G, GAP.julia_to_gap(sort(o), Val(true)), GAP.Globals.OnSetsSets)[1] for o=O]...)[1]
  end

  si = G(si)

  return G, F, si
end

@doc Markdown.doc"""
    find_prime(f::fmpq_poly, extra::Int = 5; pStart::Int = 2*degree(f)) -> p, c

Tries to find a useful prime for the computation of Galois group. Useful means
  - degree of the splitting field not too small
  - degree of the splitting field not too large
Starts searching at `pStart`, returns the prime as well as the cycle types identified.

If `pStart` is negative, `-pStart` is used and no search performed.
"""
function find_prime(f::fmpq_poly, extra::Int = 5; pStart::Int = 2*degree(f))
  if pStart < 0
    p = -pStart
    lf = factor(f, GF(p))
    return p, [sort(map(degree, collect(keys(lf.fac))))]
  end

  d_max = degree(f)^2
  d_min = min(2, d_max)
  p_best = 1
  cnt = 5
  ct = Set{Array{Int, 1}}()
  ps = Array{Tuple{Int, Int}, 1}()

  #find a q-adic splitting field of "good degree":
  # - too small, then the Frobenius automorphisms is not containing lots of
  #   information
  # - too large, all is slow.
  # careful: group could be (C_2)^n hence d_min might be small...
  no_p = 0
  for p = Hecke.PrimesSet(pStart, -1)
    lf = factor(f, GF(p))
    if any(x->x>1, values(lf.fac))
      continue
    end
    no_p +=1 
    push!(ct, sort(map(degree, collect(keys(lf.fac)))))
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
    if length(ps) > degree(f)
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
  return p_best, ct
end


#TODO
# - split: find primes (done) and abort on Sn/An
# - subfield lattice directly
# - subfield to block system outsource (also for lattice)
# - more invars
# - re-arrange to make cheap this first, expensive later
# - proof, live and later
# - "isInt", ie. the test for being in Z: outsource... (done(ish))
# - tschirnhausen transformation: slowly increasing degree and coeffs.
# - for larger degrees: improve starting group by doing Galois groups of subfields
# - make the descent a seperate function (done)
# - for reducibles...: write for irr. poly and for red. poly (done)
# - more base rings
# - applications: subfields of splitting field (done), towers, solvability by radicals
@doc Markdown.doc"""
    galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true, pStart::Int = 2*degree(K)) -> PermGroup, GaloisCtx

Computes the Galois group of the splitting field of the defining polynomial of `K`.
Currently the polynomial needs to be monic.

The group is returned as an explicit permutation group permuting the roots as contained
in the contex object (the 2nd return value). The roots live in a suitable unramifed
extension of the p-adics.
"""
function galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true, pStart::Int = 2*degree(K))

  p, ct = find_prime(K.pol, pStart = pStart)

  GC = GaloisCtx(Hecke.Globals.Zx(K.pol), p)

  G, F, si = starting_group(GC, K, useSubfields = useSubfields)

  return descent(GC, G, F, si, extra = extra)

end

@doc Markdown.doc"""
    descent(GC::GaloisCtx, G::PermGroup, F::GroupFilter, si::PermGroupElem; grp_id = transitive_group_identification, extra::Int = 5)

Performs a generic Stauduhar descent: starting with the group `G` that needs to be a 
supergroup of the Galois group, operating on the roots in `GC`, the context object.

The groups are filtered by `F` and the result needs to contain the permutation `si`.
For verbose output, the groups are printed through `grp_id`.
"""
function descent(GC::GaloisCtx, G::PermGroup, F::GroupFilter, si::PermGroupElem; grp_id = transitive_group_identification, extra::Int = 5)
  @vprint :GaloisGroup 2 "Have starting group with id $(grp_id(G))\n"

  p = GC.C.p
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

      c = roots(GC, clog(2*value(B), p)+extra)
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
        if e.length<2
          l = coeff(e, 0)
          lz = lift(l)
          l = Hecke.mod_sym(lz, fmpz(p)^precision(l))
          if abs(l) < value(B)
            @vprint :GaloisGroup 2 "found descent at $t and value $l\n"
            push!(fd, t)
          else
            @vprint :GaloisGroup 3 "no int: wrong size $l, $(value(B))\n"
          end
        else
          @vprint :GaloisGroup 3 "no int: wrong degree\n"
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

#=
DNW (does not work)

function subdir_invars(G, H)
  GH, emb, pro = inner_direct_product(G, H, morphisms = true)
  m = maximal_subgroup_reps(GH)
  m = [x for x = m if pro[1](x)[1] == G && pro[2](x)[1] == H]
  if length(m) == 0
    return m, []
  end
  II = SLPoly[]
  S, g = slpoly_ring(ZZ, degree(G)+degree(H))
  for x = m
    _, mx = Oscar._as_subgroup(GH, x.X)
    A = pro[1](intersect(x, emb[1](G)[1])[1])[1]
    B = pro[2](intersect(x, emb[2](H)[1])[1])[1]
    AB = inner_direct_product(A, B)
    #so G/A iso H/B iso x/(A x B)
    I = invariant(G, A)
    J = invariant(H, B)
    IJ = evaluate(I, g[1:degree(G)])+evaluate(J, rand(H), 2 .* g[degree(G)+1:end])
    push!(II, IJ) #sum(probable_orbit(x, IJ)))
  end
  return m, II
  s = maximal_subgroup_reps(m[1])
  s = [x for x = s if pro[1](x)[1] == G && pro[2](x)[1] == H]
  _, mx = Oscar._as_subgroup(m[1], s[1].X)  

    A = pro[1](intersect(s[1], emb[1](G)[1])[1])[1]
    B = pro[2](intersect(s[1], emb[2](H)[1])[1])[1]
    AB = inner_direct_product(A, B)
    I = invariant(G, A)
    J = invariant(H, B)
    IJ = evaluate(I, g[1:degree(G)])+evaluate(J, g[degree(G)+1:end])
    push!(II, sum(probable_orbit(m[1], IJ)))
    push!(m, s[1])

  return m, II
end

=#

#TODO: use above as well.
function isinteger(GC::GaloisCtx, B, e)
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

@doc Markdown.doc"""
Finds a Tschirnhausen transformation, ie a polynomial in `Zx` s.th.
`I` evaluated at the (roots in ) `r` does not have repetitions.

  ``|\{ I^s(t(r_1), ..., t(r_n)) | s in T\}| = |T|``
"""
function find_transformation(r, I::SLPoly, T::Array{PermGroupElem, 1})
  Zx = Hecke.Globals.Zx
  ts = gen(Zx)
  while true
    rt = map(ts, r)
    conj = [evaluate(I^t, rt) for t = T]
    if length(Set(conj)) == length(T)
      return ts
    end
    while true
      ts = rand(Zx, 2:length(r), -4:4)
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
  c = reverse(maximal_subgroup_chain(G, U))
  @vprint :GaloisGroup 2 "using a subgroup chain with orders $(map(order, c))\n"

  I = [invariant(c[i], c[i+1]) for i=1:length(c)-1]
  #the I[i] is a relative primitive element - it may be absolute or not...
  # right_transversal: U*g, thus G>U>V -> V a b where a runs (U/V), b (G/U):
  # G = cup U b, U = cup V a, so G = V a b
  ts = [gen(Hecke.Globals.Zx) for i = I]
  tv  = [right_transversal(c[i], c[i+1]) for i=1:length(c)-1]
  r = roots(GC, 5)
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
    a = mu[j]*a+b
  end
  @vprint :GaloisGroup 2  "have primitive element via $mu \n$a\n"

  B = upper_bound(GC, a)
  m = length(T)

  @vtime :GaloisGroup 2 r = roots(GC, clog(2*value(m*B^m), GC.C.p)+extra)
  compile!(a)
  @vtime :GaloisGroup 2 conj = [evaluate(a, t, r) for t = T]

  B = 2*m*B^m

  ps = fmpq[isinteger(GC, B, sum(conj))[2]]
  d = copy(conj)
  while length(ps) < m
#    @show length(ps)
    @vtime :GaloisGroup 2 d .*= conj
    push!(ps, isinteger(GC, B, sum(d))[2])
  end

  return number_field(Hecke.power_sums_to_polynomial(ps), check = false, cached = false)[1]
end

#based on 
#  doi:10.1006/jsco.2000.0376
# Using Galois Ideals for Computing Relative Resolvents
# Aubrey & Valibouze
# JSC (2000) 30, 635-651

function cauchy_ideal(f::fmpz_poly; parent::MPolyRing = PolynomialRing(QQ, degree(f), cached = false)[1])
  return cauchy_ideal(f(gen(Hecke.Globals.Qx)), parent=parent)
end

function cauchy_ideal(f::fmpq_poly; parent::MPolyRing = PolynomialRing(QQ, degree(f), cached = false)[1])
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

function galois_ideal(C::GaloisCtx, extra::Int = 5)
  f = C.f
  id = gens(cauchy_ideal(f))
  R = parent(id[1])
  x = gens(R)
  n = length(x)
  #we need to go down (by hand) to the starting group.
  #This is either C.G (if there was no descent in Stauduhar)
  #or the 1st group in the descent chain (C.chn)...
  if length(C.chn) == 0
    c = maximal_subgroup_chain(symmetric_group(n), C.G)
  else
    c = maximal_subgroup_chain(symmetric_group(n), C.chn[1][1])
  end

  r = roots(C, 5)
  k, mk = ResidueField(parent(r[1]))
  r = map(mk, r)

  for i=1:length(c)-1
    I = invariant(c[i+1], c[i]) 
    T = right_transversal(c[i+1], c[i])
    ts = find_transformation(r, I, T)
    B = upper_bound(C, I, ts)
    r = roots(C, clog(2*value(B), C.C.p)+extra)
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
    r = roots(C, clog(2*value(B), C.C.p)+extra)
    r = map(ts, r)
    compile!(I)
    for t = T
      e = evaluate(I, t, r)
      fl, v = isinteger(C, B, e)
      @assert fl
      push!(id, v-evaluate(I, t, x))
    end
  end
  return id
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

#TODO: do not move from fmpz_poly to fmpq_poly to fmpz_poly...
function galois_group(f::fmpz_poly; pStart::Int = 2*degree(f))
  return galois_group(f(gen(Hecke.Globals.Qx)), pStart = pStart)
end

function galois_group(f::fmpq_poly; pStart::Int = 2*degree(f))
  lf = factor(f)
  if any(x-> x > 1, values(lf.fac)) #think about this: do we drop the multiplicity or not
    error("polynomial must be squarefree")
  end

  if length(keys(lf.fac)) == 1
    @vprint :GaloisGroup 1 "poly irreducible\n"
    return galois_group(number_field(f, cached = false)[1])
  end

  @vprint :GaloisGroup 1 "factoring input...gives\n$lf\n"

  lg = sort(collect(keys(lf.fac)), lt = (a,b) -> isless(degree(b), degree(a)))
  #problem: inner_direct_product is "dropping" trivial factors.
  #trivial factors at the end are harmless, so I sort...
  #Max and/or Thomas are putting an option into Gap to not remove the
  #trivial factors (actually in general it is more complicated)
  g = prod(lg)
  p, ct = find_prime(g, pStart = pStart)

  C = [galois_group(number_field(x, cached = false)[1], pStart = -p)[2] for x = lg]
  G, emb, pro = inner_direct_product([x.G for x = C], morphisms = true)

  CC = GaloisCtx(Hecke.Globals.Zx(g), p)
  rr = roots(CC, 5)
  @assert length(Set(rr)) == length(rr)

  d = map(frobenius, rr)
  si = [findfirst(y->y==x, rr) for x = d]

  @vprint :GaloisGroup 1 "found Frobenius element: $si\n"

  k, mk = ResidueField(parent(rr[1]))
  rr = map(mk, rr)
  po = Int[]
  for GC = C
    r = roots(GC, 5)
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
  push!(F, x->all(y->pro[y](x^inv(con))[1] == C[y].G, 1:length(C)))
  return descent(CC, G, F, G(si), grp_id = x->(:-,:-))
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
include("SeriesEval.jl")
include("Qt.jl")

end

using .GaloisGrp
export galois_group, transitive_group_identification, slpoly_ring, elementary_symmetric,
       power_sum, to_elementary_symmetric, cauchy_ideal, galois_ideal, fixed_field
