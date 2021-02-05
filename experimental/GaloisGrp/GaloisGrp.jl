module GaloisGrp

using Oscar, Markdown
import Base: ^, +, -, *
import Oscar: Hecke, AbstractAlgebra, GAP
using Oscar: SLPolyRing, SLPoly, SLPolynomialRing

export galois_group, transitive_group_identification, slpoly_ring, elementary_symmetric,
       power_sum, to_elementary_symmetric

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
"""
function max_ring()
  return BoundRing( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x, x->x, "max-ring")
end

function add_ring()
  return BoundRing( (x,y) -> x*y, (x,y) -> x+y, (x,y) -> x^y, x->x, "add-ring")
end

"""
An slpoly evaluated at `cost_ring` elements `0` will count the number
of multiplications involved. A measure of the cost of evaluation at more
interesting scalars.
"""
function cost_ring()
  return BoundRing( (x,y) -> x+y+1, (x,y) -> x+y, (x,y) -> x+2*nbits(y), x->0, "cost-ring")
end

"""
An slpoly evaluated at `degree_ring` elements `1` will bound the total degree
from above.
"""
function degree_ring()
  return BoundRing( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x, x->0, "degree-ring")
end

@doc Markdown.doc"""
    cost(I::SLPoly)

Counts the number of multiplcations to evaluate `I`
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


function slpoly_ring(R::AbstractAlgebra.Ring, n::Int)
  return SLPolynomialRing(R, [ Symbol("x_$i") for i=1:n])
end

function slpoly_ring(R::AbstractAlgebra.Ring, p::Pair{Symbol, UnitRange{Int}}...)
  return SLPolynomialRing(R, p...)
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
"""
mutable struct GaloisCtx
  f::PolyElem
  C
  B::BoundRingElem
  G::PermGroup
  function GaloisCtx(f::fmpz_poly, p::Int)
    r = new()
    r.f = f
    r.C = Hecke.qAdicRootCtx(f, p, splitting_field = true)
    r.B = add_ring()(roots_upper_bound(f))
    return r
  end
end

function Base.show(io::IO, GC::GaloisCtx)
  print(io, "Galois Context for $(GC.f) and prime $(GC.C.p)")
end

@doc Markdown.doc"""
    roots(G::GaloisCtx, pr::Int)

The roots of the polynomial used to define the Galois-context in the fixed order
used in the algorithm. The roots are returned up to a precision of `pr`
p-adic digits, thus they are correct modulo ``p^pr``
"""
function Hecke.roots(G::GaloisCtx, pr::Int)
  a = Hecke.roots(G.C, pr)
  return Hecke.expand(a, all = true, flat = false, degs = Hecke.degrees(G.C.H))
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
#I assume algo can be imporved (TODO: Max?)
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
#TODO: think about the order!

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
  S, g = slpoly_ring(ZZ, degree(G))

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

  S, g = slpoly_ring(ZZ, degree(G))

  if istransitive(G) && !istransitive(H)
    @vprint :GaloisInvariant 2 "top group transitive, bottom not\n"
    return sum(probable_orbit(H, g[1]))
  end

  if !istransitive(G) 
    @vprint :GaloisInvariant 2 "both groups are intransitive\n"
    OG = orbits(G)
    OH = orbits(H)
    d = setdiff(OH, OG)
    if length(d) > 0
      @vprint :GaloisInvariant 2 "groups have different orbits\n"
      return sum(probable_orbit(H, d[1][1]))
    end
    #OH == OG
    for o = OH
      h = action_homomorphism(G, Int[x for x = o])
      hG = image(h)[1]
      hH = image(h, H)[1]
      if order(hG) > order(hH)
        @vprint :GaloisInvariant 2 "differ on action on $o, recursing\n"
        @hassert :GaloisInvariant 0 ismaximal(hG, hH)
        I = invariant(hG, hH)
        return evaluate(I, g[o])
      end
    end
    if order(H) < 100
      @vprint :GaloisInvariant 2 "group small enough, using orbit sum\n"
      return sum(probable_orbit(H, prod(g[i]^i for i=1:length(g)-1)))
    end
    error("give up")
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
      if all(p->isprobably_invariant(I, p), gens(H)) &&
         any(p->!isprobably_invariant(I, p), gens(G))
        @vprint :GaloisInvariant 3 "using D-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, 1)
      if all(p->isprobably_invariant(I, p), gens(H)) &&
         any(p->!isprobably_invariant(I, p), gens(G))
        @vprint :GaloisInvariant 3 "using s1-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, m)
      if all(p->isprobably_invariant(I, p), gens(H)) &&
         any(p->!isprobably_invariant(I, p), gens(G))
        @vprint :GaloisInvariant 3 "using sm-invar for $BB\n"
        return I
      end
      I = elementary_symmetric(d, 2)
      if all(p->isprobably_invariant(I, p), gens(H)) &&
         any(p->!isprobably_invariant(I, p), gens(G))
        @vprint :GaloisInvariant 3 "using s2-invar for $BB\n"
        return I
      end
      I = D*elementary_symmetric(d, m)
      if all(p->isprobably_invariant(I, p), gens(H)) &&
         any(p->!isprobably_invariant(I, p), gens(G))
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
      if length(sH) < length(sH)
        J = invar(sG, sH)
        C = left_transversal(H, sH)
        gg = g[BB]
        F = sum(evaluate(J, [gg[t(i)] for i = BB]) for t = C)
        @vprint :GaloisInvariant 3 "using F-invar for $BB (4.1.4)\n"
        return F
      end
    end
  end

  @vprint :GaloisInvariant 3 "no special invar found, resorting to generic\n"

  m = prod(gen(S, i)^i for i=1:degree(G)-1)
  return sum(m^s for s = H)


  @vprint :GaloisInvariant 2 "bad case... all failed, going generic\n"
  R, _ = PolynomialRing(ZZ, degree(G))
  m = prod(gen(R, i)^i for i=1:degree(G)-1)
  I = sum(orbit(H, m))
  return evaluate(I, g)
end

@doc Markdown.doc"""
    resolvent(C::GaloisCtx, G::PermGroup, U::PermGroup)

Find a `G`-relative `H`-invariant `I` and form the correspondig resolvent polynomial
```\prod (x-I^t)```
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
  ts = gen(Hecke.Globals.Zx)
  while length(Set([evaluate(I^s, k_rt, ts) for s = t])) < length(t)
    ts = rand(Hecke.Globals.Zx, 2:n, -4:4)
  end

  B = 2*n*evaluate(I, [C.B for i = 1:ngens(parent(I))], ts)^n
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

  @vprint :GaloisGroup 1 "computing starting group (upper bound for Galois group\n"

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
  

#TODO
# - split: find primes and abort on Sn/An
# - subfield lattice directly
# - subfield to block system outsource (also for lattice)
# - more invars
# - re-arrange to make cheap this first, expensive later
# - proof, live and later
# - "isInt", ie. the test for being in Z: outsource...
# - tschirnhausen transformation: slowly increasing degree and coeffs.
# - for larger degrees: improve starting group by doing Galois groups of subfields
# - make the descent a seperate function
# - for reducibles...: write for irr. poly and for red. poly
# - more base rings
# - applications: subfields of splitting field, towers, solvability by radicals
function galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true)
  d_min = 2
  d_max = typemax(Int)
  p_best = 1
  cnt = 5
  ct = Set{Array{Int, 1}}()
  #find a q-adic splitting field of "good degree":
  # - too small, then the Frobenius automorphisms is not comtaining lots of
  #   information
  # - too large, all is slow.
  # careful: group could be (C_2)^n hence d_min might be small...
  for p = Hecke.PrimesSet(2*degree(K), -1)
    lf = factor(K.pol, GF(p))
    if any(x->x>1, values(lf.fac))
      continue
    end
    push!(ct, sort(map(degree, collect(keys(lf.fac)))))
    d = lcm([degree(x) for x = keys(lf.fac)])
    if d < d_max && d >= d_min
      d_max = d
      p_best = p
      cnt = 5
    else
      cnt -= 1
    end
    if cnt < 1 && p_best > 1
      break
    end
  end

  p = p_best

  @vprint :GaloisGroup 1 "using prime $p_best with degree $d_max\n"
  @vprint :GaloisGroup 2 "and cycle types $ct\n"

  GC = GaloisCtx(Hecke.Globals.Zx(K.pol), p_best)

  G, F, si = starting_group(GC, K, useSubfields = useSubfields)

  @vprint :GaloisGroup 2 "Have starting group with id $(transitive_group_identification(G))\n"

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

      @vprint :GaloisGroup 2 "testing descent $(transitive_group_identification(G)) -> $(transitive_group_identification(s)) of index $(index(G, s))\n"

      @hassert :GaloisInvariant 1 all(x->isprobably_invariant(I, x), gens(s))
      @hassert :GaloisInvariant 1 any(x->!isprobably_invariant(I, x), gens(G))

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
          G = intersect([s^x for x = fd]...)[1]
          @assert length(fd) > 1 || order(G) == order(s)
          @assert order(G) < nG
          @vprint :GaloisGroup 2 "descending to $(transitive_group_identification(G))\n"
          if order(G) == degree(K)
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
      error("asd")
      return false, lz
    end
  else
    error("asdasd")
    return false, fmpz(0)
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
  c = reverse(maximal_subgroup_chain(G, U))
  @vprint :GaloisGroup 2 "using a subgroup cahin with orders $(map(order, c))\n"

  I = [invariant(c[i], c[i+1]) for i=1:length(c)-1]
  #the I[i] is a relative primitive element - it may be absolute or not...
  # right_transversal: U*g, thus G>U>V -> V a b where a runs (U/V), b (G/U):
  # G = cup U b, U = cup V a, so G = V a b
  ts = [gen(Hecke.Globals.Zx) for i = I]
  tv  = [right_transversal(c[i], c[i+1]) for i=1:length(c)-1]
  r = roots(GC, 5)
  mu = ones(Int, length(I))
  #need tschirni per invar
  local conj
  while true
    rt = map(ts[1], r)
    conj = [evaluate(I[1]^t, rt) for t = tv[1]]
    if length(Set(conj)) == order(c[1])//order(c[2])
      break
    end
    while true
      ts[1] = rand(parent(ts[1]), 2:degree(GC.f), -4:4)
      if degree(ts[1]) > 0
        break
      end
    end
  end
  T = tv[1]
  if ts[1] == gen(Hecke.Globals.Zx)
    a = I[1]
  else
    a = evaluate(I[1], map(ts[1], gens(parent(I[1]))))
  end
  for j=2:length(I)
    local nc, rt
    while true
      rt = map(ts[j], r)
      nc = [evaluate(I[j]^t, rt) for t = tv[j]]
      if length(Set(nc)) == order(c[j])//order(c[j+1])
        break
      end
      while true
        ts[j] = rand(parent(ts[1]), 2:degree(GC.f), -4:4)
        if degree(ts[j]) >= 1
          break
        end
      end
    end
    
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

  r = roots(GC, clog(2*value(m*B^m), GC.C.p)+extra)
  conj = [evaluate(a^t, r) for t = T] #stupid: use the previous powers!

  ps = fmpq[isinteger(GC, 2*m*B^m, power_sum(conj, i))[2] for i=1:m]
  return number_field(Hecke.power_sums_to_polynomial(ps), check = false, cached = false)[1]
end

include("Group.jl")
include("POSet.jl")
include("SeriesEval.jl")

end

using .GaloisGrp
export galois_group, transitive_group_identification, slpoly_ring, elementary_symmetric,
       power_sum, to_elementary_symmetric
