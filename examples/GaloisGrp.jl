Hecke.add_verbose_scope(:GaloisGroup)
Hecke.add_verbose_scope(:GaloisInvariant)
Hecke.add_assert_scope(:GaloisInvariant)

module GaloisGrp

using Oscar
import Base: ^, +, -, *
import Oscar: Hecke, AbstractAlgebra, GAP
using Oscar.StraightLinePrograms

export galois_group, isprimitive, istransitive, transitive_group_identification


function __init__()
  GAP.Packages.install("ferret")
  GAP.Packages.load("ferret")
end

struct BoundRing  <: AbstractAlgebra.Ring 
  mul
  add
  pow
end

struct BoundRingElem <: AbstractAlgebra.RingElem
  val::fmpz
  p::BoundRing
end

+(a::BoundRingElem, b::BoundRingElem) = BoundRingElem(a.p.add(a.val, b.val), a.p)
-(a::BoundRingElem, b::BoundRingElem) = BoundRingElem(a.p.add(a.val, b.val), a.p)
*(a::BoundRingElem, b::BoundRingElem) = BoundRingElem(a.p.mul(a.val, b.val), a.p)
^(a::BoundRingElem, b::Int) = BoundRingElem(a.p.pow(a.val, b), a.p)
-(a::BoundRingElem) = a

Oscar.parent(a::BoundRingElem) = a.p
value(a::BoundRingElem) = a.val

(R::BoundRing)(a::fmpz) = BoundRingElem(abs(a), R)
(R::BoundRing)(a::Integer) = BoundRingElem(abs(fmpz(a)), R)
(R::BoundRing)(a::BoundRingElem) = a
Oscar.one(R::BoundRing) = R(1)
Oscar.zero(R::BoundRing) = R(1)

function max_ring()
  return BoundRing( (x,y) -> x+y, (x,y) -> max(x, y), (x,y) -> y*x)
end

function add_ring()
  return BoundRing( (x,y) -> x*y, (x,y) -> x+y, (x,y) -> x^y)
end

Oscar.mul!(a::BoundRingElem, b::BoundRingElem, c::BoundRingElem) = b*c
Oscar.addeq!(a::BoundRingElem, b::BoundRingElem) = a+b

Oscar.parent_type(::BoundRingElem) = BoundRing
Oscar.elem_type(::BoundRing) = BoundRingElem

#my 1st invariant!!!
function sqrt_disc(a::Array{<:Any, 1})
  return prod([a[i] - a[j] for i = 1:length(a)-1 for j = i+1:length(a)])
end

function slpoly_ring(R::AbstractAlgebra.Ring, n::Int)
  return SLPolyRing(R, [ Symbol("x_$i") for i=1:n])
end

Base.one(R::StraightLinePrograms.SLPolyRing) = R(one(base_ring(R)))
(R::StraightLinePrograms.SLPolyRing)(a::StraightLinePrograms.SLPoly) = a
Oscar.ngens(S::SLPolyRing) = length(gens(S))

function root_bound(f::fmpz_poly)
  a = coeff(f, degree(f))
  return max(fmpz(1), maximum([ceil(fmpz, abs(coeff(f, i)//a)) for i=0:degree(f)]))
end

#TODO: check where this should be done properly
function (f::RelSeriesElem)(a::RingElem)
  y = a
  v = valuation(f)
  p = precision(f)
  z = zero(y)
  for i = Int(p):-1:Int(v)
      z *= y
      c = coeff(f, i)
      if !iszero(c)
          z += c
      end
  end
  z *= y^Int(v)
  return z
end
#roots are sums of m distinct roots of f
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

mutable struct GaloisCtx
  f::PolyElem
  C
  B::BoundRingElem
  function GaloisCtx(f::fmpz_poly, p::Int)
    r = new()
    r.f = f
    r.C = Hecke.qAdicRootCtx(f, p, splitting_field = true)
    r.B = add_ring()(root_bound(f))
    return r
  end
end

function Base.show(io::IO, GC::GaloisCtx)
  println(io, "Galois Context for $(GC.f) and prime $(GC.C.p)\n")
end

function Hecke.roots(G::GaloisCtx, pr::Int)
  a = Hecke.roots(G.C, pr)
  return Hecke.expand(a, all = true, flat = false, degs = Hecke.degrees(G.C.H))
end

function bound(G::GaloisCtx, f)
  return Oscar.evaluate(f, [G.B for i=1:degree(G.f)])
end

function bound(G::GaloisCtx, f, ts::fmpz_poly)
  B = ts(G.B)
  return Oscar.evaluate(f, [B for i=1:degree(G.f)])
end

function ^(f::MPolyElem, s::Oscar.GAPGroupElem{PermGroup})
  G = parent(s)
  @assert ngens(parent(f)) == degree(G)

  g = Generic.MPolyBuildCtx(parent(f))
  for (c, e) = Base.Iterators.zip(Generic.MPolyCoeffs(f), Generic.MPolyExponentVectors(f))
    s_e = zeros(Int, degree(G))
    for i=1:degree(G)
      s_e[s(i)] = e[i]
    end
    push_term!(g, c, s_e)
  end
  return finish(g)
end

function orbit(G::Oscar.PermGroup, f::MPolyElem)
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

function maximal_subgroup_reps(G::PermGroup)
  return Oscar._as_subgroups(G, GAP.Globals.MaximalSubgroupClassReps(G.X))
end

function transitive_group_identification(G::PermGroup)
  if degree(G) > 20
    return -1
  end
  return GAP.Globals.TransitiveIdentification(G.X)
end

#TODO:
# with Max, sit down and discus right/left action, transversals and conjugation
# invariants via BlockSys
# BlockSys from subfields

#TODO:
#- ansehen der ZykelTypen um Sn/An zu erkennen
#- Kranz-Produkte fuer die BlockSysteme (und deren Schnitt)
#- BlockSysteme fuer die Gruppen
#- Bessere Abstraktion um mehr Grundkoerper/ Ringe zu erlauben
#- Bessere Teilkpoerper: ich brauche "nur" maximale
#- sanity-checks
#- "datenbank" fuer Beispiele

# TODO: add a GSet Julia type which does something similar Magma's,
# or also to GAP's ExternalSet (but don't us ExternalSet, to avoid the overhead)


function all_blocks(G::PermGroup)
  # TODO: this returns an array of arrays;
  # TODO: AllBlocks assumes that we act on MovedPoints(G), which
  # may NOT be what we want...
  return GAP.gap_to_julia(Vector{Vector{Int}}, GAP.Globals.AllBlocks(G.X))
end

# TODO: update stabilizer to use GSet
# TODO: allow specifying actions other than the default
function stabilizer(G::Oscar.GAPGroup, seed, act)
    return Oscar._as_subgroup(G, GAP.Globals.Stabilizer(G.X, GAP.julia_to_gap(seed), act))
end

# TODO: add type BlockSystem

# TODO: perhaps get rid of set_stabilizer again, once we have proper Gsets
function set_stabilizer(G::Oscar.GAPGroup, seed::Vector{Int})
    return stabilizer(G, GAP.julia_to_gap(seed), GAP.Globals.OnSets)
end

# TODO: add lots of more orbit related stuff

function block_system(G::PermGroup, B::Vector{Int})
  orb = GAP.Globals.Orbit(G.X, GAP.julia_to_gap(B), GAP.Globals.OnSets)
  GAP.gap_to_julia(Vector{Vector{Int}}, orb)
end

# given a perm group G and a block B, compute a homomorphism into Sym(B^G)
function action_on_blocks(G::PermGroup, B::Vector{Int})
  orb = GAP.Globals.Orbit(G.X, GAP.julia_to_gap(B), GAP.Globals.OnSets)
  act = GAP.Globals.ActionHomomorphism(G.X, orb, GAP.Globals.OnSets)
  H = GAP.Globals.Image(act)
  T = Oscar._get_type(H)
  H = T(H)
  return Oscar._hom_from_gap_map(G, H, act)
end

istransitive(G::PermGroup) = GAP.Globals.IsTransitive(G.X, GAP.julia_to_gap(1:degree(G)))
isprimitive(G::PermGroup) = GAP.Globals.IsPrimitive(G.X, GAP.julia_to_gap(1:degree(G)))

Base.sign(G::PermGroup) = GAP.Globals.SignPermGroup(G.X)

Base.isodd(G::PermGroup) = sign(G) == -1
Base.iseven(n::PermGroup) = !isodd(n)

function elementary_symmetric(g::Array{<:Any, 1}, i::Int)
  return sum(prod(g[i] for i = s) for s = Hecke.subsets(Set(1:length(g)), i))
end

function power_sum(g::Array{<:Any, 1}, i::Int)
  return sum(a^i for a = g)
end

function Oscar.discriminant(g::Array{<:Any, 1})
  return prod(a-b for a = g for b = g if a!=b)
end

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

function ^(f::StraightLinePrograms.SLPoly, p::Oscar.GAPGroupElem{PermGroup})
  g = gens(parent(f))
  h = typeof(f)[]
  for i=1:ngens(parent(f))
    push!(h, g[p(i)])
  end
  return evaluate(f, h)
end

function isprobably_invariant(g, p)
  R = parent(g)
  k = GF(next_prime(2^20))
  n = ngens(R)
  lp = [rand(k) for i=1:n]
  return evaluate(g, lp) == evaluate(g^p, lp)
end

function short_right_transversal(G::PermGroup, H::PermGroup, s)
  C = GAP.Globals.ConjugacyClasses(H.X)
  cs = GAP.Globals.CycleStructurePerm(s.X)
  can = []
  for i=1:length(C)
    c = C[i]
    r = GAP.Globals.Representative(c)
    if cs == GAP.Globals.CycleStructurePerm(r)
      push!(can, r)
    end
  end

  R = []
  for c = can
    d = GAP.Globals.RepresentativeAction(G.X, c, s.X)
    if d != GAP.Globals.fail
      push!(R, Oscar.group_element(G, d))
      @assert Oscar.group_element(G, c)^R[end] == s
    end
  end

  S = []
  C = centralizer(G, s)[1]
  for r = R
    CH = centralizer(H^r, s)[1]
    for t = right_transversal(C, CH)
      push!(S, r*t)
    end
  end

  return S
end

function invariant(G::PermGroup, H::PermGroup)
  @vprint :GaloisInvariant 2 "Searching $G-relative $H-invariant\n"

  S = slpoly_ring(ZZ, degree(G))
  g = gens(S)

  if isprimitive(G) && isprimitive(H)
    if isodd(G) && iseven(H)
      @vprint :GaloisInvariant 3 "using sqrt_disc\n"
      return sqrt_disc(g)
    end
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


  R, _ = PolynomialRing(ZZ, degree(G))
  m = prod(gen(R, i)^i for i=1:degree(G)-1)
  I = sum(orbit(H, m))
  return evaluate(I, g)
end

#TODO
# - split: find primes and abort on Sn/An
# - sqrt discrimiant to check even/odd
# - subfield lattice directly
# - subfield to block outsource (also for lattice)
# - more invars
# - re-arrange to make cheap this first, expensive later
# - short coests - when available
# - proof, live and later
# - "isInt", ie. the test for being in Z: outsource...
# - tschirnhausen transformation: slowly increasing degree and coeffs.
# - for larger degrees: improve starting group by doing Galois groups of subfields
# - make the descent a seperate function
# - for reducibles...: write for irr. poly and for red. poly
# - more base rings
# - applications: subfields of splitting field, towers, solvability by radicals
function galois_group(K::AnticNumberField, extra::Int = 5)
  d_min = 2
  d_max = typemax(Int)
  p_best = 1
  cnt = 5
  ct = Set{Array{Int, 1}}()
  for p = Hecke.PrimesSet(2^20, -1)
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
    if cnt < 1
      break
    end
  end

  p = p_best

  @vprint :GaloisGroup 1 "using $p_best with degree $d_max\n"
  @vprint :GaloisGroup 2 "and cycle types $ct\n"

  GC = GaloisCtx(Hecke.Globals.Zx(K.pol), p_best)
  c = roots(GC, 5)

  @vprint :GaloisGroup 1 "computing subfields ...\n"
  @vtime :GaloisGroup 2 S = subfields(K)

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

  @vprint :GaloisGroup 2 "group will have (all) block systems: $([x[1] for x = bs])\n"

  #selecting maximals only...
  if length(bs) > 1
    L = POSet(bs, (x,y) -> issubset(x[1], y[1]) || issubset(y[1], x[1]), 
                  (x,y) -> issubset(x[1], y[1]) - issubset(y[1], x[1]))
    bs = minimal_elements(L)              
  end
  @vprint :GaloisGroup 1 "group will have (maximal) block systems: $([x[1] for x = bs])\n"

  d = map(frobenius, c)
  si = [findfirst(y->y==x, c) for x = d]

  @vprint :GaloisGroup 1 "found Frobenius element: $si\n"

  G = symmetric_group(degree(K))
  S = G
  for b = bs
    W = wreath_product(symmetric_group(length(b[1])), symmetric_group(length(b)))
    s = S(vcat(b...))
    G = intersect(G, W^s)[1]
  end

  if length(bs) == 0 #primitive case
    @vprint :GaloisGroup 1 "group is primitive (no subfields), trying operations on pairs\n"
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
    @vprint :GaloisGroup 2 "partitions: $O\n"
    G = Oscar._as_subgroup(G, GAP.Globals.Solve(GAP.julia_to_gap(vcat(GAP.Globals.ConInGroup(G.X), [GAP.Globals.ConStabilize(GAP.julia_to_gap(sort(o), Val(true)), GAP.Globals.OnSetsSets) for o in O]))))[1]
    #@show G = intersect([stabilizer(G, GAP.julia_to_gap(sort(o), Val(true)), GAP.Globals.OnSetsSets)[1] for o=O]...)[1]
  end

  @vprint :GaloisGroup 2 "Have starting group with id $(transitive_group_identification(G))\n"
  si = G(si)

  if issquare(discriminant(K))
    isev = true
  else
    isev = false
  end

  nG = length(G)
  while true
    S = maximal_subgroup_reps(G)
    @vprint :GaloisGroup 2 "found $(length(S)) many maximal subgroups\n"

    for s = S
      istransitive(s) || continue
      isev != iseven(s) && continue 

      @vprint :GaloisGroup 2 "testing descent $(transitive_group_identification(G)) -> $(transitive_group_identification(s)) of index $(index(G, s))\n"

      I = invariant(G, s)

      @hassert :GaloisInvariant 1 all(x->isprobably_invariant(I, x), gens(s))
      @hassert :GaloisInvariant 1 any(x->!isprobably_invariant(I, x), gens(G))

      B = bound(GC, I)
      @vprint :GaloisGroup 2 "root bound: $(value(B))\n"

      c = roots(GC, clog(2*value(B), p)+extra)
      local fd
      cnt = 0
      while true
        cs = Set{typeof(c[1])}()
        fd = []

        local lt
        if index(G, s) < 100
          lt = right_transversal(G, s)
        elseif isnormal(G, s)
          lt = [one(G)] # I don't know how to get the identity
        else
          lt = short_right_transversal(G, s, si)
          #if index(G, s) < 1000 && get_assert_level(:GaloisInvariant) > 2
          #  tt = [x for x = right_transversal(G, s) if si in s^x]
          #  if Set(right_coset(s, x) for x = tt) != Set(right_coset(s, x) for x = lt)
          #    return G, s, si
          #  end
          #end
          @hassert :GaloisInvariant 2 all(x->si in s^x, lt)
        end
        #x^7-2 triggers wrong descent, occaisonly, if lt is only 10 elems long
        while length(lt) < 20 && length(lt) < index(G, s)
          r = rand(G)
          if all(x->right_coset(s, r) != right_coset(s, x), lt)
            push!(lt, r)
          end
        end

        for t = lt
          e = evaluate(I^t, c)
          if e in cs
            @vprint :GaloisGroup 2 " evaluation found duplicate, transforming...\n"
            local ts
            while true
              ts = rand(Hecke.Globals.Zx, 2:degree(K), -4:4)
              if degree(ts) > 1
                break
              end
            end
            cnt += 1
            if cnt > 5
#              @show I
#              return GC, I, G, s
              error("bad")
            end
            B = bound(GC, I, ts)
            @vprint :GaloisGroup 2 "using transformation by $ts, and new bound of $(value(B))\n"
            c = roots(GC, clog(2*value(B), p)+extra)
            c = map(ts, c)
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
          break
        end
      end
      if length(fd) > 0
        @vprint :GaloisGroup 2 "descending via $fd\n"
      else
        @vprint :GaloisGroup 2 "no descent\n"
      end
      if length(fd)>0
        G = intersect([s^x for x = fd]...)[1]
        @assert length(fd) > 1 || length(G) == length(s)
        @assert length(G) < nG
        @vprint :GaloisGroup 2 "descending to $(transitive_group_identification(G))\n"
        if length(G) == degree(K)
          return G, GC
        end
        break
      end
    end
    if length(G) == nG
      return G, GC
    else
      nG = length(G)
    end
  end
end

#TODO: do this properly, with indexed types and such
mutable struct POSet
  elem::Array{Any, 1}
  cmp::Function
  can_cmp::Function

  function POSet(elem::Array{<:Any, 1}, can_cmp::Function, cmp::Function)
    r = new()
    r.elem = elem
    r.cmp = cmp
    r.can_cmp = can_cmp
    return r
  end
end

struct POSetElem
  e::Int
  p::POSet
  function POSetElem(L::POSet, i::Int)
    return new(i, p)
  end
  function POSetElem(L::POSet, e::Any)
    return new(findfirst(x->L>can_cmp(e, L.elem[i]) && L.cmp(e, L.elem[i]) == 0, 1:length(L.elem)), L)
  end
end

function maximal_elements(L::POSet)
  d = Dict{typeof(L.elem[1]), Array{Int, 1}}()
  for i = 1:length(L.elem)
    e = L.elem[i]
    d[e] = Int[]
    for j=1:length(L.elem)
      if i == j
        continue
      end
      if L.can_cmp(e, L.elem[j]) && L.cmp(e, L.elem[j]) > 0
        push!(d[e], j)
      end
    end
  end
  return [x for (x,v) = d if length(v) == 0]
end

function minimal_elements(L::POSet)
  d = Dict{typeof(L.elem[1]), Array{Int, 1}}()
  for i = 1:length(L.elem)
    e = L.elem[i]
    d[e] = Int[]
    for j=1:length(L.elem)
      if i == j
        continue
      end
      if L.can_cmp(e, L.elem[j]) && (L.cmp(e, L.elem[j]) < 0)
        push!(d[e], j)
      end
    end
  end
  return [x for (x,v) = d if length(v) == 0]
end

end

using .GaloisGrp
export galois_group, isprimitive, istransitive, transitive_group_identification,
       GaloisGrp
