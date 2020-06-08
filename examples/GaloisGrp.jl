module GaloisGrp

using Oscar
import Base: ^, +, -, *
import Oscar: Hecke, AbstractAlgebra
using StraightLinePrograms

export galois_group, isprimitive, istransitive, transitive_group_identification

Hecke.add_verbose_scope(:GaloisGroup)
Hecke.add_verbose_scope(:GaloisInvariant)
Hecke.add_assert_scope(:GaloisInvariant)


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
function sqrt_disc(a::Array{<:AbstractAlgebra.RingElem, 1})
  return prod([a[i] - a[j] for i = 1:length(a)-1 for j = i+1:length(a)])
end

function slpoly_ring(R::AbstractAlgebra.Ring, n::Int)
  return SLPolyRing(R, [ Symbol("x_$i") for i=1:n])
end

function root_bound(f::fmpz_poly)
  a = coeff(f, 0)
  return max(fmpz(1), maximum([ceil(fmpz, abs(coeff(f, i)//a)) for i=1:degree(f)]))
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
  return Oscar._as_subgroups(GAP.Globals.MaximalSubgroupClassReps(G.X), G)
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
    return Oscar._as_subgroup(GAP.Globals.Stabilizer(G.X, GAP.julia_to_gap(seed), act), G)
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

Base.sign(g::PermGroupElem) = GAP.Globals.Sign(g.X)

Base.sign(G::PermGroup) = GAP.Globals.SignPermGroup(G.X)
#Base.sign(G::PermGroup) = all(x -> sign(x) == 1, G)

Base.isodd(G::PermGroup) = sign(G) == -1
Base.iseven(n::PermGroup) = !isodd(n)

function elementary_symmetric(g::Array{<:Any, 1}, i::Int)
  return sum(prod(g[i] for i = s) for s = Hecke.subsets(Set(1:length(g)), i))
end

Oscar.ngens(S::SLPolyRing) = length(gens(S))

function isprobably_invariant(g, p)
  R = parent(g)
  k = GF(next_prime(2^20))
  n = ngens(R)
  lp = [rand(k) for i=1:n]
  return evaluate(g, lp) == evaluate(g, [lp[p(i)] for i=1:n])
end

function invariant(G, H)
  @vprint :GaloisInvariant 2 "Searching $G-relative $H-invariant\n"

  S = slpoly_ring(ZZ, degree(G))
  g = gens(S)

  if isprimitive(G) && isprimitive(H)
    if isodd(G) && iseven(H)
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
    end
  end

  @vprint :GaloisInvariant 3 "no special invar found, resorting to generic\n"
  R, _ = PolynomialRing(ZZ, degree(G))
  m = prod(gen(R, i)^i for i=1:degree(G)-1)
  I = sum(orbit(H, m))
  return I
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
function galois_group(K::AnticNumberField)
  d_min = div(degree(K), 4)
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
    if d < d_max && d > d_min
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
    G = intersection(G, W^s)[1]
  end

  @vprint :GaloisGroup 2 "Have starting group with id $(transitive_group_identification(G))\n"
  nG = length(G)
  while true
    S = maximal_subgroup_reps(G)
    @vprint :GaloisGroup 2 "found $(length(S)) many maximal subgroups\n"

    for s = S
      istransitive(s) || continue

      @vprint :GaloisGroup 2 "testing descent $(transitive_group_identification(G)) -> $(transitive_group_identification(s))\n"

      I = invariant(G, s)

      @hassert :GaloisInvariant 1 all(x->isprobably_invariant(I, x), gens(s))
      @hassert :GaloisInvariant 1 any(x->!isprobably_invariant(I, x), gens(G))

      B = bound(GC, I)
      @vprint :GaloisGroup 2 "root bound: $(value(B))\n"

      c = roots(GC, clog(2*value(B), p)+5)
      local fd
      cnt = 0
      while true
        cs = Set{typeof(c[1])}()
        fd = []
        for t = left_transversal(G, s)
          e = evaluate(I, [c[inv(t)(i)] for i=1:degree(K)])
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
              error("bad")
            end
            B = bound(GC, I, ts)
            @vprint :GaloisGroup 2 "using transformation by $ts, and new bound of $(value(B))\n"
            c = roots(GC, clog(2*value(B), p)+5)
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
            end
          end
        end
        if length(cs) == index(G, s)
          break
        end
      end
      if length(fd) > 0
        @vprint :GaloisGroup 2 "descending via $fd\n"
      else
        @vprint :GaloisGroup 2 "no descent\n"
      end
      if length(fd)>0
        G = intersection([s^inv(x) for x = fd]...)[1]
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
