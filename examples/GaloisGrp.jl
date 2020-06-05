module GaloisGrp

using Oscar
import Base: ^, +, -, *
import Hecke
import AbstractAlgebra
using StraightLinePrograms

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
function set_stabilizer(G::Oscar.GAPGroup, seed::Vector{Int}, act)
    return stabilizer(G, GAP.julia_to_gap(seed), GAP.Globals.OnSets)
end

# TODO: add lots of more orbit related stuff

# given a perm group G and a block B, compute a homomorphism into Sym(B^G)
function action_on_blocks(G::PermGroup, B::Vector{Int})
  orb = GAP.Globals.Orbit(G.X, GAP.julia_to_gap(seed), GAP.Globals.OnSets)
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

function galois_group(f::fmpz_poly, p::Int = 31)
  GC = GaloisCtx(f, p)
  G = symmetric_group(degree(f))
  Zxx, x = PolynomialRing(ZZ, :x=>1:degree(f))
  m = prod(gen(Zxx, i)^i for i = 1:degree(f)-1)

  Zx, x = PolynomialRing(ZZ)

  nG = length(G)

  while true
    @show S = maximal_subgroup_reps(G)
    for s = S
      @show s
      @show I = sum(orbit(s, m))
      @assert all(I^x == I for x = s)
      @assert !any(I^x == I for x = right_transversal(G, s) if !isone(x))
      B = bound(GC, I)
      @show value(B)
      c = roots(GC, clog(2*value(B), p)+5)
      local fd
      cnt = 0
      while true
        cs = Set{typeof(c[1])}()
        fd = []
        for t = right_transversal(G, s)
          e = evaluate(I^t, c)
          if e in cs
            @show "darn, need tschirni", e, cs
            local ts
            while true
              @show ts = rand(Zx, 2:degree(f), -2:2)
              if degree(ts) > 1
                break
              end
            end
            cnt += 1
            if cnt > 5
              error("bad")
            end
            B = bound(GC, I, ts)
            @show value(B)
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
              println(t, " => ", l)
              push!(fd, t)
            end
          end
        end
        @show length(cs), index(G, s)
        if length(cs) == index(G, s)
          break
        end
      end
      @show "descent via", fd, s
      if length(fd)>0
        G = intersection([s^inv(x) for x = fd]...)[1]
        break
      end
    end
    if length(G) == nG
      return G
    else
      nG = length(G)
    end
  end
  return G, GC
end

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

  @show "using $p_best with degree $d_max"
  @show "cycle types ", ct

  GC = GaloisCtx(Hecke.Globals.Zx(K.pol), p_best)

  S = subfields(K)
  bs = Array{Array{Int, 1}, 1}[]
  c = roots(GC, 10)
  for (s, ms) = S
    if degree(s) == degree(K) || degree(s) == 1
      continue
    end

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
    push!(bs, collect(values(b)))
  end

  @show "Group needs to have block systems: ", bs

  d = map(frobenius, c)
  si = [findfirst(y->y==x, c) for x = d]

  @show "Group needs ", si
  return bs, si
end

end
