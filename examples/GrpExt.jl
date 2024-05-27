module GrpExt_module
using Oscar

import Base: *, ==, one, rand, show, iterate

struct GrpExt{S, T} <: AbstractAlgebra.Group
  c :: Oscar.GrpCoh.CoChain{2, S, T} 
  function GrpExt(c::Oscar.GrpCoh.CoChain{2, A, B}) where A where B
    return new{A, B}(c)
  end
end

function show(io::IO, G::GrpExt)
  print(io, "Extension via $(gmodule(G))")
end

struct GrpExtElem{S, T} <: AbstractAlgebra.GroupElem
  P::GrpExt{S, T}
  g::S
  m::T
  function GrpExtElem(P::GrpExt{A, B}, g::A, m::B) where A where B 
    return new{A, B}(P, g, m)
  end
end

function show(io::IO, g::GrpExtElem)
  print(io, "($(g.g), $(g.m))")
end

Oscar.gmodule(P::GrpExt) = P.c.C
Oscar.gmodule(P::GrpExtElem) = gmodule(parent(P))
_module(P::GrpExt) = gmodule(P).M
_group(P::GrpExt) = gmodule(P).G

Oscar.parent(g::GrpExtElem) = g.P

function one(G::GrpExt)
  return GrpExtElem(G, one(_group(G)), zero(_module(G)))
end

function Oscar.gens(G::GrpExt)
  gr = _group(G)
  mo = _module(G)
  g = gens(gr)
  m = gens(mo)
  return vcat([GrpExtElem(G, x, zero(mo)) for x = g], [GrpExtElem(G, one(gr), x) for x = m])
end

function *(g::GrpExtElem, h::GrpExtElem)
  @assert parent(g) === parent(h)
  return GrpExtElem(g.P, g.g*h.g, action(gmodule(g), h.g)(g.m) + h.m + g.P.c(g.g, h.g))
end

function Oscar.inv(g::GrpExtElem)
  gi = inv(g.g)
  return GrpExtElem(g.P, gi, -action(gmodule(g), gi)(g.m)-g.P.c(g.g, gi))
end

function rand(G::GrpExt)
  return GrpExtElem(G, rand(_group(G)), rand(_module(G)))
end

function ==(g::GrpExtElem, h::GrpExtElem)
  @assert g.P === h.P
  return g.g == h.g && g.m == h.m
end

function iterate(G::GrpExt)
  I = Base.Iterators.ProductIterator((_group(G), _module(G)))
  X = iterate(I)

  return GrpExtElem(G, X[1][1], X[1][2]), (I, X[2], G)
end

function iterate(G::GrpExt, X::Tuple)
  n = iterate(X[1], X[2])
  if n === nothing
    return n
  end
  return GrpExtElem(X[3], n[1][1], n[1][2]), (X[1], n[2], X[3])
end

Oscar.order(G::GrpExt) = order(_group(G))*order(_module(G))
Base.length(G::GrpExt) = Int(order(G))

function Oscar.order(g::GrpExtElem)
  o = order(g.g)
  g = g^Int(o)
  return o*order(g.m)
end

function(P::GrpExt{S, T})(g::S, a::T) where S where T
  return GrpExtElem(P, g, a)
end

function is_stem_extension(P::GrpExt)
  G = _group(P)
  #TODO: test central first (operation trivial)
  #  algo is explained in Cohomology.jl
  E, mE = derived_subgroup(G)
  i = isomorphism(FPGroup, E) # E -> Fp
  R = relators(codomain(i))
  g = gens(codomain(i)) # in Fp
  g = map(pseudo_inv(i), g) # in E
  g = map(mE, g) # in G
  g = [P(x, zero(_module(P))) for x = g] # in P
  s, ms = sub(_module(P), [map_word(r, g).m for r = R])
  return is_surjective(ms)
end

Oscar.copy(g::GrpExtElem) = GrpExtElem(g.P, g.g, g.m)

end #module
