module GrpExt_module
using Oscar

import Base: *, ==, one, rand, show, iterate
export GrpExt, GrpExtElem
export commutator_decomposition_map
import Oscar.GAPWrap

"""
A type representing the group extension by a 2-co-cycle.
Let ``h in H^2(G, M)`` for some group `G` and module `M`.
Then `S` is the type of elements in `G` while `T` is the
type of the module elements.

The chains is then of type `CoChain{2, S, T}`, elements
of this group are essentially pairs of elements of `G` and
`M`.
"""
struct GrpExt{S, T} <: Group
  c :: Oscar.GrpCoh.CoChain{2, S, T} 
  function GrpExt(c::Oscar.GrpCoh.CoChain{2, A, B}) where A where B
    return new{A, B}(c)
  end
end

function show(io::IO, G::GrpExt)
  print(io, "Extension via $(gmodule(G))")
end

"""
The elements of the extension, given via a 2-co-chain.
"""
struct GrpExtElem{S, T} <: AbstractAlgebra.GroupElem
  P::GrpExt{S, T} #parent
  g::S            #the group bit
  m::T            #the module bit
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
Oscar.elem_type(::Type{GrpExt{S, T}}) where {S, T} = GrpExtElem{S, T}

function Oscar.extension(c::Oscar.GrpCoh.CoChain{2,<:Oscar.GAPGroupElem})
  return GrpExt(c)
end

#g in H 
# -> g in G where the gens to be used in G are different
function shiftgens(g::FPGroupElem, G, offset)
  w = GapObj(g)
  famG = GAPWrap.ElementsFamily(GAPWrap.FamilyObj(GapObj(G)))
  if GAP.Globals.IsLetterAssocWordRep(w)
    l = copy(GAP.Globals.LetterRepAssocWord(w))
    for i in 1:length(l)
      if l[i] > 0
        l[i] = l[i] + offset
      else
        l[i] = l[i] - offset
      end
    end
    ll = GAP.Globals.AssocWordByLetterRep(famG, l)
  else
    l = copy(GAPWrap.ExtRepOfObj(w))
    for i in 1:2:length(l)
      l[i] = l[i] + offset
    end
    ll = GAPWrap.ObjByExtRep(famG, l)
  end
  return FPGroupElem(G, ll)
end

function Oscar.isomorphism(::Type{FPGroup}, E::GrpExt)
  c = E.c
  C = c.C
  G = C.G
  mGF = Oscar.isomorphism(FPGroup, G, on_gens=true) #G -> F
  F = codomain(mGF)
  M = C.M
  ac = action(C)
  mfM = inv(Oscar.isomorphism(FPGroup, M))
  fM = domain(mfM)

  N = free_group(ngens(G) + ngens(fM))

  s = map(x->shiftgens(x, N, ngens(G)), relators(fM))
  for R = relators(F)
    t = map_word(R, gens(E)[1:ngens(G)])
    push!(s, shiftgens(R, N, 0)*(shiftgens(preimage(mfM, t.m), N, ngens(G))))
  end
  for i=1:ngens(G)
    for j=1:ngens(fM)
      #g[i]*t = m[j]*g[i] = g[i] m[j]^g[i] = m[j] g[i] (cancellation in conj)
      t = preimage(mfM, ac[i](gen(M, j)))
      push!(s, gen(N, ngens(G)+j)*gen(N, i)*inv(shiftgens(t, N, ngens(G))) * inv(gen(N, i)))
    end
  end
  Q, mQ = quo(N, s)
  @assert ngens(Q) == ngens(N)
  function EtoQ(x::GrpExtElem)
    w = mGF(x.g)
    ww = map_word(w, gens(E)[1:ngens(G)]) #this performs a collection
                                          #and will transform x.g into
                                          #canonical form
    return mQ(shiftgens(w, N, 0)*shiftgens(preimage(mfM, x.m-ww.m), N, ngens(G)))

  end
  return MapFromFunc(E, Q, EtoQ, y->map_word(y, gens(E)))
  #the projection will be   hom(Q, G, vcat(gens(G), ones(G, ???)
  #the injection should be  hom(M, Q, gens(Q)[ngens(G)+1:end])
  #both can/ should be handled by shiftgens (or sylable/ create)
end

# convert syllables in canonical form into exponent vector
#Thomas
function exponent_vector(sylls::Vector{Pair{Int64, ZZRingElem}}, n)
  res = zeros(ZZRingElem, n)
  for pair in sylls
    @assert res[pair.first] == 0 #just to make sure 
    res[pair.first] = pair.second
  end
  return res
end

function Oscar.isomorphism(::Type{PcGroup}, E::GrpExt)
  c = E.c
  C = c.C
  G = C.G
  @assert isa(G, PcGroup)
  M = C.M
  mMf = Oscar.isomorphism(PcGroup, M) # M -> PcGroup
  #TODO: this group is internal only. Should it be made available?
  #TODO: this should be using the new syllable/ expo vector interface and
  #      not do any group arithmetic
  #it has 
  fM = codomain(mMf)

  cM = collector(fM)
  cG = collector(G)
  nG = ngens(G)
  cE = collector(nG + ngens(fM))

  set_relative_orders!(cE, vcat(get_relative_orders(cG), get_relative_orders(cM)))

  for i=1:ngens(fM)
    set_power!(cE, i + nG, [g[1] + nG => g[2] for g = get_power(cM, i)])
  end

  #TODO: if all is well, then the tails used to compute H^2 will match exactly
  #the data required here - thus no computation will be required  
  #well for the conjugates the action on the module
  for i=1:nG
    #XXX: only since we don't have GrpEltElem^ZZRingElem (yet)
    x = E[i]^Int(get_relative_order(cG, i))
    t = vcat(Oscar.syllables(x.g), [ g[1] + nG => g[2] for g = Oscar.syllables(mMf(x.m))])
    set_power!(cE, i, t)
  end

  for i=1:nG
    for j=i+1:nG
      x = E[j]^E[i]
      t = vcat(Oscar.syllables(x.g), [ g[1] + nG => g[2] for g = Oscar.syllables(mMf(x.m))])
      set_conjugate!(cE, j, i, t)
    end
    for j=nG+1:nG+ngens(fM)
      x = E(one(G), preimage(mMf, fM[j-nG]))^E[i]
      t = vcat(Oscar.syllables(x.g), [ g[1] + nG => g[2] for g = Oscar.syllables(mMf(x.m))])
      set_conjugate!(cE, j, i, t)
    end
  end

  Q = pc_group(cE)

  function EtoQ(x::GrpExtElem)
    @assert parent(x) == E
    t = vcat(Oscar.syllables(x.g), [ g[1] + nG => g[2] for g = Oscar.syllables(mMf(x.m))])
    return Q(t)
  end

  function QtoE(x::PcGroupElem)
    #implicitly assumes the syllables to be a normalized expo vector
    #in particular the module part needs to be in the end
    @assert parent(x) == Q
    s = Oscar.syllables(x)
    a = Pair{Int, ZZRingElem}[]
    b = Pair{Int, ZZRingElem}[]
    for i = s
      if i[1] <= nG
        push!(a, i)
      else
        push!(b, i[1]-nG => i[2])
      end
    end
    return E(G(a), preimage(mMf, fM(b)))
  end
  #XXX: this is not yet good:
  # - we use elt -> syllable
  # - but create elt from exponent vector (well this is hidden)
  # why is the collector using syllables
  #
  #the collector has a type parameter for the exponents. This is "ignored"
  #here

  return MapFromFunc(E, Q, EtoQ, QtoE) #y->map_word(y, gens(E)))
end
#=
Thomas:
Fuer die Umkehrabbildung kann man statt y -> map_word(y, gens(E))
dann y in Syllables auspacken, an der richtigen Stelle in zwei Teile
aufteilen und die mit Hilfe der Kollektoren cM und cG in Elemente
verpacken, aus denen dann ein GrpExtElem gebaut wird.

=#
function one(G::GrpExt)
  return GrpExtElem(G, one(_group(G)), zero(_module(G)))
end

function Oscar.gen(G::GrpExt, i::Int)
  i == 0 && return one(G)
  i < 0 && return inv(gen(G, -i))
  i > ngens(G) && error("index out of range")
  gr = _group(G)
  mo = _module(G)
  if i <= ngens(gr)
    return GrpExtElem(G, gen(gr, i), zero(mo))
  end
  return GrpExtElem(G, one(gr), gen(mo, i-ngens(gr)))
end

function Oscar.ngens(G::GrpExt)
  return ngens(_group(G)) + ngens(_module(G))
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
  #XXX: will this work for FPGroupElems? They may not be normalized
  #     if the element is obtained through E(g, m)
  #     if it comes from arithmetic, I think it will be fine
  return g.g == h.g && g.m == h.m
end

function Base.hash(g::GrpExtElem, h::UInt)
  return hash(g.g, hash(g.m, hash(g.P, h)))
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

function Oscar.is_one(m::Map{FinGenAbGroup, FinGenAbGroup})
  domain(m) === codomain(m) || return false
  return all(x->x == m(x), domain(m))
end

function Oscar.is_central(P::GrpExt)
  return all(is_one, action(gmodule(P)))
end

function Oscar.is_stem_extension(P::GrpExt)
  is_central(P) || return false
  G = _group(P)
  #  part of algo is explained in Cohomology.jl
  #
  # a commutator (comm(X, Y)) does not depend on X.m and Y.m
  # We need to check M <= P'
  # in P, M looks like (1, M)
  # any word evaluated in P will have the G part as the evaluation
  # only on the G part (the G part is the projection, hence a hom)
  # any word evaluated in P to get a (1, M) will be a relator of G'
  # Writing such a relator as a word in commutators make the evaluation
  # depend only on the co-chain, not on other choices.
  E, mE = derived_subgroup(G)
  h = commutator_decomposition_map(mE)
  i = isomorphism(FPGroup, E) # E -> Fp
  R = relators(codomain(i))

  g1 = gens(codomain(i)) # in Fp
  g2 = map(pseudo_inv(i), g1) # in E
  g = [map_word(h(x), [P(x, zero(_module((P)))) for x = gens(G)]) for x = g2]
  ss = [map_word(r, g) for r = R]

  @assert all(x->isone(x.g), ss)
  s, ms = sub(_module(P), [x.m for x = ss])
  return is_surjective(ms)
end

@doc raw"""
    sub(G::Oscar.GAPGroup, U::Oscar.GAPGroup, g::Oscar.GAPGroupElem) -> GAPGroup, Map

Compute the subgroup of `G` generated by `U` and `g`.
"""
function Oscar.sub(G::Oscar.GAPGroup, U::Oscar.GAPGroup, g::Oscar.GAPGroupElem)
  H = GAP.Globals.ClosureGroup(U.X, g.X)
  return Oscar._as_subgroup(G, H)
end

function commutator_decomposition_map(G::Oscar.GAPGroup)
  return commutator_decomposition_map(derived_subgroup(G)[2])
end

@doc raw"""
    commutator_decomposition_map(md::Map{<:Oscar.GAPGroup, <:Oscar.GAPGroup})
  
Given
```G' \to G```
compute a map 
```\phi: G' \to F```
into the free group on the number of generators of ``G``, such that the
image is composed of commutators and conjugates only.
"""
function commutator_decomposition_map(md::Map{<:Oscar.GAPGroup, <:Oscar.GAPGroup})
  d = domain(md)
  G = codomain(md)
  g = gens(G)
  F = free_group(length(g))
  s = elem_type(G)[]
  t = elem_type(F)[]
  u, mu = sub(G, [one(G)])
  for i=1:ngens(G), j=i+1:ngens(G)
    gij = comm(G[i], G[j])
    if is_one(gij) || gij in s || gij in u
      continue
    end
    push!(s, gij)
    push!(t, comm(F[i], F[j]))
    u, mu = sub(G, u, gij)
    if order(u) == order(d)
      break
    end
  end
  
  j = 1
  while order(u) < order(d)
    n = length(s)
    for i=j:n, h = 1:ngens(G)
      hg = s[i]^g[h]
      if hg in s || hg in u
        continue
      end
      push!(s, hg)
      push!(t, t[i]^F[h])
      u, mu = sub(G, u, hg)
      if order(u) == order(d)
        break
      end
    end
    j = n
  end
  h = hom(free_group(length(s)), u, s)

  return MapFromFunc(d, F, x->map_word(preimage(h, x), t))
end  
    
Oscar.copy(g::GrpExtElem) = GrpExtElem(g.P, g.g, g.m)

function Oscar.permutation_group(G::GrpExt)
  if isa(group(gmodule(G)), PcGroup)
    return permutation_group(codomain(isomorphism(PcGroup, G)))
  else
    return permutation_group(codomain(isomorphism(FPGroup, G)))
  end
end

end #module

using .GrpExt_module

export commutator_decomposition_map
