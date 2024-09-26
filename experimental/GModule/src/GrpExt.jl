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
struct GrpExt{S, T} <: AbstractAlgebra.Group
  c :: Oscar.GrpCoh.CoChain{2, S, T} 
  function GrpExt(c::Oscar.GrpCoh.CoChain{2, A, B}) where A where B
    return new{A, B}(c)
  end
end

function show(io::IO, G::GrpExt)
  print(io, "Extension via $(gmodule(G))")
end

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
  mGF = Oscar.isomorphism(FPGroup, G, on_gens=true)
  F = codomain(mGF)
  M = C.M
  ac = action(C)
  mfM = inv(Oscar.isomorphism(FPGroup, M))
  fM = domain(mfM)

  N = free_group(ngens(G) + ngens(fM))

  s = map(x->shiftgens(x, N, ngens(G)), relators(fM))
  for R = relators(F)
    t = map_word(R, gens(E)[1:ngens(G)])
    push!(s, shiftgens(R, N, 0)*(shiftgens(preimage(mfM, t.m), N, ngens(fM))))
  end
  for i=1:ngens(G)
    for j=1:ngens(fM)
      #g[i]*t = m[j]*g[i] = g[i] m[j]^g[i] = m[j] g[i] (cancellation in conj)
      t = preimage(mfM, ac[i](gen(M, j)))
      push!(s, gen(N, ngens(G)+j)*gen(N, i)*inv(shiftgens(t, N, ngens(fM))) * inv(gen(N, i)))
    end
  end
  Q, mQ = quo(N, s)
  @assert ngens(Q) == ngens(N)
  function EtoQ(x::GrpExtElem)
    w = mGF(x.g)
    ww = map_word(w, gens(E)[1:ngens(G)]) #this performs a collection
                                          #and will transform x.g into
                                          #canonical form
    return mQ(shiftgens(w, N, 0)*shiftgens(preimage(mfM, x.m-ww.m), N, ngens(fM)))

  end
  return MapFromFunc(E, Q, EtoQ, y->map_word(y, gens(E)))
  #the projection will be   hom(Q, G, vcat(gens(G), ones(G, ???)
  #the injection should be  hom(M, Q, gens(Q)[ngens(G)+1:end])
  #both can/ should be handled by shiftgens (or sylable/ create)
end

#=
function isomorphism(::Type{PcGroup}, E::GrpExt)
  c = E.c
  C = c.C
  G = Group(C)
  @assert isa(G, PcGroup)
  M = Module(C)
  ac = action(C)
  iac = inv_action(C)
  fM, mfM = pc_group_with_isomorphism(M)

  N = free_group(ngens(G) + ngens(fM))
  Gp = GAP.Globals.Pcgs(GapObj(G))
  @assert length(Gp) == ngens(G)
#  @assert all(x->Gp[x] == GapObj(gen(G, x)), 1:ngens(G))
  Go = GAP.Globals.RelativeOrders(Gp)

  Mp = GAP.Globals.Pcgs(GapObj(fM))
  @assert length(Mp) == ngens(fM) == ngens(M)
#  @assert all(x->Mp[x] == GapObj(gen(fM, x)), 1:ngens(M))
  #problem/ TODO: Z/100Z has a useful GAP-pc-group has 4 gens (of
  #order 2, 2, 5, 5
  #so need to switch GAP to the other Pc-Groups and/or drop this 
  #assert
  Mo = GAP.Globals.RelativeOrders(Mp)

  CN = GAP.Globals.SingleCollector(GapObj(N), GAP.Globals.Concatenation(Go, Mo))
  FN = GAP.Globals.FamilyObj(GapObj(N[1]))

  for i=1:ngens(fM)
    lp = deepcopy(GAPWrap.ExtRepOfObj(Mp[i]^Mo[i]))
    for k=1:2:length(lp)
      lp[k] += ngens(G)
    end
    m = GAP.Globals.ObjByExtRep(FN, lp)
    GAP.Globals.SetPower(CN, i+ngens(G), m)
    for j=i+1:ngens(fM)
      p = Mp[j]^Mp[i]
      @assert p == Mp[j]
      lp = deepcopy(GAPWrap.ExtRepOfObj(p))
      for k=1:2:length(lp)
        lp[k] += ngens(G)
      end
      GAP.Globals.SetConjugate(CN, j+ngens(G), i+ngens(G), GAP.Globals.ObjByExtRep(FN, lp))
    end
  end

  fMtoN = function(x)
    lp = deepcopy(GAPWrap.ExtRepOfObj(GapObj(x)))
    for k=1:2:length(lp)
      @assert lp[k] > 0
      lp[k] += ngens(G)
    end
    return GAP.Globals.ObjByExtRep(FN, lp)
  end

  word = function(y)
    z = GAPWrap.UnderlyingElement(y)
    return map(Int, GAP.Globals.LetterRepAssocWord(z))
  end

  #for W = (w1, ... w_n) compute ((w1, 0), ..., (wn, 0))
  #and return the tail only.
  word_to_elem = function(W)
    t = zero(M)
    g = one(G)
    r = one(N)
    for w = W
      if w > 0
        t = ac[w](t) + c(g, gen(G, w))
        g = g*gen(G, w)
        r = r*gen(N, w)
      else
        t = iac[-w](t) + c(g, inv(gen(G, -w))) - c(gen(G, -w), inv(gen(G, -w)))
        g = g*inv(gen(G, -w))
        r = r*inv(gen(N, -w))
      end
    end
    return t
    return fMtoN(mfM(t))
  end

  #to lift the pc-relations:
  # F^p = w (order relation)
  #  compute (F, 0)^p = (?, t) = (?, 0)(1, t)
  #  compute (w, 0)   = (?, s) = (?, 0)(1, s)
  #  so (?, 0) = (w, 0)(1,s)^-1= (w, 0)(1,-s) if chain is normalized
  #  thus (F, 0)^p = (?, 0)(1, t) = (w, 0)(1,-s)(1, t)
  #  the ? should be identical, namely the collected version of w
  #  then (F, 0)^p = (w, t-s) might be the answer
  # F^G = w (conjugate relation): same
  #  (F, 0)^(G, 0) = (?, t) = (?, 0)(1, t)
  #  (w, 0)        = (?, s) = (?, 0)(1, s)
  #  thus (F, 0)^(G, 0) = (w, t-s)
  for i=1:ngens(G)
    p = Gp[i]^Go[i]
    pp = GAP.Globals.ObjByExtRep(FN, GAPWrap.ExtRepOfObj(p))
    m = fMtoN(mfM(word_to_elem([i for k=1:Go[i]])-word_to_elem(word(p))))
    GAP.Globals.SetPower(CN, i, pp*m)
    for j=i+1:ngens(G)
      p = Gp[j]^Gp[i]
      m = fMtoN(mfM(word_to_elem([-i, j, i])-word_to_elem(word(p))))
      pp = GAP.Globals.ObjByExtRep(FN, GAPWrap.ExtRepOfObj(p))
      GAP.Globals.SetConjugate(CN, j, i, pp*m)
    end
    for j=1:ngens(fM)
      m = fMtoN(mfM(action(C, gen(G, i), preimage(mfM, gen(fM, j)))))
      GAP.Globals.SetConjugate(CN, j+ngens(G), i, m)
    end
  end

#  l = GAP.Obj([])
#  GAP.Globals.FinitePolycyclicCollector_IsConfluent(CN, l)
#  @show l

#  z = GAP.Globals.GroupByRwsNC(CN)
#  s = GAP.Globals.GapInputPcGroup(z, GAP.Obj("Z"))
#  @show GAP.gap_to_julia(s)
  Q = PcGroup(GAP.Globals.GroupByRws(CN))
  fQ = GAP.Globals.FamilyObj(GapObj(one(Q)))
  mQ = hom(N, Q, gens(N), gens(Q); check = false)

  @assert ngens(Q) == ngens(N)
  MtoQ = hom(fM, Q, gens(fM), gens(Q)[ngens(G)+1:end]; check = false)
  QtoG = hom(Q, G, gens(Q), vcat(gens(G), [one(G) for i=1:ngens(fM)]); check = false)
  @assert domain(mfM) == M
  @assert codomain(mfM) == fM
#  @assert is_surjective(QtoG)
#  @assert is_injective(MtoQ)

  mfG = epimorphism_from_free_group(G)
  mffM = epimorphism_from_free_group(fM)

  function GMtoQ(wg, m)
    wm = GAP.gap_to_julia(GAPWrap.ExtRepOfObj(GapObj(preimage(mffM, mfM(m)))))
    for i=1:2:length(wm)
      push!(wg, wm[i]+ngens(G))
      push!(wg, wm[i+1])
    end
    return mQ(FPGroupElem(N, GAP.Globals.ObjByExtRep(FN, GAP.Obj(wg))))
  end

  return Q, mfM*MtoQ, QtoG, GMtoQ
end
=#

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
  # a commutator (comm(X, Y)) does not depend on X,m and Y.m
  # We need to check M <= P'
  # in P, M ;ooks like (1, M)
  # any word evaluated in P will have the G part as the evaluation
  # only on the G part (the G part is the projection, hence a hom)
  # any word evaluated in P to get a (1, M) will be a relator of G'
  # Writig such a relator as a word in commutators make the evaluation
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

Computes the subgroup of `G` generated by `U` and `g`.
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

end #module

using .GrpExt_module

export commutator_decomposition_map
