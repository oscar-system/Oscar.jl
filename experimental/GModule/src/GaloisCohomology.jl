module GaloisCohomology_Mod
using Oscar
import Oscar: gmodule
import Oscar: GrpCoh
import Oscar.GrpCoh: CoChain, MultGrpElem, MultGrp, GModule, is_consistent, 
                     Group
import Base: parent
import Oscar: direct_sum
import Oscar: pretty, Lowercase, Indent, Dedent, terse

export is_coboundary, idele_class_gmodule, relative_brauer_group
export local_invariants, global_fundamental_class, shrink
export local_index, units_mod_ideal, brauer_group

import AbstractAlgebra: @show_name, @show_special, extra_name


_can_cache_aut(::PadicField) = nothing
function _can_cache_aut(k) 
  a = get_attribute(k, :AutGroup)
  if a === nothing
    a = Dict{Tuple, Map}()
    set_attribute!(k, :AutGroup => a)
  end
  return a
end
#not type restricted: used for number fields as well as
#LocalFields
#could be extended to allow more group types
function Oscar.automorphism_group(::Type{PermGroup}, k)
  a = _can_cache_aut(k)
  if a !== nothing && haskey(a, (PermGroup, ))
    mG = a[(PermGroup, )]
    return domain(mG), mG
  end
  G, mG = automorphism_group(k)
  mH = isomorphism(PermGroup, G)
  mmH = inv(mH)*mG
  if a !== nothing 
    a[(PermGroup, )] = mmH
  end
  return codomain(mH), mmH
end

function Oscar.automorphism_group(::Type{PermGroup}, K, k)
  a = _can_cache_aut(k)
  if a !== nothing && haskey(a, (PermGroup, k))
    mG = a[(PermGroup, k)]
    return domain(mG), mG
  end

  G, mG = automorphism_group(K, k)
  mH = isomorphism(PermGroup, G)
  mmH = inv(mH)*mG
  if a !== nothing 
    a[(PermGroup, k)] = mmH
  end
  return codomain(mH), mmH
end

function Oscar.absolute_automorphism_group(::Type{PermGroup}, k)
  G, mG = absolute_automorphism_group(k)
  mH = isomorphism(PermGroup, G)
  return codomain(mH), inv(mH)*mG
end

"""
    units_mod_ideal(I::AbsSimpleNumFieldOrderIdeal; n_quo::Int = 0) -> FinGenAbGroup, Map{Grp, AbsSimpleNumFieldOrder}

Compute the unit group of the order modulo `I`. If `n_quo` is non-zero, the quotient
modulo `n_quo` is computed.
"""
function units_mod_ideal(I::AbsSimpleNumFieldOrderIdeal; n_quo::Int = 0)
  #TODO: add places for sign condition (RatResidueRing in Magma)
  #      use the n_quo already in the creation.
  R, mR = quo(order(I), I)
  U, mU = unit_group(R)
  if n_quo != 0
    u, mu = quo(U, n_quo)
    U = u
    mU = pseudo_inv(mu)*mU
  end
  return U, mU*pseudo_inv(mR)
end

"""
    gmodule(H::PermGroup, mR::MapRayClassGrp, mG = automorphism_group(PermGroup, Hecke.nf(order(codomain((mR)))))[2])

The natural `ZZ[H]` module where `H`, a subgroup of the
  automorphism group acts on the ray class group.
"""
function Oscar.gmodule(H::PermGroup, mR::MapRayClassGrp, mG = automorphism_group(PermGroup, Hecke.nf(order(codomain((mR)))))[2])
  k = Hecke.nf(order(codomain(mR)))
  G = domain(mG)

  ac = Hecke.induce_action(mR, [image(mG, G(g)) for g = gens(H)])
  return GModule(H, ac)
end

"""
    gmodule(R::ClassField, mG = automorphism_group(PermGroup, base_field(R))[2]; check::Bool = true)

The natural `ZZ[G]` module where `G`, the
  automorphism group, acts on the ideal group defining the class field.
"""
function Oscar.gmodule(R::ClassField, mG = automorphism_group(PermGroup, base_field(R))[2]; check::Bool = true)
  k = base_field(R)
  G = domain(mG)
  mR = R.rayclassgroupmap
  mq = R.quotientmap
  if check
    c, mc = conductor(R)
    @req all(x -> c == Hecke.induce_image(mG(x), c), gens(G)) "field is not normal"
    s1 = Set(mc)
    @req all(x -> s1 == Set(Hecke.induce_image(mG(x), y) for y = s1), gens(G)) "field is not normal"
  end

  ac = Hecke.induce_action(mR, [image(mG, g) for g = gens(G)], mq)
  return GModule(G, ac)
end

"""
    gmodule(H::PermGroup, R::ClassField, mG = automorphism_group(PermGroup, k))

The natural `ZZ[H]` module where `H`, a subgroup of the 
  automorphism group, acts on the ideal group defining the class field.
"""
function Oscar.gmodule(H::PermGroup, R::ClassField, mG = automorphism_group(PermGroup, k))
  k = base_field(R)
  G = domain(mG)
  mR = R.rayclassgroupmap
  mq = R.quotientmap

  ac = Hecke.induce_action(mR, [image(mG, G(g)) for g = gens(H)], mq)
  #TODO: think about adding a restriction map?
  return GModule(G, ac)
end

#TODO: think: this should probably all use MultGrpElem???
#      NO, the "module" is the abstract abelian group
function _gmodule(k::AbsSimpleNumField, H::PermGroup, mu::Map{FinGenAbGroup, <:Union{AbsSimpleNumField, FacElemMon{AbsSimpleNumField}}}, mG = automorphism_group(PermGroup, k)[2])
  u = domain(mu)
  U = [mu(g) for g = gens(u)]
  G = domain(mG)
  ac = [hom(u, u, [preimage(mu, mG(G(g))(x)) for x = U]) for g = gens(H)]
  return gmodule(H, ac)
end

function Oscar.gmodule(H::PermGroup, mu::Map{FinGenAbGroup, FacElemMon{AbsSimpleNumField}}, mG = automorphism_group(PermGroup, base_ring(codomain(mu)))[2])
  return _gmodule(base_ring(codomain(mu)), H, mu, mG)
end

function Oscar.gmodule(H::PermGroup, mu::Map{FinGenAbGroup, AbsSimpleNumFieldOrder}, mG = automorphism_group(PermGroup, Hecke.nf(codomain(mu)))[2])
  #TODO: preimage for sunits can fail (inf. loop) if
  # (experimentally) the ideals in S are not coprime or include 1
  # or if the s-unit is not in the image (eg. action and not closed set S)
  u = domain(mu)
  U = [mu(g) for g = gens(u)]
  zk = codomain(mu)
  k = Hecke.nf(zk)
  G = domain(mG)
  ac = [hom(u, u, [preimage(mu, zk(mG(G(g))(k(x)))) for x = U]) for g = gens(H)]
  return gmodule(H, ac)
end

function Oscar.gmodule(H::PermGroup, mu::Map{FinGenAbGroup, AbsSimpleNumField}, mG =  automorphism_group(PermGroup, codomain(mu))[2])
  return _gmodule(codomain(mu), H, mu, mG)
end

function is_coboundary(c::CoChain{2,PermGroupElem,MultGrpElem{AbsSimpleNumFieldOrderFractionalIdeal}})

  zk = order(first(values(c.d)).data)
  
  I = AbsSimpleNumFieldOrderIdeal[]
  for v = values(c.d)
    if !isone(v.data)
      for x = c.C.G
        J = action(c.C, x, v)
        n, d = integral_split(J.data)
        n in I || push!(I, n)
        d in I || push!(I, d)
      end
    end
  end
  cp = coprime_base(I)

  G = c.C.G
  k = number_field(zk)
  MI = c.C.M
  if length(cp) == 0
    z = zero(MI)
    return true, CoChain{1, PermGroupElem, elem_type(MI)}(c.C, Dict((g,) => z for g = G))
  end
  O = orbits(gset(G, (i, g) -> numerator(action(c.C, g, MI(i//1)).data), cp, closed = true))
  local res
  frst = true
  for _o = O
    for i = 1:2
      if i == 1
        o = _o
      else
        o = vcat([collect(keys(factor(x))) for x = _o]...)
      end
      M = abelian_group([0 for x = o])
      h = MapFromFunc(c.C.M, M, x->M([valuation(x.data, y) for y = o]))
      D = gmodule(G, [hom(M, M, [h(action(c.C, g, c.C.M(i//1))) for i = o]) for g = gens(c.C.G)])
#      @assert all(is_bijective, D.ac)

      cc = map_entries(h, c, parent = D)
      h2, mh2, icb = Oscar.GrpCoh.H_two(D)
      fl, x = icb(cc)
      if !fl
        if i == 2
          return false, nothing
        else
          continue
        end
      end
      h = MapFromFunc(M, MI, x->MI(prod((o[i]//1)^Int(x[i]) for i=1:length(o))))
      if frst
        frst = false
        res = map_entries(h, x, parent = c.C)
      else
        res += map_entries(h, x, parent = c.C)
      end
      break
    end
  end

  #problem:
  #if G = C_2 and I = P1 P2 for 2 conjugate ideals P1, P2, then
  #H^2(G, I) = C_2 (trivial operation), but H^2(G, <P1, P2}) = C_1
  #thus this test will fail
  #project down onto Galois orbits of ideals in cp
  # if it is split on that orbit, remove it
  #idea would be H^2(<cp>) = sum H^2(<o>) where the sum runs over the
  #Galois orbits in cp
  #if an orbit can be removed, fine. If not
  #either the ideals should factor (and can be removed)
  #or not possible
  #so: probably better to split this into Galois orbits...
  return true, res
end

function Oscar.map_entries(MI::Oscar.GrpCoh.MultGrp{AbsSimpleNumFieldOrderFractionalIdeal}, c::Oscar.GrpCoh.CoChain{2, <:GAPGroupElem, Oscar.GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}})
  k = c.C.M.data
  zk = maximal_order(k)
  D = gmodule(c.C.G, [MapFromFunc(MI, MI, x->MI(hom(k, k, action(c.C, g, c.C.M(gen(k))).data)(x.data))) for g = gens(c.C.G)])
  h = MapFromFunc(c.C.M, MI, x->MI(x.data*zk))
  return map_entries(h, c, parent = D)
end

Oscar.isfinite(M::Generic.FreeModule{ZZRingElem}) = rank(M) == 0

"""
    is_coboundary(c::CoChain{2,PermGroupElem,MultGrpElem{AbsSimpleNumFieldElem}})

For a 2-cochain with values in the multiplicative group of a number field,
decide if it is a coboundary (trivial in the Brauer group). If so also
return a splitting 1-cochain.
"""
function is_coboundary(c::CoChain{2,PermGroupElem,MultGrpElem{AbsSimpleNumFieldElem}})
  #TODO/think: test locally. Does not need the normalisation
  @vprint :GaloisCohomology 1 "testing if 2-chain is a boundary\n"

  zk = maximal_order(parent(first(values(c.d)).data))
  @vprint :GaloisCohomology 2 ".. reducing via ideals first ..\n"
  MI = Oscar.GrpCoh.MultGrp(Hecke.FracIdealSet(zk))
  fl, s = is_coboundary(map_entries(MI, c))
  fl || return fl, nothing
  ss = map_entries(MapFromFunc(s.C.M, c.C.M, x->c.C.M(Hecke.short_elem(inv(x.data)))), s, parent = c.C)
  c += Oscar.GrpCoh.differential(ss)

  @vprint :GaloisCohomology 2 ".. gathering primes in the support ..\n"

  cp = coprime_base(vcat([numerator(norm(x.data*denominator(x.data))) for x = values(c.d)],
                         map(x->denominator(x.data), values(c.d))))
  @vprint :GaloisCohomology 2 ".. coprime done, now factoring ..\n"
  s = Set(reduce(vcat, [prime_divisors(x) for x = cp], init = [1]))
  while 1 in s
    pop!(s, 1)
  end

  @vprint :GaloisCohomology 2 ".. class group ..\n"
  Cl, mCl = class_group(zk)
  if length(s) == 0
    S = Set{AbsSimpleNumFieldOrderIdeal}()
  else
    S = Set(reduce(vcat, [[x[1] for x = prime_decomposition(zk, p)] for p = s]))
  end
  
  @vprint :GaloisCohomology 2 ".. enlarge primes ..\n"
  q, mq = quo(Cl, [preimage(mCl, x) for x = S])
  p = 2
  while order(q) > 1
    p = next_prime(p)
    if p in s
      continue
    end
    lp = prime_decomposition(zk, p)
    cP = [mq(preimage(mCl, x[1])) for x= lp]
    if all(iszero, cP)
      continue
    end
    S = union(S, Set([x[1] for x = lp]))
    q, mmq = quo(q, cP)
    mq = mq*mmq
  end

  @vprint :GaloisCohomology 2 ".. S-units ..\n"
  if length(S) == 0
    u, mu = Hecke.unit_group_fac_elem(zk)
  else
    u, mu = Hecke.sunit_group_fac_elem(collect(S))
  end
  C = gmodule(Group(c.C), mu)

  @vprint :GaloisCohomology 2 ".. cohomology ..\n"
  H2, _, z = cohomology_group(C, 2)
  @vprint :GaloisCohomology 2 ".. map to abstract chain ..\n"
  cc = CoChain{2,PermGroupElem,FinGenAbGroupElem}(C, Dict((h, preimage(mu, FacElem(v.data))) for (h,v) = c.d))
  @vprint :GaloisCohomology 2 ".. test for boundary ..\n"
  fl, d = z(cc; reduce = true)
  if !fl
    @vprint :GaloisCohomology 2 ".. no boundary\n"
    return fl, d
  end
  @vprint :GaloisCohomology 2 ".. explicit boundary\n"
  return fl, CoChain{1,elem_type(c.C.G),elem_type(c.C.M)}(c.C, Dict((h, c.C.M(evaluate(mu(v)))) for (h,v) = d.d)) - ss
end

function isunramified(p::AbsSimpleNumFieldOrderIdeal)
  return ramification_index(p) == 1
end


"""
    decomposition_group(K::AbsSimpleNumField, mK::Map, mG::Map, mGp; _sub::Bool = false)
For a completion C of a number field K, implicitly given as the map
    mK:  K -> C
And the automorphism group G of K given via
    mG:  G -> Aut(K)
`mG` defaults to the full automorphism group of `K` as a permutation group.
And the automorphism group Gp of Kp, given via
    mGp: Gp -> Aut(Kp)
defaulting to the full automorphism group over the prime field.

Find the embedding of Gp -> G, realizing the local automorphism group
as a subgroup of the global one.
"""
function Oscar.decomposition_group(K::AbsSimpleNumField, mK::Map, mG::Map = automorphism_group(K)[2], mGp::Map = automorphism_group(codomain(mK), absolute_base_field(codomain(mK))); _sub::Bool = false)
  Kp = codomain(mK)
  @assert domain(mK) == K

  Gp = domain(mGp)
  G = domain(mG)

  im = elem_type(G)[]
  pm = elem_type(Gp)[]
  elG = [g for g = G]
  imK = [mK(mG(g)(gen(K))) for g = elG]
  if _sub
    gn = Gp
  else
    gn = gens(Gp)
  end
  for s = gn
    h = mGp(s)(mK(gen(K)))
    z = findall(isequal(h), imK)
    if length(z) == 0
      _sub && continue
      z = argmax([valuation(h-x) for x = imK], dims = 1)
    end
    @assert length(z) == 1
    push!(im, elG[z[1]])
    push!(pm, s)
  end
  if _sub 
    Gp, _mG = sub(Gp, pm)
    s = [preimage(_mG, x) for x = pm]
    return hom(Gp, G, s, im), _mG
  end
  return hom(Gp, G, pm, im)
end

"""
    decomposition_group(K::AbsSimpleNumField, emb::Hecke.NumFieldEmb, mG::Map)

For a real or complex embedding `emb`, find the unique automorphism
that acts on this embedding as complex conjugation.
`mG`, when given,  should be the map from the automorphism group. Otherwise
it will be automatically supplied.
"""
function Oscar.decomposition_group(K::AbsSimpleNumField, emb::Hecke.NumFieldEmb, mG::Map = automorphism_group(PermGroup, K)[2])
  G = domain(mG)
  if is_real(emb)
    return sub(G, [one(G)])[2]
  end
  g = gen(K)
  lG = [g for g  = G]
  l = findall(x->overlaps(conj(emb(g)), emb(mG(x)(g))), lG)
  if length(l) == 0
    return sub(G, [one(G)])[2]
  end
  @assert length(l) == 1
  sigma = lG[l[1]]
  return sub(G, [sigma])[2]
end

#= TODO
 - (DONE) induce a gmodule into a larger group
 - (DONE) direct sum/prod of gmodules
 - maps (a pair of G->H and N -> M or so)?
 - quotient?
 - the local/ global fund class, ie. normalize the cochain
 - map a local chain into a ray class group
=#

"""
    gmodule(K::Hecke.LocalField, k::Union{Hecke.LocalField, PadicField, QadicField} = base_field(K); Sylow::Int = 0, full::Bool = false)

For a local field extension K/k, return the multiplicative
group of K as a Gal(K/k) module. Strictly, it returns a quotient
of `K^*` that is cohomologically equivalent. For unramified fields,
the quotient is just `ZZ`, for tame fields, it will be isomorphic to
  `ZZ times F_q^*`
for the residue field `F_q`. In the wild case a suitable induced
module is factored out.

Returns: 
 - the gmodule
 - the map from G = Gal(K/k) -> Set of actual automorphisms
 - the map from the module into K
"""
function Oscar.gmodule(K::Hecke.LocalField, k::Union{Hecke.LocalField, PadicField, QadicField} = base_field(K); Sylow::Int = 0, full::Bool = false)

  #if K/k is unramified, then the units are cohomological trivial,
  #   so Z (with trivial action) is correct for the gmodule
  #if K/k is tame, then the 1-units are cohomological trivial, hence
  #   Z time k^* is enough...

  e = divexact(absolute_ramification_index(K), absolute_ramification_index(k))
  f = divexact(absolute_degree(K), e)
  @vprint :GaloisCohomology 1 "the local mult. group as a Z[G] module for e=$e and f = $f\n"
  @vprint :GaloisCohomology 2 " .. the automorphism group ..\n"

  G, mG = automorphism_group(PermGroup, K, k)

  if e == 1 && !full
    @vprint :GaloisCohomology 2 " .. unramified, only the free part ..\n"
#    @show :unram
    A = abelian_group([0])
    Hecke.assure_has_hnf(A)
    pi = uniformizer(K)
    return gmodule(G, [hom(A, A, [A[1]]) for g = gens(G)]),
      mG,
      MapFromFunc(A, K, x->pi^x[1], y->Int(e*valuation(y))*A[1])
  end

  if e % prime(K) != 0 && !full #tame!
    @vprint :GaloisCohomology 2 " .. tame, no 1-units ..\n"
#    @show :tame
    k, mk = residue_field(K)
    u, mu = unit_group(k)
    pi = uniformizer(K)
    # move to a Teichmueller lift?
    gk = preimage(mk, mu(u[1]))
    pr = precision(gk)
    gkk = setprecision(gk^order(k), pr)
    while !iszero(gkk - gk)
      gk = gkk
      gkk = setprecision(gk^order(k), pr)
    end
    A = abelian_group([0, order(u)])
    Hecke.assure_has_hnf(A)
    h = Map[]
    for g = gens(G)
      im = [A[1]+preimage(mu, mk(mG(g)(pi)*inv(pi)))[1]*A[2], preimage(mu, mk(mG(g)(gk)))[1]*A[2]]
      push!(h, hom(A, A, im))
    end
    return gmodule(G, h),
      mG,
      MapFromFunc(A, K, x->pi^x[1] * gk^x[2],
        function(y)
          v = Int(e*valuation(y))
          y *= pi^-v
          return v*A[1] + preimage(mu, mk(y))[1]*A[2]
        end)
  end
 
#  @show :wild
  @vprint :GaloisCohomology 2 " .. wild case (or requested), unit group ..\n"
  U, mU = unit_group(K)
  n = divexact(absolute_degree(K), absolute_degree(k))
  @assert order(G) == n

  @vprint :GaloisCohomology 2 " .. find lattice (normal basis) ..\n"
  b = absolute_basis(K)
  # need a normal basis for K/k, so the elements need to be k-lin. indep
  local o, best_o
  cnt = 0
  while true
    a = sum(b[i]*rand(-5:5) for i=1:length(b))
    o = [mG(g)(a) for g = G]
    m = matrix(k, n, n, vcat([coordinates(x, k) for x = o]...))
    dm = det(m)
    cnt += 1
    if cnt > 10
      error("dnw")
    end
    if iszero(dm) #|| valuation(dm) > 5
      continue
    else
      break
    end
  end

  #o needs to be expanded to be an absolute basis
  b = absolute_basis(k)
  o = [x*y for x = b for y = o]


  @vprint :GaloisCohomology 2 " .. quotient ..\n"
  if prime(k) == 2
    #we need val(p^k) > 1/(p-1)
    #val(p) = 1 and the only critical one is p=2, where k>1 is
    #necessary
    ex = 2
  else
    ex = 1
  end
  #x -> 1+pi*x is in general, not injective, not even for a basis
  # if valuation(dm) == 0, then by Lorenz Alg II, 26.F10 it should
  # be, but we're not using it. This was used to avoid exp
  Q, mQ = quo(U, [preimage(mU, exp(prime(k)^ex*x)) for x = o])
  S, mS = snf(Q)
  Q = S
  mQ = mQ*inv(mS)

  if Sylow > 0
    @assert is_prime(Sylow)
    G, mS = sylow_subgroup(G, Sylow)
    mG = mS*mG
  end

  @vprint :GaloisCohomology 2 " .. the module ..\n"
  hh = [hom(Q, Q, [mQ(preimage(mU, mG(i)(mU(preimage(mQ, g))))) for g = gens(Q)]; check = false) for i=gens(G)]
  Hecke.assure_has_hnf(Q)
  return gmodule(G, hh), mG, pseudo_inv(mQ)*mU
end

#=  Not used
function one_unit_cohomology(K::Hecke.LocalField, k::Union{Hecke.LocalField, PadicField, QadicField} = base_field(K))

  U, mU = Hecke.one_unit_group(K)
  G, mG = automorphism_group(PermGroup, K, k)

  b = absolute_basis(K)
  local o
  while true
    a = uniformizer(K)^30*sum(b[i]*rand(-5:5) for i=1:length(b))
    o = [mG(g)(a) for g = G]
    if length(Set(o)) == order(G)
      break
    end
  end

  S, mS = sub(U, [preimage(mU, 1+x) for x = o])
  Q, mQ = quo(U, S)
  hh = [hom(Q, Q, [mQ(preimage(mU, mG(i)(mU(preimage(mQ, g))))) for g = gens(Q)]) for i=gens(G)]
  return gmodule(G, hh)
end

=#
#= TODO
 - (DONE) induce a gmodule into a larger group
 - (DONE) direct sum/prod of gmodules
 - maps (a pair of G->H and N -> M or so)?
 - quotient?
 - the local/ global fund class, ie. normalize the cochain
 - map a local chain into a ray class group
=#

#= TODO
  for Z, Z/nZ, F_p and F_q moduln -> find Fp-presentation
  for finite Z, Z/nZ, F_p and F_q moduln -> find pc-presentation
  #done: for FinGenAbGroup          -> find Fp-presentation
  #done: for FinGenAbGroup          -> find pc-presentation
  #done: for a in H^2 find Fp-presentation
  #done: for a in H^2 find pc-presentation
  for a in H^2 find (low degree) perm group using the perm group we have?
  Magma's DistinctExtensions
  probably aut(GrpAb), ...

Sort: 
 - move the additional FinGenAbGroupHom stuff elsewhere
 - move (and fix) the ModuleHom stuff
 - add proper quo for Modules (done)
 - split generic coho/ gmodule and number theory  (partly done)

  features   
   - inflation (done), restriction (done), long exact sequence  
   - induction (done)/ coinduction ...
   - restriction (of gmodules to Sylow subgroups)
   - think about Debeerst: if P, Q are above the some prime then
     Ind_G_P^G L_P = Ing_G_Q^G L_Q??? (no - aprently yes)
   - use prod_Q|P L_Q rather than prod Ind...  

  dreams
   - we we extend to H^-1, ^-2?
   - H^3 (in some cases)
   - cup products
   - the relative cohomology
     https://arxiv.org/pdf/1809.01209.pdf
     https://doi.org/10.1017/S2040618500033050

  GModule for 
    - (done for mult grp) local field (add (trivial) and mult)
    - (done) (S-)units
    - (done) Ali's stuff.... (in progress: see Hecke/src/LocalField/neq.jl)
    - local field (add (trivial) and mult)
    - (S-)units
=#    

#TODO: what do we need to return?
# - mG (if we cache this in the field, not necessary)
# - the local stuff?
# - the S-Units?
# - ???
# - a different type containing all this drivel? (probably)
#    YES - need to have maps to and from local stuff
#    use Klueners/ Acciaro to map arbitrary local into idele
#    use ...               to project to ray class
# - a magic(?) function to get idele approximations in and out?

"""
M has to be a torsion free Z module with a C_2 action by sigma.
Return data for the decomposition into indecomposables.
They will be of type
 - Z with trivial and non-trivial action
 - Z[C_2]

Two arrays are returned:
 - generators for the 1-dim modules
 - C_2 generators for the 2-dim ones

Follows Debeerst rather closely...

(Helper for the idele class stuff)
"""
function debeerst(M::FinGenAbGroup, sigma::Map{FinGenAbGroup, FinGenAbGroup})
  @assert domain(sigma) == codomain(sigma) == M
  @assert all(x->sigma(sigma(x)) == x, gens(M))
  @assert is_free(M) && torsion_free_rank(M) == ngens(M)

  K, mK = kernel(id_hom(M)+sigma)
  fl, mX = has_complement(mK)
  @assert fl
  X = domain(mX)
  _X, _mX = snf(X)

  S, mS = image(sigma -id_hom(M))
  fl, mSK = is_subgroup(S, K)
  @assert fl


  _K, _mK = snf(K)
  _S, _mS = snf(S)
  @assert is_trivial(_S) || torsion_free_rank(_S) == ngens(_S) 
  @assert torsion_free_rank(_K) == ngens(_K) 

  m = matrix(FinGenAbGroupHom(_mS * mSK * inv((_mK))))
  # elt in S * m = elt in K
  # so
  # elt in S * U^-1 U m V V^-1 = elt_in K
  # elt in S * U^-1 snf = elt_in * V
  s, U, V = snf_with_transform(m)
  if is_trivial(S)
    r = 0
  else
    r = maximum(findall(x->isone(s[x,x]), 1:ngens(_S)))
  end

  mu = hom(_S, _S, inv(U))
  mv = hom(_K, _K, V)
  @assert is_trivial(S) || all(i-> M(_mS(mu(gen(_S, i)))) == s[i,i] * M(_mK(mv(gen(_K, i)))), 1:ngens(S))
  b = [_mK(mv(x)) for x = gens(_K)]

  Q, mQ = quo(S, image(sigma -id_hom(M), K)[1])
  B, mB = sub(Q,  [mQ(preimage(mSK, x)) for x = b[1:r]])
  @assert order(B) == order(Q)

  phi = FinGenAbGroupHom(_mX*mX*(sigma -id_hom(M))*pseudo_inv(mS)*mQ)
  @assert is_surjective(phi)
  A = vcat([preimage(mB, phi(k)).coeff for k = gens(_X)]...)
  h, t = hnf_with_transform(A)
  #t*A = h = diag(1) vcat 0
  x = [sum(t[i,j]*_X[j] for j=1:ngens(_X)) for i=1:ngens(_X)]
  sm1 = sigma - id_hom(M)
  sm1_K = hom(K, M, [sm1(mK(x)) for x= gens(K)])
  lambda = vcat([preimage(sm1_K, sm1(mX(_mX(x[i]))) - mK(b[i])) for i=1:r],
                [preimage(sm1_K, sm1(mX(_mX(x[i])))) for i=r+1:length(x)])
  x = map(_mX*mX, x)
  lambda = map(mK, lambda)
  y = x .- lambda
  b = map(mK, b)

  #just checking the action on the 2-dim stuff. magic.
  #= (s-1)x = b + (s-1)l  1..r
     y = x-l
     s(y) = s(x) - s(l) = s(x)-x + x - s(l) + l - l
           = (s-1)x + x -(s-1)l -l
           = b + (s-1)l + x -(s-1)l - l = b + x-l = b + y
  =#        
  #=
  h2 = []
  for i=1:r
    a = abelian_group([0,0])
    push!(h2, (hom(a, M, [b[i], y[i]]), -y[i]-b[i]))
  end
  h_minus = []
  for i=r+1:length(b)
    a = abelian_group([0])
    push!(h_minus, (hom(a, M, [b[i]]), b[i]))
  end
  h_plus = []
  for i=r+1:length(y)
    a = abelian_group([0])
    push!(h_plus, (hom(a, M, [y[i]]), y[i]))
  end
  =#

  return vcat(b[r+1:end], y[r+1:end]), [-y[i] - b[i] for i=1:r]
end

function Hecke.extend_easy(m::Hecke.CompletionMap, L::FacElemMon{AbsSimpleNumField})
  k = base_ring(L)
  @assert k == domain(m)

  #want a map: L-> codomain(m)
  function to(a::FacElem{AbsSimpleNumFieldElem})
    return prod(m(k)^v for (k,v) in a)
  end
  function from(a::Hecke.LocalFieldElem)
    return FacElem(preimage(m, a))
  end
  return MapFromFunc(L, codomain(m), to, from)
end

function Hecke.extend_easy(m::Hecke.CompletionMap, mu::Map, L::FacElemMon{AbsSimpleNumField})
  k = base_ring(L)
  @assert k == domain(m)
  @assert codomain(mu) == codomain(m)

  cache = Dict{AbsSimpleNumFieldElem, FinGenAbGroupElem}()
  #want a map: L-> codomain(m) -> domain(mu)
  function to(a::FacElem{AbsSimpleNumFieldElem})
    s = domain(mu)[0]
    for (k,v) = a
      if haskey(cache, k)
        s += v*cache[k]
      else
        kk = preimage(mu, m(k))
        cache[k] = kk
        s += v*kk
      end
    end
    return s
  end
  function from(a::Hecke.LocalFieldElem)
    return FacElem(preimage(m, mu(a)))
  end
  return MapFromFunc(L, domain(mu), to, from)
end


"""
    IdeleParent

A container structure for computations around idele class cohomology.
Describes a gmodule that is cohomologically equivalent to the 
idele class group, following Chinberg/Debeerst/Aslam.
"""
mutable struct IdeleParent
  k::AbsSimpleNumField
  mG::Map # AutGrp -> Automorphisms
  S::Vector{AbsNumFieldOrderIdeal} # for each prime number ONE ideal above
  C::Vector{Map} # the completions at S
  D::Vector{Map} # Gp -> Aut
  L::Vector{Map} # the mult. group map at S

  #for P in S the modules used actually is
  #    Ind_G_p^G L[P]
  #        = sum L[P] otimes s_i
  # (for s_i a fixed system of coset reps G//G_P)
  # L[P] otimes s_i "should be" the completion data at P^s_i - one of the other ideals
  # should be L[P] ni l -> C[P] -> k -> inv(s_i)(..) to get a proper rep in k
  # completion at P^s is C[P] but with the map twisted by s

  mU::Map #S-unit group map
  M::FinGenAbGroup  # the big module, direct product from
    # infinite gmodule x finite ones
  mq::Map # "projection" of M -> the actual module in the end

  data

  function IdeleParent()
    return new()
  end
end

"""
    idele_class_gmodule(k::AbsSimpleNumField, s::Vector{Int} = Int[])

Following Debeerst:
  Algorithms for Tamagawa Number Conjectures. Dissertation, University of Kassel, June 2011.
or Ali, 

Find a gmodule C s.th. C is cohomology-equivalent to the cohomology
of the idele class group. The primes in `s` will always be used.
"""
function idele_class_gmodule(k::AbsSimpleNumField, s::Vector{Int} = Int[]; redo::Bool=false)
  @vprint :GaloisCohomology 2 "Ideal class group cohomology for $k\n"
  I = get_attribute(k, :IdeleClassGmodule)
  if !redo && I !== nothing
    return I
  end

  I = IdeleParent()
  I.k = k
  G, mG = automorphism_group(PermGroup, k)
  I.mG = mG

  zk = maximal_order(k)

  sf = subfields(k)
  sf = [x[1] for x = sf if degree(x[1]) > 1]
  zf = map(maximal_order, sf)
  cf = map(class_group, zf)
  cf = Tuple{FinGenAbGroup, <:Map}[x for x = cf]

  @vprint :GaloisCohomology 2 " .. gathering primes ..\n"
  s = push!(Set{ZZRingElem}(s), Set{ZZRingElem}(prime_divisors(discriminant(zk)))...)
  for i=1:length(sf)
    l = factor(prod(s)*zf[i])
    q, mq = quo(cf[i][1], [preimage(cf[i][2], P) for P = keys(l)])
    cf[i] = (q, pseudo_inv(mq)*cf[i][2])
  end

  #think: does the quotient have to be trivial - or coprime to |G|?
  #coprime should be enough

  for p = PrimesSet(2, -1)
    p in s && continue
    all(x->order(x[1]) == 1, cf) && break
    new = false
    for i=1:length(sf)
      l = factor(p*zf[i])
      q, mq = quo(cf[i][1], [preimage(cf[i][2], P) for P = keys(l)])
      if order(q) != order(cf[i][1])
        new = true
      end
      cf[i] = (q, pseudo_inv(mq)*cf[i][2])
    end
    if new
      push!(s, p)
    end
  end

  S = collect(keys(factor(prod(s)*zk)))
  @vprint :GaloisCohomology 2 " .. need $(length(S)) prime ideals ..\n"

  s = [findfirst(x->minimum(x) == t, S) for t = s]
  @vprint :GaloisCohomology 2 " .. split into $(length(s)) G-orbits ..\n"

  @vprint :GaloisCohomology 2 " .. S-units (for all) ..\n"
  U, mU = sunit_group_fac_elem(S)
  I.mU = mU
  z = MapFromFunc(codomain(mU), k, evaluate, FacElem)
  E = gmodule(G, mU, mG)
  Hecke.assure_has_hnf(E.M)
  @hassert :GaloisCohomology -1 is_consistent(E)

  if is_totally_real(k)
    @vprint :GaloisCohomology 2 " .. real field, easy case ..\n"
    mG_inf = Oscar.decomposition_group(k, real_embeddings(k)[1], mG)
    G_inf = domain(mG_inf)
    Et = gmodule(G_inf, mU, mG)
    @hassert :GaloisCohomology 1 is_consistent(Et)
    iEt = Oscar.GrpCoh.induce(Et, mG_inf, E, id_hom(U))
  else
    #TODO: failing on x^8 + 70*x^4 + 15625
    @vprint :GaloisCohomology 2 " .. complex field, hard case ..\n"
    mG_inf = Oscar.decomposition_group(k, complex_embeddings(k)[1], mG)
    G_inf = domain(mG_inf)
    sigma = action(E, mG_inf(G_inf[1]))
    @assert order(G_inf[1]) == 2 == order(G_inf)

    @assert order(U[1]) >0
    q, mq = quo(U, [U[1]]) 
    q, _mq = snf(q)
    mq = mq*pseudo_inv(_mq)
    sigma_q = hom(q, q, [mq(sigma(preimage(mq, x))) for x = gens(q)])
    x, y = debeerst(q, sigma_q)
    # just to verify... Gunter Malle: the C_2 modules are visible over GF(2)...
    _M = gmodule(GF(2), gmodule(G_inf, [sigma_q]))
    _i = indecomposition(_M)
    @assert length(findall(x->dim(x[1]) == 2, _i)) == length(y)
    @assert length(findall(x->dim(x[1]) == 1, _i)) == length(x)
      #possibly: now the H^2 is correct, but the H^1 is not...
      # x^8 - 12*x^7 + 44*x^6 - 24*x^5 - 132*x^4 + 120*x^3 + 208*x^2 - 528*x + 724

    #theta:
    theta = U[1] #should be a generator for torsion, torsion is even,
                 #hence this elem cannot be a square

    x = [preimage(mq, i) for i = x]
    y = [preimage(mq, i) for i = y]

    z, mz = sub(U, [sigma(U[1]) - U[1]])
    theta_i = [sigma(t)-t for t = x]
    inv = Int[]
    not_inv = Int[]
    for i=1:length(x)
      w = theta_i[i]
      fl, pe = has_preimage_with_preimage(mz, w)
      if fl
        push!(inv, i)
        zz = mq(x[i])
        x[i] -= pe[1]*U[1]
        @assert zz == mq(x[i])
        theta_i[i] = sigma(x[i]) - x[i]
        @assert iszero(theta_i[i])
      else
        push!(not_inv, i)
      end
    end

    @assert length(not_inv) > 0
    @assert length(not_inv) + length(inv) == length(x)
    x = vcat(x[not_inv], x[inv]) #reordering
    theta_i = vcat(theta_i[not_inv], theta_i[inv])
    
    U_t, mU_t = sub(U, [U[1]])
    
    sm1 = hom(U_t, U, [sigma(mU_t(g)) - mU_t(g) for g = gens(U_t)])
    eta_i = [preimage(sm1, theta - theta_i[i]) for i=1:length(not_inv)]

    eta_i = map(mU_t, eta_i)
    V = abelian_group(elementary_divisors(U))
    im_psi = [U[1], x[1]+ eta_i[1]]
    for i=2:length(not_inv)
      push!(im_psi, x[i] - x[1] + eta_i[i] - eta_i[1])
      @assert sigma(im_psi[end]) == im_psi[end]
      #should be chosen to be pos. at place, flip signs...
    end
    for i=length(not_inv)+1:length(x)
      push!(im_psi, x[i])
      @assert sigma(im_psi[end]) == im_psi[end]
      #should be chosen to be pos. at place, flip signs...
    end
    for i=1:length(y)
      push!(im_psi, y[i])
      push!(im_psi, sigma(y[i]))
    end
    psi = hom(V, U, im_psi)
    @assert all(i->psi(V[i]) == im_psi[i], 1:length(im_psi))
    @assert is_bijective(psi)
    F = abelian_group([0 for i=2:length(x)])
    Hecke.assure_has_hnf(F)
    W, pro, inj = direct_product(V, F, task = :both)
    @assert isdefined(W, :hnf)

    ac = FinGenAbGroupHom(pro[1]*psi*sigma*pseudo_inv(psi)*inj[1]) +
         FinGenAbGroupHom(pro[2]*hom(F, W, FinGenAbGroupElem[inj[1](V[i+1]) - inj[2](F[i-1]) for i=2:length(x)]))

    Et = gmodule(G_inf, [ac])
    @assert is_consistent(Et)
    mq = pseudo_inv(psi)*inj[1]
    iEt = Oscar.GrpCoh.induce(Et, mG_inf, E, mq)
  end
  #test if the G-action is the same:
  # induce returns a map U -> E that should be a Z[G]-hom
  function is_G_lin(U, E, mUE, acU)
    G = E.G
    for g = gens(G)
      for u = gens(U)
        a = mUE(u)
        b = mUE(acU(g)(u))
        @assert b == action(E, g, a)
      end
    end
    return true
  end
  @hassert :GaloisCohomology 1 is_G_lin(U, iEt[1], iEt[2], g->action(E, g))
  @hassert :GaloisCohomology 1 is_consistent(iEt[1])
  
  S = S[s]
  I.S = S

  #TODO: precision: for some examples the default is too small
  @vprint :GaloisCohomology 2 " .. gathering the completions ..\n"
  Hecke.pushindent()
  L = [completion(k, x, 40*ramification_index(x)) for x = S]
  I.C = [x[2] for x = L]
  Hecke.popindent()
  @vprint :GaloisCohomology 2 " .. gathering the local modules ..\n"
  Hecke.pushindent()
  C = [gmodule(x[1], absolute_base_field(x[1])) for x = L];
  I.D = [x[2] for x = C]
  I.L = [x[3] for x = C]
  @hassert :GaloisCohomology 1 all(x->is_consistent(x[1]), C)
  D = [Oscar.GrpCoh.induce(C[i][1], Oscar.decomposition_group(k, L[i][2], mG, C[i][2]), E, (mU*Hecke.extend_easy(L[i][2], C[i][3], codomain(mU)))) for i=1:length(S)]
  @hassert :GaloisCohomology 1 all(x->is_consistent(x[1]), D)
  @hassert :GaloisCohomology 1 all(x->is_G_lin(U, D[x][1], D[x][2], g->action(E, g)), 1:length(D))
  Hecke.popindent()
  @vprint :GaloisCohomology 2 " .. the big product and the quotient\n"
  @assert isdefined(iEt[1].M, :hnf)
  @assert all(x->isdefined(x[1].M, :hnf), D)

  F = direct_product(iEt[1], [x[1] for x = D]..., task = :both)
  I.M = F[1].M

  @hassert :GaloisCohomology 1 is_consistent(F[1])

  h = iEt[2]*F[3][1]+sum(D[i][2]*F[3][i+1] for i=1:length(S));
  @vtime :GaloisCohomology 2 q, mq = quo(F[1], h)
  @hassert :GaloisCohomology 1 is_consistent(q)
  @vtime :GaloisCohomology 2 q, _mq = simplify(q)
  @vtime :GaloisCohomology 2 mq = FinGenAbGroupHom(mq * pseudo_inv(_mq))
  @hassert :GaloisCohomology 1 is_consistent(q)
  I.mq = mq
  I.data = (q, F[1])
  set_attribute!(k, :IdeleClassGmodule=>I)
  return I
end

"""
returns a `Dict` containing the prime ideals where the idele
represented by `a` as an element of `I` can be non-trivial, 
ie. at the places in `I`.
The values are elements in the corresponding completion.
"""
function idele_finite(I::IdeleParent, a::FinGenAbGroupElem)
  a = preimage(I.mq, a)
  #a lives in a direct product
  #  1st the units - induced, thus a direct product
  #  2nd ... the data at primes in S:
  #   induced from the completion at P
  d = Dict{AbsSimpleNumFieldOrderIdeal, Tuple{AbsSimpleNumFieldElem, Hecke.LocalFieldElem}}()
  for i=1:length(I.S)
    p = I.S[i]
    lp = prime_decomposition(order(p), minimum(p))
    for P = lp
      Kp, nKp, mGp, mUp, pro, inj = completion(I, P[1])
      x = mUp(pro(a))
      d[P[1]] = (preimage(nKp, x), x)
    end
  end
  return d
end

#from is_local_norm in Hecke
#based on Klueners/ Acciaro
function _local_norm(m0::AbsSimpleNumFieldOrderIdeal, a::AbsNumFieldOrderElem, p::AbsNumFieldOrderIdeal)
  v1 = valuation(a, p)
  v2 = valuation(m0, p)
  n0 = divexact(m0, p^v2)
  o0 = p^(v1 + v2)
  y = crt(order(p)(1), n0, a, o0)
  Y = y*order(p)
  Y = divexact(Y, p^v1)
  return Y
end


#maybe we need Idele's as independent objects?
#realizes C -> Cl (or the coprime version into a ray class group:
#for the idele `a` in `I` find an "equivalent" ideal.
function Oscar.ideal(I::IdeleParent, a::FinGenAbGroupElem; coprime::Union{AbsSimpleNumFieldOrderIdeal, Nothing})
  a = preimage(I.mq, a)
  zk = maximal_order(I.k)
  o_zk = zk
  if coprime !== nothing && order(coprime) !== zk
    o_zk = order(coprime)
    coprime = zk*coprime
    @assert order(coprime) === zk
  end
  id = 1*zk
  for p = I.S
    lp = prime_decomposition(zk, minimum(p))
    for P = lp
      Kp, nKp, mGp, mUp, pro, inj = completion(I, P[1])
      x = mUp(pro(a))
      if coprime === nothing
        id *= fractional_ideal(zk, P[1])^Int(valuation(x)*absolute_ramification_index(parent(x)))
      else
        if valuation(x) > 0
          id *= _local_norm(coprime, zk(preimage(nKp, x)), P[1])
        else
          id *= inv(_local_norm(coprime, zk(preimage(nKp, inv(x))), P[1]))
        end
      end
    end
  end
  return o_zk*id
end

function Oscar.galois_group(A::ClassField)
  return permutation_group(codomain(A.quotientmap))
end

"""
    galois_group(A::ClassField, ::QQField; idele_parent::IdeleParent = idele_class_gmodule(base_field(A)))

Determine the Galois group (over `QQ`) of the extension parametrized by `A`.
`A` needs to be normal over `QQ`.

# Examples
```jldoctest; setup = :(using Oscar), filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> QQx, x = QQ[:x];

julia> k, a = number_field(x^4 - 12*x^3 + 36*x^2 - 36*x + 9);

julia> a = ray_class_field(8*3*maximal_order(k), n_quo = 2);

julia> a = filter(is_normal, subfields(a, degree = 2));

julia> I = idele_class_gmodule(k);

julia> b = [galois_group(x, QQ, idele_parent = I) for x = a];

julia> [describe(x[1]) for x = b]
3-element Vector{String}:
 "C4 x C2"
 "Q8"
 "D8"

```
"""    
function Oscar.galois_group(A::ClassField, ::QQField; idele_parent::Union{IdeleParent,Nothing} = nothing) 

  m0, m_inf = defining_modulus(A)
  @assert length(m_inf) == 0
  if !is_normal(A)
    A = normal_closure(A)
  end
  mR = A.rayclassgroupmap
  mQ = A.quotientmap
  zk = order(m0)
  @req order(automorphism_group(Hecke.nf(zk))[1]) == degree(zk) "base field must be normal"
  if gcd(degree(A), degree(base_field(A))) == 1
    s, ms = split_extension(FPGroup, gmodule(A))
    return permutation_group(s), ms
  end
  if idele_parent === nothing
    idele_parent = idele_class_gmodule(base_field(A))
  end

  qI = cohomology_group(idele_parent, 2)
  q, mq = snf(qI[1])
  a = qI[2](image(mq, q[1])) # should be a 2-cycle in q
  @assert Oscar.GrpCoh.istwo_cocycle(a)
  gA = gmodule(A, idele_parent.mG)
  qA = cohomology_group(gA, 2)
  n = degree(Hecke.nf(zk))
  aa = map_entries(a, parent = gA) do x
    x = parent(x)(Hecke.mod_sym(x.coeff, gcd(n, degree(A))))
    J = ideal(idele_parent, x, coprime = m0)
    mQ(preimage(mR, numerator(J)) - preimage(mR, denominator(J)*zk))
  end
  @assert Oscar.GrpCoh.istwo_cocycle(aa)
  return permutation_group(aa), (aa, gA)
end

function Oscar.components(A::FinGenAbGroup)
  return get_attribute(A, :direct_product)
end

function Oscar.completion(I::IdeleParent, P::AbsNumFieldOrderIdeal)
  s = [minimum(x) for x = I.S]
  p = findfirst(isequal(minimum(P)), s)
  @assert p !== nothing

  mKp = I.C[p]
  Kp = codomain(mKp)
  mUp = I.L[p]
  mGp = I.D[p]

  inj = canonical_injection(I.M, p+1) #units are first
  pro = canonical_projection(I.M, p+1)


  @assert domain(inj) == codomain(pro)

  J = components(I.M)[p+1]
  if mKp.P == P #easy case
    return Kp, mKp, mGp, mUp, pro * canonical_projection(J, 1), canonical_injection(J, 1)*inj
  end

  prm = get_attribute(J, :induce)[2]
  mG = I.mG

  z = findall(pr -> mG(pr)(mKp.P) == P, prm)
  pr = inv(prm[z[1]])
  
  nKp = MapFromFunc(I.k, Kp, x->mKp(mG(pr)(x)), y->mG(inv(pr))(preimage(mKp, y)))

  return Kp, nKp, mGp, mUp, pro * canonical_projection(J, z[1]), canonical_injection(J, z[1])*inj 
end

function Oscar.map_entries(mp::Union{Map, Function}, C::GrpCoh.CoChain{N, G, M}; parent::GModule) where {N, G, M}
  d = Dict( k=> mp(v) for (k,v) = C.d)
  if isdefined(C, :D)
    return GrpCoh.CoChain{N, G, elem_type(parent.M)}(parent, d, x->mp(C.D(x)))
  else
    return GrpCoh.CoChain{N, G, elem_type(parent.M)}(parent, d)
  end
end

#=
 Currently automorphism_group maps are not cached and thus will return
 different groups each time. Thus to be consistent bewtween calls one 
 should compute the group (map) once and pass it around.

 Maybe we introduce the caching - however, there might be some data
 dependency problems adding "random" fields into the attributes
=#

"""
For a 2-cochain with with values in the multiplicative group of a number field,
compute the local invariant of the algebra at the completion at the prime
ideal or embeddings given.
"""
function local_index(CC::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}, P::AbsSimpleNumFieldOrderIdeal, mG::Map = automorphism_group(PermGroup, Hecke.nf(order(P)))[2]; B = nothing, index_only::Bool = false)
  return local_index([CC], P, mG, B = B, index_only = index_only)[1]
end

function local_index(C::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}, emb::Hecke.NumFieldEmb, mG::Map = automorphism_group(PermGroup, emb.K)[2]; index_only::Bool = false)
  return local_index([C], emb, mG)[1]
end

function local_index(CC::Vector{GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}}, emb::Hecke.NumFieldEmb, mG::Map = automorphism_group(PermGroup, Hecke.nf(order(P)))[2]; index_only::Bool = false)
  if is_real(emb)
    return [Hecke.QmodnZ()(0) for x = CC]
  end
  k = number_field(emb)
  _m = decomposition_group(k, emb, mG)

  Q = Hecke.QmodnZ()
  Gp = domain(_m)
  if order(Gp) == 1
    return [Q(0) for i = CC]
  end
#  @show [collect(values(x.d)) for x = CC]
  g = _m(gens(Gp)[1])
  @assert order(Gp) == 2 && order(g) == 2


  r = elem_type(Q)[]
  for C = CC
    x = C(g, g)
    y = emb(x.data)
#    @assert is_real(y)
    if is_positive(real(y))
      push!(r, Q(0))
    else
      push!(r, Q(1//2))
    end
  end
  return r
end

function cyclic_generator(G::Oscar.GAPGroup)
  while true
    g = rand(G)
    order(g) == order(G) && return g
  end
end

"""
Let `mkK: k -> K` be an embedding map of number fields,
`ml: k -> l` the completion map at some prime `p` in `k` and
`mL: K -> L` the completion map at some prime `P` in `K` above `p`.
This function returns the embedding `l -> L` induced by this data.
"""
function induce_hom(ml::Hecke.CompletionMap, mL::Hecke.CompletionMap, mkK::NumFieldHom)
  @assert domain(ml) == domain(mkK)
  @assert domain(mL) == codomain(mkK)
  l = codomain(ml)
  L = codomain(mL)
  #to make sure domain/codomain have the same precision...
  #otherwise, we'll get cool errors if the data does not 
  #define a homomorphism - due to not enough prec.
  rel_e = divexact(absolute_ramification_index(L), absolute_ramification_index(l))
  pr_ml = ml.precision
  pr_mL = mL.precision
  setprecision!(ml, max(pr_ml*rel_e, pr_mL))
  setprecision!(mL, max(pr_ml*rel_e, pr_mL))
  im_data = mL(mkK(preimage(ml, gen(l))))
  #CompletionMap is always Eisenstein/Unram
  #so need embedding of the unram parts:
  bl = base_field(l)
  bL = base_field(L)
  @assert isa(bl, QadicField)
  @assert isa(bL, QadicField)
  @assert degree(bL) % degree(bl) == 0
  rL, mrL = residue_field(bL)
  rl, mrl = residue_field(bl)
  rt = coeff(mL(mkK(preimage(ml, l(preimage(mrl, gen(rl)))))), 0)
  f = map_coefficients(x->preimage(mrL, rL(x)), defining_polynomial(rl))
  rt = Hecke.refine_roots1(f, [rt])[1]
#  setprecision!(ml, pr_ml)
#  setprecision!(mL, pr_mL)
  return hom(l, L, hom(bl, bL, rt), im_data)
end

function local_index(CC::Vector{GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}}, P::AbsSimpleNumFieldOrderIdeal, mG::Map = automorphism_group(PermGroup, Hecke.nf(order(P)))[2]; B::Any = nothing, index_only::Bool = false)
  
  k = Hecke.nf(order(P))

  if B !== nothing && haskey(B.lp, P)
    data = B.lp[P]
    mL = data[1]
    L = domain(mL)
    mU = data[2]
    mGp = data[3]
    emb, _m = data[4]
    C = data[5]
    c, mq = data[6]
    q = domain(mq)
    cn = data[7]
  else
    L, mL = completion(k, P)#, 40*ramification_index(P))
    C, mGp, mU = gmodule(L, absolute_base_field(L))

    G = domain(mG)
    emb, _m = decomposition_group(k, mL, mG, mGp, _sub = true)
    # _m : domain(emb) -> Gp
    #emb :             -> G 
    @assert domain(emb) === domain(_m)
    @assert codomain(emb) === G
    @assert codomain(_m) === domain(mGp)

    C = restrict(C, _m)

    c = cohomology_group(C, 2)
    q, mq = snf(c[1])
    @assert order(q) == order(domain(emb))


    if order(q) == 1
      if B !== nothing
        B.lp[P] = (mL, mU, mGp, (emb, _m), C, (c, mq), 1)
      end
      return [Hecke.QmodnZ()(0) for x = CC]
    end

    if index_only
      @assert B === nothing
      cn = ZZRingElem(1)
    else
      pp = minimum(B.mkK, P)   
      if is_cyclic(domain(_m)) &&
         ramification_index(pp) == ramification_index(P)

        Gp = domain(_m)
        #need THE Frobenius
        k, mk = residue_field(L)
        gk = gen(k)
        im_gk = gk^norm(pp)
        fr = [x for x = Gp if mk(mGp(_m(x))(preimage(mk, gk))) == im_gk]
        @assert length(fr) == 1
        g = fr[1]
        #and a uniformizer of the small field.
        x = preimage(mU, mL(B.mkK(B.k(uniformizer(pp)))))
        #should be a non-norm...
        can = CoChain{2, PermGroupElem, FinGenAbGroupElem}(C, Dict{NTuple{2, PermGroupElem}, FinGenAbGroupElem}((g^i, g^j) => i+j<order(q) ? zero(parent(x)) : x for i=0:order(q)-1 for j=0:order(q)-1))
      else
        l, ml = completion(B.k, pp)
        setprecision!(ml, precision(codomain(mL)))
        mlL = induce_hom(ml, mL, B.mkK)
        s = Hecke.Hecke.local_fundamental_class_serre(mlL)
        can = CoChain{2, PermGroupElem, FinGenAbGroupElem}(C, Dict{NTuple{2, PermGroupElem}, FinGenAbGroupElem}((g, h) => preimage(mU, s(mGp(_m(g)), mGp(_m(h)))) for g = domain(_m) for h = domain(_m)))
      end
      @assert Oscar.GrpCoh.istwo_cocycle(can)

      cn = preimage(mq, preimage(c[2], can))[1]
    end

    if B !== nothing
      B.lp[P] = (mL, mU, mGp, (emb, _m), C, (c, mq), cn)
    end
  end


  D = [GrpCoh.CoChain{2, PermGroupElem, FinGenAbGroupElem}(C, 
    Dict((g, h) => preimage(mU, mL(x(emb(g), emb(h)).data)) 
       for g = domain(emb) 
         for h = domain(emb))) for x = CC]

  if !isone(cn)
    cn = invmod(cn, order(q))
  end
  vals = [preimage(mq, cn * preimage(c[2], x)) for x = D]
  return [Hecke.QmodnZ()(iszero(x) ? 0 : x[1]//order(q)) for x = vals]
end

"""
    RelativeBrauerGroup

For `K/k` a number field, where `k` has to be Antic or QQField and `K` Antic
return a container for the relative Brauer group parametrizing central
simple algebras with center `k` that are split by `K` (thus can be realized
as a 2-cochain with values in `K`)
"""
@attributes mutable struct RelativeBrauerGroup
  K::AbsSimpleNumField
  k::Union{AbsSimpleNumField, QQField}
  mG::Map #PermGroup -> Aut(K/k)
  lp::Dict{AbsNumFieldOrderIdeal, Any} #for each ideal: (in K)
                             # completion map
                             # mult. group map
                             # aut map
                             # decomposition group map
                             # gmodule
                             # cohomology_group, 2
  li::Dict{<:Hecke.NumFieldEmb, Any} # for each infinite place in K
                               # decomposition group map
  map::Map         # Brauer -> CoCycle
  mkK::Map         #embedding k -> K
  function RelativeBrauerGroup(mkK::Map{<:Union{QQField, AbsSimpleNumField}, AbsSimpleNumField})
    B = RelativeBrauerGroup(codomain(mkK), domain(mkK))
    B.mkK = mkK
    return B
  end
  function RelativeBrauerGroup(K::AbsSimpleNumField, k)
    B = new()
    B.K = K
    B.k = k
    B.lp = Dict{NumFieldOrderIdeal, Any}() #ideals in k
    B.li = Dict{Hecke.NumFieldEmb, Any}()
    return B
  end
end

function extra_name(B::RelativeBrauerGroup)
  sK = get_name(B.K)
  sk = get_name(B.k)
  if sK !== nothing && sk !== nothing
    return "Br($sK, $sk)"
  end
  return nothing
end

function Base.show(io::IO, B::RelativeBrauerGroup)
  @show_name(io, B)
  @show_special(io, B)
  io = pretty(io)
  print(io, "Relative Brauer group for ", Lowercase(), B.K, " over ", Lowercase(), B.k)
end

"""
    RelativeBrauerGroupElem

Elements of the Brauer group can be in 2 forms:
 - via their local invariants: for places in `k` associate
   the invariant as an element in QQ/ZZ.
 - as a 2-cochain with values in `K`
"""
mutable struct RelativeBrauerGroupElem
  parent :: RelativeBrauerGroup
  data :: Dict{Union{NumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}
  chain:: Oscar.GrpCoh.CoChain{2, PermGroupElem, Oscar.GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}} # values in K
  function RelativeBrauerGroupElem(B::RelativeBrauerGroup, d::Dict)
    r = new()
    r.parent = B
    r.data = d
    return r
  end
end


#TODO:
# - maps between Br(k) and Br(K) for k -> K
# - tensor with field extensions?
# - map into the relative/ lift into the non-relative?
# - convert to cyclic case (easier cocycle) (Thomas Preu?)
#   Cent. Eur. J. Math. • 11(12) • 2013 • 2138-2149
#   DOI: 10.2478/s11533-013-0319-4
#   sum(c(al^j, al) for j=0:2) G = <al>, sum os over 0:order-1
#   should be the generator as a cyclic algebra.
#   as G is abelian (cyclic) I think right and left cocycles are the same
"""
  The Brauer group of k. The map translates between 2-cycles in extensions
  of k and the local data used to represent elements.
"""
@attributes mutable struct BrauerGroup
  k::AbsSimpleNumField
  map::Map
  function BrauerGroup(k::AbsSimpleNumField)
    B = new()
    B.k = k
    return B
  end
end

function extra_name(B::BrauerGroup)
  sK = get_name(B.k)
  if sK !== nothing 
    return "Br($sK)"
  end
  return nothing
end


function Base.show(io::IO, B::BrauerGroup)
  @show_name(io, B)
  @show_special(io, B)
  io = pretty(io)
  print(io, "Brauer group for ", Lowercase(), B.k)
end

mutable struct BrauerGroupElem
  parent :: BrauerGroup
  data :: Dict{Union{NumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}
  chain:: Oscar.GrpCoh.CoChain{2, PermGroupElem, Oscar.GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}} # values in K
  function BrauerGroupElem(B::BrauerGroup, d::Dict)
    r = new()
    r.parent = B
    r.data = d
    return r
  end
end

Oscar.elem_type(::Type{RelativeBrauerGroup}) = RelativeBrauerGroupElem
Oscar.parent_type(::Type{RelativeBrauerGroupElem}) = RelativeBrauerGroup
Oscar.parent(a::RelativeBrauerGroupElem) = a.parent

Oscar.elem_type(::Type{BrauerGroup}) = BrauerGroupElem
Oscar.parent_type(::Type{BrauerGroupElem}) = BrauerGroup
Oscar.parent(a::BrauerGroupElem) = a.parent

(B::RelativeBrauerGroup)(d::Dict{Union{AbsNumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}) = RelativeBrauerGroupElem(B, d)
(B::BrauerGroup)(d::Dict{Union{AbsNumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}) = BrauerGroupElem(B, d)

function Base.:+(a::T, b::T) where T <: Union{BrauerGroupElem, RelativeBrauerGroupElem}
  @assert parent(a) == parent(b)
  d = copy(a.data)
  for (k,v) = b.data
    if haskey(d, k)
      d[k] += v
    else
      d[k] = v
    end
  end
  return T(parent(a), d)
end

function Base.:-(a::T, b::T) where T <: Union{BrauerGroupElem, RelativeBrauerGroupElem}
  @assert parent(a) == parent(b)
  d = copy(a.data)
  for (k,v) = b.data
    if haskey(d, k)
      d[k] -= v
    else
      d[k] = -v
    end
  end
  return T(parent(a), d)
end

function Base.:*(a::Union{Integer, ZZRingElem}, b::T)  where T <: Union{BrauerGroupElem, RelativeBrauerGroupElem}
  d = copy(b.data)
  for (k,v) = d
    d[k] = a*v
  end
  return T(parent(b), d)
end

function Base.show(io::IO, a::RelativeBrauerGroupElem)
  #online
  print(io, "Element of relative Brauer group of $(parent(a).k)")
end

function Base.show(io::IO, m::MIME"text/plain", a::RelativeBrauerGroupElem)
  io = pretty(io)
  print(io, "Element of relative Brauer group of ", Lowercase(), parent(a).k)
  io = terse(io)
  io = IOContext(io, :compact => true) # FIXME: for now also enable compact printing
  print(io, Indent())
  data = sort(collect(a.data); by =(x -> first(x) isa AbsSimpleNumFieldEmbedding ? Inf : minimum(first(x))))
  for (p,v) in data
    print(io, "\n", p, " -> ", v)
  end
  print(io, Dedent())
end

function local_invariants(B::RelativeBrauerGroup, CC::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}})
  return B(CC).data
end

function Oscar.order(b::RelativeBrauerGroupElem)
  return lcm([order(v) for v = values(b.data)]...)
end

function (B::RelativeBrauerGroup)(d::Dict{<:Any, Hecke.QmodnZElem}; check::Bool = true)
  d = Dict{Union{NumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}((k,v) for (k,v) = d)
  if check
    for (k,v) = d  
      if isa(k, NumFieldOrderIdeal)
        P = B.mkK(k)
        lP = factor(P) #K/k should be normal or Galois Coho does not work...
                       #so any prime will do
        P = first(keys(lP))               
        loc_deg = divexact(ramification_index(P)*inertia_degree(P), ramification_index(k)*inertia_degree(k))
        @req loc_deg % order(v) == 0 "wrong index at $k - will not be split by larger field"
      else
        #need places above ...
      end
    end
  end
  return RelativeBrauerGroupElem(B, d)
end


function (B::BrauerGroup)(d::Dict{<:Any, Hecke.QmodnZElem})
  d = Dict{Union{NumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}((k,v) for (k,v) = d)
  return BrauerGroupElem(B, d)
end

"""
Given a 2-cochain with values in `K`, represent the algebra by they local
invariants in the relative Brauer group 
    `Br(K/k) = H^2(Aut(K/k), K^*)`
"""
function (B::RelativeBrauerGroup)(CC::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}})
  k = CC.C.M.data
  @assert k == B.K
  zk = maximal_order(k)
  S = [discriminant(zk)]
  mG = B.mG
  for v = values(CC.d)
    n = norm(v.data)
    if !isone(denominator(n))
      push!(S, denominator(n))
    end
    if !isone(numerator(n))
      push!(S, numerator(n))
    end
  end
  S = coprime_base(S)
  T = []
  for s = S
    isone(s) && continue
    lp = collect(keys(factor(s*zk)))
    append!(T, lp)
  end
  #=
    G = the group involved in CC
    Then we're talking about Br(k/k^G) and need invariants only
    for places over k^G (the field fixed by G)
    Hence, we need to sort into G-orbits and take one/orbit only

    Same for the infinite places
  =#
  O = []
  G = domain(mG)
  while length(T) > 0
    I = pop!(T)
    push!(O, I)
    for g = G
      gI = mG(g)(I)
      p = findfirst(isequal(gI), T)
      if p !== nothing
        deleteat!(T, p)
      end
    end
  end
  d = Dict{Any, Hecke.QmodnZElem}()
  for p = O
    d[minimum(B.mkK, p)] = local_index(CC, p, mG, B = B)
  end

  T = copy(complex_embeddings(k))
  O = []
  while length(T) > 0
    I = pop!(T)
    push!(O, I)
    for g = G
      gI = compose(mG(g), I)
      p = findfirst(isequal(gI), T)
      if p !== nothing
        deleteat!(T, p)
      end
    end
  end

  for p = O
    d[restrict(p, B.mkK)] = local_index(CC, p, mG)
  end
  b = RelativeBrauerGroupElem(B, d)
  b.chain = CC
  return b
end


"""
    relative_brauer_group(K::AbsSimpleNumField, k)

For `k` a subfield or `QQ`, create the relative Brauer group as
an infinite direct sum of the local Brauer groups.

The second return value is a map translating between the local data
and explicit 2-cochains.

```jldoctest
julia> G = SL(2,5)
SL(2,5)

julia> T = character_table(G);

julia> R = gmodule(T[9])
G-module for G acting on vector space of dimension 6 over abelian closure of QQ

julia> S = gmodule(CyclotomicField, R)
G-module for G acting on vector space of dimension 6 over cyclotomic field of order 5

julia> B, mB = relative_brauer_group(base_ring(S), character_field(S));

julia> B
Relative Brauer group for cyclotomic field of order 5 over number field of degree 1 over QQ

julia> b = B(S)
Element of relative Brauer group of number field of degree 1 over QQ
  <2> -> 1//2 + Z
  <5> -> 0 + Z
  Real embedding of number field -> 1//2 + Z
```
"""
function relative_brauer_group(K::AbsSimpleNumField, k::Union{QQField, AbsSimpleNumField} = QQ)
  G, mG = automorphism_group(PermGroup, K)
  if k != QQ 
    fl, mp = is_subfield(k, K)
    p = mp(gen(k))
    @assert fl
    s, ms = sub(G, [x for x = G if mG(x)(p) == p])
    G = s
    mG = ms*mG
  else
    mp = MapFromFunc(QQ, K, K, QQ)
  end
  B = RelativeBrauerGroup(mp)
  B.mG = mG

  function elem_to_cocycle(b::RelativeBrauerGroupElem)
    B = parent(b)
    K = B.K
    lp = Set([minimum(p) for p = keys(b.data) if isa(p, NumFieldOrderIdeal)])
    lp = union!(lp, Set(ramified_primes(maximal_order(K))))
    ZK = maximal_order(K)
    C, mC = class_group(ZK)
    lP = vcat([collect(keys(factor(p*ZK))) for p = lp]...)
    q, mq = quo(C, [preimage(mC, x) for x = lP])
    p = 2
    while order(q) > 1
      while p in lp
        p = next_prime(p)
      end
      P = collect(keys(factor(p*ZK)))
      cP = [preimage(mq, preimage(mC, x)) for x = P]
      if all(iszero, cP)
        continue
      end
      push!(lp, p)
      append!(lP, P)
      q, _mq = quo(q, cP)
      mq = _mq * mq
    end

    S, mS = sunit_group(lP)

    MC = Oscar.GrpCoh.MultGrp(K)
    mMC = MapFromFunc(K, MC, MC, y->y.data)

    mG = B.mG
    G = domain(mG)
    mu = gmodule(MC, G, [hom(MC, MC, mG(g)) for g = gens(G)])
    M = gmodule(G, mS, mG)
    @assert Oscar.GrpCoh.is_consistent(M)
    z = cohomology_group(M, 2);

    for x = gens(z[1])
      @assert Oscar.GrpCoh.istwo_cocycle(z[2](x))
    end
    q, mq = snf(z[1])


    p = []
    for x = lp
      push!(p, prime_decomposition(ZK, x)[1][1])
    end

    zz = [map_entries(mS*mMC, z[2](image(mq, x)), parent = mu) for x = gens(q)]
    k = collect(keys(zz[1].d))

    @assert all(Oscar.GrpCoh.istwo_cocycle, zz)

    em = complex_embeddings(B.k)
    EM = [extend(x, mp)[1] for x = em]
    lb = RelativeBrauerGroupElem[]
    for x = zz
      d = Dict{Union{NumFieldOrderIdeal, Hecke.NumFieldEmb}, Hecke.QmodnZElem}()
      for P = lP
        d[minimum(mp, P)] = local_index(x, P, mG, B = B)
      end
      for i = length(em)
        d[em[i]] = local_index(x, EM[i], mG)
      end
      push!(lb, RelativeBrauerGroupElem(B, d))
    end
  
    fl, x = can_solve_with_solution(lb, b)
    @assert fl
 
    return map_entries(mS*mMC, z[2](image(mq, q(x.coeff))), parent = mu)
  end

  function cocycle_to_elem(C::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}})
    b = B(C)
    return b
  end

  B.map =  MapFromFunc(B, Oscar.GrpCoh.AllCoChains{2, PermGroupElem, Oscar.GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}(),
                       elem_to_cocycle, 
                       cocycle_to_elem)
  return B, B.map
end

"""
    brauer_group(K::AbsSimpleNumField)

The Brauer group as
an infinite direct sum of the local Brauer groups.

The second return value is a map translating between the local data
and explicit 2-cochains.

# EXAMPLES

```jldoctest
julia> k = rationals_as_number_field()[1];

julia> lp = collect(keys(factor(30*maximal_order(k))));

julia> qz = Hecke.QmodnZ();

julia> B, mB = brauer_group(k);

julia> b = B(Dict(lp[1]=>qz(1//3), lp[2]=>qz(2//3)));

julia> c = mB(b);

julia> C = structure_constant_algebra(c);

julia> Hecke.local_schur_indices(C);

julia> c = mB(b+b);

```

"""
function brauer_group(K::AbsSimpleNumField)
  B = BrauerGroup(K)

  function elem_to_cocycle(b::BrauerGroupElem)
    L = number_field(grunwald_wang(b)) # a splitting field
        # TODO: cache splitting fields or relative_brauer_groups?
    L, _ = simple_extension(L)
    La, mpLa_L = absolute_simple_field(L)
    _ = lll(maximal_order(La))
    C, mC = relative_brauer_group(La, K)
    return mC(C(b.data))
  end

  function cocycle_to_elem(C::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}})
    L = C.C.M.data
    D, mD = relative_brauer_group(L, K)
    b = D(C)
    return B(b.data)
  end

  B.map =  MapFromFunc(B, Oscar.GrpCoh.AllCoChains{2, PermGroupElem, Oscar.GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}(),
                       elem_to_cocycle, 
                       cocycle_to_elem)
  return B, B.map
end

function (B::RelativeBrauerGroup)(M::GModule{<:Group, AbstractAlgebra.Generic.FreeModule{AbsSimpleNumFieldElem}})
  @assert B.K == base_ring(M)
  c = factor_set(M, B.mG)
  return preimage(B.map, c)
end

"""
Return the local component at `p`.
"""
function (a::RelativeBrauerGroupElem)(p::Union{NumFieldOrderIdeal, Hecke.NumFieldEmb})
  if haskey(a.data, p)
    return a.data[p]
  else
    return Hecke.QmodnZ()(1)
  end
end

#write (or try to write) `b` as a ZZ-linear combination of the elements in `A`
function Oscar.can_solve_with_solution(A::Vector{RelativeBrauerGroupElem}, b::RelativeBrauerGroupElem)
  @assert all(x->parent(x) == parent(b), A)
  lp = Set(collect(keys(b.data)))
  for a = A
    for p = keys(a.data)
      push!(lp, p)
    end
  end
  lp = collect(lp)
  push!(A, b)

  li = [x for x = lp if isa(x, Hecke.NumFieldEmb)]
  lp = [x for x = lp if isa(x, NumFieldOrderIdeal)]

  d = [lcm([order(x(k)) for x = A]...) for k = vcat(lp, li)]
  F = free_abelian_group(length(A)-1)         
  G = abelian_group(d)
  function q_to_a(x, d)
    return divexact(d, denominator(x.elt)) * numerator(x.elt)
  end
  
  a = [G(vcat([q_to_a(x(lp[i]), d[i]) for i = 1:length(lp)],
               [x(li[i]) == 0 ? 0 : 1 for i = 1:length(li)])) for x = A] 
 
  h = hom(F, G, a[1:end-1])
  pop!(A)
  return has_preimage_with_preimage(h, a[end])
end  

"""
    structure_constant_algebra(a::RelativeBrauerGroupElem)

Compute an explicit matrix algebra with the local invariants given
by the element
"""
function Hecke.structure_constant_algebra(a::RelativeBrauerGroupElem)
  return A = structure_constant_algebra(parent(a).map(a), parent(a).mG, parent(a).mkK)
end

function Hecke.grunwald_wang(b::Union{BrauerGroupElem, RelativeBrauerGroupElem})
  d1 = Dict((x=>Int(order(y))) for (x,y) = b.data if isa(x, NumFieldOrderIdeal))
  d2 = Dict((x=>Int(order(y))) for (x,y) = b.data if !isa(x, NumFieldOrderIdeal))

  if length(d2) == 0
    return Hecke.grunwald_wang(d1)
  else
    return Hecke.grunwald_wang(d1, d2)
  end
end

function local_index(C::GModule{<:Oscar.GAPGroup, Generic.FreeModule{AbsSimpleNumFieldElem}}, p::Union{Integer, ZZRingElem})
  K = base_ring(C)
  k, mkK = Oscar.GModuleFromGap._character_field(C)
  A, mA = automorphism_group(PermGroup, K)
  gl = mkK(gen(k))
  s, ms = sub(A, [a for a = A if mA(a)(gl) == gl])
  lp = factor(p*maximal_order(K))
  r = []
  mp = ms*mA
  c = Oscar.GModuleFromGap._two_cocycle(mp, C, two_cycle = true)
  for P = keys(lp)
    push!(r, P => local_index(c, P, mp))
  end
  return Dict(r)
end

"""
    structure_constant_algebra(CC::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}, mG::Map = automorphism_group(PermGroup, CC.C.M.data)[2], mkK::Union{<:Map, Nothing} = nothing)

Return the cross-product algebra defined by the factor system given
as a 2-cochain.
"""
function Hecke.structure_constant_algebra(CC::GrpCoh.CoChain{2, PermGroupElem, GrpCoh.MultGrpElem{AbsSimpleNumFieldElem}}, mG::Map = automorphism_group(PermGroup, CC.C.M.data)[2], mkK::Union{<:Map, Nothing} = nothing)

  k = CC.C.M.data
  G = domain(mG)
  all_G = collect(G)
  n = order(Int, G)
  if mkK !== nothing && domain(mkK) != QQ
    k, mp = relative_simple_extension(mkK)
  else
    mp = id_hom(k)
  end
  M = zeros(base_field(k), n*n, n*n, n*n)
  #basis is prod. basis of basis(k) and e_sigma sigma in G
  # e_sigma * e_tau is via cocycle
  # bas * bas is direct
  # e_sigma * bas = bas^sigma * e_sigma
  # bas * e_sigma = bas * e_sigma
  B = [(b, s) for b = basis(k) for s = all_G]
  one = elem_type(k)[ isone(b[1]) && isone(b[2]) ? Base.one(k) : zero(k) for b = B]
  function mul(a::Tuple, b::Tuple)
    # a*A * b * B
    # -> a* b^A * A*B
    # -> a*b^A * (sigma(A, B)*AB)
    c = a[1] * preimage(mp, mG(a[2])(mp(b[1]))) * preimage(mp, CC(a[2], b[2]).data)
    C = a[2]*b[2]
    z = zeros(base_field(k), n*n)
    bk = basis(k)
    for i=1:length(bk)
      p = findfirst(isequal((bk[i], C)), B)
      z[p] = coeff(c, i-1)
    end
    return z
  end
  for i = 1:n*n
    for j = 1:n*n
      M[i, j, : ] = mul(B[i], B[j])
    end
  end
  return Hecke.structure_constant_algebra(base_field(k), M; check = false, one = one)
end

function serre(A::IdeleParent, P::AbsNumFieldOrderIdeal)
  C = A.data[1]
  Kp, mKp, mGp, mUp, pro, inj = completion(A, P)
  mp = decomposition_group(A.k, mKp, A.mG, mGp)
  qr = restrict(C, mp)
  s = Hecke.Hecke.local_fundamental_class_serre(Kp, absolute_base_field(Kp))
#  Oscar.GModuleFromGap.istwo_cocycle(Dict( (g, h) => s(mGp(g), mGp(h)) for g = domain(mGp) for h = domain(mGp)), mGp)

  z = gmodule(domain(mGp), [hom(domain(mUp), domain(mUp), [preimage(mUp, mGp(g)(mUp(u))) for u = gens(domain(mUp))]) for g = gens(domain(mGp))])

  c = CoChain{2, PermGroupElem, FinGenAbGroupElem}(z, Dict{NTuple{2, PermGroupElem}, FinGenAbGroupElem}((g, h) => preimage(mUp, s(mGp(g), mGp(h))) for g = domain(mGp) for h = domain(mGp)))

  @assert Oscar.GrpCoh.istwo_cocycle(c)

  return c
end

function serre(A::IdeleParent, P::Union{Integer, ZZRingElem})
  C = A.data[1]
  t = findfirst(isequal(ZZ(P)), [minimum(x) for x = A.S])
  Inj = canonical_injection(A.M, t+1)
  Pro = canonical_projection(A.M, t+1)

  inj = canonical_injection(domain(Inj), 1)
  pro = canonical_projection(domain(Inj), 1)

  Kp, mKp, mGp, mUp, _, _ = completion(A, A.S[t])
  @assert domain(inj) == domain(mUp) 
  mp = decomposition_group(A.k, mKp, A.mG, mGp)
 
  tt = serre(A, A.S[t])
  @assert tt.C.G == domain(mGp)

  I = domain(Inj)    
  zz = gmodule(C.G, [Inj * action(A.data[2], g) * Pro for g = gens(C.G)])
  #the induced module
  mu = cohomology_group(zz, 2)
  q, mq = snf(mu[1])
  g = mu[2](mq(q[1]))
  #g is the (non-canonical) generator of H^2 of the induced module
  hg = map_entries(Inj*A.mq, g, parent = C)
  #hg is the (non-canonical) generator of H^2 of the induced module
  #in the global idele class cohomology
  #             
  #Z[G_p] K_p   <-> Ind_G_p^G K_p = Z[G] prod K_p over all conjugate primes
  #              -> Z[G] C the idele class group
  #the image should be the restriction I think
  gg = map_entries(pro, g, parent = tt.C)
  #gg is the non-canomical generator in Z[G_p] K_p
  gg = Oscar.GrpCoh.CoChain{2, PermGroupElem, FinGenAbGroupElem}(tt.C, Dict( (g, h) => gg(mp(g), mp(h)) for g = tt.C.G for h = tt.C.G))

  nu = cohomology_group(tt.C, 2)
  ga = preimage(nu[2], gg)
  ta = preimage(nu[2], tt)
  #here we can compare: tt is canonical, gg is non-canonical
  #we return the multiple that makes gg canonical.
  #and the non-canonical gg in Z[G] C, which should be the
  #restriction, hence a multiple of the global generator
  return findfirst(x->x*ga == ta, 1:order(tt.C.G)), hg
  #so i*hg should restrict to the local fund class...
end

"""
    global_fundamental_class(A::IdeleParent)

Tries to find the canonical generator of `H^2(C)` the  2nd
cohomology group of the idele class group.

Currently only works if this can be inferred from the local data in
the field, ie. if the `lcm` of the degrees of the completions is the
full field degree.
"""    
function global_fundamental_class(A::IdeleParent)
  C = A.data[1]
  d = lcm([ramification_index(P) * inertia_degree(P) for P = A.S])
  G = C.G
  if d != order(G)
    error("sorry - no can do(yet)")
  end

  z = cohomology_group(C, 2)

  q, mq = snf(z[1])
  @assert ngens(q) == 1
  g = z[2](mq(gen(q, 1))) # to get a 2-CoCycle
  #g is the non-canonical generator
  @assert Oscar.GrpCoh.istwo_cocycle(g)
  n = order(q)
  g = mq(gen(q, 1))

  scale = []

  for P = A.S
    s = serre(A, minimum(P))
    d = inertia_degree(P) * ramification_index(P) #local degree
    #s[1] = j s.th. canonical at P = j * actual gen at P
    #s[2] should be the restriction to the local component
    #     hence a factor times the global gen
    # we need:
    #  mu sth 
    #   o_p (mu H^2(C)[1]) = j_P* s_P[2] for all P
    #  where o_p = deg(K)/order(s_p[2) = deg(K)/deg(K_p)

    k = findfirst(k->k * preimage(z[2], s[2]) == divexact(n, d)*g, 1:d)
    #so (n/d) * k * g = res(g, G_p) = s[2]
    #   s[1] * k * (n/d) * g = canonical at P
    push!(scale, ((s[1]*k) % d, d))
  end
  #put together..
  #want x s.th. (n/d[2]) * x = d[1] mod d[2]
  #cannot use CRT as the d[2] are not necc. coprime
  a = abelian_group([x[2] for x = scale])
  b = abelian_group([n])
  h = hom(b, a, [sum(a[i] * divexact(n, scale[i][2]) for i=1:length(scale))])
  p = preimage(h, sum(a[i] * scale[i][1] for i=1:length(scale)))
  @assert gcd(p[1], n) == 1  
  #so p[1]* g is canonical!!!
  return p[1]
end

function Oscar.cohomology_group(A::IdeleParent, i::Int)
  return Oscar.cohomology_group(A.data[1], i)
end

"""
For a C a G-module with operation on M and an element o of M,
compute the G-orbit (not the Z[G] orbit!) of o.
"""
function Oscar.orbit(C::GModule, o)
  @assert parent(o) == C.M
  return orbit(C.G, (x,y) -> action(C, y, x), o)
end

#TODO: reduce torsion: the part coprime to |G| can go...
"""
    shrink(C::GModule{PermGroup, FinGenAbGroup}, attempts::Int = 10)

Tries to find cohomologically trivial submodules to factor out.
Return a cohomologically equivalent module with fewer generators and
the quotient map.
"""
function shrink(C::GModule{PermGroup, FinGenAbGroup}, attempts::Int = 10)
  mq = hom(C.M, C.M, gens(C.M))
  q = C
  first = true
  while true
    prog = false
    for i=1:attempts
      o = Oscar.orbit(q, rand(gens(q.M)))
      if length(o) == order(group(q))
        s, ms = sub(q.M, collect(o), false)
        if torsion_free_rank(s) == length(o)
          q, _mq = quo(q, ms, false)
          if first
            mq = _mq
            first = false
          else
            mq = mq*_mq
          end
          q, _mq = simplify(q)
          mq = mq*inv(_mq)
          prog = true
          break
        end
      end
    end
    prog || return q, mq
  end
end

function Oscar.simplify(C::GModule{PermGroup, FinGenAbGroup})
  s, ms = snf(C.M)
  S = GModule(s, C.G, [FinGenAbGroupHom(ms*x*pseudo_inv(ms)) for x = C.ac])
  if isdefined(C, :iac)
    S.iac = [FinGenAbGroupHom(ms*x*pseudo_inv(ms)) for x = C.iac]
  end
  return S, ms
end

#deliberately at the end as it messes up the syntax-highlighting for me
#confusion about " and ()...
function Base.show(io::IO, I::IdeleParent)
  io = pretty(io)
  #print(io, "Relative Brauer group for ", Lowercase(), B.K, " over ", Lowercase(), B.k)
  println(io, "Idele group of")
  print(io, Indent())
  println(io, Lowercase(), I.k)
  plcs = sort(collect(Set(minimum(x) for x = I.S)))
  print(io, "using prime ideals over [", join(plcs, ", "), "] as places")
  print(io, Dedent())
end

end # module GrpCoh

using .GaloisCohomology_Mod
export is_coboundary, 
       idele_class_gmodule,
       relative_brauer_group,
       units_mod_ideal,
       brauer_group

