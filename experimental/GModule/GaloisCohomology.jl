module GaloisCohomology_Mod
using Oscar
import Oscar: GrpCoh
import Oscar.GrpCoh: CoChain, MultGrpElem, MultGrp, GModule, is_consistent, 
                     Group
import Base: parent
import Oscar: direct_sum
export is_coboundary, idel_class_gmodule


Oscar.elem_type(::Type{Hecke.NfMorSet{T}}) where {T <: Hecke.LocalField} = Hecke.LocalFieldMor{T, T}
parent(f::Hecke.LocalFieldMor) = Hecke.NfMorSet(domain(f))

function Oscar.automorphism_group(::Type{PermGroup}, k)
  G, mG = automorphism_group(k)
  H = symmetric_group(degree(k))
  gens(G) #to make sure gens are actually there...
  H = sub(H, [H(G.mult_table[:, i]) for i=G.gens])[1]

  function HtoG(p::PermGroupElem)
    m = [i^p for i=1:degree(k)]
    i = Base.findfirst(x->G.mult_table[:, x] == m, 1:degree(k))
    return mG(GrpGenElem(G, i))
  end

  function GtoH(a::NfToNfMor)
    g = preimage(mG, a)
    return H(G.mult_table[:, g.i])
  end

  return H, MapFromFunc(HtoG, GtoH, H, codomain(mG))
end

function Oscar.automorphism_group(::Type{PermGroup}, K, k)
  G, mG = automorphism_group(K, k)
  H = symmetric_group(length(G))
  gens(G) #to make sure gens are actually there...
  H = sub(H, [H(G.mult_table[:, i]) for i=G.gens])[1]

  function HtoG(p::PermGroupElem)
    m = [i^p for i=1:length(G)]
    i = Base.findfirst(x->G.mult_table[:, x] == m, 1:length(G))
    return mG(GrpGenElem(G, i))
  end

  function GtoH(a::NfToNfMor)
    g = preimage(mG, a)
    return H(G.mult_table[:, g.i])
  end

  return H, MapFromFunc(HtoG, GtoH, H, codomain(mG))
end


"""
The natural `ZZ[H]` module where `H`, a subgroup of the
  automorphism group acts on the ray class group.
"""
function Oscar.gmodule(H::PermGroup, mR::MapRayClassGrp, mG = automorphism_group(PermGroup, k)[2])
  k = nf(order(codomain(mR)))
  G = domain(mG)

  ac = Hecke.induce_action(mR, [image(mG, G(g)) for g = gens(H)])
  return GModule(H, ac)
end

"""
The natural `ZZ[G]` module where `G`, the
  automorphism group, acts on the ideal group defining the class field.
"""
function Oscar.gmodule(R::ClassField, mG = automorphism_group(PermGroup, k)[2])
  k = base_field(R)
  G = domain(mG)
  mR = R.rayclassgroupmap
  mq = R.quotientmap

  ac = Hecke.induce_action(mR, [image(mG, g) for g = gens(G)], mq)
  return GModule(G, ac)
end

"""
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
function _gmodule(k::AnticNumberField, H::PermGroup, mu::Map{GrpAbFinGen, FacElemMon{AnticNumberField}}, mG = automorphism_group(PermGroup, k)[2])
  u = domain(mu)
  U = [mu(g) for g = gens(u)]
  G = domain(mG)
  ac = [hom(u, u, [preimage(mu, mG(G(g))(x)) for x = U]) for g = gens(H)]
  return gmodule(H, ac)
end

function Oscar.gmodule(H::PermGroup, mu::Map{GrpAbFinGen, FacElemMon{AnticNumberField}}, mG = automorphism_group(PermGroup, base_ring(codomain(mu)))[2])
  return _gmodule(base_ring(codomain(mu)), H, mu, mG)
end

function Oscar.gmodule(H::PermGroup, mu::Hecke.MapUnitGrp{NfOrd}, mG = automorphism_group(PermGroup, k)[2])
  #TODO: preimage for sunits can fail (inf. loop) if
  # (experimentally) the ideals in S are not coprime or include 1
  # or if the s-unit is not in the image (eg. action and not closed set S)
  u = domain(mu)
  U = [mu(g) for g = gens(u)]
  zk = codomain(mu)
  k = nf(zk)
  G = domain(mG)
  ac = [hom(u, u, [preimage(mu, zk(mG(G(g))(k(x)))) for x = U]) for g = gens(H)]
  return gmodule(H, ac)
end

function Oscar.gmodule(H::PermGroup, mu::Map{GrpAbFinGen, AnticNumberField})
  return _gmodule(codomain(mu), H, mu)
end

function is_coboundary(c::CoChain{2,PermGroupElem,MultGrpElem{nf_elem}})
  @vprint :GaloisCohomology 1 "testing if 2-chain is a boundary\n"

  zk = maximal_order(parent(first(values(c.d)).data))
  @vprint :GaloisCohomology 2 ".. gathering primes in the support ..\n"
  cp = coprime_base(vcat([numerator(norm(x.data*denominator(x.data))) for x = values(c.d)],
                         map(x->denominator(x.data), values(c.d))))
  s = Set(reduce(vcat, [collect(keys(factor(x).fac)) for x = cp], init = [1]))
  while 1 in s
    pop!(s, 1)
  end

  @vprint :GaloisCohomology 2 ".. class group ..\n"
  Cl, mCl = class_group(zk)
  if length(s) == 0
    S = Set{NfOrdIdl}()
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
  cc = CoChain{2,PermGroupElem,GrpAbFinGenElem}(C, Dict((h, preimage(mu, FacElem(v.data))) for (h,v) = c.d))
  @vprint :GaloisCohomology 2 ".. test for boundary ..\n"
  fl, d = z(cc)
  if !fl
    @vprint :GaloisCohomology 2 ".. no boundary\n"
    return fl, d
  end
  @vprint :GaloisCohomology 2 ".. explicit boundary\n"
  MK = MultGrp(number_field(zk))
  return fl, CoChain{1,PermGroupElem,elem_type(MK)}(c.C, Dict((h, MK(evaluate(mu(v)))) for (h,v) = d.d))
end

function isunramified(p::NfOrdIdl)
  return ramification_index(p) == 1
end


"""
For a completion C of a number field K, implicitly given as the map
    mK:  K -> C
And the automorphism group G of K given via
    mG:  G -> aut(K)
and the automorphism group Gp of Kp, given via
    mGp: Gp -> Aut(Kp)
Find the embedding of Gp -> G, realizing the local automorphism group
as a subgroup of the global one.
"""
function Oscar.decomposition_group(K::AnticNumberField, mK::Map, mG::Map = automorphism_group(K)[2], mGp::Map = automorphism_group(codomain(mK), prime_field(codomain(mK))))
  Kp = codomain(mK)
  @assert domain(mK) == K

  Gp = domain(mGp)
  G = domain(mG)

  im = elem_type(G)[]
  elG = [g for g = G]
  imK = [mK(mG(g)(gen(K))) for g = elG]
  for s = gens(Gp)
    h = mGp(s)(mK(gen(K)))
    z = findall(isequal(h), imK)
    if length(z) == 0
      z = argmax([valuation(h-x) for x = imK], dims = 1)
    end
    @assert length(z) == 1
    push!(im, elG[z[1]])
  end
  return hom(Gp, G, im)
end

"""
  For a real or complex embedding `emb`, find the unique automorphism
  that acts on this embedding as complex conjugation.
"""
function Oscar.decomposition_group(K::AnticNumberField, emb::Hecke.NumFieldEmb, mG::Map = automorphism_group(K)[2])
  G = domain(mG)
  if is_real(emb)
    return sub(G, [one(G)])[2]
  end
  g = gen(K)
  lG = [g for g  = G]
  l = findall(x->overlaps(conj(emb(g)), emb(mG(x)(g))), lG)
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
For a local field extension K/k, return a gmodule for the multiplicative
group of K as a Gal(K/k) module.

Returns: 
 - the gmodule
 - the map from G = Gal(K/k) -> Set of actual automorphisms
 - the map from the module into K
"""
function Oscar.gmodule(K::Hecke.LocalField, k::Union{Hecke.LocalField, FlintPadicField, FlintQadicField} = base_field(K); Sylow::Int = 0, full::Bool = false)

  #if K/k is unramified, then the units are cohomological trivial,
  #   so Z (with trivial action) is correct for the gmodule
  #if K/k is tame, then the 1-units are cohomologycal trivial, hence
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
      MapFromFunc(x->pi^x[1], y->Int(e*valuation(y))*A[1], A, K)
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
      MapFromFunc(x->pi^x[1] * gk^x[2],
        function(y)
          v = Int(e*valuation(y))
          y *= pi^-v
          return v*A[1] + preimage(mu, mk(y))[1]*A[2]
        end, A, K)
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
  Q, mQ = quo(U, [preimage(mU, 1+prime(k)^4*x) for x = o])
  S, mS = snf(Q)
  Q = S
  mQ = mQ*inv(mS)

  if Sylow > 0
    @assert isprime(Sylow)
    G, mS = sylow_subgroup(G, Sylow)
    mG = mS*mG
  end

  @vprint :GaloisCohomology 2 " .. the module ..\n"
  hh = [hom(Q, Q, [mQ(preimage(mU, mG(i)(mU(preimage(mQ, g))))) for g = gens(Q)]) for i=gens(G)]
  Hecke.assure_has_hnf(Q)
  return gmodule(G, hh), mG, pseudo_inv(mQ)*mU
end

#=  Not used
function one_unit_cohomology(K::Hecke.LocalField, k::Union{Hecke.LocalField, FlintPadicField, FlintQadicField} = base_field(K))

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

export GModule
export action
export cohomology_group
export confluent_fp_group
export extension
export fp_group
export gmodule
export induce
export is_coboundary
export pc_group
export word

#= TODO
  for Z, Z/nZ, F_p and F_q moduln -> find Fp-presentation
  for finite Z, Z/nZ, F_p and F_q moduln -> find pc-presentation
  #done: for GrpAbFinGen          -> find Fp-presentation
  #done: for GrpAbFinGen          -> find pc-presentation
  #done: for a in H^2 find Fp-presentation
  #done: for a in H^2 find pc-presentation
  for a in H^2 find (low degree) perm group using the perm group we have?
  Magma's DistinctExtensions
  probably aut(GrpAb), ...

Sort: 
 - move the additional GrpAbFinGenMap stuff elsewhere
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
#    use Klueners/ Acciaro to map arbitrary local into idel
#    use ...               to project to ray class
# - a magic(?) function to get idel-approximations in and out?

"""
M has to be a torsion free Z module with a C_2 action by sigma.
Returns data for the decomposition into indecomposables.
They will be of type
 - Z with trivial and non-trivial action
 - Z[C_2]

Two arrays are returned:
 - generators for the 1-dim modules
 - C_2 generators for the 2-dim ones

Follows Debeerst rather closely...

(Helper for the idel-class stuff)
"""
function debeerst(M::GrpAbFinGen, sigma::Map{GrpAbFinGen, GrpAbFinGen})
  @assert domain(sigma) == codomain(sigma) == M
  @assert all(x->sigma(sigma(x)) == x, gens(M))
  @assert is_free(M) && rank(M) == ngens(M)

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
  @assert istrivial(_S) || rank(_S) == ngens(_S) 
  @assert rank(_K) == ngens(_K) 

  m = matrix(GrpAbFinGenMap(_mS * mSK * inv((_mK))))
  # elt in S * m = elt in K
  # so
  # elt in S * U^-1 U m V V^-1 = elt_in K
  # elt in S * U^-1 snf = elt_in * V
  s, U, V = snf_with_transform(m)
  if istrivial(S)
    r = 0
  else
    r = maximum(findall(x->isone(s[x,x]), 1:ngens(_S)))
  end

  mu = hom(_S, _S, inv(U))
  mv = hom(_K, _K, V)
  @assert istrivial(S) || all(i-> M(_mS(mu(gen(_S, i)))) == s[i,i] * M(_mK(mv(gen(_K, i)))), 1:ngens(S))
  b = [_mK(mv(x)) for x = gens(_K)]

  Q, mQ = quo(S, image(sigma -id_hom(M), K)[1])
  B, mB = sub(Q,  [mQ(preimage(mSK, x)) for x = b[1:r]])
  @assert order(B) == order(Q)

  phi = GrpAbFinGenMap(_mX*mX*(sigma -id_hom(M))*pseudo_inv(mS)*mQ)
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

function (G::GrpAbFinGen)(x::GrpAbFinGenElem)
  fl, m = is_subgroup(parent(x), G)
  @assert fl
  return m(x)
end

function Hecke.extend_easy(m::Hecke.CompletionMap, L::FacElemMon{AnticNumberField})
  k = base_ring(L)
  @assert k == domain(m)

  #want a map: L-> codomain(m)
  function to(a::FacElem{nf_elem})
    return prod(m(k)^v for (k,v) = a.fac)
  end
  function from(a::Hecke.LocalFieldElem)
    return FacElem(preimage(m, a))
  end
  return MapFromFunc(to, from, L, codomain(m))
end

function Hecke.extend_easy(m::Hecke.CompletionMap, mu::Map, L::FacElemMon{AnticNumberField})
  k = base_ring(L)
  @assert k == domain(m)
  @assert codomain(mu) == codomain(m)

  cache = Dict{nf_elem, GrpAbFinGenElem}()
  #want a map: L-> codomain(m) -> domain(mu)
  function to(a::FacElem{nf_elem})
    s = domain(mu)[0]
    for (k,v) = a.fac
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
  return MapFromFunc(to, from, L, domain(mu))
end


mutable struct IdelParent
  k::AnticNumberField
  mG::Map # AutGrp -> Automorohisms
  S::Vector{NfAbsOrdIdl} # for each prime number ONE ideal above
  C::Vector{Map} # the completions at S
  D::Vector{Map} # Gp -> Aut
  L::Vector{Map} # the mult. group map at C

  #for P in S the modules used actually is
  #    Ind_G_p^G L[P]
  #        = sum L[P] otimes s_i
  # (for s_i a fixed system of coset reps G//G_P)
  # L[P] otimes s_i "should be" the completion data at P^s_i - one of the other ideals
  # should be L[P] ni l -> C[P] -> k -> inv(s_i)(..) to get a proper rep in k
  # completion at P^s is C[P] but with the map twisted by s

  mU::Map #S-unit group map
  M::GrpAbFinGen  # the big module, direct product from
    # infinite gmodule x finite ones
  mq::Map # "projection" of M -> the acutal module in the end

  data

  function IdelParent()
    return new()
  end
end

"""
Following Debeerst:
  Algorithms for Tamagawa Number Conjectures. Dissertation, University of Kassel, June 2011.
or Ali, 

Find a gmodule C s.th. C is cohomology-equivalent to the cohomology
of the idel-class group.
"""
function idel_class_gmodule(k::AnticNumberField, s::Vector{Int} = Int[])
  @vprint :GaloisCohomology 1 "Ideal class group cohomology for $k\n"
  I = IdelParent()
  I.k = k
  G, mG = automorphism_group(PermGroup, k)
  I.mG = mG

  zk = maximal_order(k)

  sf = subfields(k)
  sf = [x[1] for x = sf if degree(x[1]) > 1]
  zf = map(maximal_order, sf)
  cf = map(class_group, zf)
  cf = Tuple{GrpAbFinGen, <:Map}[x for x = cf]

  @vprint :GaloisCohomology 2 " .. gathering primes ..\n"
  s = push!(Set{ZZRingElem}(s), Set{ZZRingElem}(keys(factor(discriminant(zk)).fac))...)
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
  z = MapFromFunc(x->evaluate(x), y->FacElem(y), codomain(mU), k)
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
    T = abelian_group([order(U[1]), 0])             
    ac_T = hom(T, T, [sigma(U[1])[1]*T[1], T[1]+T[2]])

    x = [preimage(mq, i) for i = x]
    y = [preimage(mq, i) for i = y]

    z, mz = sub(U, [sigma(U[1]) - U[1]])
    theta_i = [sigma(t)-t for t = x]
    inv = Int[]
    not_inv = Int[]
    for i=1:length(x)
      w = theta_i[i]
      fl, pe = haspreimage(mz, w)
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
      #should be chosen to be pos. at place, flip signs...
    end
    for i=length(not_inv)+1:length(x)
      push!(im_psi, x[i])
      #should be chosen to be pos. at place, flip signs...
    end
    for i=1:length(y)
      push!(im_psi, y[i])
      push!(im_psi, sigma(y[i]))
    end
    psi = hom(V, U, im_psi)
    @assert is_bijective(psi)
    F = abelian_group([0 for i=2:length(x)])
    Hecke.assure_has_hnf(F)
    W, pro, inj = direct_product(V, F, task = :both)
    @assert isdefined(W, :hnf)

    ac = GrpAbFinGenMap(pro[1]*psi*sigma*pseudo_inv(psi)*inj[1])+ GrpAbFinGenMap(pro[2]*hom(F, W, [inj[1](preimage(psi, x[i])) - inj[2](F[i-1]) for i=2:length(x)]))
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
  C = [gmodule(x[1], prime_field(x[1])) for x = L];
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
  I.data = F[1]

  @hassert :GaloisCohomology 1 is_consistent(F[1])

  h = iEt[2]*F[3][1]+sum(D[i][2]*F[3][i+1] for i=1:length(S));
  @vtime :GaloisCohomology 2 q, mq = quo(F[1], h)
  @hassert :GaloisCohomology 1 is_consistent(q)
  @vtime :GaloisCohomology 2 q, _mq = simplify(q)
  @vtime :GaloisCohomology 2 mq = GrpAbFinGenMap(mq * pseudo_inv(_mq))
  @hassert :GaloisCohomology 1 is_consistent(q)
  I.mq = mq
  function idel(a::GrpAbFinGenElem)
    a = preimage(mq, a) # in F
    u = F[2][1](a) #in iEt need to get to the S-Unit somehow, maybe
    v = [m(a) for m = F[2][2:end]] #in the induced GModules
    #= TODO
     - the induced stuff is equivalent to doing all completions:
       for infinite places, same as finite ones over the same prime
       Galois operates transitive, and the :induce has coset reps
       -> sort places against coset reps
       -> sort primes against coset reps
       return a dictionary where keys are the 
         places, prime ideals
       and values are
         s.th. for the infinite places, for the reals we need the sign?
         elements in the completion for the prime ideals
         (think: nf_elem as the completions do not exist? only in 
           spirit via the cosets?)

     - we need also the inverse operation...
     =#
    return u, v
  end

  return q, idel, I
end

function Oscar.components(A::GrpAbFinGen)
  return get_attribute(A, :direct_product)
end

function Oscar.completion(I::IdelParent, P::NfAbsOrdIdl)
  s = [minimum(x) for x = I.S]
  p = findfirst(isequal(minimum(P)), s)
  @assert p !== nothing

  mKp = I.C[p]
  Kp = codomain(mKp)
  mUp = I.L[p]
  mGp = I.D[p]

  inj = Hecke.canonical_injection(I.M, p+1) #units are first
  pro = Hecke.canonical_projection(I.M, p+1)


  @assert domain(inj) == codomain(pro)

  J = components(I.M)[p+1]
  if mKp.P == P #easy case
    return Kp, mKp, mGp, mUp,  pro * Hecke.canonical_projection(J, 1) ,  Hecke.canonical_injection(J, 1)*inj
  end

  prm = get_attribute(J, :induce)[2]
  mG = I.mG

  z = findall(pr -> mG(pr)(mKp.P) == P, prm)
  pr = inv(prm[z[1]])
  
  nKp = MapFromFunc(x->mKp(mG(pr)(x)), y->mG(inv(pr))(preimage(mKp, y)), I.k, Kp)

  return Kp, nKp, mGp, mUp, pro * Hecke.canonical_projection(J, z[1]), Hecke.canonical_injection(J, z[1])*inj 
end

function Oscar.map_entries(mp::Map, C::GrpCoh.CoChain{N, G, M}; parent::GModule) where {N, G, M}
  d = Dict( k=> mp(v) for (k,v) = C.d)
  return GrpCoh.CoChain{N, G, elem_type(codomain(mp))}(parent, d)
end

function serre(C::GModule, A::IdelParent, P::NfAbsOrdIdl)
  Kp, mKp, mGp, mUp, pro, inj = completion(A, P)
  mp = decomposition_group(A.k, mKp, A.mG, mGp)
  qr = restrict(C, mp)
  s = Hecke.Hecke.local_fundamental_class_serre(Kp, prime_field(Kp))
#  Oscar.GModuleFromGap.istwo_cocycle(Dict( (g, h) => s(mGp(g), mGp(h)) for g = domain(mGp) for h = domain(mGp)), mGp)

  z = gmodule(domain(mGp), [hom(domain(mUp), domain(mUp), [preimage(mUp, mGp(g)(mUp(u))) for u = gens(domain(mUp))]) for g = gens(domain(mGp))])

  c = CoChain{2, PermGroupElem, GrpAbFinGenElem}(z, Dict{NTuple{2, PermGroupElem}, GrpAbFinGenElem}((g, h) => preimage(mUp, s(mGp(g), mGp(h))) for g = domain(mGp) for h = domain(mGp)))

  @assert Oscar.GrpCoh.istwo_cocycle(c)

  return c
end

function serre(C::GModule, A::IdelParent, P::Union{Integer, ZZRingElem})
  t = findfirst(isequal(ZZ(P)), [minimum(x) for x = A.S])
  Inj = Hecke.canonical_injection(A.M, t+1)
  Pro = Hecke.canonical_projection(A.M, t+1)

  inj = Hecke.canonical_injection(domain(Inj), 1)
  pro = Hecke.canonical_projection(domain(Inj), 1)

  Kp, mKp, mGp, mUp, _, _ = completion(A, A.S[t])
  @assert domain(inj) == domain(mUp) 
  mp = decomposition_group(A.k, mKp, A.mG, mGp)
 
  tt = serre(C, A, A.S[t])
  @assert tt.C.G == domain(mGp)

  I = domain(Inj)    
  zz = gmodule(C.G, [Inj * action(A.data, g) * Pro for g = gens(C.G)])
  mu = cohomology_group(zz, 2)
  q, mq = snf(mu[1])
  g = mu[2](mq(q[1]))
  hg = map_entries(Inj*A.mq, g, parent = C)
  gg = map_entries(pro, g, parent = tt.C)
  gg = Oscar.GrpCoh.CoChain{2, PermGroupElem, GrpAbFinGenElem}(tt.C, Dict( (g, h) => gg.d[mp(g), mp(h)] for g = tt.C.G for h = tt.C.G))

  nu = cohomology_group(tt.C, 2)
  ga = preimage(nu[2], gg)
  ta = preimage(nu[2], tt)
  return findfirst(x->x*ga == ta, 1:order(tt.C.G)), hg
  #so i*hg should restrict to the local fund class...
end


function global_fundamental_class(C::GModule, A::IdelParent)
  d = lcm([ramification_index(P) * inertia_degree(P) for P = A.S])
  G = C.G
  if d != order(G)
    error("sorry - no can do(yet)")
  end

  z = cohomology_group(C, 2)

  q, mq = snf(z[1])
  @assert ngens(q) == 1
  g = z[2](mq(gen(q, 1))) # to get a 2-CoCycle
  @assert Oscar.GrpCoh.istwo_cocycle(g)

  scale = []

  for P = A.S
    s = serre(C, A, minimum(P))
    push!(scale, s)
  end
  #put to gether..
  return scale, z, mq 
end

function Oscar.orbit(C::GModule{PermGroup, GrpAbFinGen}, o::GrpAbFinGenElem)
  or = Set([o])
  done = false
  while !done
    sz = length(or)
    done = true
    for f = C.ac
      while true
        or = union(or, [f(x) for x = or])
        if length(or) == sz
          break
        end
        done = false
        sz = length(or)
      end
    end
  end
  return collect(or)
end

"""
    shrink(C::GModule{PermGroup, GrpAbFinGen}, attempts::Int = 10)

Tries to find cohomologically trivial submodules to factor out.
Returns a cohomologically equivalent module with fewer generators and
the quotient map.
"""
function shrink(C::GModule{PermGroup, GrpAbFinGen}, attempts::Int = 10)
  local mq
  q = C
  first = true
  while true
    prog = false
    for i=1:attempts
      o = Oscar.orbit(q, rand(gens(q.M)))
      if length(o) == order(group(q))
        s, ms = sub(q.M, o)
        if rank(s) == length(o)
          q, _mq = quo(q, ms)
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

function Oscar.direct_sum(G::GrpAbFinGen, H::GrpAbFinGen, V::Vector{<:Map{GrpAbFinGen, GrpAbFinGen}})
  dG = get_attribute(G, :direct_product)
  dH = get_attribute(H, :direct_product)

  if dG === nothing || dH === nothing
    error("both groups need to be direct products")
  end
  @assert length(V) == length(dG) == length(dH)

  @assert all(i -> domain(V[i]) == dG[i] && codomain(V[i]) == dH[i], 1:length(V))
  h = hom(G, H, cat([matrix(V[i]) for i=1:length(V)]..., dims=(1,2)), check = !true)
  return h

end

function Oscar.simplify(C::GModule{PermGroup, GrpAbFinGen})
  s, ms = snf(C.M)
  S = GModule(s, C.G, [GrpAbFinGenMap(ms*x*pseudo_inv(ms)) for x = C.ac])
  if isdefined(C, :iac)
    S.iac = [GrpAbFinGenMap(ms*x*pseudo_inv(ms)) for x = C.iac]
  end
  return S, ms
end

function Base.show(io::IO, I::IdelParent)
  print(io, "Idel-group for $(I.k) using $(sort(collect(Set(minimum(x) for x = I.S)))) as places")
end

end # module GrpCoh

using .GaloisCohomology_Mod
export is_coboundary, idel_class_gmodule


#=
x^4 - 60*x^2 + 16

=#
