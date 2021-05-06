module GrpCoh

using Oscar
import AbstractAlgebra: Group, Module

mutable struct CohomologyModule
  G::Group
  M::Union{Module, GrpAbFinGen}
  ac::Array{Map, 1}
end

AbstractAlgebra.Group(C::CohomologyModule) = C.G
AbstractAlgebra.Module(C::CohomologyModule) = C.M
action(C::CohomologyModule) = C.ac

function H_zero(C::CohomologyModule)
  G = Group(C)
  M = Module(C)
  id = hom(M, M, gens(M))
  ac = action(C)
  k = kernel(id - ac[1])[1]
  for i=2:length(ac)
    k = intersect(k, kernel(id - ac[i])[1])
  end
  return k
end

function H_one(C::CohomologyModule)
  #= idea, after Holt:
  H^1 = crossed homs. due to actio on the right
  f(ab) = f(a)^b + f(b)
  if G=<g_1, ..., g_r | r_1, ..., r_l>
  then X in H^1 iff X(r_i) = 0 for all i
  X:G->M is given as X in M^r, where X(g_i) = X[i]
  X(r_i) corresponds to some map phi_i : M^r -> M
  phi_i = oplus h_j M for some homs h_j coming from the word in r
  so, a kernel computation again

  =#
  G = Group(C)
  n = ngens(G)
  M = Module(C)
  D, pro, inj = direct_product([M for i=1:n]..., task = :both)

  X = GAP.Globals.IsomorphismFpGroupByGenerators(G.X, GAP.Globals.GeneratorsOfGroup(G.X))
  X = GAP.Globals.Range(X)
  R = GAP.Globals.RelatorsOfFpGroup(X)
  K = D
  ac = action(C)
  idM = hom(M, M, gens(M)) #identity map to start with
  for r = R
    W = GAP.Globals.LetterRepAssocWord(r)
    g = idM
    P = hom(D, M, [id(M) for i=1:ngens(D)])
    for w in GAP.Globals.Reversed(W)
      if w < 0
        g = inv(ac[-w])*g
        P -= pro[-w]*g
      else
        P += pro[w]*g
        g = ac[w]*g
      end
    end
    K = intersect(kernel(P)[1], K)
  end
  #K is Z[1]  - the co-cycles
  #TODO: is kernel(g) directly faster than the method above (H_zero)
  #      where kernel(g) is computed slice by slice?
  #TODO: cache the expensive objects!!!
  g = sum((ac[i] - idM)*inj[i] for i=1:n)
  return quo(K, image(g)[1])[1]
end

function confluent_rws(G)
  C = GAP.Globals.ConfluentMonoidPresentationForGroup(G.X)
  #has different generators than G! So the action will have to
  #be adjusted to those words. I do not know if a RWS (Confluent) can
  #just be changed...
  k = C[2] #hopefully the monhom entry in 4.12 it will be the name
  M = GAP.Globals.Range(k)
  g = [GAP.Globals.PreImageElm(k, x) for x = GAP.Globals.GeneratorsOfMonoid(M)]
  g = map(GAP.Globals.UnderlyingElement, g)
  g = map(GAP.Globals.LetterRepAssocWord, g)
  @assert all(x->length(x) == 1, g)
  g = [Int(x[1]) for x = g]
  #the Gap code uses RelationsOfFpMonoid rather than
  # ReducedConfluentRewritingSystem...
  # no idea why
#  R = GAP.Globals.Rules(GAP.Globals.ReducedConfluentRewritingSystem(M))
  R = GAP.Globals.RelationsOfFpMonoid(M)

  ru = Array{Tuple{Array{Int, 1}, Array{Int, 1}}, 1}()
  for r = R
    push!(ru, (map(x->g[Int(x)], GAP.Globals.LetterRepAssocWord(r[1])), 
               map(x->g[Int(x)], GAP.Globals.LetterRepAssocWord(r[2]))))
  end

  #now to express the new gens as words in the old ones:
  
  mFp = GAP.Globals.IsomorphismFpGroupByGenerators(G.X, GAP.Globals.GeneratorsOfGroup(G.X))
  ge = Array{Vector{Int}, 1}()
  Fp = GAP.Globals.Source(k)
  gFp = GAP.Globals.GeneratorsOfGroup(Fp)
  for i=1:length(gFp)
    i = GAP.Globals.ImageElm(mFp, GAP.Globals.PreImageElm(C[1], gFp[i]))
    push!(ge, map(Int, GAP.Globals.LetterRepAssocWord(GAP.Globals.UnderlyingElement(i))))
  end

  #actually, need the map (isomorphism) between G <-> Fp group with the
  #RWS in ru
  return ru, ge
end

mutable struct CollectCtx
  r::Array{Tuple{Vector{Int}, Vector{Int}}, 1} #the rules, RWS

  d1::Dict{Int, Int} #rules where lhs has length 1

  d2::Dict{Tuple{Int, Int}, Array{Int, 1}} # length 2 prefixes

  f::Function #w::Array{Int, 1}, r::Int, p::Int
              #to be called in addition (to play with the tail(s))
              #w the word, to be "reduced" using rule no r at pos p
  T::Any #data used in f

  function CollectCtx(R::Array{Tuple{Vector{Int}, Vector{Int}}, 1})
    n = new()
    n.r = R
    n.d1 = Dict{Int, Int}()
    n.d2 = Dict{Tuple{Int, Int}, Array{Int, 1}}()
    for i = 1:length(R)
      r = R[i]
      if length(r[1]) == 1
#        @assert length(r[2]) == 1
        n.d1[r[1][1]] = r[2][1]
        continue
      end
      @assert length(r[1]) > 1
      p = (r[1][1], r[1][2])
      if Base.haskey(n.d2, p)
        push!(n.d2[p], i)
      else
        n.d2[p] = [i]
      end
    end
    for p = keys(n.d2)
      sort!(n.d2[p], lt = (a,b) -> isless(R[a], R[b]))
    end
    return n
  end
end

function collect(w::Array{Int, 1}, C::CollectCtx)
  d1 = C.d1
  d2 = C.d2
  R = C.r
  do_f = isdefined(C, :f)
  for i=1:length(w)
    if haskey(d1, w[i])
      w[i] = d1[w[i]]
    end
  end

  nc = 0
  i = 1
  while true
    nc += 1
    if i>=length(w)
      return w
    end
    if haskey(d2, (w[i], w[i+1]))
      for r = d2[(w[i], w[i+1])]
        if length(R[r][1]) + i-1 <= length(w) &&
           R[r][1] == w[i:i+length(R[r][1])-1]
          if do_f
            C.f(w, r, i)
          end
          w = vcat(w[1:i-1], R[r][2], w[i+length(R[r][1]):end])
          i = 0
          break
        end
      end
    end
    i += 1
  end
  return w
end

function Oscar.hom(V::Module, W::Module, v::Vector{<:ModuleElem}; check::Bool = true)
  return Generic.ModuleHomomorphism(V, W, vcat([x.v for x = v]))
end
function Oscar.hom(V::Module, W::Module, v::MatElem; check::Bool = true)
  return Generic.ModuleHomomorphism(V, W, v)
end
function Oscar.inv(M::Generic.ModuleHomomorphism)
  return hom(codomain(M), domain(M), inv(mat(M)))
end

function Oscar.direct_product(M::Module...; task::Symbol = :none)
  D, inj, pro = direct_sum(M...)
  if task == :none
    return D
  elseif task == :both
    return D, pro, inj
  elseif task == :sum
    return D, inj
  elseif task == :prod
    return D, pro
  end
  error("illegal task")
end

Base.:+(a::Generic.ModuleHomomorphism, b::Generic.ModuleHomomorphism) = hom(domain(a), codomain(a), mat(a) + mat(b))
Base.:-(a::Generic.ModuleHomomorphism, b::Generic.ModuleHomomorphism) = hom(domain(a), codomain(a), mat(a) - mat(b))
Base.:-(a::Generic.ModuleHomomorphism) = hom(domain(a), codomain(a), -mat(a))

Base.zero(G::GrpAbFinGen) = G[0]

function H_two(C::CohomologyModule)
  G = Group(C)
  M = Module(C)
  id = hom(M, M, gens(M), check = false)
  Ac = action(C)
  iAc = map(inv, Ac)

  R, ge = confluent_rws(G)
  #now map the action generators (for gens(G)) to the gens for the RWS
  ac = []
  iac = []
  for g = ge
    f = id
    for i = g
      if i < 0
        f = f*iAc[-i]
      else
        f = f*Ac[i]
      end
    end
    push!(ac, f)
    push!(iac, inv(f))
  end

  c = CollectCtx(R)

  #rules with length(LHS) == 1 and rules of the form
  # [a a^-1] -> [], [a^-1 1] -> [] do not get tails
  pos = Array{Int, 1}()
  n = 0
  for i = 1:length(R)
    r = R[i]
    if length(r[1]) == 1
      push!(pos, 0)
      continue
    end
    if length(r[1]) == 2 && length(r[2]) == 0 && r[1][1] == -r[1][2]
      push!(pos, 0)
      continue
    end
    n += 1
    push!(pos, n)
  end

  D, pro, inj = direct_product([M for i=1:n]..., task = :both)

  #when collecting (i.e. applying the RWS we need to also
  #use the tails:  g v h -> gh h(v) 
  #and if [gh] -> [x] with tail t, then
  #       gh -> x t, so 
  #       g v h -> gh h(v) -> x t+h(v)
  # Hulpke calls this the normal version: reduced group word
  # at the beginning, module at the end, the tail.
  # collect will call the extra function c.f if set in the
  # CollectCtx
  c.f = function(w::Array{Int, 1}, r::Int, p::Int)
    #w = ABC and B == r[1], B -> r[2] * tail[r]
    # -> A r[2] C C(tail)
    # C = c1 c2 ... C(tail):
    @assert w[p:p+length(R[r][1])-1] == R[r][1]

    if pos[r] == 0
      return
    end
    T = pro[pos[r]]
    for i=w[p+length(R[r][1]):end]
      if i < 0
        T = T*iac[-i]
      else
        T = T*ac[i]
      end
    end
    c.T += T
  end

  E = D
  all_T = []
  Z = hom(D, M, [M[0] for i=1:ngens(D)], check = false)
  for i = 1:length(R)
    r = R[i]
    #rules of LHS length 1 do not generate equations
    if length(r[1]) == 1
      continue
    end
    for j=1:length(R)
      s = R[j]
      if length(s[1]) == 1
        continue
      end
      #we want overlaps, all of them:
      #r[1] = AB, s[1] = BC this is what we need to find...
      #(then we call collect on r[2]C and As[2] they should agree)
      for l=1:min(length(s[1]), length(r[1]))
        if r[1][end-l+1:end] == s[1][1:l]
          #TODO  AB    -> Ss  s,t are tails
          #       BC   -> Tt
          #      (AB)C -> SsC -> SC C(s)
          #      A(BC) -> ATt -> AT t
          if pos[i] > 0 
            c.T = pro[pos[i]]
            for h = s[1][l+1:end]
              if h < 0
                c.T = c.T * iac[-h]
              else
                c.T = c.T * ac[h]
              end
            end
          else
            c.T = Z
          end
          z1 = collect(vcat(r[2], s[1][l+1:end]), c)
          T = c.T
          c.T = Z
          z2 = collect(vcat(r[1][1:end-l], s[2]), c)
          if pos[j] > 0
            c.T += pro[pos[j]]
          end
          @assert z1 == z2
          push!(all_T, T-c.T)
        end
      end
    end
  end

  Q, jinj = direct_product([M for i in all_T]..., task = :sum)
  mm = sum(all_T[i]*jinj[i] for i = 1:length(all_T))
  E = kernel(mm)[1]


  B, B_pro = direct_product([M for i=1:length(ac)]..., task = :prod)
  C = hom(B, D, [zero(D) for i=1:ngens(B)], check = false)
  for i=1:length(R)
    if pos[i] == 0
      continue
    end
    r = R[i]
    if length(r[1]) == 1
      continue
    end
    #we have words r[1] and r[2] or shape g_1 g_2 .... 
    #they need to be replaced by g_1 pro[1] g_2 pro[2]
    #and then sorted: g_1 pro[1] g_2 pro[2] ... ->
    #                 g_1 g_2 (pro[1] * g_2 + pro[2]) ...
    if r[1][1] < 0
      T = -B_pro[-r[1][1]]*iac[-r[1][1]]
    else
      T = B_pro[r[1][1]]
    end
    for j=2:length(r[1])
      if r[1][j] < 0
        T = (T-B_pro[-r[1][j]])*iac[-r[1][j]] 
      else
        T = T*ac[r[1][j]] + B_pro[r[1][j]]
      end
    end

    if length(r[2]) == 0
      S = hom(B, M, [M[0] for g = gens(B)], check = false)
    elseif r[2][1] < 0
      S = -B_pro[-r[2][1]]*iac[-r[2][1]]
    else
      S = B_pro[r[2][1]]
    end
    for j=2:length(r[2])
      if r[2][j] < 0
        S = (S-B_pro[-r[2][j]])*iac[-r[2][j]]
      else
        S = S*ac[r[2][j]] + B_pro[r[2][j]]
      end
    end

#    @assert issubset(image((T-S)*inj[pos[i]])[1], E)

    C += (T-S)*inj[pos[i]]
  end
  i = image(C)[1]
  H2, mH2 = quo(E, image(C)[1])
  return H2, mH2
  #now the rest...
  #(g, m)*(h, n) = (gh, m^h+n+gamma(g, h)) where gamma is "the" 2-cocycle
  #using tails:
  # gmhn -> gh h(m)+n -> x t+h(m) + n where x is the reduced
  #                                   word under collection and t is the 
  #                                   "tail"
  # so gamma(g, h) = t
  # given gamma need the tails:
  # need to implement the group operation for the extension
  # (g, v)(h, u) -> (gh, v^h + u + gamma(g, h))
  # then the rules with tails need to be evaluated at
  # the group generators (g_i, 0) 
  # r -> s gives a relation r s^-1 which should evaluate, using gamma
  # to (0, t) where t is the tail for this rule
end

function gamma(f::PermGroup, g::PermGroup, mH2, pro, cocycle)
end

Base.:-(M::GrpAbFinGenMap) = hom(domain(M), codomain(M), [-M(g) for g = gens(domain(M))], check = false)

function Oscar.automorphism_group(::Type{PermGroup}, k::AnticNumberField)
  G, mG = automorphism_group(k)
  H = symmetric_group(degree(k))
  H = sub(H, [H(G.mult_table[i, :]) for i=G.gens])[1]

  function HtoG(p::PermGroupElem)
    m = [i^p for i=1:degree(k)]
    i = Base.findfirst(x->G.mult_table[x, :] == m, 1:degree(k))
    return mG(GrpGenElem(G, i))
  end

  function GtoH(a::NfToNfMor)
    g = preimage(mG, a)
    return H(G.mult_table[g.i, :])
  end

  return H, MapFromFunc(HtoG, GtoH, H, codomain(mG))
end

function cohomology_module(H::PermGroup, mR::MapRayClassGrp)
  k = nf(order(codomain(mR)))
  G, mG = automorphism_group(PermGroup, k)

  ac = Hecke.induce_action(mR, [image(mG, G(g)) for g = gens(H)])
  return CohomologyModule(H, domain(mR), ac)
end


#= TODO
  for Z, Z/nZ, F_p and F_q moduln -> find Fp-presentation
  for GrpAbFinGen            -> find Fp-presentation
  for a in H^2 find Fp-presentation
  for a in H^2 find (low degree) perm group using the perm group we have?
  Magma's DistinctExtensions
  probably aut(GrpAb), ...

Sort: 
 - move the additional GrpAbFinGenMap stuff elsewhere
 - move (and fix) the ModuleHom stuff
 - add proper quo for Modules

  group better: the CohomologyModule should
   - cache the H^i's that are computed
   - we need to sort caching for the IsomorphismFpGroupByGenerators as it is
     used frequently

  features   
   - need to provide translations to other reps for cocycles, coboundaries
     (gamma in H^i -> function from G^i -> M)
     (for cyclic G special?)
   - test for triviality and return the boundary!
   - inflation, restriction, long exact sequence  

  dreams
   - we we extend to H^-1, ^-2?
   - H^3 (in some cases)
   - cup products
   - the relative cohomology
     https://arxiv.org/pdf/1809.01209.pdf
     https://doi.org/10.1017/S2040618500033050
   - understand Derek Holt and use BSGS for large perm groups
     rather than the RWS (or use BSGS to get an RWS?)


  CohomologyModule for 
    - abelian_extension (do for subgroups of the aut of the base field!)
    - local field (add (trivial) and mult)
    - (S-)units
    - Ali's stuff....
=#    

end # module GrpCoh
