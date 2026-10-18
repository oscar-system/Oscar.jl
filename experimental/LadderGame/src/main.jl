# Compute A\G/U
# TODO indices should start at 1...
# we need: a subgroup ladder L = (A_i)_(0<=i<=a)
# with A_0 = G and A_a = A
# Transversals: For each (A_i+1, A_i)
#           T_i in R(A_i\A_i+1) or R(A_i+1\A_i) as appropriate
# Test functions: For each i with A_i+1 <= A_i,
#           tau_i : A_i+1\A_i >-> NN (injection)
# Need a transversal chain for U:
#           subgroup chain (U_i)_(0<=i<=k)
#           transversals (S_i)_(0<i<=k) of R(U_i\U_i-1)

# O(1) zero-cost parent caster: skips GAP Schreier-Sims membership sifting
@inline function reparent(G::PermGroup, x::PermGroupElem)
  parent(x) === G ? x : Oscar.PermGroupElem(G, x.X)
end


function ladder_game!(L::SubgroupLadder, G::PermGroup, U::PermGroup)
  sprev = L[1]
  for s in L
    s.A == s.Aprev && (_ladder_start!(s, U), continue)
    s.is_up_step ? up_step(s, sprev, U) : down_step(s, sprev, U)
    sprev = s
  end
  A = L[end].A
  res = GroupDoubleCoset{PermGroup, PermGroupElem}[]
  for (g, chain) in L[end].DSt
    C = GroupDoubleCoset(G, A, U, g)
    stab_order = order(data(chain)[1][1])
    n = divexact(order(A), stab_order) * order(U)
    C.size[] = ZZRingElem(n)
    push!(res, C)
  end
  return res
end

function _ladder_start!(S::LadderStep, U::PermGroup)
  # @assert S.A == S.Aprev "Initial ladder step must be trivial"
  if !isdefined(S, :test)
    S.test = FunctionTest([one(S.A)], x -> one(S.A); ambient = _ambient_group(S))
  end
  # isdefined(S, :T) || (S.T = [one(S.A)])
  # isdefined(S, :Tmap) || (S.Tmap = map_from_func(S.A, S.A, x -> one(S.A)))
  S.DSt = Dict(one(S.A) => get_transversal_chain(U))
  S.Im = Dict(one(S.A) => (one(S.A), one(U))) # should both be one(S.A) for the sake of parents?
  return S
end

function down_step(S::LadderStep, Sprev::LadderStep, U::PermGroup)
  A = S.A
  Aprev = S.Aprev
  G = _ambient_group(S)

  if !isdefined(S, :test)
    S.test = GAPTransversalTest(S.Aprev, S.A; ambient = _ambient_group(S))
  end
  T = elements(S.test)

  DSt = Dict{PermGroupElem, TransversalChain}()
  Im = Dict{PermGroupElem, Tuple{PermGroupElem, PermGroupElem}}()

  for (a, st) in Sprev.DSt
    inv_a = inv(a)
    seen = IndexedSet{PermGroupElem}()
    S2 = st[1][1]
    for t in T
      t in seen && continue
      d = t*a
      S1 = _my_intersect(S2, A^d)

      # not used in calculation
      # S2 = _my_intersect(U, Aprev^a)
      # @assert is_subset(S1, S2)
      # @assert st[1][1] == S2
      # @assert S2 == _my_intersect(U, Aprev^d)

      st2, Sta = _induce_chain(S1, S.test, st ; conj=d)
      map!(x -> reparent(U, x), Sta, Sta)

      # not used in calculation
      # @assert length(Sta) == index(S2, S1)
      # Tss = right_transversal(S2, S1)
      # ss = map_from_func(S2, S2, x -> Tss[index_of_coset(Tss, x)])
      # @assert allunique(map(x -> ss(S2(x)), Sta))

      DSt[d] = st2
      for u in Sta
        tt = coset_rep(S.test, d*reparent(G, u)*inv_a)
        push!(seen, tt)
        Im[tt*a] = (d, u)
        # @assert tt*a*inv(u)*inv(d) in A
      end
    end
    # @assert Set(collect(T)) == Set(tmp)
  end
  S.DSt = DSt
  S.Im = Im

  return S
end

function up_step(S::LadderStep, Sprev::LadderStep, U::PermGroup)
  A = S.A
  Aprev = S.Aprev

  if !isdefined(S, :test)
    S.test = GAPTransversalTest(S.A, S.Aprev; ambient = _ambient_group(S))
  end
  T = elements(S.test)

  DSt = Dict{PermGroupElem, TransversalChain}()
  Im = Dict{PermGroupElem, Tuple{PermGroupElem, PermGroupElem}}()
  cls_store = Dict{PermGroupElem, Tuple{PermGroupElem, PermGroupElem}}()

  for (a, st) in Sprev.DSt
    if haskey(cls_store, a)
      (rep, rep_u) = cls_store[a]
      new = false
    else
      rep = nothing
      new = true
    end

    tt = PermGroupElem[]
    for t in T
      (~new && ~is_one(t)) && continue
      ta = t*a
      (at, ut) = _get_rep(ta, Sprev)

      # @assert ta*inv(ut)*inv(at) in Aprev
      # @assert at in keys(Sprev.DSt)

      rep === nothing && ((rep, rep_u) = (at, ut))
      tr_elem = rep_u*inv(ut)

      haskey(cls_store, at) || (cls_store[at] = (rep, tr_elem))
      at == rep && push!(tt, ut*inv(rep_u))
      is_one(t) && (Im[a] = (rep, tr_elem))
    end
    if new
      rep_ = rep::PermGroupElem
      StTC = TransversalChain([[ (_my_intersect(A^rep_, U), tt) ]; data(st)])
      DSt[rep_] = StTC
      # @assert index(StTC[1][1], StTC[2][1]) == length(tt)
    end
  end

  S.DSt = DSt
  S.Im = Im

  return S
end

# This gives a transversal chain
function _induce_chain(V::PermGroup, tf::AbstractTestFunction, C::TransversalChain ; conj=one(C[1][1]))
  # Tm is a transversal map for C[1][1]/V
  # @assert is_subset(V, C[1][1])

  c = Tuple{PermGroup, Vector{PermGroupElem}}[ (trivial_subgroup(V)[1], [one(V)]) ]

  last_tUU = PermGroupElem[ one(C[end][1]) ]
  tU = IndexedSet([ one(_ambient_group(tf)) ]) # magma version: { Tm(one(parent(codomain(Tm)))^inv(conj)) } but why???
  conj_inv = inv(conj)
  for i in (length(C)-1):-1:1
    map!(x -> reparent(C[i][1],x), last_tUU, last_tUU)
    tUU = copy(last_tUU)

    U = _my_intersect(C[i][1], V)
    tV = PermGroupElem[]
    sizehint!(tV, length(C[i][2]))
    target_count = length(last_tUU) * length(C[i][2])
    for t in last_tUU, s in C[i][2]
      x = t*s
      if x in V
        push!(tV, x)
        length(tU)*length(tV) == target_count && break
      else
        xx = coset_rep(tf, reparent(parent(conj), x)^conj_inv)
        if ~(xx in tU)
          push!(tU, xx)
          push!(tUU, x)
          length(tU)*length(tV) == length(tUU)*length(C[i][2]) && break
        end
      end
    end
    # @assert length(tV) == index(U, c[end][1])
    # @assert length(tU) == index(C[i][1], U)
    # @assert length(tU)*length(tV) == length(last_tUU)*length(C[i][2])

    last_tUU = tUU
    length(tV) != 1 && push!(c, (U, tV) )
  end

  return TransversalChain(reverse(c)), last_tUU
end

# Recursive recognition method
function _get_rep(g::PermGroupElem, S::LadderStep)
  F = _get_previous_steps(S)
  ambient = _ambient_group(S)

  dprev, uprev = F[1].Im[one(ambient)]
  for s in F
    p = s.is_up_step ? dprev : coset_rep(s.test, g*inv(dprev*uprev)) * dprev
    (d, u) = s.Im[p]
    dprev = d
    uprev = u*uprev
  end
  return dprev, uprev
end



################################################################################
#
#  Young subgroup ladders
#
################################################################################

function young_subgroup(p::Vector{T} ; full::T=sum(p)) where T<:Integer
  s = full - sum(p)
  # should error if s < 0
  s == 0 && return inner_direct_product([symmetric_group(i) for i in p])
  return inner_direct_product([symmetric_group(i) for i in vcat(p, ones(T, s))])
end

# works - matches Magma
function _young_subgroup_ladder( p::Vector{T} ; full::T=sum(p)) where T<:Integer
  n = sum(p)
  L = PermGroup[]
  pp = copy(p)
  pushfirst!(L, young_subgroup(pp ; full=full))
  s = n - pp[1]
  while s!=0
    if pp[2]==1
      deleteat!(pp, 2)
      pp[1]+=1
    else
      pp[2]-=1
      insert!(pp,2,1)
      pushfirst!(L, young_subgroup(pp ; full=full))
      deleteat!(pp, 2)
      pp[1]+=1
    end
    pushfirst!(L, young_subgroup(pp ; full=full))
    s = n-pp[1]
  end

  return L
end


# TODO make transversal maps, integrate into each LadderStep
function young_subgroup_ladder( p::Vector{T} ; full::T=sum(p)) where T<:Integer
  n = sum(p)
  L = LadderStep[]
  pp = copy(p)
  Hprev = young_subgroup(pp ; full=full)
  s = n - pp[1]
  while s!=0
    if pp[2]==1
      deleteat!(pp, 2)
      pp[1]+=1
    else
      pp[2]-=1
      insert!(pp,2,1)
      H = Hprev
      Hprev = young_subgroup(pp ; full=full)
      pushfirst!(L, LadderStep(Hprev, H))
      deleteat!(pp, 2)
      pp[1]+=1
    end
    H = Hprev
    Hprev = young_subgroup(pp ; full=full)
    pushfirst!(L, LadderStep(Hprev, H))
    s = n-pp[1]
  end

  pushfirst!(L, LadderStep(Hprev,Hprev))

  return SubgroupLadder(L)
end


function _my_intersect(G1::PermGroup, V::Oscar.GAPGroup...)
  return _my_intersect([G1, V...])
end

function _my_intersect(V::AbstractVector{PermGroup})
  L = GapObj(V; recursive = true)
  K = Oscar.GAPWrap.Intersection(L)
  return Oscar._as_subgroup(V[1], K)[1]
end






# NOTE don't use this as-is
# TODO make this our own, optimize, needs to use up-steps when relevant



# Find a subgroup Anext with Aprev > Anext >= A, reducing the index as much as
# convenient, by finding block systems of Aprev acting on cosets of A.
# Falls back to A itself if the action is primitive (no further splitting possible).
function find_good_descent_step(Aprev::PermGroup, A::PermGroup; max_step_index::Int = 10000)
  idx = index(Aprev, A)
  idx <= max_step_index && return A

  c = right_coset(A, one(Aprev))
  Omega = orbit(Aprev, c)

  B = blocks(Omega)
  block = representative(B)

  # primitive edge case: no nontrivial block system exists,
  # no further splitting possible via this method
  length(block) == 1 && return A

  Anext, _ = stabilizer(B, block)
  return Anext
end



function build_strong_ladder(G::PermGroup, A::PermGroup; max_step_index::Int = 10_000)
  groups = PermGroup[G]
  Aprev = G
  while Aprev != A
    Anext = find_good_descent_step(Aprev, A; max_step_index)
    push!(groups, Anext)
    Aprev = Anext
  end
  return SubgroupLadder(groups)
end
