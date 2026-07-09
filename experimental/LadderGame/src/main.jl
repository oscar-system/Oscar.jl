
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

function ladder_game!(L::SubgroupLadder, G::PermGroup, U::PermGroup)
  s = _ladder_start(L[1], U)
  # H = L[1]
  for i in 2:length(L)
    if is_subgroup(H, L[i])
      r = up_step(L[i],U,r)
    else
      r = down_step(L[i], U, r)
    end
    H = L[i]
  end
  return r
end

function _ladder_start!(S::LadderStep, U::PermGroup)
  @assert S.A == S.Aprev "Initial ladder step must be trivial"
  isdefined(S, :T) || (S.T = [one(S.A)])
  isdefined(S, :Tmap) || (S.Tmap = map_from_func(S.A, S.A, x -> one(S.A)))
  S.D = [one(S.A)]
  S.St = [get_transversal_chain(U)]
  S.I = [one(S.A)]
  S.m = map_from_func(S.I, (S.D, U), x -> (one(S.A), one(U)))
  return S
end

function down_step(S::LadderStep, Sprev::LadderStep, U::PermGroup)
  A = S.A
  Aprev = S.Aprev

  isdefined(S, :T) || (S.T = right_transversal(Aprev, A))
  isdefined(S, :Tmap) || (S.Tmap = map_from_func(Aprev, Aprev, x -> S.T[index_of_coset(T, x)]))

  T = S.T

  D = []
  St = []
  I = []
  S3 = []

  for (i, a) in enumerate(Sprev.D)
    tmp = PermGroupElem[]
    for t in T
      if ~(t in tmp)
        d = t*a
        st = Sprev.St[i]
        S1, _ = intersect(U, A^d)

        # not used in calculation
        S2, _ = intersect(U, Aprev^a)
        @assert is_subset(S1, S2)
        @assert st[1][1] == S2
        @assert S2 == intersect(U, Aprev^d)[1]

        st, Sta = _induce_chain(S1, S.Tmap, st ; conj=d)

        # not used in calculation
        @assert length(Sta) == index(S2, S1)
        Tss = right_transversal(S2, S1)
        ss = map_from_func(S2, S2, x -> Tss[index_of_coset(Tss, x)])
        @assert allunique(map(x -> ss(S2(x)), Sta))

        push!(D,d)
        push!(St, st)
        for u in Sta
          tt = S.Tmap(Aprev(d*u*inv(a)))
          push!(tmp, tt)
          push!(I, tt*a)
          push!(S3, (d,u) )
          @assert tt*a*inv(u)*inv(d) in A
        end
      end
    end
    @assert Set(collect(T)) == Set(tmp)
  end
  S.D = D
  S.St = St
  S.I = I
  S.m = map_from_func(S.I, (S.D, U), x -> (one(S.A), one(U)))
  # S.m = map_from_func(S.I, (S.D, U), x -> )
  # S3 # TODO this should be a map
  # could we use a Dict I -> (S.D, U) ???
  return S
end

function up_step(S::LadderStep, Sprev::LadderStep, U::PermGroup)
  # A::PermGroup, U::PermGroup, rec::LadderStep)
  # r = LadderStep(
  #       # Ai = A
  #       # Aim1 = rec.Ai
  #       # last = rec
  #       # g,s for ev variable "F" set to false
  # )

  A = S.A
  Aprev = S.Aprev

  isdefined(S, :T) || (S.T = right_transversal(A, Aprev))
  isdefined(S, :Tmap) || (S.Tmap = map_from_func(A, A, x -> S.T[index_of_coset(T, x)]))

  T = S.T

  D = []
  S3 = []
  St = []
  I = []

  seen = []
  data = []
  # TODO cls/cls_data should be Dict
  cls = []
  cls_data = []

  for a in Aprev.D
    if a in cls
      p = findfirst(a, cls)
      rep = cls_data[p][1]
      rep_u = cls_data[p][2]
      new = false
    else
      rep = nothing
      new = true
    end

    tt = []
    for t in T
      (~new && ~is_one(t)) && continue
      ta = t*a
      at, ut = _get_rep(ta, Aprev)

      @assert ta*inv(ut)*inv(at) in Aprev
      @assert at in Sprev.D

      rep === nothing && (rep = at, rep_u = ut)
      at in cls || (push!(cls, at), push!(cls_data, (rep, rep_u*inv(ut)) ))
      at == rep && push!(tt, ut*inv(rep_u))
      is_one(t) && (push!(I, a), push!(S3, (rep, rep_u*inv(ut)) ))
    end
    if new
      push!(D, rep)
      push!(St, [[ (intersect(A^rep, U)[1], tt) ] Sprev.St[findfirst(rep, Aprev.D)]] ) # in magma code: rec.St[Position(R.Di, rep)]
      # assert Index(St[end][1][1], St[end][2][1]) == length(tt)
    end
  end

  S.D = D
  S.St = St

  S.I = I
  # S.m = map_from_func(S.I, (S.D, U), x -> (one(S.A), one(U))) TODO fix
  # S.S3 = S3

  return S
end

# I think this gives a transversal chain ???
# can streamline wrt tUU / last_tUU
function _induce_chain(V::PermGroup, Tm::Map, C::Vector ; conj=one(C[1][1]))
  # Tm is a MAP: transversal map for C[1][1]/V
  @assert is_subset(V, C[1][1])

  c = Tuple{PermGroup, Vector{PermGroupElem}}[ (sub(V, [one(V)])[1], [one(V)]) ]

  last_tUU = [ one(C[end][1]) ]
  # tUU = change_universe(last_tUU, C[1][1]) # - unnecessary???
  tUU = copy(last_tUU)
  tU = [ Tm( one(codomain(Tm))^inv(conj) ) ] # isn't this just Tm(1)?

  for i in (length(C)-1):-1:1
    # tUU = change_universe(last_tUU, C[i][1])
    # change_universe!(last_tUU, C[i][1])

    U, _ = intersect(C[i][1], V)
    tV = []

    for t in last_tUU, s in C[i][2]
      x = t*s
      if x in V
        push!(tV, x)
        length(tU)*length(tV) == length(last_tUU)*length(C[i][2]) && break
      else
        xx = Tm(domain(Tm)(x^inv(conj)))
        if ~(xx in tU)
          push!(tU, xx)
          push!(tUU, x)
          length(tU)*length(tV) == length(tUU)*length(C[i][2]) && break
        end
      end
    end
    # assert length(tV) == index(U, c[end][1])
    # assert length(tU) == index(C[i][1], U)
    # assert length(tU)*length(tV) == length(last_tUU)*length(C[i][2])

    last_tUU = tUU
    length(tV) != 1 && push!(c, (U, tV) )
  end

  return reverse(c), tUU
end

function get_transversal_chain(U::PermGroup)
  # c = PermGroup[U]
  U0 = U
  C = Tuple{PermGroup, Vector{PermGroupElem}}[]
  while length(U0) != 1
    s = first(moved_points(U0))
    U = stabilizer(U0, s)[1]
    T = right_transversal(U0, U)
    push!(C, (U0, T) )
    U0 = U
  end
  push!(C, (U0, [one(U0)]))
  return C
end


function _get_rep(g::PermGroupElem, S::LadderStep)
  F = S.F

  dprev = one(F[1].A)
  uprev = dprev
  for s in F
    p = s.is_up_step ? dprev : (dprev * s.Tmap(g*inv(dprev*uprev)))
    (d, u) = s.m[p]
    dprev = d
    uprev = u*uprev
  end

  return dprev, uprev
end

function young_subgroup(p::Vector{T} ; full::T=sum(p)) where T<:Integer
  s = full - sum(p)
  # should error if s < 0
  s == 0 && return inner_direct_product([symmetric_group(i) for i in p])
  return inner_direct_product([symmetric_group(i) for i in vcat(p, ones(T, s))])
end

# Test method - groups only
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







function young_subgroup_ladder( p::Vector{T} ; full::T=sum(p)) where T<:Integer
  n = sum(p)
  L = LadderStep[]
  pp = copy(p)
  H = young_subgroup(pp ; full=full)
  s = n - pp[1]
  while s!=0
    if pp[2]==1

      Hprev = H
      deleteat!(pp, 2)
      pp[1]+=1
      H = young_subgroup(pp ; full=full)
      pushfirst!(L, young_subgroup(pp ; full=full))

    else

      Hprev = H
      pp[2]-=1
      insert!(pp,2,1)
      H = young_subgroup(pp ; full=full)
      pushfirst!(L, young_subgroup(pp ; full=full))


      Hprev = H
      deleteat!(pp, 2)
      pp[1]+=1
      H = young_subgroup(pp ; full=full)
      pushfirst!(L, young_subgroup(pp ; full=full))

    end
    s = n-pp[1]
  end

  return L
end
