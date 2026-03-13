################################################################################
# Ladder Game
#
#
#
################################################################################

module LadderGame

using Oscar


# Compute A\G/U
# we need: a subgroup ladder L = (A_i)_(0<=i<=a)
# Need a transversal chain for U:
#           subgroup chain (U_i)_(0<=i<=k)
#           transversals (S_i)_(0<i<=k) of R(U_i\U_i-1)
# Transversals: For each (A_i+1, A_i)
#           T_i in R(A_i\A_i+1) or R(A_i+1\A_i) as appropriate
# Test functions: For each i with A_i+1 <= A_i,
#           tau_i : A_i+1\A_i >-> NN (injection)
function ladder_game(L, G, U)
  r = ladder_start(G, U)
  H = L[1]
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

function ladder_start(G, U)
  st = [get_transversal_chain(U)]
  r = LadderData(
        # Ai = G
        # T = [ one(G) ]
        # mT =  ** map G ->G, x :-> one(G) **
        # I = [ one(G) ]
        # S3 = [ (one(G), one(G)) ]
        # Di = [ one(G) ]
        # St = st
        # g,s for env variable "F" set to false
  )
  # what about this environment variable "F" set to false?
  return r
end

function down_step(A, U, rec)
  r = LadderData(
        # Ai = A
        # Aim1 = rec.Ai
        # last = rec
        # g,s for ev variable "F" set to false
  )
  T, mT = transversal(rec.Ai, A)
  r.T = T
  r.mT = mT
  Di = []
  S3 = []
  St = []
  I = []
  for a in rec.Di
    tmp = []
    for t in T
      if ~(t in tmp)
        d = t*a
        st = rec.St[a]
        S1 = intersect(U, A^d)
        S2 = intersect(U, R.Ai^a)
        # assert S1 subset S2
        # assert st[1][1] == S2
        # assert S2 == intersect(U, R.Ai^d)
        st, Sta = _induce_chain(S1, mT, st ; conj=d)
        # assert length(Sta) == Index(S2, S1)
        _, ss = transversal(S2, S1)
        # assert length({ss(x) : x in Sta}) == #Sta
        append!(Di,d)
        append!(St, st)
        for u in Sta
          tt = mT(d*u*inv(a))
          include!(tmp, tt)
          include!(I, tt*a)
          append!(S3, (d,u) )
          # assert tt*a*inv(u)*inv(d) in A
        end
      end
    end
    # assert T == tmp
  end
  r.Di = Di
  r.I = I
  r.S3 = S3
  r.St = St
  return r
end

function up_step(A, U, rec)
  r = LadderData(
        # Ai = A
        # Aim1 = rec.Ai
        # last = rec
        # g,s for ev variable "F" set to false
  )

  T, mT = transversal(A, R.Ai)
  r.T = T
  r.mT = mT
  Di = []
  S3 = []
  St = []
  I = []

  seen = []
  data = []
  cls = []
  cls_data = []

  for a in rec.Di
    if a in cls
      p = index(cls, a)
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
      at, ut = _get_rep(ta, rec)
      # assert ta*inv(ut)*inv(at) in r.Aim1
      # assert at in rec.Di
      rep===nothing && (rep=at, rep_u = ut)
      at in cls || (include!(cls, at), append!(cls_data, (rep, rep_u*inv(ut)) ))
      at == rep && append!(tt, ut*inv(rep_u))
      is_one(t) && (include!(I, a), append!(S3, (rep, rep_u*inv(ut)) ))
    end
    if new
      append!(Di, rep)
      append!(St, [[ (intersect(A^rep, U), tt) ] rec.St[rep]] ) # in magma code: rec.St[Position(R.Di, rep)]
      # assert Index(St[end][1][1], St[end][2][1]) == length(tt)
    end
  end

  r.Di = Di
  r.I = I
  r.St = St
  r.S3 = S3

  return r
end

# I think this gives a transversal chain ???
# can streamline wrt tUU / last_tUU
function _induce_chain(V, Tm, C ; conj=nothing)
  # Tm is a MAP: transversal map for C[1][1]/V
  # assert V subgroup C[1][1]
  c = [ (subgroup(V, []), [one(V)]) ]
  conj===nothing && (conj = one(C[1][1]))
  last_tUU = [ one(C[end][1]) ]
  tUU = change_universe(last_tUU, C[1][1]) # - unnecessary???
  tU = [ Tm( one(universe(codomain(Tm)))^inv(conj) ) ] # isn't this just Tm(1)?
  for i in (length(C)-1):-1:1
    U = intersect(C[i][1], V)
    tUU = change_universe(last_tUU, C[i][1])
    tV = []
    change_universe!(last_tUU, C[i][1])
    for t in last_tUU, s in C[i][2]
      x = t*s
      if x in V
        append!(tV, x)
        length(tU)*length(tV) == length(last_tUU)*length(C[i][2]) && break
      else
        xx = Tm(x^inv(conj))
        if ~(xx in tU)
          include!(tU, xx)
          append!(tUU, x)
          length(tU)*length(tV) == length(tUU)*length(C[i][2]) && break
        end
      end
    end
    # assert length(tV) == index(U, c[end][1])
    # assert length(tU) == index(C[i][1], U)
    # assert length(tU)*length(tV) == length(last_tUU)*length(C[i][2])

    last_tUU = tUU
    length(tV) != 1 && append!(c, (U, tV) )
  end

  return reverse(c), tUU
end

function get_transversal_chain(U)
  c = [U]
  while length(U) != 1
    I = support(U)
    U = stabilizer(U, rep(I))
    append!(c, U)
  end
  C = []
  for i in 1:(length(c)-1)
    append!(C, (c[i], transversal(c[i], c[i+1])) )
  end
  append!(C, (c[end], [one(c[end])]) )

  return C
end

function get_rep(g, rec)
  F = rec.F
  if F!==nothing
    rec = F[end]
  else
    F = []
    RR = rec
    while is_assigned(rec.last)
      append!(F, rec)
      rec = rec.last
    end
    # check fuckery with "F"
  end

  dim1 = one(R.Ai)
  uim1 = dim1
  for i in length(F):-1:1
    F = F[i]
    if order(R.Ai) < order(R.Aim1)
      t = R.mT(g*inv(dim1*uim1))
      p = position(R.I, t*dim1)
    else
      p = position(R.I, dim1)
    end
    d = R.S3[p][1]
    u = R.S3[p][2]
    dim1 = d
    uim1 = u*uim1
  end

  return dim1, uim1
end


end
