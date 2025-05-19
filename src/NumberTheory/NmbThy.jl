function norm_equation_fac_elem_non_max(R::AbsNumFieldOrder, k::ZZRingElem; abs::Bool = false)
  #Idea
  #every solution in R is also one in S, the maximal order, hence
  #associated to some solution in S
  #Step 1: solve completely in S
  #Step 2: compute X = S^*/R^* (finite)
  #Step 3: for each solution b find all x in X s.th. bx in R
  #
  # Improvement, sometimes: f = conductor(R, S)
  # then 1 -> R^* -> S^* -> (S/f)^*/(R/f)^* 
  # is exact (and used to compute R^*)
  # if k is coprime to f, then every solution b in S is also coprime to
  # f, hence in (S/f)^*/(R/f)^*
  # so instead of all X, just use all preimages in S^* of b
  #
  # finally, adjust for signs...
  #
  S = maximal_order(R)
  U, mU = unit_group_fac_elem(S)
  q, mqS, mSq = Hecke.OO_mod_F_mod_O_mod_F(R)
  s = norm_equation_fac_elem(S, k; abs)
  f = conductor(R, S)
  mu = hom(U, q, [preimage(mqS, mSq(mU(x))) for x = gens(U)])

  p, phi = quo(U, kernel(mu)[1])
  t = typeof(s)()
  have_unit = false
  if isodd(degree(R)) 
    u_m1 = FacElem(number_field(R)(-1))
    have_unit = true
  end
  for x = s
    if is_one(gcd(minimum(f), k))
      y = preimage(mqS, mSq(x))
      fl, d = has_preimage_with_preimage(mu, y)
      fl || continue
      u = mU(-d)
      if !abs && !have_unit && norm(u) == -1
        u_m1 = u
        have_unit = true
      end
      push!(t, x*u)
      @assert is_zero(preimage(mqS, mSq(t[end])))
    else
      for z = p
        u = mU(preimage(phi, z))
        if !abs && !have_unit && norm(u) == -1
          u_m1 = u
          have_unit = true
        end
        xx = u*x
        xxx = evaluate(xx)
        if denominator(xxx, R) == 1
          push!(t, xx)
        end
      end
    end
  end
  if !abs && !have_unit
    U, mU = unit_group(R)
    for g = gens(U)
      u = mU(g)
      if norm(u) == -1
        u_m1 = u
        have_unit = true
        break
      end
    end
  end
  tt = typeof(t)()
  for x = t
    if !abs && norm(x) != k
      if have_unit
        push!(tt, u_m1*x)
      end
    else
      push!(tt, x)
    end
  end

  return tt
end

function norm_equation_fac_elem(R::AbsNumFieldOrder, k::ZZRingElem; abs::Bool = false)
  if !Hecke.is_maximal(R)
    return norm_equation_fac_elem_non_max(R, k; abs)
  end
  lp = factor(k)
  S = Tuple{Vector{Tuple{Hecke.ideal_type(R), Int}}, Vector{ZZMatrix}}[]
  for (p, k) = lp.fac
    P = prime_decomposition(R, p)
    s = solve_non_negative(matrix(ZZ, 1, length(P), [degree(x[1]) for x = P]), matrix(ZZ, 1, 1, [k]))
    push!(S, (P, ZZMatrix[view(s, i:i, 1:ncols(s)) for i=1:nrows(s)]))
  end
  sol = FacElem{AbsSimpleNumFieldElem, AbsSimpleNumField}[]
  for x in Base.Iterators.ProductIterator(Tuple(t[2] for t = S))
    I = ideal(R, 1)
    for i = 1:length(S)
      I *= prod(Hecke.ideal_type(R)[S[i][1][j][1]^Int(x[i][j]) for j=1:length(S[i][1])])
    end
    fl, g = Hecke.is_principal_fac_elem(I::Hecke.ideal_type(R))
    if fl
      push!(sol, g)
    end
  end
  if !abs
    u, mu = unit_group_fac_elem(R)
    i = findfirst(x -> norm(mu(x)) == -1, gens(u))
    ns = QQFieldElem[norm(x) for x = sol]
    if i === nothing
      return [sol[i] for i in 1:length(sol) if ns[i] == k]
    end
    U = mu(u[i])
    return FacElem{AbsSimpleNumFieldElem, AbsSimpleNumField}[ ns[i] == k ? sol[i] : sol[i] * U for i = 1:length(sol)]
  end
  return sol
end

norm_equation_fac_elem(R::AbsNumFieldOrder, k::Base.Integer; abs::Bool = false) = 
                            norm_equation_fac_elem(R, ZZRingElem(k), abs = abs)

function norm_equation(R::AbsNumFieldOrder, k::ZZRingElem; abs::Bool = false)
  s = norm_equation_fac_elem(R, k, abs = abs)
  return elem_type(R)[R(evaluate(x)) for x = s]
end

norm_equation(R::AbsNumFieldOrder, k::Base.Integer; abs::Bool = false) = norm_equation(R, ZZRingElem(k), abs = abs)

function norm_equation_fac_elem(R::Hecke.RelNumFieldOrder{AbsSimpleNumFieldElem,Hecke.AbsSimpleNumFieldOrderFractionalIdeal}, a::AbsSimpleNumFieldOrderElem)

  @assert Hecke.is_maximal(R)
  Ka, mKa, mkK = collapse_top_layer(Hecke.nf(R))
  Ra = maximal_order(Ka)
  class_group(Ra)
  k = Hecke.nf(parent(a))
  class_group(parent(a))

  lp = factor(Ra(mkK(k(a)))*Ra)
  la = factor(a*parent(a))
  S, mS = Hecke.sunit_mod_units_group_fac_elem(collect(keys(lp)))
  s, ms = Hecke.sunit_mod_units_group_fac_elem(collect(keys(la)))
  No = hom(S, s, elem_type(s)[preimage(ms, norm(mkK, mS(g))) for g = gens(S)])
  q, mq = quo(S, kernel(No)[1])
  q, mms = snf(q)
  mq = mq*inv(mms)

  C = reduce(vcat, (matrix(ZZ, 1, ngens(q), [valuation(mS(preimage(mq, q[i])), p) for i=1:ngens(q)]) for p = keys(lp)))
  
  A = reduce(vcat, (matrix(ZZ, 1, ngens(q), [valuation(norm(mkK, mS(preimage(mq, g))), p) for g in gens(q)]) for p = keys(la)))
  b = matrix(ZZ, length(la), 1, [valuation(a, p) for p = keys(la)])

  so = solve_mixed(A, b, C)
  u, mu = Hecke.unit_group_fac_elem(parent(a))
  U, mU = Hecke.unit_group_fac_elem(Ra)
  No = hom(U, u, elem_type(u)[preimage(mu, norm(mkK, mU(g))) for g = gens(U)])
  sol = []
  for i = 1:nrows(so)
    aa = mS(preimage(mq, q(sub(so, i:i, 1:ncols(so)))))
    b = norm(mkK, aa)
    c = b*inv(FacElem(k(a)))
    d = preimage(mu, c)
    fl, p = has_preimage_with_preimage(No, d)
    if fl
      push!(sol, FacElem(Dict(mKa(x) => v for (x, v) = (aa*inv(mU(p))).fac)))
    end
  end

  return sol
end

norm_equation(R::Hecke.RelNumFieldOrder{AbsSimpleNumFieldElem,Hecke.AbsSimpleNumFieldOrderFractionalIdeal}, a::AbsSimpleNumFieldOrderElem) = map(x -> R(evaluate(x)), norm_equation_fac_elem(R, a))

function is_irreducible(a::AbsSimpleNumFieldOrderElem)
  if iszero(a)
    return false
  end
  O = parent(a)
  S = collect(keys(factor(a*O)))
  if length(S) == 0
    return false
  end
  s, ms = Hecke.sunit_mod_units_group_fac_elem(S)
  V = matrix(ZZ, [ZZRingElem[valuation(ms(x), y) for y = S] for x = gens(s)])
  b = matrix(ZZ, [ZZRingElem[valuation(a, y) for y = S]])
  sol = transpose(solve(V, b))

  #want to write sol = x+y where
  # Cx, Cy > 0
  # if this is possible, then a is not irreducible as a
  # is then ms(Ax) * ms(Ay) and neither is trivial.

  I = identity_matrix(ZZ, length(S))
  A = hcat(I, I)
  #so A*(x|y) = x+y = sol is the  1. condition
  C = cat(V, V, dims=(1,2))
  # C(x|y) >=0 iff Cx >=0 and Cy >=0
  #Cx <> 0 iff (1,..1)*Cx >= 1
  one = matrix(ZZ, 1, length(S), [1 for p = S])
  zer = matrix(ZZ, 1, length(S), [0 for p = S])
  C = vcat(C, hcat(one, zer), hcat(zer, one))
  d = matrix(ZZ, 2*length(S)+2, 1, [0 for i = 1:2*length(S) + 2])
  d[end-1, 1] = 1
  d[end, 1] = 1
  pt = solve_mixed(A, sol, C, d)
  return nrows(pt) == 0
end

@doc raw"""
    irreducibles(S::Vector{AbsSimpleNumFieldOrderIdeal}) -> Vector{AbsNumFieldOrderElem}

Return all irreducibles whose support is contained in $S$.
"""
function irreducibles(S::Vector{AbsSimpleNumFieldOrderIdeal})
  if length(S) == 0
    return []
  end
  @assert all(is_prime, S)
  #TODO: try to get a better bound - in general if S is too large
  #      it cannot work 
  
  O = order(S[1])
  @assert all(x-> order(x) == O, S)

  s, ms = Hecke.sunit_mod_units_group_fac_elem(S)
  if length(S) == 1
    return [O(evaluate(ms(s[1])))]
  end
  c, mc = class_group(O)

  V = transpose(matrix(ZZ, [[valuation(ms(x), y) for y = S] for x = gens(s)]))

  cone = cone_from_inequalities(-V)
  @assert is_pointed(cone) # otherwise the Hilbert basis is not unique
  hb = hilbert_basis(cone)
  res = [O(evaluate(ms(s(map(ZZRingElem, Array(v)))))) for v in hb]
  return res
end

#= From Lars:
using BenchmarkTools

function test(A,b, pref)
   Polymake.prefer(pref) do
      solve_non_negative(matrix(A), matrix(b))
   end
end

@benchmark test(A,b, "projection") seconds=30
@benchmark test(A,b, "libnormaliz") seconds=30
@benchmark test(A,b, "cdd") seconds=30
@benchmark test(A,b, "lrs") seconds=30
@benchmark test(A,b, "ppl") seconds=30
@benchmark test(A,b, "bbox") seconds=30
@benchmark test(A,b, "to") seconds=30

@benchmark test(AA,b, "projection") seconds=30
@benchmark test(AA,b, "libnormaliz") seconds=30
@benchmark test(AA,b, "cdd") seconds=30
@benchmark test(AA,b, "lrs") seconds=30
@benchmark test(AA,b, "ppl") seconds=30
@benchmark test(AA,b, "bbox") seconds=30
@benchmark test(AA,b, "to") seconds=30

Polymake.prefer("to") do
   solve_non_negative(matrix(A), matrix(b))
end

=#

@doc raw"""
    factorizations(a::AbsSimpleNumFieldOrderElem) -> Vector{Fac{OrdElem}}

Return all factorizations of $a$ into irreducibles.
"""
function factorizations(a::AbsSimpleNumFieldOrderElem)
  O = parent(a)
  S = collect(keys(factor(a*O)))
  if length(S) == 0
    return []
  end
  irr = irreducibles(S)
  b = transpose(matrix(ZZ, [ZZRingElem[valuation(a, y) for y = S]]))
  A = transpose(matrix(ZZ, [ZZRingElem[valuation(x, y) for y = S] for x = irr]))
  #solving Ax = b for x >=0  and A >=0 implies that columns of A with
  #entries > those in b can never be part of a solution
  #thus we prune them
  #of course, this could/ should be done in solve_non_negative...
  #hence in polymake
  i = 1
  while i<= ncols(A)
    if any(j->A[j, i] > b[j], 1:nrows(A))
      deleteat!(irr, i)
      if i==1
        A = A[:, 2:ncols(A)]
      elseif i == ncols(A)
        A = A[:, 1:ncols(A)-1]
      else
        A = hcat(A[:, 1:i-1], A[:, i+1:ncols(A)])
      end
    else
      i += 1
    end
  end
  sol = solve_non_negative(A, b)
  res = Fac{AbsSimpleNumFieldOrderElem}[]
  for j=1:nrows(sol)
    x = Dict{typeof(a), Int}()
    y = a
    for i=1:length(irr)
      if sol[j,i] != 0
        x[irr[i]] = sol[j,i]
        y = divexact(y, irr[i]^sol[j,i])
      end
    end
    push!(res, Fac(y, x))
  end
  return res
end
