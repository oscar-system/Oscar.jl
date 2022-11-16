export norm_equation_fac_elem, norm_equation, irreducibles, factorisations

function norm_equation_fac_elem(R::NfAbsOrd, k::fmpz; abs::Bool = false)
  @assert Hecke.is_maximal(R)
  lp = factor(k)
  S = Tuple{Vector{Tuple{Hecke.ideal_type(R), Int}}, Vector{fmpz_mat}}[]
  for (p, k) = lp.fac
    P = prime_decomposition(R, p)
    s = solve_non_negative(matrix(FlintZZ, 1, length(P), [degree(x[1]) for x = P]), matrix(FlintZZ, 1, 1, [k]))
    push!(S, (P, fmpz_mat[view(s, i:i, 1:ncols(s)) for i=1:nrows(s)]))
  end
  sol = FacElem{nf_elem, AnticNumberField}[]
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
    ns = fmpq[norm(x) for x = sol]
    if i === nothing
      return [sol[i] for i in 1:length(sol) if ns[i] == k]
    end
    U = mu(u[i])
    return FacElem{nf_elem, AnticNumberField}[ ns[i] == k ? sol[i] : sol[i] * U for i = 1:length(sol)]
  end
  return sol
end

norm_equation_fac_elem(R::NfAbsOrd, k::Base.Integer; abs::Bool = false) = 
                            norm_equation_fac_elem(R, fmpz(k), abs = abs)

function norm_equation(R::NfAbsOrd, k::fmpz; abs::Bool = false)
  s = norm_equation_fac_elem(R, k, abs = abs)
  return elem_type(R)[R(evaluate(x)) for x = s]
end

norm_equation(R::NfAbsOrd, k::Base.Integer; abs::Bool = false) = norm_equation(R, fmpz(k), abs = abs)

function norm_equation_fac_elem(R::Hecke.NfRelOrd{nf_elem,Hecke.NfOrdFracIdl}, a::NfAbsOrdElem{AnticNumberField,nf_elem})

  @assert Hecke.is_maximal(R)
  Ka, mKa, mkK = absolute_field(nf(R))
  Ra = maximal_order(Ka)
  class_group(Ra)
  k = nf(parent(a))
  class_group(parent(a))

  lp = factor(Ra(mkK(k(a)))*Ra)
  la = factor(a*parent(a))
  S, mS = Hecke.sunit_mod_units_group_fac_elem(collect(keys(lp)))
  s, ms = Hecke.sunit_mod_units_group_fac_elem(collect(keys(la)))
  No = hom(S, s, elem_type(s)[preimage(ms, norm(mkK, mS(g))) for g = gens(S)])
  q, mq = quo(S, kernel(No)[1])
  q, mms = snf(q)
  mq = mq*inv(mms)

  C = vcat([matrix(FlintZZ, 1, ngens(q), [valuation(mS(preimage(mq, q[i])), p) for i=1:ngens(q)]) for p = keys(lp)])
  
  A = vcat([matrix(FlintZZ, 1, ngens(q), [valuation(norm(mkK, mS(preimage(mq, g))), p) for g in gens(q)]) for p = keys(la)])
  b = matrix(FlintZZ, length(la), 1, [valuation(a, p) for p = keys(la)])

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
    fl, p = haspreimage(No, d)
    if fl
      push!(sol, FacElem(Dict(mKa(x) => v for (x, v) = (aa*inv(mU(p))).fac)))
    end
  end

  return sol
end

norm_equation(R::Hecke.NfRelOrd{nf_elem,Hecke.NfOrdFracIdl}, a::NfAbsOrdElem{AnticNumberField,nf_elem}) = map(x -> R(evaluate(x)), norm_equation_fac_elem(R, a))

function is_irreducible(a::NfAbsOrdElem{AnticNumberField,nf_elem})
  if iszero(a)
    return false
  end
  O = parent(a)
  S = collect(keys(factor(a*O)))
  if length(S) == 0
    return false
  end
  s, ms = Hecke.sunit_mod_units_group_fac_elem(S)
  V = transpose(matrix(ZZ, [fmpz[valuation(ms(x), y) for y = S] for x = gens(s)]))
  b = transpose(matrix(ZZ, [fmpz[valuation(a, y) for y = S]]))
  sol = solve(V, b)

  #want to write sol = x+y where
  # Cx, Cy > 0
  # if this is possible, then a is not irreducible as a
  # is then ms(Ax) * ms(Ay) and neither is trivial.

  I = identity_matrix(FlintZZ, length(S))
  A = hcat(I, I)
  #so A*(x|y) = x+y = sol is the  1. condition
  C = cat(V, V, dims=(1,2))
  # C(x|y) >=0 iff Cx >=0 and Cy >=0
  #Cx <> 0 iff (1,..1)*Cx >= 1
  one = matrix(FlintZZ, 1, length(S), [1 for p = S])
  zer = matrix(FlintZZ, 1, length(S), [0 for p = S])
  C = vcat(C, hcat(one, zer), hcat(zer, one))
  d = matrix(FlintZZ, 2*length(S)+2, 1, [0 for i = 1:2*length(S) + 2])
  d[end-1, 1] = 1
  d[end, 1] = 1
  pt = solve_mixed(A, sol, C, d)
  return nrows(pt) == 0
end

@doc Markdown.doc"""
    irreducibles(S::Vector{NfAbsOrdIdl{AnticNumberField,nf_elem}}) -> Vector{NfAbsOrdElem}

Return all irreducibles whose support is contained in $S$.
"""
function irreducibles(S::Vector{NfAbsOrdIdl{AnticNumberField,nf_elem}})
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
  res = [O(evaluate(ms(s(map(fmpz, Array(v)))))) for v in hb]
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

@doc Markdown.doc"""
    factorisations(a::NfAbsOrdElem{AnticNumberField,nf_elem}) -> Vector{Fac{OrdElem}}

Return all factorisations of $a$ into irreducibles.
"""
function factorisations(a::NfAbsOrdElem{AnticNumberField,nf_elem})
  O = parent(a)
  S = collect(keys(factor(a*O)))
  if length(S) == 0
    return []
  end
  irr = irreducibles(S)
  b = transpose(matrix(ZZ, [fmpz[valuation(a, y) for y = S]]))
  A = transpose(matrix(ZZ, [fmpz[valuation(x, y) for y = S] for x = irr]))
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
  res = Fac{NfAbsOrdElem{AnticNumberField,nf_elem}}[]
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


################################################################################
#
#   disc_log
#
@doc Markdown.doc"""
    disc_log(b::T, x::T) where {T <: FinFieldElem}

Return an integer `s` such that $b^s = x$.
If no such `x` exists, an exception is thrown.

# Examples
```jldoctest
julia> F = GF(3,4); a = gen(F)^21;

julia> disc_log(gen(F), a)
21
```
"""
function disc_log(b::T, x::T) where {T <: FinFieldElem}
  @assert parent(b) === parent(x)
  return disc_log_bs_gs(b, x, order(parent(b)))
end
