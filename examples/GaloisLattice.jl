module GaloisLattice

using Oscar
import Oscar: AbstractAlgebra, Hecke, GaloisGrp.GaloisCtx

mutable struct SubField
  coeff_field::Union{Nothing, SubField}
  fld::AbstractAlgebra.Field # Q or a NumField

  grp::PermGroup #the fix group

  pe::SLPoly{fmpz} #, FlintIntegerRing} # the invariant evaluating to the PE
  conj::Vector{PermGroupElem}        # the absolute conjugates
  ts::fmpz_poly # tschirnhaus

  exact_den::Union{fmpq, NumFieldElem} # in this field! f'(alpha)
  dual_basis::Vector # symbolic: coeffs of f/t-pe
  basis::Vector # in fld, symbolic: [pe^i//exact_den for i=0:n-1]

  #Caches:
  num_basis::MatElem{qadic}
  num_dual_basis::Vector{Vector{qadic}}

  function SubField()
    return new()
  end
end

function Base.show(io::IO, S::SubField)
  print(io, "subfield for $(S.grp) by $(S.fld)")
end

Oscar.number_field(S::SubField) = S.fld, gen(S.fld)

#let G, C = galois_group(..)
#starting with QQ:
#  QQ = fixed_field(C, G)
#  pe = 1 (constant)
#  conj = [()]
#
# now suppose SubField is given (fixed_field of U) and V is a (maximal) subgroup
# the fixed_field(V) is a (minimal) extension of U, given via some invariant I
# (if I is U-relative V-inv, then I is "a" generic primtive element)
# goal is to write minpoly(I) as exact elements of SubField.
# conjugates (relative) of I are U//V operating on I
#  we might need a tschirni to make them different...
#  then we compute all conjugates (as slpoly), then the coeffs of the poly
#  are the elem. symm. in the roots.
#  those should be in the SubField, hence recursively written exactly:
#    in the subfields, we have the dual (kronnecker) basis as slpoly
#    conjugate, and sum -> trace are the coeff in the next subfield down
# the dual basis 
#I don't think this should be done symbolically.
#
#= next attempt
  K = k(b)         V
  |                |
  k = Q(a)         U
  |                |
  Q                G
  and hope recursion kicks in
  a = I() for a G-relative U-invariant
  (so minpoly a is easy: conjugates are I^t for t = G//U, eval and done)

  b = J() for a U-relative V-invar
  conjugates are J^s for s = U//V (the abs. conjugates should be G//V = {ts for t = G//U and s = U//V})
             (check right vs left)
  
  k has basis 1, a, ... a^n-1 and dual basis (from the poly) b_1, ... b_n however we need to 
  work with basis 1/f`(a) (1, a, ... a^n-1) and the scaled dual basis (the coeffs of f/(t-a) in the 
  correct order)
  call this basis c_i and keep b_i for the dual basis
  Then sum c_i Z contains the maximal order and for gamma = sum g_i c_i, g_i = trace(gamma b_i)

  => all conjugates of b are via the recursive cosets
  => if minpoly(a) is known and the conjugates with low precision, then all precision is possible
     if b_i, c_i are given as polys in a, theior conjugates are also known
     the abs. T2 (or all conjugates over Q) are also known => better bounds, no need for complicated
       slpolys do describe basis

  g(b) = 0 and C_i = 1/g`(b)(1, b, ..., b^l-1) with dual basis B_i
  f(a) = 0     c_i = 1/f`(a)(1, a, ..., a^n-1) with dual basis b_i
  => abs. basis is c_iC_j with dual basis b_iB_j:
    trace(c_i C_j b_k B_l) = trace(c_i b_k trace(C_j B_l))
      inner trace is 0,1, then outer trace is also 0,1
  this might work!!!


  Given the "primitive" element I as an SLPoly (or similar)
  t = G//U a transversal, then I^t are the "conjugates"
  find tschirni to make them pairwise different after eval
  c = U a transversal - this indexes the conjugates of the subfield
  (thus I^(tc) is all conjugates)
  compute bound on power sums (I^t)^i for 0<i<#t
    mult by bound on dual basis of subfield
  Compute all conjugates of I up to this precision
  Compute dual basis of subfield to this precision
  Compute the power sums (over t) for all c
  Mult be dual basis
  add (get the trace down to Z)
  isInt
  use to represent power sums as subfield elements via basis
  get poly
  compute basis, dual basis, bound on dual basis

  Sub-functions
  given SubField (via I, U)
    find poly
  given SubField with poly
    find basis, dual basis and bounds
  access basis (symbolic)
  access any elem at precision
  access dual_basis numerically at precision
    Think: precision can be increased via
      - eval I at precision
        cost: (#mult in I)*abs degree field
      - lifting the roots, given the poly
        need to lift abs deg subfield many polys
          each poly O^~(rel deg field) ops.
          => O^~(abs deg field) (or less)
        Plus: numerical poly: rel deg * abs deg subfield ^2 (using matrix)
                                      * O^~ abs deg subfield (using tree eval)
        if poly is over smaller field savings...  
      - given basis numerically and dual basis numerically at low
        precision, dual basis can be lifted as well
  access dual_basis bound
  given conjugate vector (and bound), find exact element


  K = k[t]/f, then basis 1, t, ..., t^n-1
  trace-dual: d_i/f' in k[t] wheren f/(t-alpha) = sum d_i(alpha) t^i
  thus t^(j-1)/f' is dual to d_i AND Z_K subset sum Z_k t^(j-1)/f'
=#

function rationals_as_subfield(C::GaloisCtx)
  S = SubField()
  S.grp = C.G # Q is fixed by the entire group
  S.exact_den = fmpq(1)
  I = SLPolyRing(ZZ, degree(C.f))
  S.pe = I(1)
  S.conj = [one(C.G)]
  S.fld = QQ
  S.dual_basis = [fmpq(1)]
  S.basis = [fmpq(1)]
  return S
end

function fixed_field(C::GaloisCtx, S::SubField, U::PermGroup; invar=nothing, max_prec::Int = typemax(Int))
  @assert is_subgroup(S.grp, U)[1]
  t = right_transversal(S.grp, U)
  @assert isone(t[1])
  if invar !== nothing
    PE = invar
  else
    PE = Oscar.GaloisGrp.invariant(S.grp, U)
  end

  rt = roots(C, 5)
  ts = Oscar.GaloisGrp.find_transformation(rt, PE, t)
  B1 = length(t)*Oscar.GaloisGrp.upper_bound(C, PE^(1+length(t)), ts)
  B2 = dual_basis_bound(S)
  B = B2*B1 #maybe dual_basis_bound should do bound_ring stuff?
  pr = Oscar.GaloisGrp.bound_to_precision(C, B)
  @show pr = min(pr, max_prec)
  rt = roots(C, pr)
  if ts != gen(parent(ts))
    rt = map(ts, rt)
  end
  ps = [[] for i=t]
  con = []
  dbc = dual_basis_conj(C, S, pr) #dbc[i] = array of the i-th conjugate of all dual basis elts
  num_basis = zero_matrix(parent(rt[1]), length(t), length(t)*length(S.conj))
  for i = 1:length(t)*length(S.conj)
    num_basis[1, i] = one(parent(rt[1]))
  end
  for i=1:length(S.conj)
    c = S.conj[i]
    p = [evaluate(PE, x*c, rt) for x = t]
    num_basis[2, (i-1)*length(t)+1:i*length(t)] = p
    pp = deepcopy(p)
    append!(con, [x*c for x = t])
    push!(ps[1], sum(p) .* dbc[i])
    for j=2:length(t)
      p .*= pp
      if j<length(t)
        num_basis[j+1, (i-1)*length(t)+1:i*length(t)] = p
      end
      push!(ps[j], sum(p) .* dbc[i])
    end
  end
  #so ps[i] is the i-th power sum
  #   ps[i][j]
  b = basis_abs(S)
  pp = [sum(b[j]*Oscar.GaloisGrp.isinteger(C, B, sum(y[i][j] for i=1:length(S.conj)))[2] for j=1:length(S.conj)) for y = ps]
  f = power_sums_to_polynomial(pp)

  SS = SubField()
  SS.fld, a = number_field(f, cached = false, check = false)
  SS.exact_den = derivative(f)(a)
  SS.basis = basis(SS.fld)
  KT, T = PolynomialRing(SS.fld, cached = false)
  SS.dual_basis = collect(coefficients(divexact(map_coefficients(SS.fld, f, parent = KT), T-a)))
  SS.coeff_field = S
  SS.conj = con
  SS.num_basis = num_basis
  SS.pe = PE
  SS.ts = ts
  SS.grp = U

  return SS
end

function refined_derived_series(G::PermGroup)
  s = GAP.Globals.PcSeries(GAP.Globals.Pcgs(G.X))
  return  Oscar._as_subgroups(G,s)
end

function fixed_field(C::GaloisCtx, s::Vector{PermGroup}; invar=nothing, max_prec::Int = typemax(Int))
  k = rationals_as_subfield(C)
  if order(s[1]) != order(C.G)
    st = 1
  else
    st = 2
  end
  j = 1
  for i=st:length(s)
    if invar !== nothing && i <= length(invar)
      k = fixed_field(C, k, s[i], invar = invar[i], max_prec = max_prec)
    else
      k = fixed_field(C, k, s[i], max_prec = max_prec)
    end
    k.fld.S = Symbol("a$j")
    j += 1
  end
  return k
end

function dual_basis_bound(S::SubField)
  if S.fld == QQ
    return fmpz(1)
  end
  return upper_bound(fmpz, maximum(x->maximum(abs, Oscar.conjugates(x)), S.dual_basis))*dual_basis_bound(S.coeff_field)
end

function Hecke.length(x::NumFieldElem,  abs_tol::Int = 32, T = arb)
  return sum(x^2 for x = Oscar.conjugates(x, abs_tol, T))
end

function conjugates(C::GaloisCtx, S::SubField, a::fmpq, pr::Int = 10)
  rt = roots(C, pr)
  return [parent(rt[1])(a)]
end

function recognise(C::GaloisCtx, S::SubField, I::SLPoly)
  return recognise(C, S, [I])[1]
end

function recognise(C::GaloisCtx, S::SubField, J::Vector{<:SLPoly})
  B = dual_basis_bound(S) * length(S.conj) *
            maximum(I->Oscar.GaloisGrp.upper_bound(C, I, S.ts), J)
  pr = Oscar.GaloisGrp.bound_to_precision(C, B)
  r = roots(C, pr)
  if S.ts != gen(parent(S.ts))
    r = map(S.ts, r)
  end
  b = basis_abs(S)
  db = dual_basis_conj(C, S, pr)

  D = []
  for I = J
    c = [evaluate(I, t, r) for t = S.conj] 
    d = zero(S.fld)
    for j=1:length(b)
      fl, v = Oscar.GaloisGrp.isinteger(C, B, sum(db[i][j] * c[i] for i=1:length(S.conj)))
      fl || return nothing
      d += v*b[j]
    end
    push!(D, d)
  end
  return D    
end

function as_radical_extension(K::NumField, aut::Map, zeta::NumFieldElem)
  g = gen(K)
  d = degree(K)
  #assumes K is cyclic, aut generates K/ceoff(K), zeta has order d
  #best to assume d is prime
  local r
  while true
    r = g
    z = one(K)
    for i=2:d
      z *= zeta
      g = aut(g)
      r += z*g
    end
    iszero(r) || break
    g = rand(K, -10:10)
  end
  s = coeff(r^d, 0)
  @assert s == r^d
  L, b = number_field(gen(parent(defining_polynomial(K)))^d-s, cached = false)
  @assert base_field(L) == base_field(K)
  return L, hom(L, K, r)
end

function solve(f::fmpz_poly; max_prec::Int=typemax(Int))
  @vprint :GaloisGroup 0 "computing initial galois group...\n"
  @vtime :GaloisGroup 0 G, C = galois_group(f)
  lp = [p for p = keys(factor(order(G)).fac) if p > 2]
  if length(lp) > 0
    @vprint :GaloisGroup 0 "need to add roots-of-one: $lp\n"
    @vtime :GaloisGroup 0 G, C = galois_group(f*prod(cyclotomic(Int(p), gen(parent(f))) for p = lp))
  end
  r = roots(C, 2)
  #the indices of zeta
  pp = [findfirst(isone, [x^p for x = r]) for p = lp]
  #and the indices of the roots of f
  rt = findall(iszero, map(f, r))
  s = PermGroup[]
  for i=1:length(lp)
    ap = findall(isone, [x^lp[i] for x = r])
    @assert length(ap) == lp[i]-1    
    push!(s, stabilizer(G, pp[1:i])[1])
  end

  s = vcat(s[1:end-1], refined_derived_series(s[end]))
  S = slpoly_ring(ZZ, degree(G))[1]
  @vprint :GaloisGroup 0 "computing tower...\n"
  @vtime :GaloisGroup 0 All = fixed_field(C, s, invar = gens(S)[pp], max_prec = max_prec)
                    #here one could actually specify the invariant
                          #at least for the cyclos

  fld_arr = [All]
  while fld_arr[1].fld !== QQ
    pushfirst!(fld_arr, fld_arr[1].coeff_field)
  end
  
  cyclo = fld_arr[length(pp)+1]
  @vprint :GaloisGroup 0 "finding roots-of-1...\n"
  @vprint :GaloisGroup 0 zeta = [recognise(C, cyclo, gens(parent(cyclo.pe))[i]) for i=pp]
  @assert all(i->isone(zeta[i]^lp[i]), 1:length(pp))
  aut = []
  @vprint :GaloisGroup 0 "finding automorphisms...\n"
  for i=length(pp)+2:length(fld_arr)
    @vprint :GaloisGroup 0 "..on level $(i-length(pp)-1)...\n"
    K = fld_arr[i]
    @vtime :GaloisGroup 0 push!(aut, hom(K.fld, K.fld, recognise(C, K, K.pe^K.conj[2])))
  end
  for i=1:length(pp)
    fld_arr[i+1].fld.S = Symbol("z_$(lp[i])")
  end
  @vprint :GaloisGroup 0 "find roots...\n"
  @vtime :GaloisGroup 0 R = recognise(C, All, gens(S)[rt])
  #now, rewrite as radicals..
  #the cyclos are fine:
  K = number_field(fld_arr[length(pp)+1])[1]
  i = length(pp)+2
  h = hom(K, K, gen(K))
  h_data = Any[gen(K)]
  h_data = []
  @vprint :GaloisGroup 0 "transforming to radical...\n"
  while i <= length(fld_arr)
    @vprint :GaloisGroup 0 "level $(length(h_data)+1)\n"
    L = number_field(fld_arr[i])[1]
    @assert domain(h) === base_field(L)
    @assert codomain(h) === K
    if degree(L) == 2
      f = defining_polynomial(L)
      f = map_coefficients(h, f)
      t = coeff(f, 1)
      if !iszero(t)
        x = gen(parent(f))
        t = divexact(t, 2)
        f = f(x-t)
        @assert iszero(coeff(f, 1))
        K = number_field(f, cached = false, check = false)[1]
        push!(h_data, gen(K)-t)
        h = hom(L, K, h_data...)
      else
        K_ = number_field(f, cached = false, check = false)[1]
        @assert base_field(K_) === K
        K = K_
        @assert base_field(L) === domain(h)
        @assert base_field(K) === codomain(h)
        push!(h_data, gen(K))
        h = hom(L, K, h_data...)
      end
    else
      @vtime :GaloisGroup 0 Ra, hh = as_radical_extension(L, aut[i-length(pp)-1], zeta[findfirst(isequal(degree(L)), lp)])
      #hh: new -> old

      @vtime :GaloisGroup 0 g = map_coefficients(h, parent(defining_polynomial(L))(preimage(hh, gen(L))))
      @vtime :GaloisGroup 0 K = number_field(map_coefficients(h, defining_polynomial(Ra)), cached = false, check = false)[1]
      push!(h_data, K(g))
      h = hom(L, K, h_data...)
    end
    K.S = L.S
    i += 1
  end
  
  return K, map(h, R)
end

function conjugates(C::GaloisCtx, S::SubField, a::NumFieldElem, pr::Int = 10)
  @assert parent(a) == S.fld
  if !isdefined(S, :num_basis) || precision(S.num_basis[1,1]) < pr
    rt = roots(C, pr)
    if isdefined(S, :ts)
      rt = map(S.ts, rt)
    end
    for i=1:length(S.conj)
      c = S.conj[i]
      p = evaluate(S.pe, c, rt)
      S.num_basis[2, i] = p
      for j=3:degree(S.fld)
        S.num_basis[j, i] = p*S.num_basis[2, i]
      end
    end
  end
  nb = S.num_basis
  K = base_ring(nb)
  coef = zero_matrix(K, 1, degree(S.fld))
  res = zero_matrix(K, 1, ncols(nb))
  tmp = zero_matrix(K, 1, ncols(nb))
  for i=0:degree(S.fld)-1
    d = conjugates(C, S.coeff_field, coeff(a, i), pr)
    for j=1:length(d)
      tmp[1, (j-1)*degree(S.fld)+1:j*degree(S.fld)] = d[j]*nb[i+1, (j-1)*degree(S.fld)+1:j*degree(S.fld)]
    end
    res += tmp
  end
  return res
end


function dual_basis_conj(C::GaloisCtx, S::SubField, pr::Int = 10)
  r = roots(C, pr)
  if S.fld == QQ
    return [[parent(r[1])(1)]]
  end
  dbc = []
  for b = dual_basis_abs(S)
    c = conjugates(C, S, b, pr)
    push!(dbc, c)
  end
  return [[dbc[i][j] for i = 1:length(dbc[1])] for j=1:length(dbc)]

end

function dual_basis_abs(S::SubField)
  if S.fld == QQ
    return [fmpq(1)]
  end
  d = dual_basis_abs(S.coeff_field)
  b = S.dual_basis
  return [i*j for j = d for i = b]
end

function basis_abs(S::SubField)
  if S.fld == QQ
    return [fmpq(1)]
  end
  d = basis_abs(S.coeff_field)
  b = S.basis .* inv(S.exact_den)
  return [i*j for j = d for i = b]
end

end # GaloisLattice
