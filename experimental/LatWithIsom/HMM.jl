

###############################################################################
#
#  Local determinants morphism
#
###############################################################################

function _get_quotient_split(P::Hecke.NfRelOrdIdl, i::Int)
  OE = order(P)
  E = nf(OE)
  Eabs, EabstoE = Hecke.absolute_simple_field(E)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)

  function dlog(x::Hecke.NfRelElem)
    @assert parent(x) == E
    d = denominator(x)
    xabs = OEabs(d*EabstoE\(x))
    dabs = OEabs(d)
    F = prime_decomposition(OEabs, EabstoE\(minimum(P)))
    for PP in F
      @assert valuation(EasbtoE\(x), PP[1]) >= 0
      api = OEabs(anti_uniformizer(PP[1]))
      exp = valuation(OEabs(d), PP[1])
      dabs *= api^exp
      xabs *= api^exp
    end

    xabs_image = mURPabs\(mRPabs(xabs))
    dabs_image = mURPabs\(mRPabs(dabs))

    ret = xabs_image - dabs_image

    return ret
  end

  function exp(k::GrpAbFinGenElem)
    @assert parent(k) === URPabs
    x = EabstoE(Eabs(mRPabs\(mURPabs(k))))
    @assert dlog(x) == k
    return x
  end

  return URPabs, exp, dlog
end

function _get_quotient_inert(P::Hecke.NfRelOrdIdl, i::Int)
  OE = order(P)
  OK = base_ring(OE)
  E = nf(OE)
  K = base_field(E)
  p = minimum(P)

  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  Pabs = EabstoE\P
  OEabs = order(Pabs)
  Rp, mRp = quo(OK, p^i)
  URp, mURp = unit_group(Rp)

  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)
  f = hom(URPabs, URp, [mURp\(mRp(OK(norm(EabstoE(elem_in_nf(mRPabs\(mURPabs(URPabs[i])))))))) for i in 1:ngens(URPabs)])

  K, mK = kernel(f)

  S, mS = snf(K)

  function exp(k::GrpAbFinGenElem)
    @assert parent(k) === S
    return EabstoE(elem_in_nf(mRPabs\(mURPabs(mK(mS(k))))))
  end

  function dlog(x::Hecke.NfRelElem)
    @assert parent(x) === E
    d = denominator(x)
    xabs = OEabs(EabstoE\(d*x))
    dabs = OEabs(d)
    F = prime_decomposition(OEabs, EabstoE\(minimum(P)))
    for PP in F
      @assert valuation(EabstoE\(x), PP[1]) >= 0
      api = OEabs(anti_uniformizer(PP[1]))
      exp = valuation(OEabs(d), PP[1])
      dabs *= api^exp
      xabs *= api^exp
    end

    xabs_image = mURPabs\(mRPabs(xabs))
    dabs_image = mURPabs\(mRPabs(dabs))

    ret = mS\(mK\(xabs_image - dabs_image))

    return ret
  end

  return S, exp, dlog
end

function _get_quotient_ramified(P::Hecke.NfRelOrdIdl, i::Int)
  OE = order(P)
  E = nf(OE)
  p = minimum(P)
  e = valuation(different(OE), P)

  if i < e
    S = abelian_group(Int[])
    return S, x -> one(E), x -> id(S)
  end

  t = e-1

  psi(x) = x <= t ? x : t + 2*(x-t)

  pi = uniformizer(P)

  jj = t+1//2
  while ceil(psi(jj)) != i
    jj += 1//2
  end
  j = Int(ceil(jj))

  Pi = P^i

  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  OK = order(p)
  Rp, mRp = quo(OK, p^j)
  URp, mURp = unit_group(Rp)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)

  f = hom(URPabs, URp, [mURp\(mRp(OK(norm(EabstoE(elem_in_nf(mRPabs\(mURPabs(URPabs[i])))))))) for i in 1:ngens(URPabs)])

  K, mK = kernel(f)

  S, mS = snf(K)

  function exp(k::GrpAbFinGenElem)
    @assert parent(k) === S
    return EabstoE(elem_in_nf(mRPabs\(mURPabs(mK(mS(k))))))
  end

  function dlog(x::Hecke.NfRelElem)
    @assert parent(x) === E
    d = denomiator(x)
    xabs = OEabs(EasbtoE\(d*x))
    dabs = OEabs(d)
    F = prime_decomposition(OEabs, EabstoE\(minimum(P)))
    for PP in F
      @assert valuation(EabstoE\(x), PP[1]) >= 0
      api = OEabs(anti_uniformizer(PP[1]))
      exp = valuation(OEabs(d), PP[1])
      dabs *= api^exp
      xabs *= api^exp
    end

    xabs_image = mURPabs\(mRPabs(xabs))
    dabs_image = mURPabs\(mRPabs(dabs))

    ret = mS\(mK\ (xabs_image - dabs_image))
    return ret
  end

  return S, exp, dlog
end

function _get_quotient(O::Hecke.NfRelOrd, p::Hecke.NfOrdIdl, i::Int)
  @assert is_prime(p)
  @assert is_maximal(order(p))
  @assert order(p) === base_ring(O)
  F = prime_decomposition(O, p)
  P = F[1][1]
  if length(F) == 2
    S, dexp, dlog = _get_quotient_split(P, i)
  elseif F[1][2] == 1
    S, dexp, dlog = _get_quotient_inert(P, i)
  else
    S, dexp, dlog = _get_quotient_ramified(P, i)
  end
  return S, dexp, dlog, P
end

function _get_product_quotient(E::Hecke.NfRel, Fac)
  groups = []
  exps = []
  dlogs = []
  Ps = []

  if length(Fac) == 0
    A = abelian_group()
    function dlog_0(x::Hecke.NfRelElem); return one(A); end;
    function exp_0(x::GrpAbFinGenElem); return one(E); end;
    return A, dlog_0, exp_0
  end

  for i in 1:length(Fac)
    p, e = Fac[i]
    KK, f, g, P = _get_quotient(OE, p, e)
    push!(groups, KK)
    push!(exps, f)
    push!(dlogs, g)
    push!(Ps, P)
  end

  if length(groups) == 1
    return groups[1], dlogs[1], exps[1]
  end

  G, inj, proj = biproduct(groups)

  function dlog(x::Vector{Hecke.NfRelElem})
    return sum([inj[i](dlogs[i](x)) for i in 1:length(Fac)])
  end

  function exp(x::GrpAbFinGenElem)
    v = Hecke.NfRelOrdElem[OE(exps[i](proj[i](x))) for i in 1:length(Fac)]
    a = crt(v, [Ps[i]^(3*Fac[i][2]) for i in 1:length(Ps)])
    @assert dlog(a) == x
    return a
  end

  return G, dlog, exp
end

function _local_determinant_morphism(Lf::LatWithIsom)
  @assert is_of_hermitian_type(Lf)
  qL, fqL = discriminant_group(Lf)
  OqL = orthogonal_group(qL)
  G, _ = centralizer(OqL, fqL)
  H = hermitian_structure(Lf)
  l = get_attribute(Lf, :transfert_data)
  E = base_ring(H)
  OE = maximal_order(E)
  DKQ = different(base_ring(OE))*OE
  DEK = different(OE)
  DEQ = DEK*DKQ
  DinvH = inv(DEQ)*H
  res = Hecke.SpaceRes(ambient_space(Lf), ambient_space(H))
  Lv = trace_lattice(DinvH, res, l = l)
  @assert cover(q) === lattice(Lv)
  @assert ambient_isometry(Lv) == ambient_isometry(Lf)
  gene_herm = [_transfer_discriminant_isometries(res, g) for g in gens(G)]
  @assert all(m -> m*gram_matrix_of_rational_span(DinvH)*map_entries(involution(E), transpose(m)) == gram_matrix_of_rational_span(DinvH), gene_herm)

  S = elementary_divisors(DinvH, H)

  N = norm(L)

  Fsharpdata = Tuple{NfOrdIdl, Int}[]
  for p in S
    lp = prime_decomposition(OE, p)
    P = lp[1][1]
    k = valuation(N, p)
    a = valuation(DEQ, P)
    push!(Fsharpdata, (p, 2*k+a))
  end

  RmodFsharp, Fsharplog, Fsharpexp = _get_product_quotient(E, Fsharpdata)

  Fdata = Tuple{NfOrdIdl, Int}[]
  for p in S
    if !_is_special(H, p)
      continue
    end
    lp = prime_decomposition(OE, p)
    P = lp[1][1]
    e = valuation(DEK, P)
    push!(Fdata, (p, e))
  end

  RmodF, Flog, Fexp = _get_product_quotient(E, Fdata)

  A = [Flog(Fsharpexp(g)) for g in gens(RmodFsharp)]
  f = hom(gens(RmodFsharp), A)
  FmodFsharp, j = kernel(f)

  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  OEabs = maximal_order(Eabs)
  UOEabs, mUOEabs = unit_group(OEabs)
  UOK, mUOK = unit_group(base_ring(OE))

  fU = hom(UOEabs, UOK, [mUOK\(norm(OE(EabstoE(Eabs(mUOEabs(k)))))) for k in gens(UOEabs)])
  KU, jU = kernel(fU)

  gene_norm_one = Hecke.NfRelOrdElem[OE(EabstoE(Eabs(mUOEabs(jU(k))))) for k in gens(KU)]

  FOEmodFsharp, m = sub(RmodFsharp, )

end

function _is_special(L::Hecke.HermLat, p::Hecke.NfOrdIdl)
  OE = base_ring(L)
  E = nf(OE)
  lp = prime_decomposition(OE, p)
  if lp[1][2] != 2 || !is_even(rank(L))
    return false
  end
  _, _R, S = jordan_decomposition(L, p)
  R = [nrows(m) for m in _R]
  for r in R
    if r != 2
      return false
    end
  end
  P = lp[1][1]
  e = valuation(different(OE), P)
  for i in 1:length(S)
    if !iseven(e- S[i])
      return false
    end
  end
  u = uniformizer(P)
  s = involution(L)
  su = s(u)
  H = block_diagonal_matrix([matrix(E, 2, 2, [0 u^(S[i]); su^(S[i]) 0]) for i in 1:length(S)])
  return is_locally_isometric(L, hermitian_lattice(E, gram = H), p)
end


###############################################################################
#
#  Local hermitian lifting
#
###############################################################################

function _transfer_discriminant_isometries(res::Hecke.SpaceRes, g::AutomorphismGroupElem{TorQuadModule})
  E = base_ring(codomain(res))
  OE = maximal_order(OE)
  q = domain(g)
  @assert ambient_space(cover(q)) === domain(res)
  gE = zero_matrix(OE, 0, degree(H))
  vE = zeros(E, 1, degree(H))
  for i in 1:degree(H)
    vE[i] = one(E)
    vQ = res\vE
    vq = q(vQ)
    gvq = g(vq)
    gvQ = lift(gvq)
    gvE = res(gvQ)
    gE = vcat(gE, gvE)
    vE[i] = zero(E)
  end
  return gE
end

function _get_piecewise_actions_modular(H::Hecke.HermLat, g::MatrixElem{Hecke.NfRelOrdElem}, p::Hecke.NfOrdIdl, a::Int)
  @assert g*gram_matrix_of_rational_span(H)*map_entries(involution(H), transpose(g)) == gram_matrix_of_rational_span(H)
  E = base_field(H)
  OE = base_ring(H)
  @assert base_ring(g) === OE
  gE = map_entries(E, g)
  B, _, exp = jordan_decomposition(H, p)
  j = findfirst(i -> i == -a, exp)
  if j !== nothing
    popat!(B, j)
    popat!(exp, j)
  end
  local_act = typeof(g)[]
  for b in B
    gbE = solve(b, b*gE)
    gb = map_entries(OE, gb)
    push!(local_act, gb)
  end
  return exp, local_act
end

function _local_hermitian_lifting(G::MatrixElem{Hecke.NfRelElem}, F::MatrixElem{Hecke.NfRelOrdElem}, rho::Hecke.NfRelElem, l::Int, P::Hecke.NfRelOrdIdl; check = true)
  @assert trace(rho) == 1
  E = base_ring(G)
  @assert G == map_entries(involution(E), transpose(G))
  OE = base_ring(F)
  @assert nf(OE) === E
  OK = base_ring(OE)
  DEK = different(OE)
  DKQ = different(OK)*OE
  DEQ = DKQ*DEK
  e = valuation(rho, P)
  e = 1 - e
  @assert valuation(DEK, P) == e
  a = valuation(DEQ, P)
  FE = map_entries(E, F)
  RE = G - FE*G*map_entries(involution(E), transpose(FE))

  if check
    @assert min([valuation(v, P) for v in collect(inv(G)) if !iszero(v)]) >= l+a
    @assert min([valuation(v, P) for v in diagonal(inv(G)) if !iszero(v)]) >= e+a
    @assert min([valuation(v, P) for v in collect(R) if !iszero(v)]) >= l-a
    @assert min([valuation(v, P) for v in diagonal(R) if !iszero(v)]) >= l+e-1-a
  end

  UE = zero_matrix(E, nrows(RE), ncols(RE))
  for i in 1:nrows(RE)
    for j in 1:i-1
      UE[i, j] = RE[i, j]
    end
  end

  diagE = RE - UE - map_entries(involution(E), transpose(UE))

  newFE = FE + (UE + rho*diagE)*map_entries(involution(E), inv(transpose(FE)))*inv(G)
  newF = map_entries(OE, newFE)

  if check
    l2 = 2*l+1
    @assert min([valuation(v, P) for v in collect(F-newF) if !is_zero(v)]) >= l+1
    R2E = G-newFE*G*map_entries(involution(E), transpose(newFE))
    @assert min([valuation(v, P) for v in collect(R2E) if !iszero(v)]) >= l2-a
    @assert min([valuation(v, P) for v in diagonal(R2e) if !iszero(v)]) >= l2+e-1-a
  end

  return newF
end

function _find_rho(P::Hecke.NfRelOrdIdl)
  OE = order(P)
  E = nf(OE)
  dya = valuation(2, norm(P)) > 0
  !dya && return E(1//2)
  K = base_field(E)
  while true
    Et, t = E["t"]
    g = E(-rand(K, -5:5)^2-1)
    nu = 2*valuation(E(2), P)-2
    nug = valuation(g, P)
    if nu == nug
      d = denominator(g+1, OE)
      rt = roots(t^2 - (g+1)*d^2, max_roots = 1, ispure = true, is_normal=true)
      if !is_empty(rt)
        break
      end
    end
  end
  rho = (1+rt[1])//2
  @assert valuation(rho, P) == -1
  @assert trace(rho) == 1
  return rho
end

