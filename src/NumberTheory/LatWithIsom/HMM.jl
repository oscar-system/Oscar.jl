function _get_quotient_split(P::Hecke.NfRelOrdIdl, i::Int)
  OE = order(P)
  E = nf(OE)
  Eabs, EabstoE = Hecke.absolute_simple_field(E, cached = false)

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

  Eabs, EabstoE = Hecke.absolute_simple_field(E, cached=false)
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

  Eabs, EabstoE = Hecke.absolute_simple_field(E, cached=false)
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

function _get_big_quotient(E::Hecke.NfRel, Fac::Vector{Tuple{Hecke.NfOrdIdl, Int}})
  OE = maximal_order(e)
  groups = []
  exps = []
  dlogs = []
  Ps = []

  if length(Fac) == 0
    A = abelian_group()
    function dlog(x::Hecke.NfRelElem); return one(A); end;
    function exp(x::GrpAbFinGenElem); return one(E); end;
    return A, dlog, exp
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

  function dlog(x::Hecke.NfRelElem)
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

