function _get_quotient_split(P::Hecke.NfRelOrdIdl, i::Int)
  OE = order(P)
  E = nf(OE)
  Eabs, EabstoE = Hecke.absolute_simple_field(E, cached = false)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)

  function dlog(x::Hecke.NfRelElem)
    @assert parent(x) === E
    @assert is_integral(x)
    return muRPabs\(mRPabs(OEabs(EabstoE\x)))
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
    @assert is_integral(x)
    return mS\(mK\(mURPabs\(mRPabs(OEabs(EabstoE\x)))))
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
    @assert is_integral(x)
    return mS\(mK\(mURPabs\(mRPabs(OEabs(EabstoE\x)))))
  end

  return S, exp, dlog
end

function _get_quotient(P::Hecke.NfRelOrdIdl, i::Int)
  @assert is_prime(P)
  @assert is_maximal(order(P))
  OE = order(P)
  F = prime_decomposition(OE, minimum(P))
  if length(F) == 2
    S, dexp, dlog = _get_quotient_split(P, i)
  elseif F[1][2] == 1
    S, dexp, dlog = _get_quotient_inert(P, i)
  else
    S, dexp, dlog = _get_quotient_ramified(P, i)
  end
  return S, dexp, dlog
end

#function _get_big_quotient(E::)

