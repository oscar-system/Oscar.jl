
###############################################################################
#
#  Computations of the finite quotients E_0/E^i: subsection 6.8. of BH23
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
    d = denominator(x, OE)
    xabs = d*(EabstoE\(x))
    dabs = copy(d)
    F = prime_decomposition(OEabs, minimum(Pabs))
    for PP in F
      @assert valuation(EabstoE\(x), PP[1]) >= 0
      api = anti_uniformizer(PP[1])
      exp = valuation(OEabs(d), PP[1])
      dabs *= api^exp
      xabs *= api^exp
    end

    xabs_image = mURPabs\mRPabs(OEabs(xabs))
    dabs_image = mURPabs\mRPabs(OEabs(dabs))

    ret = xabs_image - dabs_image

    return ret
  end

  function exp(k::GrpAbFinGenElem)
    @assert parent(k) === URPabs
    x = EabstoE(Eabs(mRPabs\mURPabs(k)))
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
  f = hom(URPabs, URp, [mURp\(mRp(OK(norm(EabstoE(elem_in_nf(mRPabs\(mURPabs(x)))))))) for x in gens(URPabs)])

  K, mK = kernel(f)

  S, mS = snf(K)

  function exp(k::GrpAbFinGenElem)
    @assert parent(k) === S
    return EabstoE(elem_in_nf(mRPabs\(mURPabs(mK(mS(k))))))
  end

  function dlog(x::Hecke.NfRelElem)
    @assert parent(x) === E
    d = denominator(x, OE)
    xabs = EabstoE\(d*x)
    dabs = copy(d)
    F = prime_decomposition(OEabs, minimum(Pabs))
    for PP in F
      @assert valuation(EabstoE\(x), PP[1]) >= 0
      api = anti_uniformizer(PP[1])
      exp = valuation(OEabs(d), PP[1])
      dabs *= api^exp
      xabs *= api^exp
    end

    xabs_image = mURPabs\(mRPabs(OEabs(xabs)))
    dabs_image = mURPabs\(mRPabs(OEabs(dabs)))

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

  psi(x) = t + 2*(x-t)

  pi = uniformizer(P)

  jj = t+1//2
  while ceil(psi(jj)) != i
    jj += 1//2
  end
  j = Int(ceil(jj))

  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  OK = order(p)
  Rp, mRp = quo(OK, p^j)
  URp, mURp = unit_group(Rp)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)

  f = hom(URPabs, URp, [mURp\(mRp(OK(norm(EabstoE(elem_in_nf(mRPabs\(mURPabs(x)))))))) for x in gens(URPabs)])

  K, mK = kernel(f)

  S, mS = snf(K)

  function exp(k::GrpAbFinGenElem)
    @assert parent(k) === S
    return EabstoE(elem_in_nf(mRPabs\(mURPabs(mK(mS(k))))))
  end

  function dlog(x::Hecke.NfRelElem)
    @assert parent(x) === E
    d = denominator(x, OE)
    xabs = EabstoE\(d*x)
    dabs = copy(d)
    F = prime_decomposition(OEabs, minimum(EabstoE\P))
    for PP in F
      @assert valuation(EabstoE\(x), PP[1]) >= 0
      api = anti_uniformizer(PP[1])
      exp = valuation(OEabs(d), PP[1])
      dabs *= api^exp
      xabs *= api^exp
    end

    xabs_image = mURPabs\(mRPabs(OEabs(xabs)))
    dabs_image = mURPabs\(mRPabs(OEabs(dabs)))

    ret = mS\(mK\ (xabs_image - dabs_image))
    return ret
  end

  return S, exp, dlog
end

function _get_quotient(O::Hecke.NfRelOrd, p::Hecke.NfOrdIdl, i::Int)
  @assert is_prime(p)
  @assert is_maximal(order(p))
  @assert order(p) === base_ring(O)
  E = nf(O)
  F = prime_decomposition(O, p)
  P = F[1][1]
  if i == 0
    A = abelian_group()
    function dlog_0(x::Hecke.NfRelElem); return id(A); end;
    function dexp_0(x::GrpAbFinGenElem); return one(E); end;
    return A, dexp_0, dlog_0, P
  end
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
  OE = maximal_order(E)
  groups = GrpAbFinGen[]
  exps = []
  dlogs = []
  Ps = []

  if length(Fac) == 0
    A = abelian_group()
    function dlog_0(x::Hecke.NfRelElem); return id(A); end;
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

  G, proj, inj = biproduct(groups...)

  function dlog(x::Vector)
    if length(x) == 1
      return sum([inj[i](dlogs[i](x[1])) for i in 1:length(Fac)]) 
    else
      @assert length(x) == length(Fac)
      return sum([inj[i](dlogs[i](x[i])) for i in 1:length(Fac)])
    end
  end

  function exp(x::GrpAbFinGenElem)
    v = Hecke.NfRelElem[exps[i](proj[i](x)) for i in 1:length(Fac)]
    @assert dlog(v) == x
    return v
  end

  for i in 1:10
    a = rand(G)
    @assert dlog(exp(a)) == a
  end

  return G, dlog, exp
end

###############################################################################
#
#  Local determinants morphism, alias \delta in BH23
#
###############################################################################

# We collect all the prime ideals p for which the local quotient
# (D^{-1}L^#/L)_p is not unimodular

function _elementary_divisors(L::HermLat, D::Hecke.NfRelOrdIdl)
  Ps = collect(keys(factor(D)))
  primess = NfOrdIdl[]
  for P in Ps
    p = minimum(P)
    ok, a = is_modular(L, p)
    ok && a == -valuation(D, P) && continue
    push!(primess, p)
  end
  for p in primes(genus(L))
    (p in primess) && continue
    push!(primess, p)
  end
  return primess
end

# We compute here the map delta from Theorem 6.15 of BH22. Its kernel is
# precisely the image of the map O(L, f) \to O(D_L, D_f).

function _local_determinants_morphism(Lf::LatWithIsom)
  @assert is_of_hermitian_type(Lf)

  qL, fqL = discriminant_group(Lf)
  OqL = orthogonal_group(qL)

  # Since any isometry of L centralizing f induces an isometry of qL centralising
  # fqL, G is the group where we want to compute the image of O(L, f). This
  # group G corresponds to U(D_L) in the notation of BH22.
  G, _ = centralizer(OqL, fqL)
  # This is the associated hermitian O_E-lattice to (L, f): we want to make qL
  # (aka D_L) correspond to the quotient D^{-1}H^#/H by the trace construction,
  # where D is the absolute different of the base algebra of H (a cyclotomic
  # field).
  H = hermitian_structure(Lf)

  E = base_field(H)
  OE = maximal_order(E)
  DKQ = different(base_ring(OE))*OE
  DEK = different(OE)
  DEQ = DEK*DKQ

  H2 = inv(DEQ)*dual(H)
  @assert is_sublattice(H2, H) # This should be true since the lattice in Lf is integral

  res = get_attribute(Lf, :transfer_data)

  # We want only the prime ideal in O_K which divides the quotient H2/H. For
  # this, we collect all the primes dividing DEQ or for which H is not locally
  # unimodular. Then, we check for which prime ideals p, the local quotient
  # (H2/H)_p is non trivial.
  S = _elementary_divisors(H, DEQ)

  N = norm(H)

  # We want to produce the product of the F(H)/F^#(L). For
  # this, we create the map between the alternative products R/F^#(L) \to R/F(L)
  # whose kernel is exactly what we want. Here R is just a big enough group.
  # Note that here the products can be constructed since there are only finitely
  # many primes in both cases for which the local quotients are non-trivial.

  Fsharpdata = Tuple{NfOrdIdl, Int}[]
  for p in S
    lp = prime_decomposition(OE, p)
    P = lp[1][1]
    n = valuation(N*OE, P)
    a = valuation(DEQ, P)
    push!(Fsharpdata, (p, n+a))
  end

  RmodFsharp, Fsharplog, Fsharpexp = _get_product_quotient(E, Fsharpdata)

  # Here thanks to results due to M. Kirschmer, some of the p's used for the
  # previous product of quotients might produce trivial factors. We can detect
  # them and this is the goal of the `_is_special` routine. For those particular
  # prime, we use the trivial group as factor
  #
  # Note: we do not remove the factor to be able to map the corresponding the
  # factors between the two products we construct. We do this componentwise to
  # avoid computing unncessary crt. This will hold for the rest of the code, we
  # for those particular objects, the `dlog` maps take vectors, corresponding to
  # finite adeles.
  Fdata = Tuple{NfOrdIdl, Int}[]
  for p in S
    if !_is_special(H, p)
      push!(Fdata, (p, 0))
    else
      lp = prime_decomposition(OE, p)
      P = lp[1][1]
      e = valuation(DEK, P)
      push!(Fdata, (p, e))
    end
  end

  RmodF, Flog, _ = _get_product_quotient(E, Fdata)

  A = [Flog(Fsharpexp(g)) for g in gens(RmodFsharp)]
  f = hom(gens(RmodFsharp), A)
  FmodFsharp, j = kernel(f)

  # Now according to Theorem 6.15 of BH22, it remains to quotient out the image
  # of the units in E of norm 1.
  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  OEabs = maximal_order(Eabs)
  UOEabs, mUOEabs = unit_group(OEabs)
  OK = base_ring(OE)
  UOK, mUOK = unit_group(OK)

  fU = hom(UOEabs, UOK, [mUOK\norm(OE(mUOK(m))) for m in gens(UOK)])
  KU, jU = kernel(fU)

  gene_norm_one = Hecke.NfRelElem[EabstoE(Eabs(mUOEabs(jU(k)))) for k in gens(KU)]

  # Now according to Theorem 6.15 of BH22, it remains to quotient out
  FOEmodFsharp, m = sub(RmodFsharp, elem_type(RmodFsharp)[Fsharplog([x for i in 1:length(S)]) for x in gene_norm_one])

  I = intersect(FOEmodFsharp, FmodFsharp)
  # Q is where the determinant of our lifts to good precision will live. So
  # we just need to create the map from G to Q.
  Q, mQ = quo(FmodFsharp, I)

  SQ, SQtoQ = snf(Q)

  function dlog(x::Vector)
    @assert length(x) == length(Fsharpdata)
    return SQtoQ\(mQ(j\(Fsharplog(x))))
  end

  imgs = elem_type(SQ)[]
  # For each of our matrices in gene_herm, we do successive P-adic liftings in
  # order to approximate an isometry of D^{-1}H^#, up to a certain precision
  # (given by Theorem 6.25 in BH22). We do this for all the primes we have to
  # consider up to now, and then map the corresponding determinant adeles inside
  # Q. Since our matrices were approximate lifts of the generators of G, we can
  # create the map we wanted from those data.
  for g in gens(G)
    ds = elem_type(E)[]
    for p in S
      lp = prime_decomposition(OE, p)
      P = lp[1][1]
      k = valuation(N, p)
      a = valuation(DEQ, P)
      e = valuation(DEK, P)
      g_approx = _approximate_isometry(H2, g, P, e, a, k, res)
      push!(ds, det(g_approx))
    end
    push!(imgs, dlog(ds))
  end

  GSQ, SQtoGSQ, _ = Oscar._isomorphic_gap_group(SQ)
  f = hom(G, GSQ, gens(G), SQtoGSQ.(imgs), check=false)

  return f
end

# We check whether for the prime ideal p E_O(L_p) != F(L_p).
function _is_special(L::Hecke.HermLat, p::Hecke.NfOrdIdl)
  OE = base_ring(L)
  E = nf(OE)
  lp = prime_decomposition(OE, p)
  if lp[1][2] != 2 || !Oscar.iseven(rank(L))
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
  u = elem_in_nf(uniformizer(P))
  s = involution(L)
  su = s(u)
  H = block_diagonal_matrix([matrix(E, 2, 2, [0 u^(S[i]); su^(S[i]) 0]) for i in 1:length(S)])
  return is_locally_isometric(L, hermitian_lattice(E, gram = H), p)
end


###############################################################################
#
#  Local hermitian lifting - path to algorithm 8 of BH23
#
###############################################################################

Oscar.canonical_unit(x::NfOrdQuoRingElem) = one(parent(x))

# Starting from an isometry of the torsion quadratic module `domain(g)`, for
# which we assume that the cover M has full rank, we compute a fake lift to the
# ambient space of the cover. This is not an isometry, but induces g on `domain(g)`.
#function _get_ambient_isometry(g::AutomorphismGroupElem{TorQuadModule})
#  q = domain(g)
#  L = cover(q)
#  @assert rank(L) == degree(L)
#  d = degree(L)
#  M1 = reduce(vcat, [matrix(QQ, 1, d, lift(t)) for t in gens(q)])
#  M2 = reduce(vcat, [matrix(QQ, 1, d, lift(g(t))) for t in gens(q)])
#  return solve(M1, M2)
#end

# Using the function used for the transfer construction, between the ambient
# space of the cover of our torsion module, and the ambient space of the
# corresponding hermitian structure, we transfer the fake lift computed with the
# previous function. This will be an invertible matrix, but the corresponding
# automorphism is not an isometry.
function _transfer_discriminant_isometry(res::AbstractSpaceRes, g::AutomorphismGroupElem{TorQuadModule}, Bp::T, P::Hecke.NfRelOrdIdl, i::Int) where T <: MatrixElem{Hecke.NfRelElem{nf_elem}}
  E = base_ring(codomain(res))
  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  Pabs = EabstoE\P
  OEabs = order(Pabs)
  q = domain(g)
  @assert ambient_space(cover(q)) === domain(res)
  #Bp = reduce(vcat, [matrix(E, 1, ncols(Bp), res(lift(q(res\vec(collect(Bp[i, :])))))) for i in 1:nrows(Bp)])
  B2 = zero(Bp)
  for i in 1:nrows(Bp)
    vE = vec(collect(Bp[i, :]))
    vQ = res\vE
    vq = q(vQ)
    gvq = g(vq)
    gvQ = lift(gvq)
    gvE = res(gvQ)
    B2[i, :] = gvE
  end

  Bpabs = map_entries(a -> EabstoE\a, Bp)
  B2abs = map_entries(a -> EabstoE\a, B2)
  d = lcm(lcm([denominator(a, OEabs) for a in Bpabs]), lcm([denominator(a, OEabs) for a in B2abs]))
  dBpabs = change_base_ring(OEabs, d*Bpabs)
  dB2abs = change_base_ring(OEabs, d*B2abs)
  Q, p = quo(OEabs, Pabs^(i+valuation(d, Pabs))) 
  BpQ = map_entries(p, dBpabs)
  B2Q = map_entries(p, dB2abs)
  println(BpQ, B2Q)
  KQ = solve_left(BpQ, B2Q)
  K = map_entries(a -> EabstoE(Eabs(p\a)), KQ)
  @assert _scale_valuation(K*Bp-B2, P) >= i
  println(K)
  return K
end

# the minimum P-valuation among all the non-zero entries of M
function _scale_valuation(M::T, P::Hecke.NfRelOrdIdl) where T <: MatrixElem{Hecke.NfRelElem{nf_elem}}
  @assert nf(order(P)) === base_ring(M)
  iszero(M) && return inf
  return minimum([valuation(v, P) for v in collect(M) if !iszero(v)])
end

# the minimum P-valuation among all the non-zero diagonal entries of M
function _norm_valuation(M::T, P::Hecke.NfRelOrdIdl) where T <: MatrixElem{Hecke.NfRelElem{nf_elem}}
  @assert nf(order(P)) === base_ring(M)
  iszero(diagonal(M)) && return inf
  r =  minimum([valuation(v, P) for v in diagonal(M) if !iszero(v)])
  return r
end

# This is algorithm 8 of BH22: under the good assumptions, then we can do a
# P-adic lifting of a matrix which represents an isometry up to a certain
# precision. In this way, we approximate our matrix by another matrix, to a
# given precision and the new matrix defines also an isometry up to a finer
# precision than the initial matrix.
#
# We use this method iteratively to lift isometries (along a surjective map), by
# looking at better representatives until we reach a good enough precision for
# our purpose.
function _local_hermitian_lifting(G::T, F::T, rho::Hecke.NfRelElem, l::Int, P::Hecke.NfRelOrdIdl, e::Int, a::Int; check = true) where T <: MatrixElem{Hecke.NfRelElem{nf_elem}}
  @assert trace(rho) == 1
  E = base_ring(G)
  # G here is a local gram matrix
  @assert G == map_entries(involution(E), transpose(G))
  @assert base_ring(F) === E

  # R represents the defect, how far F is to be an isometry of G
  R = G - F*G*map_entries(involution(E), transpose(F))
  # These are the necessary conditions for the input of algorithm 8 in BH22
  if check
    @assert _scale_valuation(inv(G), P) >= 1+a
    @assert _norm_valuation(inv(G), P) + valuation(rho, P) >= 1+a
    @assert _scale_valuation(R, P) >= l-a
    @assert _norm_valuation(R, P) + valuation(rho,P) >= l-a
  end

  # R is s-symmetric, where s is the canonical involution of E/K. We split R
  # into U + D + s(U), i.e we take U to be the strict upper triangular part of
  # R and D to be the diagonal.
  U = zero_matrix(E, nrows(R), ncols(R))
  for i in 1:nrows(R)
    for j in i+1:ncols(R)
      U[i, j] = R[i, j]
    end
  end

  diag = R - U - map_entries(involution(E), transpose(U))

  # this newF is suppose to be a better lift than F, i.e. it is congruent to F
  # modulo P^{l+1} and the corresponding defect R2 has a higher P-valuation (so
  # P-adic, we are close to have a proper isometry)
  newF = F + (U + rho*diag)*map_entries(involution(E), inv(transpose(F)))*inv(G)

  l2 = 2*l+1
  if check
    @assert _scale_valuation(F-newF, P) >= l+1
    R2 = G-newF*G*map_entries(involution(E), transpose(newF))
    @assert _scale_valuation(R2, P) >= l2-a
    @assert _norm_valuation(R2, P) + valuation(rho, P) >= l2-a
  end

  return newF, l2
end

# Starting from our fake isometry g on H (which will be here D^{-1}H^#), and a
# prime ideal P, we iteratively lift g to a finer fake isometry, i.e. a matrix
# defining a P-local isometry up to the precision P^{2*k+a}.
function _approximate_isometry(H::Hecke.HermLat, g::AutomorphismGroupElem{TorQuadModule}, P::Hecke.NfRelOrdIdl, e::Int, a::Int, k::Int, res::AbstractSpaceRes)
  E = base_field(H)
  @assert nf(order(P)) === E
  ok, b = is_modular(H, minimum(P))
  if ok && b == -a
    return identity_matrix(E, 1)
  end

  Bp, v = _local_basis_modular_submodules(H, minimum(P), a, res)
  Gp = Bp*gram_matrix(ambient_space(H))*map_entries(involution(E), transpose(Bp))
  Fp = _transfer_discriminant_isometry(res, g, Bp, P, v)
  # This is the local defect. By default, it should have scale P-valuations -a
  # and norm P-valuation e-1-a
  Rp = Gp - Fp*Gp*map_entries(involution(E), transpose(Fp))
  rho = _find_rho(P, e)

  l = 0
  while _scale_valuation(Rp, P) < 2*k+a
    Fp, l = _local_hermitian_lifting(Gp, Fp, rho, l, P, e, a)
    Rp = Gp - Fp*Gp*map_entries(involution(E), transpose(Fp))
  end

  return Fp
end

function _local_basis_modular_submodules(H::Hecke.HermLat, p::Hecke.NfOrdIdl, a::Int, res::AbstractSpaceRes)
  L = restrict_scalars(H, res)
  B, _ , exps = jordan_decomposition(H, p)
  if exps[end] == -a
    pop!(B)
    pop!(exps)
  end
  subs = eltype(B)[]
  for b in B
    H2 = lattice_in_same_ambient_space(H, b)
    if !is_sublattice(H, H2)
      L2 = restrict_scalars(H2, res)
      L2 = intersect(L, L2)
      B2 = basis_matrix(L2)
      gene = [res(vec(collect(B2[i, :]))) for i in 1:nrows(B2)]
      H2 = lattice(ambient_space(H), gene)
      b = local_basis_matrix(H2, p, type = :submodule)
    end
    push!(subs, b)
  end
  return reduce(vcat, subs), a-exps[end]
end

# We need a special rho for Algorithm 8 of BH23: we construct such an element
# here, which will be used to lift fake isometries up to a better precision.
function _find_rho(P::Hecke.NfRelOrdIdl, e)
  OE = order(P)
  E = nf(OE)
  lp = prime_decomposition(OE, minimum(P))
  dya = is_dyadic(P)
  !dya && return E(1//2)
  lp = prime_decomposition(OE, minimum(P))
  if lp[1][2] == 1
    return Hecke._special_unit(P, minimum(P))
  end
  K = base_field(E)
  Eabs, EabstoE = Hecke.absolute_simple_field(E)
  Pabs = EabstoE\P
  OEabs = order(Pabs)
  while true
    Eabst, t = Eabs["t"]
    g = EabstoE\(E(-rand(K, -5:5)^2-1))
    nu = 2*valuation(Eabs(2), Pabs)-2*e+2
    nug = valuation(g, Pabs)
    if nu == nug
      d = denominator(g+1, OEabs)
      rt = roots(t^2 - (g+1)*d^2, max_roots = 1, ispure = true, is_normal=true)
      if !is_empty(rt)
        rho = (1+rt[1])//2
        @assert valuation(rho, Pabs) == 1-e
        @assert trace(EabstoE(rho)) == 1
        return EabstoE(rho)
      end
    end
  end
end

