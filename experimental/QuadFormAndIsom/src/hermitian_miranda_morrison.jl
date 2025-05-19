###############################################################################
#
#  Computations of the finite quotients E_0/E^i: subsection 6.8. of [BH23]
#  ============================================================================
#
#  The 5 following functions are an import of identical functions written by
#  Tommy Hofmann on Magma.
#
###############################################################################

# Split case

function _get_quotient_split(P::Hecke.RelNumFieldOrderIdeal, i::Int)
  OE = order(P)
  E = number_field(OE)
  Eabs, EabstoE = absolute_simple_field(E)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)

  function dlog(x::Hecke.RelSimpleNumFieldElem)
    @hassert :ZZLatWithIsom 1 parent(x) == E
    d = denominator(x, OE)
    xabs = d*(EabstoE\(x))
    dabs = copy(d)
    F = prime_decomposition(OEabs, minimum(Pabs))
    for PP in F
      @hassert :ZZLatWithIsom 1 valuation(EabstoE\(x), PP[1]) >= 0
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

  function exp(k::FinGenAbGroupElem)
    @hassert :ZZLatWithIsom 1 parent(k) === URPabs
    x = EabstoE(Eabs(mRPabs\mURPabs(k)))
    @hassert :ZZLatWithIsom 1 dlog(x) == k
    return x
  end

  return URPabs, exp, dlog
end

# Inert case

function _get_quotient_inert(P::Hecke.RelNumFieldOrderIdeal, i::Int)
  OE = order(P)
  OK = base_ring(OE)
  E = number_field(OE)
  K = base_field(E)
  p = minimum(P)

  Eabs, EabstoE = absolute_simple_field(E)
  Pabs = EabstoE\P
  OEabs = order(Pabs)
  Rp, mRp = quo(OK, p^i)
  URp, mURp = unit_group(Rp)

  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)
  f = hom(URPabs, URp, elem_type(URp)[mURp\(mRp(OK(norm(EabstoE(elem_in_nf(mRPabs\(mURPabs(x)))))))) for x in gens(URPabs)])

  K, mK = kernel(f)

  S, mS = snf(K)

  function exp(k::FinGenAbGroupElem)
    @hassert :ZZLatWithIsom 1 parent(k) === S
    return EabstoE(elem_in_nf(mRPabs\(mURPabs(mK(mS(k))))))
  end

  function dlog(x::Hecke.RelSimpleNumFieldElem)
    @hassert :ZZLatWithIsom 1 parent(x) === E
    d = denominator(x, OE)
    xabs = EabstoE\(d*x)
    dabs = copy(d)
    F = prime_decomposition(OEabs, minimum(Pabs))
    for PP in F
      @hassert :ZZLatWithIsom 1 valuation(EabstoE\(x), PP[1]) >= 0
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

# Ramified case

function _get_quotient_ramified(P::Hecke.RelNumFieldOrderIdeal, i::Int)
  OE = order(P)
  E = number_field(OE)
  p = minimum(P)
  e = valuation(different(OE), P)

  if i < e
    S = abelian_group()
    return S, _ -> one(E), _ -> id(S)
  end

  t = e-1

  psi(x) = t + 2*(x-t)

  pi = uniformizer(P)

  jj = t+1//2
  while ceil(psi(jj)) != i
    jj += 1//2
  end
  j = Int(ceil(jj))

  Eabs, EabstoE = absolute_simple_field(E)
  OK = order(p)
  Rp, mRp = quo(OK, p^j)
  URp, mURp = unit_group(Rp)

  Pabs = EabstoE\P
  OEabs = order(Pabs)
  RPabs, mRPabs = quo(OEabs, Pabs^i)
  URPabs, mURPabs = unit_group(RPabs)

  f = hom(URPabs, URp, elem_type(URp)[mURp\(mRp(OK(norm(EabstoE(elem_in_nf(mRPabs\(mURPabs(x)))))))) for x in gens(URPabs)])

  K, mK = kernel(f)

  S, mS = snf(K)

  function exp(k::FinGenAbGroupElem)
    @hassert :ZZLatWithIsom 1 parent(k) === S
    return EabstoE(elem_in_nf(mRPabs\(mURPabs(mK(mS(k))))))
  end

  function dlog(x::Hecke.RelSimpleNumFieldElem)
    @hassert :ZZLatWithIsom 1 parent(x) === E
    d = denominator(x, OE)
    xabs = EabstoE\(d*x)
    dabs = copy(d)

    F = factor(EabstoE\(minimum(P)*OE))
    for PP in F
      @hassert :ZZLatWithIsom 1 valuation(EabstoE\(x), PP[1]) >= 0
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

# We obtain the quotient here: we first check what the ramification type of p
# in O is, in order to distribute to the appropriate function above.
function _get_quotient(
    O::Hecke.RelNumFieldOrder,
    p::Hecke.AbsSimpleNumFieldOrderIdeal,
    i::Int
  )
  @hassert :ZZLatWithIsom 1 is_prime(p)
  @hassert :ZZLatWithIsom 1 is_maximal(order(p))
  @hassert :ZZLatWithIsom 1 order(p) === base_ring(O)
  E = number_field(O)
  F = prime_decomposition(O, p)
  P = F[1][1]
  if i == 0
    A = abelian_group()
    function dlog_0(x::Hecke.RelSimpleNumFieldElem); return id(A); end;
    function dexp_0(x::FinGenAbGroupElem); return one(E); end;
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

# This is the product quotient from Remark 6.16 in [BH23], the finite
# abelian group where one do the local determinants computations
# for the hermitian version of Miranda-Morrison theory

function _get_product_quotient(
    E::Hecke.RelSimpleNumField, 
    Fac::Vector{Tuple{AbsSimpleNumFieldOrderIdeal, Int}}
  )
  OE = maximal_order(E)
  groups = FinGenAbGroup[]
  exps = Function[]
  dlogs = Function[]
  Ps = ideal_type(OE)[]

  if length(Fac) == 0
    A = abelian_group()
    function dlog_0(x::Vector{<:Hecke.RelSimpleNumFieldElem}); return id(A); end;
    function exp_0(x::FinGenAbGroupElem); return one(E); end;
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

  G, proj, inj = biproduct(groups...)

  function dlog(x::Vector{<:Hecke.RelSimpleNumFieldElem})
    if length(x) == 1
      return sum(inj[i](dlogs[i](x[1])) for i in 1:length(Fac))
    else
      @hassert :ZZLatWithIsom 1 length(x) == length(Fac)
      return sum(inj[i](dlogs[i](x[i])) for i in 1:length(Fac))
    end
  end

  function exp(x::FinGenAbGroupElem)
    v = elem_type(E)[exps[i](proj[i](x)) for i in 1:length(Fac)]
    @hassert :ZZLatWithIsom 1 dlog(v) == x
    return v
  end

  @v_do :ZZLatWithIsom 3 for i in 1:10
      a = rand(G)
      @hassert :ZZLatWithIsom 1 dlog(exp(a)) == a
    end

  return G, dlog, exp
end

###############################################################################
#
#  Local determinants morphism, alias \delta in [BH23]
#
###############################################################################

# Given a hermitian lattice $L$ over a degree 2 extensions of number fields
# $E/K$, and given an $O_E$-ideal $D$, return the prime $O_K$-ideals $p$ so that
# the quotient $D^{-1}L^\#_p/L_p$ is not unimodular.
#
# According to [BH23], such a quotient is unimodular if and only if
# * either $L_p$ is unimodular, and $D$ and $p$ are coprime;
# * or $L$ is $P^{-a}$-modular where $P$ is the largest prime $O_E$-ideal over
#   $p$ which is fixed by the canonical involution, and $a = \val_P(D)$.
function _elementary_divisors(L::HermLat, D::Hecke.RelNumFieldOrderIdeal)
  Ps = collect(keys(factor(D)))
  primess = AbsSimpleNumFieldOrderIdeal[]
  for P in Ps
    p = minimum(P)
    ok, a = is_modular(L, p)
    ok && a == -valuation(D, P) && continue
    push!(primess, p)
  end
  minPs = minimum.(Ps)
  for p in primes(genus(L))
    ((p in primess) || (p in minPs)) && continue
    push!(primess, p)
  end
  return unique!(primess)
end

# For each generator g, we find a good enough approximation on the ambient space
# of L restricting to g, which we transport by the trace construction on the
# hermitian structure of (L, f), and we create a better approximation of the
# new isometry of D^{-1}L^#/L on the hermitian ambient space of L, up to the
# valuation given by Theorem 6.25 in [BH23].
#
# This approximation is obtained by applying successive hermitian hensel
# lifting as described in Algorithm 8 of [BH23].
#
# Once this new approximation obtained, we compute its determinant which we map
# in the big product quotient constructed according to BH23, Section 6.
#
# The major part of the following algorithm follows the implementation of a
# similar function on Magma by Tommy Hofmann. Only the last loop about
# determinants approximations is new in this code.
#
# Given an even indefinite integer lattice with isometry $(L, f)$ of finite
# hermitian type, with $L$ of rank at least 3, compute the map $\delta$
# defined in Theorem 6.15 of [BH23]. Its kernel is precisely the image of the
# discriminant representation $O(L, f)\to O(D_L, D_f)$. The first output is
# $\delta$ seen as a map from $O(D_L, D_f)$ to a suitably chosen finite abelian
# group, and the second output is the embedding $O(D_L, D_f)\to O(D_L)$.
function _local_determinants_morphism(Lf::ZZLatWithIsom)
  @hassert :ZZLatWithIsom 1 is_of_hermitian_type(Lf)

  # We want to compute the image of the centralizer as a subgroup of OqL. For
  # this, if Lf is not full, we need to consider an isomorphic pair Lf2 of full
  # rank and then transport the generators along an appropriate map.
  qL, fqL = discriminant_group(Lf)
  OqL = orthogonal_group(qL)
  fqL = OqL(hom(fqL); check=false)

  if rank(Lf) != degree(Lf)
    Lf2 = integer_lattice_with_isometry(integer_lattice(; gram=gram_matrix(Lf)), isometry(Lf); ambient_representation=false, check=false)
    qL2, fqL2 = discriminant_group(Lf2)
    OqL2 = orthogonal_group(qL2)
    phi12 = hom(qL, qL2, identity_matrix(ZZ, ngens(qL)))
    @hassert :ZZLatWithIsom 1 is_isometry(phi12)
  else
    Lf2 = Lf
    qL2, fqL2, OqL2 = qL, fqL, OqL
    phi12 = id_hom(qL)
  end

  # Since any isometry of L centralizing f induces an isometry of qL centralising
  # fqL, G is the group where we want to compute the image of O(L, f). This
  # group G corresponds to U(D_L) in the notation of [BH23].
  G2, _ = centralizer(OqL2, fqL2)
  gensG2 = gens(G2)
  gensG = elem_type(OqL)[OqL(compose(phi12, compose(hom(g), inv(phi12))); check=false) for g in gensG2]
  G, GinOqL = sub(OqL, gensG)
  @hassert :ZZLatWithIsom 1 fqL in G
  GtoG2 = hom(G, G2, gensG, gensG2; check=false)
  @hassert :ZZLatWithIsom 1 GtoG2(fqL) == fqL2

  # This is the associated hermitian O_E-lattice to (L, f): we want to make qL
  # (aka D_L) correspond to the quotient D^{-1}H^#/H by the trace construction,
  # where D is the absolute different of the base algebra of H (a cyclotomic
  # field).
  H = hermitian_structure(Lf2)

  E = base_field(H)
  OE = maximal_order(E)
  DKQ = different(base_ring(OE))*OE
  DEK = different(OE)
  DEQ = DEK*DKQ

  H2 = inv(DEQ)*dual(H)
  @hassert :ZZLatWithIsom 1 is_sublattice(H2, H) # This should be true since the lattice in Lf is integral

  # This is the map used for the trace construction: it is stored on Lf when we
  # have constructed H. We need this map because it sets the rule of the
  # transfer between the quadratic world in which we study (L, f) and the
  # hermitian world in which H lives. In particular, the trace lattice of H2
  # with respect to res will be exactly the dual of L.
  res = get_attribute(Lf2, :transfer_data)

  # We want only the prime ideal in O_K which divides the quotient H2/H. For
  # this, we collect all the primes dividing DEQ or for which H is not locally
  # unimodular. Then, we check for which prime ideals p, the local quotient
  # (H2/H)_p is non trivial.
  S = _elementary_divisors(H, DEQ)

  N = norm(H)

  # We want to produce the product of the F(H_p)/F^#(H_p) for the p's in S. For
  # this, we create the map between the alternative products R/F^#(H_p) \to R/F(H_p)
  # whose kernel is exactly what we want. Here R is just a big enough group.
  # Note that here the products can be constructed since there are only finitely
  # many primes in both cases for which the local quotients are nontrivial.

  Fsharpdata = Tuple{AbsSimpleNumFieldOrderIdeal, Int}[]
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
  # Note: we do not remove the factor to be able to map the corresponding
  # factors between the two products we construct. We do this componentwise to
  # avoid computing unnecessary crt. This will hold for the rest of the code:
  # for those particular objects, the `dlog` maps take vectors, corresponding to
  # finite adeles.
  list_ker = eltype(S)[]
  Fdata = Tuple{AbsSimpleNumFieldOrderIdeal, Int}[]
  for _i in 1:length(S)
    p = S[_i]
    if !_is_special(H, p)
      push!(Fdata, (p, 0))
      if Fsharpdata[_i][2] == 0
        push!(list_ker, p)
      end
    else
      lp = prime_decomposition(OE, p)
      P = lp[1][1]
      e = valuation(DEK, P)
      push!(Fdata, (p, e))
      if Fsharpdata[_i][2] == e
        push!(list_ker, p)
      end
    end
  end

  RmodF, Flog, _ = _get_product_quotient(E, Fdata)

  A = elem_type(RmodF)[Flog(Fsharpexp(g)) for g in gens(RmodFsharp)]
  f = hom(RmodFsharp, RmodF, A; check=false)
  FmodFsharp, j = kernel(f)

  # Now according to Theorem 6.15 of [BH23], it remains to quotient out the image
  # of the units in E of norm 1.
  Eabs, EabstoE = absolute_simple_field(E)
  OEabs = maximal_order(Eabs)
  UOEabs, mUOEabs = unit_group(OEabs)
  OK = base_ring(OE)
  UOK, mUOK = unit_group(OK)

  fU = hom(UOEabs, UOK, elem_type(UOK)[mUOK\norm(OE(mUOK(m))) for m in gens(UOK)]; check=false)
  KU, jU = kernel(fU)

  gene_norm_one = elem_type(E)[EabstoE(Eabs(mUOEabs(jU(k)))) for k in gens(KU)]

  FOEmodFsharp, m = sub(RmodFsharp, elem_type(RmodFsharp)[Fsharplog(typeof(x)[x for i in 1:length(S)]) for x in gene_norm_one])

  I = intersect(FOEmodFsharp, FmodFsharp)
  # Q is where the determinant of our lifts to good precision will live. So
  # we just need to create the map from G to Q.
  Q, mQ = quo(FmodFsharp, I)

  SQ, SQtoQ = snf(Q)
  oSQ = order(SQ)

  function dlog(x::Vector)
    @hassert :ZZLatWithIsom 1 length(x) == length(Fsharpdata)
    return SQtoQ\(mQ(j\(Fsharplog(x))))
  end

  imgs = elem_type(SQ)[]
  # For each of our matrices in gensG2, we do successive P-adic liftings in
  # order to approximate an isometry of D^{-1}H^#, up to a certain precision
  # (given by Theorem 6.25 in [BH23]). We do this for all the primes we have to
  # consider up to now, and then map the corresponding determinant adeles inside
  # Q.
  for g in gensG2
    ds = elem_type(E)[]
    for p in S
      if isone(oSQ) || p in list_ker
        push!(ds, one(E))
      else
        lp = prime_decomposition(OE, p)
        P = lp[1][1]
        k = valuation(N*OE, P)
        a = valuation(DEQ, P)
        e = valuation(DEK, P)
        g_approx = _approximate_isometry(H, H2, g, P, e, a, k, res)
        push!(ds, det(g_approx))
      end
    end
    push!(imgs, dlog(ds))
  end

  GSQ, SQtoGSQ, _ = Oscar._isomorphic_gap_group(SQ)
  f2 = hom(G2, GSQ, gensG2, SQtoGSQ.(imgs); check=false)
  f = compose(GtoG2, f2)
  @hassert :ZZLatWithIsom 1 isone(f(fqL))
  return f, GinOqL # Needs the second map to map the kernel of f into OqL
end

# We check whether for the prime ideal p: E_O(L_p) != F(L_p).
#
# There exist some prime ideals p which are elementary divisors for the
# quotients D^{-1}L^#/L and, according to Kirschmer, for which the above
# quotient can still be trivial. We call special those which do not satisfy
# this second property: considering them, we can actually reduce the size of
# the finite abelian group in which we will perform our local approximate
# determinants computations.
#
# This function is an import of a function written by Markus Kirschmer in the
# Magma package about hermitian lattices.
function _is_special(L::HermLat, p::AbsSimpleNumFieldOrderIdeal)
  OE = base_ring(L)
  E = number_field(OE)
  lp = prime_decomposition(OE, p)
  if lp[1][2] != 2 || !iseven(rank(L))
    return false
  end
  _, R, S = jordan_decomposition(L, p)
  R = Int[nrows(m) for m in R]
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
  u = E(uniformizer(P))
  s = involution(L)
  su = s(u)
  H = block_diagonal_matrix(dense_matrix_type(E)[matrix(E, 2, 2, [0 u^(S[i]); su^(S[i]) 0]) for i in 1:length(S)])
  return is_locally_isometric(L, hermitian_lattice(E; gram=H), p)
end

###############################################################################
#
#  Local hermitian lifting --- path to algorithm 8 of [BH23]
#
###############################################################################

# Once we have fixed a good local basis matrix of D^{-1}H^#, at the prime
# ideal p below P, we transfer any g\in O(D_L, D_f) via the trace construction
# (fixed by res). The new map we obtain should be a local invertible map with
# integer (for the given local field at p) entries which defines an isometry of
# D^{-1}H^# modulo H.
function _transfer_discriminant_isometry(
    res::AbstractSpaceRes,
    g::AutomorphismGroupElem{TorQuadModule},
    Bp::T,
    Bp2::T,
    pr::T
  ) where T <: MatrixElem{Hecke.RelSimpleNumFieldElem{AbsSimpleNumFieldElem}}
  q = domain(g)
  @hassert :ZZLatWithIsom 1 ambient_space(cover(q)) === domain(res)

  # B2 will be a local basis at p of the image of D^{-1}H^# under the induced
  # action of g via res.
  B2 = zero(Bp)
  for i in 1:nrows(Bp)
    B2[i, :] = res(lift(g(q(res\(Bp[i,:])))))*pr
  end
  # We should have a global exact solution. For now it does not seem to be slow
  # If for larger (quite large) examples it happens to be slow, then one could
  # modify this a bit and try to look for an inexact approximation modulo H.
  ok, K = can_solve_with_solution(Bp2, B2; side=:left)
  @hassert :ZZLatWithIsom 1 ok

  return K
end

# the minimum P-valuation among all the non-zero entries of M
function _scale_valuation(
    M::T,
    P::Hecke.RelNumFieldOrderIdeal
  ) where T <: MatrixElem{Hecke.RelSimpleNumFieldElem{AbsSimpleNumFieldElem}}
  OE = order(P)
  E = number_field(OE)
  @hassert :ZZLatWithIsom 1 base_ring(M) === E
  iszero(M) && return inf
  O = fractional_ideal(OE, one(E))
  to_sum = [ M[i, j] * O for j in 1:nrows(M) for i in 1:j ]
  d = length(to_sum)
  for i in 1:d
    push!(to_sum, involution(E)(to_sum[i]))
  end
  s = sum(to_sum; init=fractional_ideal(maximal_order(E), zero(E)))
  return valuation(s, P)
end

# the minimum P-valuation among all the non-zero diagonal entries of M
function _norm_valuation(
    M::T,
    P::Hecke.RelNumFieldOrderIdeal
  ) where T <: MatrixElem{Hecke.RelSimpleNumFieldElem{AbsSimpleNumFieldElem}}
  E = number_field(order(P))
  @hassert :ZZLatWithIsom 1 base_ring(M) === E
  iszero(M) && return inf
  return Hecke._get_norm_valuation_from_gram_matrix(M, P)
end

# This is algorithm 8 of [BH23]: under the good assumptions, we can do a
# P-adic lifting of a matrix which represents an isometry up to a certain
# precision. In this way, we approximate our matrix by another matrix, to a
# given precision and the new matrix defines also an isometry up to a finer
# precision than the initial matrix.
#
# We use this method iteratively to lift isometries (along a surjective map),
# by looking at better representatives until we reach a good enough precision
# for our purpose (i.e. computing its determinant).
function _local_hermitian_lifting(
    G::T,
    F::T,
    rho::Hecke.RelSimpleNumFieldElem,
    l::Int,
    P::Hecke.RelNumFieldOrderIdeal,
    P2::Hecke.RelNumFieldOrderIdeal,
    split::Bool,
    e::Int,
    a::Int;
    check = true
  ) where T <: MatrixElem{Hecke.RelSimpleNumFieldElem{AbsSimpleNumFieldElem}}
  @hassert :ZZLatWithIsom 1 trace(rho) == 1
  E = base_ring(G)
  s = involution(E)
  # G here is a local gram matrix
  @hassert :ZZLatWithIsom 1 G == map_entries(s, transpose(G))
  @hassert :ZZLatWithIsom 1 base_ring(F) === E

  # R represents the defect, how far F is to be an isometry of G
  R = G - F*G*map_entries(s, transpose(F))
  # These are the necessary conditions for the input of algorithm 8 in [BH23]
  if check
    @hassert :ZZLatWithIsom 1 _scale_valuation(inv(G), P) >= 1+a
    @hassert :ZZLatWithIsom 1 _norm_valuation(inv(G), P) + valuation(rho, P) >= 1+a
    @hassert :ZZLatWithIsom 1 _scale_valuation(R, P) >= l-a
    @hassert :ZZLatWithIsom 1 _norm_valuation(R, P) + valuation(rho,P) >= l-a
    if split
      @hassert :ZZLatWithIsom 1 _scale_valuation(inv(G), P2) >= 1+a
      @hassert :ZZLatWithIsom 1 _norm_valuation(inv(G), P2) + valuation(rho, P2) >= 1+a
      @hassert :ZZLatWithIsom 1 _scale_valuation(R, P2) >= l-a
      @hassert :ZZLatWithIsom 1 _norm_valuation(R, P2) + valuation(rho, P2) >= l-a
    end
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

  diag = R - U - map_entries(s, transpose(U))

  # This newF is supposed to be a better lift than F, i.e. it is congruent to F
  # modulo P^{l+1} and the corresponding defect R2 has a higher P-valuation (so
  # P-adically, we are close to have an actual isometry)
  newF = F + (U + rho*diag)*map_entries(s, inv(transpose(F)))*inv(G)

  l2 = 2*l+1
  if check
    @hassert :ZZLatWithIsom 1 minimum(valuation(m, P) for m in (F-newF) if !iszero(m)) >= l+1
    R2 = G-newF*G*map_entries(s, transpose(newF))
    @hassert :ZZLatWithIsom 1 _scale_valuation(R2, P) >= l2-a
    @hassert :ZZLatWithIsom 1 _norm_valuation(R2, P) + valuation(rho, P) >= l2-a
    if split
      @hassert :ZZLatWithIsom 1 minimum(valuation(m, P2) for m in (F-newF) if !iszero(m)) >= l+1
      @hassert :ZZLatWithIsom 1 _scale_valuation(R2, P2) >= l2-a
      @hassert :ZZLatWithIsom 1 _norm_valuation(R2, P2) + valuation(rho, P2) >= l2-a
    end
  end

  return newF, l2
end

# Starting from an isometry g of H2/H (here H2=  D^{-1}H^#), and a
# prime ideal P, we iteratively lift g to invertible maps which define
# isometries of H2 up to a given precision specified in [BH23], locally
# at the prime ideal p below P.
#
# Here:
#  - e is the P-valuation of the relative different
#  - a is the P-valuation of the absolute different
#  - k is the P-valuation of the norm of H
#  - res is the map representing the functor used for the trace equivalence
#  - g is the torsion quadratic module automorphism, which centralizes D_f in
#    D_L (H is the hermitian structure associated to (L, f) along res) which
#    aims to approximately transfer to H2 along res at P
function _approximate_isometry(
    H::HermLat,
    H2::HermLat,
    g::AutomorphismGroupElem{TorQuadModule},
    P::Hecke.RelNumFieldOrderIdeal,
    e::Int,
    a::Int,
    k::Int,
    res::AbstractSpaceRes
  )
  E = base_field(H)
  s = involution(E)
  P2 = s(P)
  split = P2 != P
  # In the split case, we need to check valuation at both primes above p
  P2.is_prime = 1
  @hassert :ZZLatWithIsom 1 number_field(order(P)) === E
  p = minimum(P)
  ok, b = is_modular(H, p)
  if ok && b == -a
    return identity_matrix(E, 1)
  end

  Bp, Bp2, pr = _local_basis_matrix_and_projection(H2, p, a, res)
  Gp = Bp2*gram_matrix(ambient_space(H))*map_entries(s, transpose(Bp2))
  Fp = _transfer_discriminant_isometry(res, g, Bp, Bp2, pr)
  # This is the local defect. By default, it should have scale P-valuations -a
  # and norm P-valuation e-1-a
  Rp = Gp - Fp*Gp*map_entries(s, transpose(Fp))
  rho = _find_rho(P, e)

  l = 0
  while _scale_valuation(Rp, P) < k+a && (!split || _scale_valuation(Rp, P2) < k+a)
    Fp, l = _local_hermitian_lifting(Gp, Fp, rho, l, P, P2, split, e, a)
    Rp = Gp - Fp*Gp*map_entries(s, transpose(Fp))
  end
  return Fp
end

# We use this function to compute a local basis matrix of Hv at p, whose
# integral span defines a sublattice of Hv. If H has a P^{-a}-modular block
# at p, we remove the corresponding basis vector(s) from the output to only
# keep the Jordan blocks of Hv_p which are of P-valuation different from -a.
#
# In our context of use, -a is actually the biggest scale valuation for
# Hv. Since we took care earlier to discard the cases where Hv = D^{-1}H^# is
# P^{-a}=modular locally at p, we just have to remove the last Jordan block
# from the Jordan decomposition of Hv at p. From that point, we massage a bit
# the basis matrices of the other Jordan blocks to obtain local basis matrices
# which span sublattices of Hv.
#
# If Hv is locally U + R where U if a Jordan block of scale P^{-a}, the matrix
# pr in output is the projection onto R. We multiply our local basis for R by
# pr to remove all residue from U: we use this map also to compute
# approximation of local isometries later.
function _local_basis_matrix_and_projection(
    Hv::HermLat,
    p::AbsSimpleNumFieldOrderIdeal,
    a::Int,
    res::AbstractSpaceRes
  )
  Bv, Gv , expsv = jordan_decomposition(Hv, p)
  pr = identity_matrix(base_field(Hv), rank(Hv))
  # If the last Jordan block U is P^{-a}-modular, then this part will
  # vanish in the quotient Hv/H and we want to get rid of it. However
  # we need to remember a decomposition Hv_p = R+U: for this
  # we define the projection map Hv_p \to R so later we can
  # remove any trace of U in our elements (for computing the first
  # approximation of our isometries)
  if expsv[end] == -a
    Bp = reduce(vcat, Bv[1:end-1])
    Bp0 = reduce(vcat, [Bp, zero(Bv[end])])
    Bv = reduce(vcat, Bv)
    pr = solve(Bv, Bp0; side=:right)
  else
    Bp = reduce(vcat, Bv)
  end

  # Now that we have a good looking local basis at p for R, we massage it
  # a bit to make sure that globally it spans a sublattice of Hv.
  H2v = lattice_in_same_ambient_space(Hv, Bp)
  if !is_sublattice(Hv, H2v)
    Lv = restrict_scalars(Hv, res)
    L2 = restrict_scalars(H2v, res)
    L2 = intersect(Lv, L2)
    B2 = basis_matrix(L2)
    gene = Vector{elem_type(base_field(Hv))}[res(vec(collect(B2[i, :]))) for i in 1:nrows(B2)]
    H2v = lattice(ambient_space(Hv), gene)
    @hassert :ZZLatWithIsom 1 is_sublattice(Hv, H2v)
    @hassert :ZZLatWithIsom 1 is_locally_isometric(Hv, H2v, p)
    Bp = local_basis_matrix(H2v, p; type = :submodule)
    @hassert :ZZLatWithIsom 1 is_locally_isometric(Hv, lattice_in_same_ambient_space(Hv, Bp), p)
    @hassert :ZZLatWithIsom 1 is_sublattice(Hv, lattice_in_same_ambient_space(Hv, Bp))
  end
  @hassert :ZZLatWithIsom 1 rank(Bp) == rank(Bp*pr)
  return Bp, Bp*pr, pr
end

# We need a special rho for Algorithm 8 of [BH23]: we construct such an element
# here, which will be used to lift isometries up to a better precision.
#
# If the prime ideal is non dyadic, then 1//2 will always be good for us.
# Otherwise, if the prime p is dyadic, we do in two different ways
#  - p is inert or split, in which case we can use an algorithm of Markus
#    Kirschmer, implemented on Hecke, which provides us with a "special unit",
#    which is exactly what we look for (a unit at P with trace 1);
#  - p is ramified, and then we cook up a good element. The actual code from
#    that part is taken from the Sage implementation of Simon Brandhorst
function _find_rho(P::Hecke.RelNumFieldOrderIdeal, e::Int)
  OE = order(P)
  E = number_field(OE)
  p = minimum(P)
  lp = prime_decomposition(OE, p)
  dya = is_dyadic(P)
  if !dya
    rho = E(1//2)
    @hassert :ZZLatWithIsom 1 trace(rho) == 1
    return rho
  end
  lp = prime_decomposition(OE, p)
  if lp[1][2] == 1
    rho = Hecke._special_unit(P, p)
    @hassert :ZZLatWithIsom 1 trace(rho) == 1
    return rho
  end
  K = base_field(E)
  Eabs, EabstoE = absolute_simple_field(E)
  Pabs = EabstoE\P
  OEabs = order(Pabs)
  while true
    Eabst, t = Eabs[:t]
    g = EabstoE\(E(-rand(K, -5:5)^2-1))
    nu = 2*valuation(Eabs(2), Pabs)-2*e+2
    nug = valuation(g, Pabs)
    if nu == nug
      d = denominator(g+1, OEabs)
      rt = roots(t^2 - (g+1)*d^2; max_roots=1, ispure=true, is_normal=true)
      if !is_empty(rt)
        rho = (1+rt[1])//2
        @hassert :ZZLatWithIsom 1 valuation(EabstoE(rho), P) == 1-e
        @hassert :ZZLatWithIsom 1 trace(EabstoE(rho)) == 1
        return EabstoE(rho)
      end
    end
  end
end
