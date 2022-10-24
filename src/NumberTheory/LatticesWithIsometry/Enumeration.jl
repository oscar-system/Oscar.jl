export admissible_triples, is_admissible_triple

##################################################################################
#
# This is an import to Oscar of the methods written following the paper [BH22] on
# "Finite subgroups of automorphisms of K3 surfaces".
#
##################################################################################

# The tuples in output are pairs of positive integers!
function _tuples_divisors(d::fmpz)
  div = divisors(d)
  return [(dd,abs(divexact(d,dd))) for dd in div]
end

# This is line 8 of Algorithm 1, they correspond to the possible
# discriminant for the genera A and B to glue to fit in C. d is
# the determinant of C, m the maximal p-valuation of the gcd of
# d1 and dp.
function _find_D(d::fmpz, m::Int, p::Int)
  @assert is_prime(p)
  @assert d != 0

  # If m == 0, there are no conditions on the gcd of d1 and dp
  if m == 0
    return _tuples_divisors(d)
  end
  
  D = Tuple{Int,Int}[]
  # We try all the values of g possible, from 1 to p^m
  for j=0:m 
    g = p^j
    dj = _tuples_divisors(d*g^2)
    for (d1,dp) in dj
      if mod(d1,g) == mod(dp,g) == 0
        push!(D,(d1,dp))
      end
    end
  end
  return D
end

# This is line 10 of Algorithm 1. We need the condition on the even-ness of
# C since subgenera of an even genus are even too. r is the rank of
# the subgenus, d its determinant, s and l the scale and level of C
function _find_L(r::Int, d, s::fmpz, l::fmpz, p::Int, even = true)
  L = ZGenus[]
  for (s1,s2) in [(s,t) for s=0:r for t=0:r if s+t==r]
    gen = genera((s1,s2), d, even=even)
    filter!(G -> divides(scale(G), s)[1], gen)
    filter!(G -> divides(p*l, level(G))[1], gen)
    append!(L, gen)
  end
  return L
end

@doc Markdown.doc"""
    is_admissible_triple(A::ZGenus, B::ZGenus, C::ZGenus, p::Integer) -> Bool

Given a triple of $\mathbb Z$-genera `(A,B,C)` and a prime number `p`, such
that the rank of `B` is divisible by $p-1$ and the level of `C` is a power
of `p`, return whether `(A,B,C)` is `p`-admissible in the sense of
Definition 4.13. [BH22]
"""
function is_admissible_triple(A::ZGenus, B::ZGenus, C::ZGenus, p::Integer)
  zg = genus(Zlattice(gram = matrix(QQ, 0, 0, [])))
  if ((A == zg) && (B == C)) || ((B == zg) && (A == C))
    # C can be always glued with the empty genus to obtain C
    return true
  elseif (A == zg) || (B == zg)
    # If A or B is empty but the other is not C, then there is no glueing
    return false
  end

  @req divides(rank(B), p-1)[1] "p-1 must divide the rank of B"
  @req prime_divisors(level(C)) == [p] || level(C) == one(fmpq) "The level of C must be a power of p"
  AperpB = orthogonal_sum(A,B)
  # If A and B glue to C, the sum of their ranks must match the one of C
  if rank(AperpB) != rank(C)
    return false
  end

  lA = ngens(discriminant_group(A))
  lB = ngens(discriminant_group(B))
  if all(g -> abs(det(AperpB)) != p^(2*g)*abs(det(C)), 0:min(lA, lB, divexact(rank(B), p-1)))
    return false
  end

  # A+B and C must agree locally at every primes except p
  for q in filter(qq -> qq != p, union([2], primes(AperpB), primes(C)))
    if local_symbol(AperpB,q) != local_symbol(C,q)
      return false
    end
  end

  g = divexact(valuation(div(det(AperpB), det(C)), p), 2)
  # If the determinants agree, A+B = C and, equivalently, they agree locally at p too.
  # Otherwise, their localisations at p must be rationally equal.
  if g == 0
    return local_symbol(AperpB, p) == local_symbol(C, p)
  elseif excess(local_symbol(AperpB, p)) != excess(local_symbol(C, p))
    return false
  end

  if !(divides(scale(AperpB), scale(C))[1] && divides(p*level(C), level(AperpB))[1])
    return false
  end

  # At this point, if C is unimodular at p, the glueing condition is equivalent to have 
  # an anti-isometry between the p-part of the (quadratic) discriminant forms of A and B
  qA = discriminant_group(A)
  qB = discriminant_group(B)
  if !divides(det(C), p)[1]
    return is_anti_isometric_with_anti_isometry(primary_part(qA, p)[1], primary_part(qB, p)[1])[1]
  end

  l = valuation(level(C), p)
  Ap = symbol(local_symbol(A, p))
  Bp = symbol(local_symbol(B, p))
  a_max = sum(Int[s[2] for s in Ap if s[1] == l+1])
  b_max = sum(Int[s[2] for s in Bp if s[1] == l+1])
  # For the glueing, rho_{l+1}(A_p) and rho_{l+1}(B_p) are anti-isometric, so they must have the
  # same order
  if a_max != b_max
    return false
  end
  
  # Since p^l*A^\vee/A is in the glue, its order is less than the order of the glue
  if g < a_max
    return false
  end

  a1 = sum(Int[s[2] for s in Ap if s[1] == 1])
  a2 = sum(Int[s[2] for s in Ap if s[1] >= 2])
  b1 = sum(Int[s[2] for s in Bp if s[1] == 1])
  b2 = sum(Int[s[2] for s in Bp if s[1] >= 2])

  ABp = symbol(local_symbol(AperpB, p))
  Cp = symbol(local_symbol(C, p))

  if a_max == g
    if length(Ap) > 1
      Ar = ZpGenus(p, Ap[1:end-1])
    else
      Ar = genus(matrix(ZZ,0,0,[]), p)
    end
   
    if length(Bp) > 1
      Br = ZpGenus(p, Bp[1:end-1])
    else
      Br = genus(matrix(ZZ, 0, 0, []), p)
    end
  
    ABr = orthogonal_sum(Ar, Br)
    
    for i = 1:l
      s1 = symbol(ABr)[i]
      s2 = symbol(Cp)[i]
      if s1[2] != s2[2]
        return false
      end
      if p == 2 && s1[4] > s2[4]
        return false
      elseif p != 2 && s1[3] != s2[3]
        return false
      end
    end
  end

  s3 = symbol(local_symbol(C, 2))[end]
  if p == 2
    s1 = Ap[s3[1]+1]
    s2 = Bp[s3[1]+1]
    if s1[3] != s2[3]
      return false
    end
  end

  for s in Cp
    s[1] += 2
  end
  Cp = ZpGenus(p, Cp)
  
  if !represents(local_symbol(AperpB, p), Cp)
    return false
  end
  if !represents(C, AperpB)
    return false
  end 
  
  return true
end

function is_admissible_triple(A::Union{ZLat, ZGenus}, B::Union{ZLat, ZGenus}, C::Union{ZLat, ZGenus}, p::Integer)
  L = ZGenus[typeof(D) <: ZLat ? genus(D) : D for D = (A, B, C)]
  return is_admissible_triple(L[1], L[2], L[3], p)
end

@doc Markdown.doc"""
    admissible_triples(C::ZGenus, p::Integer) -> Vector{Tuple{ZGenus, ZGenus}}

Given a $\mathbb Z$-genus `C` and a prime number `p`, return all tuples of
$\mathbb Z$-genera `(A, B)` such that `(A, B, C)` is `p`-admissible and
`B` is of rank divisible by $p-1$. See Algorithm 1 [BH22].
"""
function admissible_triples(G::ZGenus, p::Int64)
  @req is_prime(p) "p must be a prime number"
  n = rank(G)
  d = det(G)
  even = iseven(G)
  L = Tuple{ZGenus, ZGenus}[]
  for ep in 0:div(n, p-1)
    rp = (p-1)*ep
    r1 = n - rp
    m = min(ep, r1) 
    D = _find_D(d, m, p)
    for (d1, dp) in D
      L1 = _find_L(r1, d1, scale(G), level(G), p, even)
      Lp = _find_L(rp, dp, scale(G), level(G), p, even)
      for (A, B) in [(A, B) for A in L1 for B in Lp]
        if is_admissible_triple(A, B, G, p)
          push!(L, (A, B))
        end
      end
    end
  end
  return L
end

admissible_triples(L::ZLat, p::Integer) = admissible_triples(genus(L), p)

function _get_V(q, f, fq, mu, p)
  L = relations(q)
  f_ambient = inv(basis_matrix(L))*f*basis_matrix(L)
  @assert f*gram_matrix(L)*transpose(f) == gram_matrix(L)
  f_mu = mu(f_ambient)
  if !is_zero(f_mu)
    @assert det(f_mu) != 0
    L_sub = intersect(lattice_in_same_ambient_space(L, inv(f_mu)), dual(L))
    B = basis_matrix(L_sub)
    V, Vinq = sub(q, q.([vec(collect(B[i,:])) for i in 1:nrows(B)]))
  else
    V, Vinq = q, id_hom(q)
  end
  pV, pVinV = primary_part(V, p)
  pV, pVinV = sub(V, pVinV.([divexact(order(g), p)*g for g in gens(pV) if !(order(g)==1)]))
  pVinq = compose(pVinV, Vinq)
  fpV = _restrict(fq, pVinq)
  return pV, pVinq, fpV
end

function _rho_functor(q::TorQuadMod, p, l::Union{Integer, fmpz})
  N = relations(q)
  if l == 0
    Gl = N
    Gm = intersect(1//p*N, dual(N))
    rholN = torsion_quadratic_module(Gl, p*Gm, modulus = QQ(1), modulus_qf = QQ(2))
    gene = [q(lift(g)) for g in gens(rholN)]
    return sub(q, gene)[2]
  end
  k = l-1
  m = l+1
  Gk = intersect(1//p^k*N, dual(N))
  Gl = intersect(1//p^l*N, dual(N))
  Gm = intersect(1//p^m*N, dual(N))
  rholN = torsion_quadratic_module(Gl, Gk+p*Gm, modulus = QQ(1), modulus_qf = QQ(2))
  gene = [q(lift(g)) for g in gens(rholN)]
  return sub(q, gene)[2]
end

function primitive_extensions(Afa::LatticeWithIsometry, Bfb::LatticeWithIsometry, Cfc::LatticeWithIsometry, p::Integer)
  # requirement for the algorithm of BH22
  @req is_prime(p) "p must be a prime number"  
  A, B, C = lattice.([Afa, Bfb, Cfc])
  @req is_admissible_triple(A, B, C, p) "(A, B, C) must be p-admissble"
  # we need to compare the type of the output with the one of (C, f_c^p)
  t = type(lattice_with_isometry(C, isometry(Cfc)^p))

  results = LatticeWithIsometry[]
  # this is the glue valuation: it is well-defined because the triple in input is admissible
  g = div(valuation(divexact(det(A)*det(B), det(C)), p), 2)
  
  fA, fB = isometry.([Afa, Bfb])
  qA, fqA = discriminant_group(Afa)
  qB, fqB = discriminant_group(Bfb)
  GA = image_centralizer_in_Oq(Afa)
  GB = image_centralizer_in_Oq(Bfb)

  # this is where we will perform the glueing
  D, qAinD, qBinD = orthogonal_sum(qA, qB)
  OD = orthogonal_group(D)
  prim_ext = Tuple{TorQuadMod, AutomorphismGroup}[]
  OqAinOD, OqBinOD = embedding_orthogonal_group(qAinD, qBinD)

  # if the glue valuation is zero, then we glue along the trivial group and we don't
  # have much more to do. Since the triple is p-admissible, A+B = C
  if g == 0
    geneA = OqAinOD.(gens(GA))
    geneB = OqBinOD.(gens(GB))
    gene = vcat(geneA, geneB)
    GC2 = sub(OD, gene)[1]
    C2 = orthogonal_sum(A, B)[1]
    fC2 = block_diagonal_matrix(fA, fB)
    C2fC2 = lattice_with_isometry(C2, fC2)
    if type(C2fC2) == t
      set_attribute!(C2fC2, :image_centralizer_in_Oq, GC2)
      push!(results, C2fC2)
    end
    @goto exit
  end

  # these are GA/GM-invariant, fA/fB-stable, and should contain the kernels of any glue map
  VA, VAinqA, fVA = _get_V(qA, fA, fqA, minpoly(Bfb), p)
  VB, VBinqB, fVB = _get_V(qB, fB, fqB, minpoly(Afa), p)
  
  # since the glue kernels must have order p^g, in this condition, we have nothing
  if min(order(VA), order(VB)) < p^g
    return results
  end

  # scale of the dual: any glue kernel must contain the multiple of l of the respective
  # discriminant groups
  l = level(genus(C))

  # We look for the GA/GB-invariant and fA/fB-stable subgroups of VA/VB which respectively
  # contained lqA/lqB. This is done by computing orbits and stabilisers of VA/lqA (resp VB/lqB)
  # seen as a F_p-vector space under the action of GA (resp. GB). Then we check which ones
  # are fA-stable (resp. fB-stable)
  subsA = _subgroups_representatives(VAinqA, GA, g, fVA, l)
  subsB = _subgroups_representatives(VBinqB, GB, g, fVB, l)
  # once we have the potential kernels, we create pairs of anti-isometric groups since glue
  # maps are anti-isometry
  R = [(H1, H2) for H1 in subsA for H2 in subsB if is_anti_isometric_with_anti_isometry(domain(H1[1]), domain(H2[1]))[1]]
  # now, for each pair of anti-isometric potential kernels, we need to see whether
  # it is (fA,fB)-equivariant, up to conjugacy. For each working pair, we compute the
  # corresponding overlattice and check whether it satisfies the type condition
  for (H1, H2) in R
    #return H1, H2
    SAinqA, stabA = H1
    SA = domain(SAinqA)
    # SA might not be in normal form
    SBinqB, stabB = H2
    SB = domain(SBinqB)
    # we already know they are but now we get the map
    ok, phi = is_anti_isometric_with_anti_isometry(SA, SB)
    @assert ok
    OSA = orthogonal_group(SA)
    # we compute the image of the stabalizer in OSA
    actA = hom(stabA, OSA, [OSA(Oscar._restrict(x, SAinqA)) for x in gens(stabA)])
    imA = image(actA, stabA)[1]
    # we keep track of the element of the stabilizer acting trivially on SA
    kerA = [OqAinOD(x) for x in gens(kernel(actA)[1])]
    fSA = _restrict(fqA, SAinqA)
    fSA = OSA(fSA)
    OSB = orthogonal_group(SB)
    actB = hom(stabB, OSB, [OSB(Oscar._restrict(x, SBinqB)) for x in gens(stabB)])
    imB = image(actB, stabB)
    kerB = [OqBinOD(x) for x in gens(kernel(actB)[1])]
    fSB = _restrict(fqB, SBinqB)
    fSB = OSB(fSB)
    # we get all the elements of qB of order exactly p^l, which are not mutiple of an
    # element of order p^{l+1}. In theory, glue maps are classified by the orbit of phi
    # under the action of O(SB, rho_l(qB), fB)
    rBinqB = _rho_functor(qB, p, valuation(l, p))
    @assert Oscar._is_invariant(stabB, rBinqB) # otherwise there is something wrong!
    #return rBinqB, SBinqB
    rBinSB = hom(domain(rBinqB), SB, [SBinqB\(rBinqB(k)) for k in gens(domain(rBinqB))])
    @assert is_injective(rBinSB) # we indeed have rho_l(qB) which is a subgroup of SB
    # We compute the generators of O(SB, rho_l(qB))
    stabrB = [g for g in collect(OSB) if Oscar._is_invariant(g, rBinSB)]
    OSBrB, _ = sub(OSB, stabrB)
    @assert fSB in OSBrB # otherwise there is something wrong
    # phi might not "send" the restriction of fA to this of fB, but at least phi(fA)phi^-1
    # should be conjugate to fB inside O(SB, rho_l(qB)) for the glueing.
    # If not, we try the next potential pair.
    fSAinOSBrB = OSBrB(compose(inv(phi), compose(hom(fSA), phi)))
    bool, g0 = representative_action(OSBrB, fSAinOSBrB, OSBrB(fSB))
    bool || continue
    phi = compose(phi, hom(OSB(g0)))
    fSAinOSB = OSB(compose(inv(phi), compose(hom(fSA), phi)))
    # Now the new phi is "sending" the restriction of fA to this of fB.
    # So we can glue SA and SB.
    @assert fSAinOSB == fSB
    # Now it is time to compute generators for O(SB, rho_l(qB), fB), and the induced
    # image of stabA for taking the double cosets next
    center, _ = centralizer(OSBrB, OSBrB(fSB))
    center, _ = sub(OSB, OSB.(gens(center)))
    stabSAphi = [OSB(compose(inv(phi), compose(hom(actA(g)), phi))) for g in gens(stabA)]
    stabSAphi, _ = sub(center, stabSAphi)
    stabSB, _ = sub(center, actB.(gens(stabB)))
    reps = double_cosets(center, stabSAphi, stabSB)
    # now we iterate over all double cosets and for each representative, we compute the
    # corresponding overlattice in the glueing. If it has the wanted type, we compute
    # the image of the centralizer in OD from the stabA ans stabB.
    for g in reps
      g = representative(g)
      phig = compose(phi, hom(g))
      # problem of lattices that are not in the same ambient space...
      # TODO: fix `overlattice` or change something somewhere...
      S = relations(domain(phig))
      R = relations(codomain(phig))
      VS = ambient_space(S)
      diag, trafo = Hecke._gram_schmidt(gram_matrix(S), identity)
      V1 = quadratic_space(QQ, diag)
      S1 = lattice(V1, basis_matrix(S)*inv(trafo))
      VR = ambient_space(R)
      b,T = isisometric_with_isometry(VR, V1)
      @assert b
      R1 = lattice(V1, T)
      glue = [lift(qAinD(SAinqA(g))) + lift(qBinD(SBinqB(phig(g)))) for g in gens(domain(phig))]
      z = zero_matrix(QQ,0,degree(S)+degree(R))
      glue = reduce(vcat, [matrix(QQ,1,degree(S)+degree(R),g) for g in glue], init=z)
      glue = vcat(block_diagonal_matrix([basis_matrix(S1), basis_matrix(R1)]), glue)
      glue = FakeFmpqMat(glue)
      B = hnf(glue)
      B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(B))
      C2 = lattice(, B[end-rank(S)-rank(R)+1:end, :])
      fC2 = block_diagonal_matrix([fA, fB])
      return C2, fC2
      C2fC2 = lattice_with_isometry(C2, fC2)
      if type(C2fC2) != t
        continue
      end
      im2_phi = sub(OSA, [compose(phig, compose(g1, inv(phig))) for g1 in gens(imB)])
      im3, _ = sub(imA, gens(intersect(imA, im2_phi)[1]))
      stab = [(x, imB(compose(inv(phig), compose(x, phig)))) for x in gens(im3)]
      stab = [OqAinOD(x[1])*OqBinOD(x[2]) for x in stab]
      stab = union(stab, kerA)
      stab = union(stab, kerB)
      stab = Oscar._orthogonal_group(discriminant_group(C2), matrix.(stab))
      set_attribute!(C2fC2, :image_centralizer_in_Oq, stab)
      push!(results, C2fC2)
    end
  end

  @label exit
  return results
end
