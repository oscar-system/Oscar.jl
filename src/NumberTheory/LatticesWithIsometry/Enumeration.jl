export admissible_triples, is_admissible_triple

##################################################################################
#
# This is an import to Oscar of the methods written following the paper [BH22] on
# "Finite subgroups of automorphisms of K3 surfaces".
#
##################################################################################

##################################################################################
#
# Admissible triples
#
##################################################################################

# The tuples in output are pairs of positive integers!
function _tuples_divisors(d::T) where T <: Union{Integer, fmpz}
  div = divisors(d)
  return Tuple{T, T}[(dd,abs(divexact(d,dd))) for dd in div]
end

# This is line 8 of Algorithm 1, they correspond to the possible
# discriminant for the genera A and B to glue to fit in C. d is
# the determinant of C, m the maximal p-valuation of the gcd of
# d1 and dp.
function _find_D(d::T, m::Int, p::Int) where T <: Union{Integer, fmpz}
  @assert is_prime(p)
  @assert d != 0

  # If m == 0, there are no conditions on the gcd of d1 and dp
  if m == 0
    return _tuples_divisors(d)
  end
  
  D = Tuple{T, T}[]
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
      s2 = Cp[i]
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

function is_admissible_triple(A::T, B::T, C::T, p::Integer) where T <: Union{ZLat, LatticeWithIsometry}
  L = ZGenus[genus(D) for D = (A, B, C)]
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

admissible_triples(L::T, p::Integer) where T <: Union{ZLat, LatticeWithIsometry} = admissible_triples(genus(L), p)

##################################################################################
#
# Primitive extensions
#
##################################################################################

function _get_V(L, q, f, fq, mu, p)
  f_mu = mu(f)
  if !is_zero(f_mu)
    L_sub = intersect(lattice_in_same_ambient_space(L, inv(f_mu)*basis_matrix(L)), dual(L))
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

  if ambient_space(Afa) === ambient_space(Bfb)
    @req rank(intersect(A, B)) == 0 "Lattice in same ambient space must have empty intersection to glue"
    @req rank(A) + rank(B) <= dim(ambient_space(A)) "Lattice cannot glue in their ambient space"
  end

  @req is_admissible_triple(Afa, Bfb, Cfc, p) "(A, B, C) must be p-admissble"
  # we need to compare the type of the output with the one of Cfc
  t = type(Cfc)

  results = LatticeWithIsometry[]
  # this is the glue valuation: it is well-defined because the triple in input is admissible
  g = div(valuation(divexact(det(A)*det(B), det(C)), p), 2)

  fA, fB = isometry.([Afa, Bfb])
  qA, fqA = discriminant_group(Afa)
  qB, fqB = discriminant_group(Bfb)
  GA = image_centralizer_in_Oq(Afa)
  GB = image_centralizer_in_Oq(Bfb)

  # this is where we will perform the glueing
  if ambient_space(Afa) === ambient_space(Bfb)
    D, qAinD, qBinD = inner_orthogonal_sum(qA, qB)
  else
    D, qAinD, qBinD = orthogonal_sum(qA, qB)
  end

  OD = orthogonal_group(D)
  OqAinOD, OqBinOD = embedding_orthogonal_group(qAinD, qBinD)
  OqA = domain(OqAinOD)
  OqB = domain(OqBinOD)

  # if the glue valuation is zero, then we glue along the trivial group and we don't
  # have much more to do. Since the triple is p-admissible, A+B = C
  if g == 0
    geneA = gens(GA)
    geneA = AutomorphismGroupElem{TorQuadMod}[OqAinOD(a) for a in geneA]
    geneB = gens(GB)
    geneB = AutomorphismGroupElem{TorQuadMod}[OqBinOD(b) for b in geneB]
    gene = vcat(geneA, geneB)
    GC2, _ = sub(OD, gene)
    if ambient_space(Afa) === ambient_space(Bfb)
      C2 = A+B
      fC2 = block_diagonal_matrix([fA, fB])
      _B = solve_left(reduce(vcat, basis_matrix.([A,B])), basis_matrix(C2))
      fC2 = _B*fC2*inv(_B)
      @assert fC2*gram_matrix(C2)*transpose(fC2) == gram_matrix(C2)
    else
      C2 = orthogonal_sum(A, B)[1]
      fC2 = block_diagonal_matrix(fA, fB)
    end
    if type(lattice_with_isometry(C2, fC2^p)) == t
      C2fC2 = lattice_with_isometry(C2, fC2, ambient_representation=false)
      set_attribute!(C2fC2, :image_centralizer_in_Oq, GC2)
      push!(results, C2fC2)
      return results
    end
  end

  # these are GA|GB-invariant, fA|fB-stable, and should contain the kernels of any glue map
  VA, VAinqA, fVA = _get_V(A, qA, isometry(Afa), fqA, minpoly(Bfb), p)
  VB, VBinqB, fVB = _get_V(B, qB, isometry(Bfb), fqB, minpoly(Afa), p)

  # since the glue kernels must have order p^g, in this condition, we have nothing
  if min(order(VA), order(VB)) < p^g
    return results
  end

  # scale of the dual: any glue kernel must contain the multiples of l of the respective
  # discriminant groups
  l = level(genus(C))

  # We look for the GA|GB-invariant and fA|fB-stable subgroups of VA|VB which respectively
  # contained lqA|lqB. This is done by computing orbits and stabilisers of VA/lqA (resp VB/lqB)
  # seen as a F_p-vector space under the action of GA (resp. GB). Then we check which ones
  # are fA-stable (resp. fB-stable)
  subsA = _subgroups_representatives(VAinqA, GA, g, fVA, l)
  subsB = _subgroups_representatives(VBinqB, GB, g, fVB, l)

  # once we have the potential kernels, we create pairs of anti-isometric groups since glue
  # maps are anti-isometry
  R = Tuple{eltype(subsA), eltype(subsB), TorQuadModMor}[]
  for H1 in subsA, H2 in subsB
    ok, phi = is_anti_isometric_with_anti_isometry(domain(H1[1]), domain(H2[1]))
    !ok && continue
    push!(R, (H1, H2, phi))
  end

  # now, for each pair of anti-isometric potential kernels, we need to see whether
  # it is (fA,fB)-equivariant, up to conjugacy. For each working pair, we compute the
  # corresponding overlattice and check whether it satisfies the type condition
  for (H1, H2, phi) in R
    SAinqA, stabA = H1
    SA = domain(SAinqA)
    OSAinOqA = embedding_orthogonal_group(SAinqA)
    OSA = domain(OSAinOqA)
    OSAinOD = compose(OSAinOqA, OqAinOD)

    SBinqB, stabB = H2
    SB = domain(SBinqB)
    OSBinOqB = embedding_orthogonal_group(SBinqB)
    OSB = domain(OSBinOqB)
    OSBinOD = compose(OSBinOqB, OqBinOD)

    # we compute the image of the stabalizers in the respective OS* and we keep track
    # of the elements of the stabilizers action trivially in the respective S*

    actA = hom(stabA, OSA, [OSA(Oscar._restrict(x, SAinqA)) for x in gens(stabA)])
    imA, _ = image(actA)
    kerA = AutomorphismGroupElem{TorQuadMod}[OqAinOD(x) for x in gens(kernel(actA)[1])]
    fSA = OSA(_restrict(fqA, SAinqA))

    actB = hom(stabB, OSB, [OSB(Oscar._restrict(x, SBinqB)) for x in gens(stabB)])
    imB, _ = image(actB)
    kerB = AutomorphismGroupElem{TorQuadMod}[OqBinOD(x) for x in gens(kernel(actB)[1])]
    fSB = OSB(_restrict(fqB, SBinqB))

    # we get all the elements of qB of order exactly p^{l+1}, which are not mutiple of an
    # element of order p^{l+2}. In theory, glue maps are classified by the orbit of phi
    # under the action of O(SB, rho_{l+1}(qB), fB)
    rBinqB = _rho_functor(qB, p, valuation(l, p)+1)
    @assert Oscar._is_invariant(stabB, rBinqB)
    rBinSB = hom(domain(rBinqB), SB, [SBinqB\(rBinqB(k)) for k in gens(domain(rBinqB))])
    @assert is_injective(rBinSB) # we indeed have rho_{l+1}(qB) which is a subgroup of SB

    # We compute the generators of O(SB, rho_l(qB))
    OrBinOSB = embedding_orthogonal_group(rBinSB)
    OSBrB, _ = image(OrBinOSB)
    @assert fSB in OSBrB

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
    # images of stabA|stabB for taking the double cosets next
    center, _ = centralizer(OSBrB, OSBrB(fSB))
    center, _ = sub(OSB, [OSB(c) for c in gens(center)])
    stabSAphi = AutomorphismGroupElem{TorQuadMod}[OSB(compose(inv(phi), compose(hom(actA(g)), phi))) for g in gens(stabA)]
    stabSAphi, _ = sub(center, stabSAphi)
    stabSB, _ = sub(center, [actB(s) for s in gens(stabB)])
    reps = double_cosets(center, stabSAphi, stabSB)

    # now we iterate over all double cosets and for each representative, we compute the
    # corresponding overlattice in the glueing. If it has the wanted type, we compute
    # the image of the centralizer in OD from the stabA ans stabB.
    for g in reps
      g = representative(g)
      phig = compose(phi, hom(g))
      S = relations(domain(phig))
      R = relations(codomain(phig))

      if ambient_space(R) === ambient_space(S)
        _glue = Vector{fmpq}[lift(g) + lift(phig(g)) for g in gens(domain(phig))]
        z = zero_matrix(QQ,0, degree(S))
        glue = reduce(vcat, [matrix(QQ, 1, degree(S), g) for g in _glue], init=z)
        glue = vcat(basis_matrix(S+R), glue)
        glue = FakeFmpqMat(glue)
        _B = hnf(glue)
        _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
        C2 = lattice(ambient_space(S), _B[end-rank(S)-rank(R)+1:end, :])
        fC2 = block_diagonal_matrix([fA, fB])
        _B = solve_left(reduce(vcat, basis_matrix.([A,B])), basis_matrix(C2))
        fC2 = _B*fC2*inv(_B)
      else
        _glue = Vector{fmpq}[lift(qAinD(SAinqA(g))) + lift(qBinD(SBinqB(phig(g)))) for g in gens(domain(phig))]
        z = zero_matrix(QQ,0,degree(S)+degree(R))
        glue = reduce(vcat, [matrix(QQ,1,degree(S)+degree(R),g) for g in _glue], init=z)
        glue = vcat(block_diagonal_matrix([basis_matrix(S), basis_matrix(R)]), glue)
        glue = FakeFmpqMat(glue)
        _B = hnf(glue)
        _B = QQ(1, denominator(glue))*change_base_ring(QQ, numerator(_B))
        C2 = lattice(ambient_space(cover(D)), _B[end-rank(S)-rank(R)+1:end, :])
        fC2 = block_diagonal_matrix([fA, fB])
        __B = solve_left(reduce(vcat, basis_matrix.([A,B])), basis_matrix(C2))
        fC2 = __B*fC2*inv(__B)
        @assert fC2*gram_matrix(C2)*transpose(fC2) == gram_matrix(C2)
      end

      if type(lattice_with_isometry(C2, fC2^p, ambient_representation=false)) != t
        continue
      end

      ext, _ = sub(D, D.(_glue))
      perp, _ = orthogonal_submodule(D, ext)
      disc = torsion_quadratic_module(cover(perp), cover(ext), modulus = modulus_bilinear_form(perp), modulus_qf = modulus_quadratic_form(perp))
      disc, discinD = sub(D, gens(disc))
      qC2 = discriminant_group(C2)
      ok, phi2 = is_isometric_with_isometry(qC2, disc)
      @assert ok

      C2fC2 = lattice_with_isometry(C2, fC2, ambient_representation=false)
      im2_phi, _ = sub(OSA, OSA.([compose(phig, compose(hom(g1), inv(phig))) for g1 in gens(imB)]))
      im3, _ = sub(imA, gens(intersect(imA, im2_phi)[1]))
      stab = Tuple{AutomorphismGroupElem{TorQuadMod}, AutomorphismGroupElem{TorQuadMod}}[(x, imB(compose(inv(phig), compose(hom(x), phig)))) for x in gens(im3)]
      stab = AutomorphismGroupElem{TorQuadMod}[OSAinOD(x[1])*OSBinOD(x[2]) for x in stab]
      stab = union(stab, kerA)
      stab = union(stab, kerB)
      stab = Oscar._orthogonal_group(discriminant_group(C2), [compose(phi2, compose(_restrict(g, discinD), inv(phi2))).map_ab.map for g in stab])
      set_attribute!(C2fC2, :image_centralizer_in_Oq, stab)
      push!(results, C2fC2)
    end
  end

  return results
end

##################################################################################
#
# Representatives of lattices with isometry
#
##################################################################################

function _ideals_of_norm(E, d)
  @assert E isa Hecke.NfRel
  K = base_field(E)
  OK = maximal_order(K)
  OE = maximal_order(E)
  DE = different(OE)
  ids = []
  primes = [] 
  for p in prime_divisors(d)
    v = valuation(d, p)
    for (P, _) in prime_decomposition(OK, p)
      if !is_coprime(DE, ideal(OE, P))
        P = prime_decomposition(OE, P)[1][1]
      else
        P = ideal(OE, P)
      end
      nv = valuation(absolute_norm(P), p)
      push!(primes, [P^e for e in 1:divrem(v, nv)[1]])
    end
  end
  for I in Hecke.cartesian_product_iterator(primes)
    I = prod(I)
    if absolute_norm(I) == d
      push!(ids, I)
    end
  end
  return ids
end

function _possible_signatures(s1, s2, E)
  @assert E isa Hecke.NfRel
  ok, q = Hecke.is_cyclotomic_type(E)
  @assert ok
  @assert iseven(s2)
  @assert divides(2*(s1+s2), euler_phi(q))[1]
  l = divexact(s2, 2)
  K = base_field(E)
  inf = real_places(K)
  s = length(inf)
  signs = Dict{typeof(inf[1]), Int}[]
  parts = Vector{Int}[]
  perm = AllPerms(s)
  for v in AllParts(l)
    if length(v) > s
      continue
    end
    while length(v) != s
      push!(v, 0)
    end
    for vv in perm
      v2 = v[vv.d]
      v2 in parts ? continue : push!(parts, v2)
    end
  end
  for v in parts
    push!(signs, Dict(a => b for (a,b) in zip(inf, v)))
  end
  return signs
end

function representatives(Lf::LatticeWithIsometry, m::Int = 1)
  @req m >= 1 "m must be a positive integer"
  @req Oscar._is_cyclotomic_polynomial(minpoly(Lf)) "Minimal polyomial must be irreducible and cyclotomic"
  L = lattice(Lf)
  rk = rank(L)
  d = det(L)
  n = order_of_isometry(Lf)
  t = type(Lf)
  s1, _, s2 = signature_tuple(L)
  reps = LatticeWithIsometry[]
  if n*m < 3
    gene = genera((s1, s2), ZZ(d), max_scale=scale(L), even=iseven(L))
    f = (-1)^(n*m+1)*identity_matrix(QQ, rk)
    for G in gene
      repre = representatives(G)
      append!(reps, [lattice_with_isometry(LL, f, check=false) for LL in repre])
    end
    return reps
  end
  ok, rk = divides(rk, euler_phi(n*m))
  if !ok
    return reps
  end
  
  gene = []
  E, b = cyclotomic_field_as_cm_extension(n*m, cached=false)
  Eabs, EabstoE = absolute_simple_field(E)
  DE = EabstoE(different(maximal_order(Eabs)))
  K = base_field(E)
  OE = maximal_order(E)
  OK = maximal_order(K)
  msE = E(scale(L)*n*m)*inv(DE)
  ndE = ZZ(d*(absolute_norm(n*m*inv(DE))^rk))
  detE = _ideals_of_norm(E, ndE)
  signatures = _possible_signatures(s1, s2, E)
  for d in detE, sign in signatures
    append!(gene, genera_hermitian(E, rk, sign, d))
  end
  gene = unique(gene)
  for g in gene
    for H in representatives(g)
      M = trace_lattice(H)
      @assert det(lattice(M)) == det
      @assert Oscar._is_cyclotomic_polynomial(minpoly(M))
      @assert order_of_isometry(M) == n*m
      iseven(lattice(M)) == iseven(L) ? push!(reps, M) : continue
      @assert type(lattice_with_isometry(lattice(M), ambient_isometry(M)^m)) == t
    end
  end
  return reps
end

function representatives(t::Dict, m::Integer = 1; check::Bool = true)
  @req m >= 1 "m must be a positive integer"
  ke = collect(keys(t))
  n = max(ke)
  if check
    @req Set(divisors(n)) == Set(ke) "t does not define a type for lattices with isometry"
    for i in ke
      @req t[i][2] isa ZGenus "t does not define a type for lattices with isometry"
      @req typeof(t[i][1]) <: Union{ZGenus, GenusHerm} "t does not define a type for lattices with isometry"
    end
    @req all(i -> rank(t[i][1]) == rank(t[i][2]) == 0, [i for i in ke if i != n]) "Minimal polynomial should be irreducible and cyclotomic"
    @req rank(t[n][2]) == rank(t[n][1])*euler_phi(n) "Top-degree types do not agree"
  end
  G = t[n]
  s1, s2 = signature_tuple(G)
  rk = s1+s2
  d = det(G)
  reps = LatticeWithIsometry[]
  if n*m < 3
    gene = genera((s1, s2), ZZ(d), max_scale=scale(G), even=iseven(G))
    f = (-1)^(n*m+1)*identity_matrix(QQ, rk)
    for g in gene
      repre = representatives(g)
      append!(reps, [lattice_with_isometry(LL, f, check=false) for LL in repre])
    end
    return reps
  end
  ok, rk = divides(rk, euler_phi(n*m))
  if !ok
    return reps
  end
  
  gene = []
  E, b = cyclotomic_field_as_cm_extension(n*m, cached=false)
  Eabs, EabstoE = absolute_simple_field(E)
  DE = EabstoE(different(maximal_order(Eabs)))
  K = base_field(E)
  OE = maximal_order(E)
  OK = maximal_order(K)
  msE = E(scale(L)*n*m)*inv(DE)
  ndE = ZZ(d*(absolute_norm(n*m*inv(DE))^rk))
  detE = _ideals_of_norm(E, ndE)
  signatures = _possible_signatures(s1, s2, E)
  for d in detE, sign in signatures
    append!(gene, genera_hermitian(E, rk, sign, d))
  end
  gene = unique(gene)
  for g in gene
    for H in representatives(g)
      M = trace_lattice(H)
      @assert det(lattice(M)) == det
      iseven(lattice(M)) == iseven(G) ? push!(reps, M) : continue
      @assert type(lattice_with_isometry(lattice(M), ambient_isometry(M)^m)) == t
    end
  end
  return reps
end

function split(Lf::LatticeWithIsometry, p::Int)
  @req is_prime(p) "p must be a prime number"
  @req Oscar._is_cyclotomic_polynomial(minpoly(Lf)) "Minimal polynomial must be irreducible and cyclotomic"
  ok, q, d = is_prime_power_with_data(order_of_isometry(Lf))
  @req ok "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"
  reps = LatticeWithIsometry[]
  atp = admissible_triples(Lf, p)
  for (A, B) in atp
    LA = lattice_with_isometry(representative(A))
    RA = representatives(LA, p*q^d)
    LB = lattice_with_isometry(representative(B))
    RB = representatives(LB, q^d)
    for (L1, L2) in Hecke.cartesian_product_iterator([RA, RB])
      E = primitive_extensions(L1, L2, Lf, p)
      @assert all(LL -> type(lattice_with_isometry(lattice(LL), ambient_isometry(LL)^p) == type(Lf)), E)
      append!(reps, E)
    end
  end
  return reps
end

function first_p(Lf::LatticeWithIsometry, p::Int, b::Int = 0)
  @req is_prime(p) "p must be a prime number"
  @req b in [0, 1] "b must be an integer equal to 0 or 1"
  ok, q, e = is_prime_power_with_data(order_of_isometry(Lf))
  @req ok "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"
  reps = LatticeWithIsometry[]
  if e == 0
    return split(Lf, p)
  end
  x = gen(Hecke.Globals.Qx)
  A0 = kernel_lattice(Lf, q^e)
  B0 = kernel_lattice(Lf, x^(q^e-1)-1)
  A = split(A0, p)
  B = first_p(B0, p)
  for (L1, L2) in Hecke.cartesian_product_iterator([A, B])
    b == 1 && !divides(order_of_isometry(L1), p)[1] && !divides(order_of_isometry(L2), p)[1] && continue
    E = primitive_extensions(L1, L2, Lf, q)
    @assert all(LL -> type(lattice_with_isometry(lattice(LL), ambient_isometry(LL)^p)) == type(Lf), E)
    @assert b == 0 || all(LL -> order_of_isometry(LL) == p*q^e, E)
    append!(reps, E)
  end
  return reps
end

function pure_up(Lf::LatticeWithIsometry, p::Int)
  @req is_prime(p) "p must be a prime number"
  @req order_of_isometry(Lf) > 0 "Isometry must be of finite order"
  n = order_of_isometry(Lf)
  pd = prime_divisors(n)
  @req length(pd) <= 2 && p in pd "Order must be divisible by p and have at most 2 prime factors"
  if length(pd) == 2
    q = pd[1] == p ? pd[2] : pd[1]
    d = valuation(n, p)
    e = valuation(n, q)
  else
    q = 1
    d = valuation(n, p)
    e = 0
  end
  phi = minpoly(Lf)
  x = gen(parent(phi))
  chi = prod([cyclotomic(x, p^d*q^i) for i=0:e])
  @req divides(chi, phi)[1] "Minimal polynomial is not of the correct form"
  reps = LatticeWithIsometry[]
  if e == 0
    return representatives(Lf, p)
  end
  A0 = kernel_lattice(Lf, p^d*q^e)
  bool, r = divides(phi, cyclotomic(x, p^d*q^e))
  @assert bool
  B0 = kernel_lattice(Lf, r)
  A = representatives(A0 ,p)
  B = pure_up(B0, p)
  for (LA, LB) in Hecke.cartesian_product_iterator([A, B])
    E = extensions(LA, LB, Lf, q)
    @assert all(LL -> type(lattice_with_isometry(lattice(LL), ambient_isometry(LL)^p)) == type(Lf), E)
    append!(reps, E)
  end
  return reps
end

function next_p(Lf::LatticeWithIsometry, p)
  n = order_of_isometry(Lf)
  @req n > 0 "Isometry must be of finite order"
  pd = prime_divisors(n)
  @req length(pd) <= 2 "Order must have at most 2 prime divisors"
  if !(p in pd)
    return first_p(Lf, p, 1)
  end
  d = valuation(n, p)
  if n != p^d
    _, q, e = is_prime_power_with_data(divexact(n, p^d))
  else
    q = 1
    e = 0
  end
  x = gen(parent(minpoly(Lf)))
  B0 = kernel_lattice(Lf, x^(divexact(n, p)) - 1)
  A0 = kernel_lattice(Lf, prod([cyclotomic(x, p^d*q^i) for i in 0:e]))
  A = pure_up(A0, p)
  B = next_p(B0, p)
  for (LA, LB) in Hecke.cartesian_product_iterator([A, B])
    E = extensions(LA, LB, Lf, p)
    @assert all(LL -> order_of_isometry(LL) == p^(d+1)*q^e, E)
    @assert all(LL -> type(lattice_with_isometry(lattice(LL), ambient_isometry(LL)^p)) == type(Lf), E)
    append!(reps, E)
  end
  return reps
end

