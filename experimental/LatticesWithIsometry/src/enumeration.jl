
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
function _tuples_divisors(d::T) where T <: Hecke.IntegerUnion
  div = divisors(d)
  return Tuple{T, T}[(dd,abs(divexact(d,dd))) for dd in div]
end

# This is line 8 of Algorithm 1, they correspond to the possible
# discriminant for the genera A and B to glue to fit in C. d is
# the determinant of C, m the maximal p-valuation of the gcd of
# d1 and dp.
function _find_D(d::T, m::Int, p::Int) where T <: Hecke.IntegerUnion
  @hassert :LatWithIsom 1 is_prime(p)
  @hassert :LatWithIsom 1 d != 0

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
function _find_L(pG::Int, nG::Int, r::Int, d::Hecke.RationalUnion, s::ZZRingElem, l::ZZRingElem, p::Hecke.IntegerUnion, even = true; pos::Int = -1)
  L = ZGenus[]
  if r == 0 && d == 1
    return ZGenus[genus(Zlattice(gram = matrix(QQ, 0, 0, [])))]
  end
  if pos >= 0
    neg = r-pos
    gen = Zgenera((pos, neg), d, even=even)
    filter!(G -> Hecke.divides(numerator(scale(G)), s)[1], gen)
    filter!(G -> Hecke.divides(p*l, numerator(level(G)))[1], gen)
    append!(L, gen)
  else
    for (s1,s2) in [(s,t) for s=0:pG for t=0:nG if s+t==r]
      gen = Zgenera((s1,s2), d, even=even)
      filter!(G -> Hecke.divides(numerator(scale(G)), s)[1], gen)
      filter!(G -> Hecke.divides(p*l, numerator(level(G)))[1], gen)
      append!(L, gen)
    end
  end
  return L
end

@doc raw"""
    is_admissible_triple(A::ZGenus, B::ZGenus, C::ZGenus, p::Integer) -> Bool

Given a triple of $\mathbb Z$-genera `(A,B,C)` and a prime number `p`, such
that the rank of `B` is divisible by $p-1$ and the level of `C` is a power
of `p`, return whether `(A,B,C)` is `p`-admissible in the sense of
Definition 4.13. [BH22]
"""
function is_admissible_triple(A::ZGenus, B::ZGenus, C::ZGenus, p::Integer)
  zg = genus(Zlattice(gram = matrix(QQ, 0, 0, [])))
  AperpB = direct_sum(A, B)
  (signature_tuple(AperpB) == signature_tuple(C)) || (return false)
  if ((A == zg) && (B == C)) || ((B == zg) && (A == C))
    # C can be always glued with the empty genus to obtain C
    return true
  elseif (A == zg) || (B == zg)
    # If A or B is empty but the other is not C, then there is no glueing
    return false
  end

  @req Hecke.divides(rank(B), p-1)[1] "p-1 must divide the rank of B"

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

  if !(Hecke.divides(scale(AperpB), scale(C))[1] && Hecke.divides(p*level(C), level(AperpB))[1])
    return false
  end

  # At this point, if C is unimodular at p, the glueing condition is equivalent to have 
  # an anti-isometry between the p-part of the (quadratic) discriminant forms of A and B
  qA = discriminant_group(A)
  qB = discriminant_group(B)
  if !Hecke.divides(numerator(det(C)), p)[1]
    return is_anti_isometric_with_anti_isometry(primary_part(qA, p)[1], primary_part(qB, p)[1])[1]
  end

  l = valuation(level(C), p)
  Ap = local_symbol(A, p)
  Bp = local_symbol(B, p)
  a_max = symbol(Ap, l+1)[2]
  b_max = symbol(Bp, l+1)[2]
  # For the glueing, rho_{l+1}(A_p) and rho_{l+1}(B_p) are anti-isometric, so they must have the
  # same order
  if a_max != b_max
    return false
  end
  
  # Since p^l*A^\vee/A is in the glue, its order is less than the order of the glue
  if g < a_max
    return false
  end

  a1 = symbol(Ap, 1)[2]
  a2 = rank(Ap) - a1 - symbol(Ap, 0)[2]
  b1 = symbol(Bp, 1)[2]
  b2 = rank(Bp) - b1 - symbol(Bp, 0)[2]

  ABp = symbol(local_symbol(AperpB, p))
  Cp = local_symbol(C, p)

  if a_max == g
    if length(symbol(Ap)) > 1
      Ar = ZpGenus(p, symbol(Ap)[1:end-1])
    else
      Ar = genus(matrix(ZZ,0,0,[]), p)
    end
   
    if length(symbol(Bp)) > 1
      Br = ZpGenus(p, symbol(Bp)[1:end-1])
    else
      Br = genus(matrix(ZZ, 0, 0, []), p)
    end
  
    ABr = direct_sum(Ar, Br)

    for i = 0:l-1
      s1 = symbol(ABr, i)
      s2 = symbol(Cp, i)
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
    s1 = symbol(Ap, s3[1]+1)
    s2 = symbol(Bp, s3[1]+1)
    if s1[4] != s2[4]
      return false
    end
  end
  Cp = symbol(Cp)
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

function is_admissible_triple(A::T, B::T, C::T, p::Integer) where T <: Union{ZLat, LatWithIsom}
  L = ZGenus[genus(D) for D = (A, B, C)]
  return is_admissible_triple(L[1], L[2], L[3], p)
end

@doc raw"""
    admissible_triples(C::ZGenus, p::Integer) -> Vector{Tuple{ZGenus, ZGenus}}

Given a $\mathbb Z$-genus `C` and a prime number `p`, return all tuples of
$\mathbb Z$-genera `(A, B)` such that `(A, B, C)` is `p`-admissible and
`B` is of rank divisible by $p-1$.

See Algorithm 1 of [BH22].
"""
function admissible_triples(G::ZGenus, p::Int64; pA::Int = -1, pB::Int = -1)
  @req is_prime(p) "p must be a prime number"
  @req is_integral(G) "G must be a genus of integral lattices"
  n = rank(G)
  pG, nG = signature_pair(G)
  if pA >= 0
    @req pA <= pG "Wrong restrictions"
    if pB >= 0
      @req pA + pB == pG "Wrong restrictions"
    else
      pB = pG - pA
    end
  elseif pB >= 0
    @req pB <= pG "Wrong restrictions"
    pA = pG - pB
  end
  d = numerator(det(G))
  even = iseven(G)
  L = Tuple{ZGenus, ZGenus}[]
  for ep in 0:div(n, p-1)
    rp = (p-1)*ep
    if pB >= 0
      rp >= pB || continue
    end
    r1 = n - rp
    if pA >= 0
      r1 >= pA || continue
    end
    m = min(ep, r1) 
    D = _find_D(d, m, p)
    for (d1, dp) in D
      L1 = _find_L(pG, nG, r1, d1, numerator(scale(G)), numerator(level(G)), p, even, pos = pA)
      Lp = _find_L(pG, nG, rp, dp, numerator(scale(G)), numerator(level(G)), p, even, pos = pB)
      for (A, B) in [(A, B) for A in L1 for B in Lp]
        if is_admissible_triple(A, B, G, p)
          push!(L, (A, B))
        end
      end
    end
  end
  return L
end

admissible_triples(L::T, p::Integer; pA::Int = -1, pB::Int = -1) where T <: Union{ZLat, LatWithIsom} = admissible_triples(genus(L), p, pA = pA, pB = pB)

##################################################################################
#
# Representatives of lattices with isometry
#
##################################################################################

# we compute ideals of E/K whose absolute norm is equal to d

function _ideals_of_norm(E, d::QQFieldElem)
  if denominator(d) == 1
    return _ideals_of_norm(E, numerator(d))
  elseif numerator(d) == 1
    return [inv(I) for I in _ideals_of_norm(E, denominator(d))]
  else
    return [I*inv(J) for (I, J) in Hecke.cartesian_product_iterator([_ideals_of_norm(E, numerator(d)), _ideals_of_norm(E, denominator(d))])]
  end
end

function _ideals_of_norm(E, d::ZZRingElem)
  isone(d) && return [fractional_ideal(maximal_order(E), one(E))]
  @hassert :LatWithIsom 1 E isa Hecke.NfRel
  K = base_field(E)
  OK = maximal_order(K)
  OE = maximal_order(E)
  DE = different(OE)
  ids = []
  primes = [] 
  for p in prime_divisors(d)
    v = valuation(d, p)
    pd = [P[1] for P in prime_decomposition(OK, p)]
    for i in 1:length(pd)
      if !is_coprime(DE, ideal(OE, pd[i]))
        P = prime_decomposition(OE, pd[i])[1][1]
      else
        P = ideal(OE, pd[i])
      end
      nv = valuation(norm(P), pd[i])
      push!(primes, [P^e for e in 0:divrem(v, nv)[1]])
    end
  end
  for I in Hecke.cartesian_product_iterator(primes)
    I = prod(I)
    if absolute_norm(I) == d
      push!(ids, fractional_ideal(OE, I))
    end
  end
  return ids
end

# given a cyclotomic field (as cm extension) E/K, return all
# the possible signatures dictionnaries of any hermitian lattice over
# E/K of rank rk, whose trace lattice has signature (s1, s2).

function _possible_signatures(s1, s2, E, rk)
  @hassert :LatWithIsom 1 E isa Hecke.NfRel
  ok, q = Hecke.is_cyclotomic_type(E)
  @hassert :LatWithIsom 1 ok
  @hassert :LatWithIsom 1 iseven(s2)
  @hassert :LatWithIsom 1 Hecke.divides(2*(s1+s2), euler_phi(q))[1]
  l = divexact(s2, 2)
  K = base_field(E)
  inf = real_places(K)
  s = length(inf)
  signs = Dict{typeof(inf[1]), Int}[]
  parts = Vector{Int}[]
  perm = AllPerms(s)
  for v in AllParts(l)
    if any(i -> i > rk, v)
      continue
    end
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

@doc raw"""
    representatives_of_hermitian_type(Lf::LatWithIsom, m::Int = 1)
                                            -> Vector{LatWithIsom}

Given a lattice with isometry $(L, f)$ of hermitian type (i.e. the minimal polynomial 
of `f` is irreducible cyclotomic), and a positive integer `m`, return a set of
representatives of isomorphism classes of lattices with isometry of hermitian
type $(M, g)$ and such that the type of $(B, g^m)$ is equal to the type of
$(L, f)$. Note that in this case, the isometries `g`'s are of order $nm$.

See Algorithm 3 of [BH22].
"""
function representatives_of_hermitian_type(Lf::LatWithIsom, m::Int = 1)
  rank(Lf) == 0 && return LatWithIsom[Lf]

  @req m >= 1 "m must be a positive integer"
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian"

  rk = rank(Lf)
  d = det(Lf)
  n = order_of_isometry(Lf)
  s1, _, s2 = signature_tuple(Lf)

  reps = LatWithIsom[]

  if n*m < 3
    @vprint :LatWithIsom 1 "Order smaller than 3\n"
    f = (-1)^(n*m+1)*identity_matrix(QQ, rk)
    G = genus(Lf)
    repre = representatives(G)
    @vprint :LatWithIsom 1 "$(length(repre)) representative(s)\n"
    for LL in repre
      is_of_same_type(Lf, lattice_with_isometry(LL, f^m, check=false)) && push!(reps, lattice_with_isometry(LL, f, check=false))
    end
    return reps
  end

  !iseven(s2) && return reps

  @vprint :LatWithIsom 1 "Order bigger than 3\n"

  ok, rk = Hecke.divides(rk, euler_phi(n*m))

  ok || return reps

  gene = Hecke.HermGenus[]
  E, b = cyclotomic_field_as_cm_extension(n*m)
  Eabs, EabstoE = absolute_simple_field(E)
  DE = EabstoE(different(maximal_order(Eabs)))

  @vprint :LatWithIsom 1 "We have the different\n"

  ndE = d*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)

  @vprint :LatWithIsom 1 "All possible ideal dets: $(length(detE))\n"

  signatures = _possible_signatures(s1, s2, E, rk)

  @vprint :LatWithIsom 1 "All possible signatures: $(length(signatures))\n"
  for dd in detE, sign in signatures
    append!(gene, genera_hermitian(E, rk, sign, dd, min_scale = inv(DE), max_scale = numerator(dd)*DE))
  end
  gene = unique(gene)

  @vprint :LatWithIsom 1 "All possible genera: $(length(gene))\n"
  for g in gene
    @vprint :LatWithIsom 1 "g = $g\n"
    H = representative(g)
    if !is_integral(DE*scale(H))
      continue
    end
    if is_even(Lf) && !is_integral(different(fixed_ring(H))*norm(H))
      continue
    end
    @vprint :LatWithIsom 1 "$H\n"
    M, fM = Hecke.trace_lattice_with_isometry(H)
    det(M) == d || continue
    M = lattice_with_isometry(M, fM)
    @hassert :LatWithIsom 1 is_of_hermitian_type(M)
    @hassert :LatWithIsom 1 order_of_isometry(M) == n*m
    if is_even(M) != is_even(Lf)
      continue
    end
    if !is_of_same_type(Lf, lattice_with_isometry(lattice(M), ambient_isometry(M)^m))
      continue
    end
    gr = genus_representatives(H)
    for HH in gr
      M, fM = Hecke.trace_lattice_with_isometry(HH)
      push!(reps, lattice_with_isometry(M, fM))
    end
  end
  return reps
end

@doc raw"""
    representatives_of_hermitian_type(t::Dict, m::Int = 1; check::Bool = true)
                                                             -> Vector{LatWithIsom}

Given a hermitian type `t` for lattices with isometry (i.e. the minimal
polymomial of the associated isometry is irreducible cyclotomic) and an intger
`m` (set to 1 by default), return a set of representatives of isomorphism
classes of lattices with isometry of hermitian type $(L, f)$ such that the
type of $(L, f^m)$ is equal to `t`.

If `check === true`, then `t` is checked to be hermitian. Note that `n` can be 1.

See Algorithm 3 of [BH22].
"""
function representatives_of_hermitian_type(t::Dict, m::Integer = 1; check::Bool = true)
  M = _representative(t, check = check)
  M === nothing && return LatWithIsom[]
  return representatives_of_hermitian_type(M, m)
end

function _representative(t::Dict; check::Bool = true)
  !check || is_hermitian(t) || error("t must be pure")
  
  ke = collect(keys(t))
  n = maximum(ke)

  G = t[n][2]
  s1, s2 = signature_tuple(G)
  rk = s1+s2
  d = det(G)

  if n < 3
    L = representative(G)
    return trace_lattice(L, order = n)
  end

  ok, rk = Hecke.divides(rk, euler_phi(n))

  ok || reps
  
  gene = HermGenus[]
  E, b = cyclotomic_field_as_cm_extension(n, cached=false)
  Eabs, EabstoE = absolute_simple_field(E)
  DE = EabstoE(different(maximal_order(Eabs)))

  ndE = d*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)

  @vprint :LatWithIsom 1 "All possible ideal dets: $(length(detE))\n"

  signatures = _possible_signatures(s1, s2, E)

  @vprint :LatWithIsom 1 "All possible signatures: $(length(signatures))\n"

  for dd in detE, sign in signatures
    append!(gene, genera_hermitian(E, rk, sign, dd, min_scale = inv(DE), max_scale = numerator(DE*dd)))
  end
  gene = unique(gene)

  for g in gene
    H = representative(g)
    if !is_integral(DE*scale(H))
      continue
    end
    if iseven(Lf) && !is_integral(different(fixed_ring(H))*norm(H))
      continue
    end
    H = H
    M = trace_lattice(H)
    det(M) == d || continue
    @hassert :LatWithIsom 1 is_of_hermitian_type(M)
    @hassert :LatWithIsom 1 order_of_isometry(M) == n
    if iseven(M) != iseven(G)
      continue
    end
    if !is_of_type(M, t)
      continue
    end
    return M
  end
  return nothing
end

@doc raw"""
    splitting_of_hermitian_prime_power(Lf::LatWithIsom, p::Int) -> Vector{LatWithIsom}

Given a lattice with isometry $(L, f)$ of hermitian type with `f` of order $q^d$
for some prime number `q`, and given another prime number $p \neq q$, return a
set of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

Note that `d` can be 0.

See Algorithm 4 of [BH22].
"""
function splitting_of_hermitian_prime_power(Lf::LatWithIsom, p::Int; pA::Int = -1, pB::Int = -1)
  rank(Lf) == 0 && return LatWithIsom[Lf]

  @req is_prime(p) "p must be a prime number"
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian type"

  ok, q, d = is_prime_power_with_data(order_of_isometry(Lf))

  @req ok || d == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = LatWithIsom[]
  @vprint :LatWithIsom 1 "Compute admissible triples\n"
  atp = admissible_triples(Lf, p, pA = pA, pB = pB)
  @vprint :LatWithIsom 1 "$(length(atp)) admissible triple(s)\n"
  for (A, B) in atp
    LB = lattice_with_isometry(representative(B))
    RB = representatives_of_hermitian_type(LB, p*q^d)
    if is_empty(RB)
      continue
    end
    LA = lattice_with_isometry(representative(A))
    RA = representatives_of_hermitian_type(LA, q^d)
    for (L1, L2) in Hecke.cartesian_product_iterator([RA, RB])
      E = try admissible_equivariant_primitive_extensions(L1, L2, Lf, p)
      catch e return L1, L2, Lf, p
      end
      GC.gc()
      append!(reps, E)
    end
  end
  return reps
end

@doc raw"""
    splitting_of_hermitian_prime_power(t::Dict, p::Int) -> Vector{LatWithIsom}

Given a hermitian type `t` of lattice with isometry $(L, f)$ with `f` of order
$q^d$ for some prime number `q`, and given another prime number $p \neq q$,
return a set of representatives of the isomorphisms classes of lattices with
isometry $(M, g)$ such that the type of $(M, g^p)$ is equal to `t`.

Note that `d` can be 0.

See Algorithm 4 of [BH22].
"""
function splitting_of_hermitian_prime_power(t::Dict, p::Int)
  @req is_prime(p) "p must be a prime number"
  @req is_hermitian(t) "t must be hermitian"
  Lf = _representative(t)
  return splitting_of_hermitian_prime_power(Lf, p)
end

@doc raw"""
    splitting_of_prime_power(Lf::LatWithIsom, p::Int, b::Int = 0) -> Vector{LatWithIsom}

Given a lattice with isometry $(L, f)$ with `f` of order $q^e$ for some prime number
`q`, a prime number $p \neq q$ and an integer $b = 0, 1$, return a set of representatives
of the isomorphism classes of lattices with isometry $(M, g)$ such that the type of
$(M, g^p)$ is equal to the type of $(L, f)$. If `b == 1`, return only the lattices
with isometry $(M, g)$ where `g` is of order $pq^e$.

Note that `e` can be 0.

See Algorithm 5 of [BH22].
"""
function splitting_of_prime_power(Lf::LatWithIsom, p::Int, b::Int = 0)
  if rank(Lf) == 0
    (b == 0) && return LatWithIsom[Lf]
    return LatWithIsom[]
  end

  @req is_prime(p) "p must be a prime number"
  @req b in [0, 1] "b must be an integer equal to 0 or 1"

  ok, q, e = is_prime_power_with_data(order_of_isometry(Lf))

  @req ok || e == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = LatWithIsom[]

  if e == 0
    reps = splitting_of_hermitian_prime_power(Lf, p)
    (b == 1) && filter!(M -> order_of_isometry(M) == p, reps)
    return reps
  end

  x = gen(Hecke.Globals.Qx)
  A0 = kernel_lattice(Lf, q^e)
  B0 = kernel_lattice(Lf, x^(q^(e-1))-1)
  A = splitting_of_hermitian_prime_power(A0, p)
  is_empty(A) && return reps
  B = splitting_of_prime_power(B0, p)
  for (L1, L2) in Hecke.cartesian_product_iterator([A, B])
    b == 1 && !Hecke.divides(order_of_isometry(L1), p)[1] && !Hecke.divides(order_of_isometry(L2), p)[1] && continue
    E = try admissible_equivariant_primitive_extensions(L2, L1, Lf, q, p)
    catch e return L2, L1, Lf, q, p
    end
    @hassert :LatWithIsom 1 b == 0 || all(LL -> order_of_isometry(LL) == p*q^e, E)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_partial_mixed_prime_power(Lf::LatWithIsom, p::Int)
                                                 -> Vector{LatWithIsom}

Given a lattice with isometry $(L, f)$ and a prime number `p`, such that
the minimal polynomial of `f` divides $\prod_{i=0}^e\Phi_{p^dq^i}(f)$ for some
$d > 0$ and $e \geq 0$, return a set of representatives of the isomorphism classes
of lattices with isometry $(M, g)$ such that the type of $(M, g^p)$ is equal to the type
of $(L, f)$.

Note that `e` can be 0, while `d` has to be positive.

See Algorithm 6 of [BH22].
"""
function splitting_of_partial_mixed_prime_power(Lf::LatWithIsom, p::Int)
  rank(Lf) == 0 && return LatWithIsom[]

  @req is_prime(p) "p must be a prime number"
  @req is_finite(order_of_isometry(Lf)) "Isometry must be of finite order"

  n = order_of_isometry(Lf)
  pd = prime_divisors(n)

  @req 1 <= length(pd) <= 2 && p in pd "Order must be divisible by p and have at most 2 prime factors"

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
  chi = prod([cyclotomic_polynomial(p^d*q^i, parent(phi)) for i=0:e])

  @req Hecke.divides(chi, phi)[1] "Minimal polynomial is not of the correct form"

  reps = LatWithIsom[]

  if e == 0
    return splitting_of_hermitian_prime_power(Lf, p)
  end

  A0 = kernel_lattice(Lf, p^d*q^e)
  bool, r = Hecke.divides(phi, cyclotomic_polynomial(p^d*q^e, parent(phi)))
  @hassert :LatWithIsom 1 bool

  B0 = kernel_lattice(Lf, r)
  A = splitting_of_prime_power(A0, p)
  is_empty(A) && return reps
  B = splitting_of_partial_mixed_prime_power(B0, p)
  for (LA, LB) in Hecke.cartesian_product_iterator([A, B])
    E = admissible_equivariant_primitive_extensions(LB, LA, Lf, q, p)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_mixed_prime_power(Lf::LatWithIsom, p::Int)
                                          -> Vector{LatWithIsom}

Given a lattice with isometry $(L, f)$ and a prime number `p` such that
`f` is of order $p^dq^e$ for some prime number $q \neq p$, return a set
of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ of order $p^{d+1}q^e$ such that the type of $(M, g^p)$ is equal
to the type of $(L, f)$.

Note that `d` and `e` can be both zero.

See Algorithm 7 of [BH22].
"""
function splitting_of_mixed_prime_power(Lf::LatWithIsom, p::Int)
  rank(Lf) == 0 && return LatWithIsom[]

  n = order_of_isometry(Lf)

  @req is_finite(n) "Isometry must be of finite order"

  pd = prime_divisors(n)

  @req length(pd) <= 2 "Order must have at most 2 prime divisors"

  if !(p in pd)
    return splitting_of_prime_power(Lf, p, 1)
  end

  d = valuation(n, p)
  if n != p^d
    _, q, e = is_prime_power_with_data(divexact(n, p^d))
  else
    q = 1
    e = 0
  end

  reps = LatWithIsom[]

  x = gen(parent(minpoly(Lf)))
  B0 = kernel_lattice(Lf, x^(divexact(n, p)) - 1)
  A0 = kernel_lattice(Lf, prod([cyclotomic_polynomial(p^d*q^i) for i in 0:e]))
  A = splitting_of_partial_mixed_prime_power(A0, p)
  isempty(A) && return reps
  B = splitting_of_mixed_prime_power(B0, p)
  for (LA, LB) in Hecke.cartesian_product_iterator([A, B])
    E = admissible_equivariant_primitive_extensions(LB, LA, Lf, p)
    @hassert :LatWithIsom 1 all(LL -> order_of_isometry(LL) == p^(d+1)*q^e, E)
    append!(reps, E)
  end
  return reps
end

