export admissible_triples,
       is_admissible_triple,
       prime_splitting,
       prime_splitting_of_prime_power,
       prime_splitting_of_pure_type_prime_power,
       prime_splitting_of_semi_pure_type,
       representatives_of_pure_type

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
function _find_L(r::Int, d::Hecke.RationalUnion, s::fmpz, l::fmpz, p::IntegerUnion, even = true)
  L = ZGenus[]
  if r == 0 && d == 1
    return ZGenus[genus(Zlattice(gram = matrix(QQ, 0, 0, [])))]
  end
  for (s1,s2) in [(s,t) for s=0:r for t=0:r if s+t==r]
    gen = genera((s1,s2), d, even=even)
    filter!(G -> divides(numerator(scale(G)), s)[1], gen)
    filter!(G -> divides(p*l, numerator(level(G)))[1], gen)
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
  AperpB = orthogonal_sum(A, B)
  (signature_tuple(AperpB) == signature_tuple(C)) || (return false)
  if ((A == zg) && (B == C)) || ((B == zg) && (A == C))
    # C can be always glued with the empty genus to obtain C
    return true
  elseif (A == zg) || (B == zg)
    # If A or B is empty but the other is not C, then there is no glueing
    return false
  end

  @req divides(rank(B), p-1)[1] "p-1 must divide the rank of B"

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
  if !divides(numerator(det(C)), p)[1]
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
  
    ABr = orthogonal_sum(Ar, Br)
    
    for i = 1:l
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

function is_admissible_triple(A::T, B::T, C::T, p::Integer) where T <: Union{ZLat, LatticeWithIsometry}
  L = ZGenus[genus(D) for D = (A, B, C)]
  return is_admissible_triple(L[1], L[2], L[3], p)
end

@doc Markdown.doc"""
    admissible_triples(C::ZGenus, p::Integer) -> Vector{Tuple{ZGenus, ZGenus}}

Given a $\mathbb Z$-genus `C` and a prime number `p`, return all tuples of
$\mathbb Z$-genera `(A, B)` such that `(A, B, C)` is `p`-admissible and
`B` is of rank divisible by $p-1$.

See Algorithm 1 of [BH22].
"""
function admissible_triples(G::ZGenus, p::Int64)
  @req is_prime(p) "p must be a prime number"
  @req is_integral(G) "G must be a genus of integral lattices"
  n = rank(G)
  d = numerator(det(G))
  even = iseven(G)
  L = Tuple{ZGenus, ZGenus}[]
  for ep in 0:div(n, p-1)
    rp = (p-1)*ep
    r1 = n - rp
    m = min(ep, r1) 
    D = _find_D(d, m, p)
    for (d1, dp) in D
      L1 = _find_L(r1, d1, numerator(scale(G)), numerator(level(G)), p, even)
      Lp = _find_L(rp, dp, numerator(scale(G)), numerator(level(G)), p, even)
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
# Representatives of lattices with isometry
#
##################################################################################

# we compute ideals of E/K whose absolute norm is equal to d

function _ideals_of_norm(E, d::fmpq)
  if denominator(d) == 1
    return _ideals_of_norm(E, numerator(d))
  elseif numerator(d) == 1
    return [inv(I) for I in _ideals_of_norm(E, denominator(d))]
  else
    return [I*inv(J) for (I, J) in Hecke.cartesian_product_iterator([_ideals_of_norm(E, numerator(d)), _ideals_of_norm(E, denominator(d))])]
  end
end

function _ideals_of_norm(E, d::fmpz)
  isone(d) && return [1*maximal_order(E)]
  @assert E isa Hecke.NfRel
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
      push!(ids, I)
    end
  end
  return ids
end

# given a cyclotomic field (as cm extension) E/K, return all
# the possible signatures dictionnaries of any hermitian lattice over
# E/K of rank rk, whose trace lattice has signature (s1, s2).

function _possible_signatures(s1, s2, E, rk)
  @assert E isa Hecke.NfRel
  ok, q = Hecke.is_cyclotomic_type(E)
  @assert ok
  @assert iseven(s2)
  @assert divides(2*(s1+s2), euler_phi(q))[1]
  n = 
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

@doc Markdown.doc"""
    representatives_of_pure_type(Lf::LatticeWithIsometry, m::Int = 1)
                                            -> Vector{LatticeWithIsometry}

Given a lattice with isometry $(L, f)$ of pure type (i.e. the minimal
polynomial of `f` is irreducible cyclotomic), and a positive integer `m` (set to
1 by default), return a set of representatives of isomorphism classes of lattices
with isometry of pure type $(M, g)$ and such that the type of $(B, g^m)$ is equal
to the type of $(L, f)$. Note that in this case, the isometries `g`'s are of order
$nm$.

See Algorithm 3 of [BH22].
"""
function representatives_of_pure_type(Lf::LatticeWithIsometry, m::Int = 1)
  rank(Lf) == 0 && return LatticeWithIsometry[]

  @req m >= 1 "m must be a positive integer"
  @req is_of_pure_type(Lf) "Minimal polyomial must be irreducible and cyclotomic"

  L = lattice(Lf)
  rk = rank(L)
  d = det(L)
  n = order_of_isometry(Lf)
  s1, _, s2 = signature_tuple(L)

  reps = LatticeWithIsometry[]

  if n*m < 3
    @info "Order smaller than 3"
    f = (-1)^(n*m+1)*identity_matrix(QQ, rk)
    G = genus(Lf)
    repre = representatives(G)
    @info "$(length(repre)) representative(s)"
    for LL in repre
      is_of_same_type(Lf, lattice_with_isometry(LL, f^m, check=false)) && push!(reps, lattice_with_isometry(LL, f, check=false))
    end
    return reps
  end

  @info "Order bigger than 3"

  ok, rk = divides(rk, euler_phi(n*m))

  ok || return reps

  gene = []
  E, b = cyclotomic_field_as_cm_extension(n*m)
  Eabs, EabstoE = absolute_simple_field(E)
  DE = EabstoE(different(maximal_order(Eabs)))

  @info "We have the different"

  ndE = d*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)

  @info "All possible ideal dets: $(length(detE))"

  signatures = _possible_signatures(s1, s2, E, rk)

  @info "All possible signatures: $(length(signatures))"

  for dd in detE, sign in signatures
    append!(gene, genera_hermitian(E, rk, sign, dd, min_scale = inv(DE), max_scale = numerator(DE*dd)))
  end
  gene = unique(gene)

  @info "All possible genera: $(length(gene))"

  for g in gene
    @info "g = $g"
    H = representative(g)
    if !is_integral(DE*scale(H))
      continue
    end
    @info "$H"
    M = trace_lattice(H)
    @assert det(lattice(M)) == d
    @assert is_cyclotomic_polynomial(minpoly(M))
    @assert order_of_isometry(M) == n*m
    if iseven(lattice(M)) != iseven(L)
      continue
    end
    if !is_of_same_type(Lf, lattice_with_isometry(lattice(M), ambient_isometry(M)^m))
      continue
    end
    append!(reps, [trace_lattice(HH) for HH in genus_representatives(H)])
  end
  return reps
end

@doc Markdown.doc"""
    representatives_of_pure_type(t::Dict, m::Int = 1; check::Bool = true)
                                          -> Vector{LatticeWithIsometry}

Given a type `t` for lattices with isometry of pure type (i.e. the minimal
polymomial of the associated isometry is irreducible cyclotomic) and an intger
`m` (set to 1 by default), return a set of representatives of isomorphism
classes of lattices with isometry of pure type $(L, f)$ such that the
type of $(L, f^m)$ is equal to `t`.

If `check === true`, then `t` is checked to be pure. Note that `n` can be 1.

See Algorithm 3 of [BH22].
"""
function representatives_of_pure_type(t::Dict, m::Integer = 1; check::Bool = true)
  M = _representative(t, check = check)
  M === nothing && return LatticeWithIsometry[]
  return representatives_of_pure_type(M, m)
end


function _representative(t::Dict; check::Bool = true)
  !check || is_pure(t) || error("t must be pure")
  
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

  ok, rk = divides(rk, euler_phi(n))

  ok || reps
  
  gene = []
  E, b = cyclotomic_field_as_cm_extension(n, cached=false)
  Eabs, EabstoE = absolute_simple_field(E)
  DE = EabstoE(different(maximal_order(Eabs)))

  ndE = d*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)

  @info "All possible ideal dets: $(length(detE))"

  signatures = _possible_signatures(s1, s2, E)

  @info "All possible signatures: $(length(signatures))"

  for dd in detE, sign in signatures
    append!(gene, genera_hermitian(E, rk, sign, dd, min_scale = inv(DE), max_scale = numerator(DE*dd)))
  end
  gene = unique(gene)

  for g in gene
    H = representative(g)
    if !is_integral(DE*scale(H))
      continue
    end
    H = H
    M = trace_lattice(H)
    @assert det(lattice(M)) == d
    @assert is_cyclotomic_polynomial(minpoly(M))
    @assert order_of_isometry(M) == n
    if iseven(lattice(M)) != iseven(G)
      continue
    end
    if !is_of_type(M, t)
      continue
    end
    return M
  end
  return nothing
end

@doc Markdown.doc"""
    prime_splitting_of_pure_type_prime_power(Lf::LatticeWithIsometry, p::Int)
                                                  -> Vector{LatticeWithIsometry}

Given a lattice with isometry $(L, f)$ of pure type (i.e. the minimal polynomial
of `f` is irreducible cyclotomic) with `f` of order $q^d$ for some prime number `q`,
and a prime number $p \neq q$, return a set of representatives of the isomorphisms
classes of lattices with isometry $(M, g)$ such that the type of $(M, g^p)$ is equal
to the type of $(L, f)$.

Note that `d` can be 0.

See Algorithm 4 of [BH22].
"""
function prime_splitting_of_pure_type_prime_power(Lf::LatticeWithIsometry, p::Int)
  rank(Lf) == 0 && return LatticeWithIsometry[]

  @req is_prime(p) "p must be a prime number"
  @req is_of_pure_type(Lf) "Minimal polynomial must be irreducible and cyclotomic"

  ok, q, d = is_prime_power_with_data(order_of_isometry(Lf))

  @req ok || d == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = LatticeWithIsometry[]
  @info "Compute admissible triples"
  atp = admissible_triples(Lf, p)
  @info "$(atp) admissible triple(s)"
  for (A, B) in atp
    LB = lattice_with_isometry(representative(B))
    RB = representatives_of_pure_type(LB, p*q^d)
    if is_empty(RB)
      continue
    end
    LA = lattice_with_isometry(representative(A))
    RA = representatives_of_pure_type(LA, q^d)
    for (L1, L2) in Hecke.cartesian_product_iterator([RA, RB])
      E = admissible_equivariant_primitive_extensions(L1, L2, Lf, p, check=false)
      GC.gc()
      append!(reps, E)
    end
  end
  return reps
end


@doc Markdown.doc"""
    prime_splitting_of_pure_type_prime_power(t::Dict, p::Int)
                                                  -> Vector{LatticeWithIsometry}

Given a type `t` of lattice with isometry of pure type $(L, f)$ (i.e. the minimal
polynomial of `f` is irreducible cyclotomic) with `f` of order $q^d$ for some
prime number `q`, and a prime number $p \neq q$, return a set of representatives
of the isomorphisms classes of lattices with isometry $(M, g)$ such that the type
of $(M, g^p)$ is equal to `t`.

Note that `d` can be 0.

See Algorithm 4 of [BH22].
"""
function prime_splitting_of_pure_type_prime_power(t::Dict, p::Int)
  @req is_prime(p) "p must be a prime number"
  @req is_pure(t) "t must be pure"
  Lf = _representative(t)
  return prime_splitting_of_pure_type_prime_power(Lf, p)
end

@doc Markdown.doc"""
    prime_splitting_of_prime_power(Lf::LatticeWithIsometry, p::Int, b::Int = 0)
                                                     -> Vector{LatticeWithIsometry}

Given a lattice with isometry $(L, f)$ with `f` of order $q^e$ for some prime number
`q`, a prime number $p \neq q$ and an integer $b = 0, 1$, return a set of representatives
of the isomorphism classes of lattices with isometry $(M, g)$ such that the type of
$(M, g^p)$ is equal to the type of $(L, f)$. If `b == 1`, return only the lattices
with isometry $(M, g)$ where `g` is of order $pq^e$.

Note that `e` can be 0.

See Algorithm 5 of [BH22].
"""
function prime_splitting_of_prime_power(Lf::LatticeWithIsometry, p::Int, b::Int = 0)
  rank(Lf) == 0 && return LatticeWithIsometry[]

  @req is_prime(p) "p must be a prime number"
  @req b in [0, 1] "b must be an integer equal to 0 or 1"

  ok, q, e = is_prime_power_with_data(order_of_isometry(Lf))

  @req ok || e == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = LatticeWithIsometry[]

  if e == 0
    reps = prime_splitting_of_pure_type_prime_power(Lf, p)
    (b == 1) && filter!(M -> order_of_isometry(M) == p, reps)
    return reps
  end

  x = gen(Hecke.Globals.Qx)
  A0 = kernel_lattice(Lf, q^e)
  B0 = kernel_lattice(Lf, x^(q^e-1)-1)
  A = prime_splitting_of_pure_type_prime_power(A0, p)
  is_empty(A) && return reps
  B = prime_splitting_of_prime_power(B0, p)
  for (L1, L2) in Hecke.cartesian_product_iterator([A, B])
    b == 1 && !divides(order_of_isometry(L1), p)[1] && !divides(order_of_isometry(L2), p)[1] && continue
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, q, check=false)
    @assert b == 0 || all(LL -> order_of_isometry(LL) == p*q^e, E)
    append!(reps, E)
  end
  return reps
end

@doc Markdown.doc"""
    prime_splitting_of_semi_pure_type(Lf::LatticeWithIsometry, p::Int)
                                                 -> Vector{LatticeWithIsometry}

Given a lattice with isometry $(L, f)$ and a prime number `p`, such that
the minimal of `f` divides $\prod_{i=0}^e\Phi_{p^dq^i}(f)$ for some
$d > 0$ and $e \geq 0$, return a set of representatives of the isomorphism classes
of lattices with isometry $(M, g)$ such that the type of $(M, g^p)$ is equal to the type
of $(L, f)$.

Note that `e` can be 0, while `d` has to be positive.

See Algorithm 6 of [BH22].
"""
function prime_splitting_of_semi_pure_type(Lf::LatticeWithIsometry, p::Int)
  rank(Lf) == 0 && return LatticeWithIsometry[]

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
  x = gen(parent(phi))
  chi = prod([cyclotomic_polynomial(p^d*q^i, parent(phi)) for i=0:e])

  @req divides(chi, phi)[1] "Minimal polynomial is not of the correct form"

  reps = LatticeWithIsometry[]

  if e == 0
    return representatives_of_pure_type(Lf, p)
  end

  A0 = kernel_lattice(Lf, p^d*q^e)
  bool, r = divides(phi, cyclotomic_polynomial(p^d*q^e, parent(phi)))
  @assert bool

  B0 = kernel_lattice(Lf, r)
  A = representatives_of_pure_type(A0, p)
  is_empty(A) && return reps
  B = prime_splitting_of_semi_pure_type(B0, p)
  for (LA, LB) in Hecke.cartesian_product_iterator([A, B])
    E = admissible_equivariant_primitive_extensions(LA, LB, Lf, q, check = false)
    append!(reps, E)
  end
  return reps
end

@doc Markdown.doc"""
    prime_splitting(Lf::LatticeWithIsometry, p::Int) -> Vector{LatticeWithIsometry}

Given a lattice with isometry $(L, f)$ and a prime number `p` such that
`f` is of order $p^dq^e$ for some prime number $q \neq p$, return a set
of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ of order $p^{d+1}q^e$ such that the type of $(M, g^p)$ is equal
to the type of $(L, f)$.

Note that `d` and `e` can be both zero.

See Algorithm 7 of [BH22].
"""
function prime_splitting(Lf::LatticeWithIsometry, p::Int)
  rank(Lf) == 0 && return LatticeWithIsometry[]

  n = order_of_isometry(Lf)

  @req is_finite(n) "Isometry must be of finite order"

  pd = prime_divisors(n)

  @req length(pd) <= 2 "Order must have at most 2 prime divisors"

  if !(p in pd)
    return prime_splitting_of_prime_power(Lf, p, 1)
  end

  d = valuation(n, p)
  if n != p^d
    _, q, e = is_prime_power_with_data(divexact(n, p^d))
  else
    q = 1
    e = 0
  end

  reps = LatticeWithIsometry[]

  x = gen(parent(minpoly(Lf)))
  B0 = kernel_lattice(Lf, x^(divexact(n, p)) - 1)
  A0 = kernel_lattice(Lf, prod([cyclotomic_polynomial(p^d*q^i) for i in 0:e]))
  A = prime_splitting_of_semi_pure_type(A0, p)
  isempty(A) && return reps
  B = prime_splitting(B0, p)
  for (LA, LB) in Hecke.cartesian_product_iterator([A, B])
    E = admissible_equivariant_primitive_extensions(LA, LB, Lf, p, check = false)
    @assert all(LL -> order_of_isometry(LL) == p^(d+1)*q^e, E)
    append!(reps, E)
  end
  return reps
end

