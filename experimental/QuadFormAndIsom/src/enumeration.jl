#################################################################################
#
# This is an import to Oscar of the methods written following the paper [BH23] on
# "Finite subgroups of automorphisms of K3 surfaces".
#
##################################################################################

##################################################################################
#
# Admissible triples
#
##################################################################################

# We collect the pairs `(d', d/d')` for all divisors `d'` of `d`.
function _tuples_divisors(d::T) where T <: IntegerUnion
  return Tuple{T, T}[(dd, abs(divexact(d, dd))) for dd in divisors(d)]
end

# This is line 8 of Algorithm 1, they correspond to the possible
# discriminant for the genera A and B to glue to fit in C. d is
# the determinant of C, m the maximal p-valuation of the gcd of
# d1 and dp.
function _find_D(d::T, m::Int, p::IntegerUnion) where T <: IntegerUnion
  @hassert :ZZLatWithIsom 1 is_prime(p)
  @hassert :ZZLatWithIsom 1 d != 0

  # If m == 0, there are no conditions on the gcd of d1 and dp
  if m == 0
    return _tuples_divisors(d)
  end

  D = Tuple{T, T}[]
  # We try all the values of g possible, from 1 to p^m
  for g in powers(p, m)
    dj = _tuples_divisors(d*g^2)
    for (d1, dp) in dj
      if mod(d1, g) == mod(dp, g) == 0
        push!(D, (d1, dp))
      end
    end
  end
  return D
end

# This is line 10 of Algorithm 1. We need the condition on the parity of
# C since subgenera of an even genus are even too. r is the rank of
# the subgenus, d its determinant, s and l the scale and level of C
function _find_L(pG::Int, nG::Int, r::Int, d::RationalUnion, s::ZZRingElem, l::ZZRingElem, p::IntegerUnion, even = true; pos::Int=-1, neg::Int=-1)
  def = ZZGenus[genus(integer_lattice(; gram=matrix(QQ, 0, 0, [])))]
  if r == 0 && d == 1
    return def
  end
  if pos >= 0
    pos > pG && return def
    if neg >= 0
      pos + neg != r && return def
    else
      neg = r-pos
    end
    neg > nG && return def
    gen = integer_genera((pos, neg), d; even)
    filter!(G -> is_divisible_by(numerator(scale(G)), s), gen)
    filter!(G -> is_divisible_by(p*l, numerator(level(G))), gen)
  elseif neg >= 0
    neg > nG && return def
    pos = r-neg
    pos > pG && return def
    gen = integer_genera((pos, neg), d; even)
    filter!(G -> is_divisible_by(numerator(scale(G)), s), gen)
    filter!(G -> is_divisible_by(p*l, numerator(level(G))), gen)
  else
    gen = ZZGenus[]
    for s1 in max(0, r-nG):min(pG, r)
      s2 = r-s1
      L = integer_genera((s1, s2), d; even)
      filter!(G -> is_divisible_by(numerator(scale(G)), s), L)
      filter!(G -> is_divisible_by(p*l, numerator(level(G))), L)
      append!(gen, L)
    end
  end
  return gen
end

@doc raw"""
    is_admissible_triple(A::ZZGenus, B::ZZGenus, C::ZZGenus, p::Integer) -> Bool

Given a triple of $\mathbb Z$-genera $(A, B, C)$ and a prime number $p$, such
that the rank of $B$ is divisible by $p-1$, return whether $(A, B, C)$ is
$p$-admissible.

See Lemma 4.15. in [BH23](@cite) for a definition of $p$-admissible.

# Examples
A standard example is the following: let $(L, f)$ be a lattice with isometry of
prime order $p$, let $F:= L^f$ and $C:= L_f$ be respectively the invariant
and coinvariant sublattices of $(L, f)$. Then, the triple of genera
$(g(F), g(C), g(L))$ is $p$-admissible.

```jldoctest
julia> L = root_lattice(:A, 5);

julia> f = matrix(QQ, 5, 5, [1  1  1  1  1;
                             0 -1 -1 -1 -1;
                             0  1  0  0  0;
                             0  0  1  0  0;
                             0  0  0  1  0]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> F = invariant_lattice(Lf);

julia> C = coinvariant_lattice(Lf);

julia> is_admissible_triple(genus(F), genus(C), genus(Lf), 5)
true
```
"""
function is_admissible_triple(A::ZZGenus, B::ZZGenus, C::ZZGenus, p::IntegerUnion)
  zg = genus(integer_lattice(; gram=matrix(QQ, 0, 0, [])))
  AperpB = direct_sum(A, B)
  (signature_tuple(AperpB) == signature_tuple(C)) || (return false)
  if ((A == zg) && (B == C)) || ((B == zg) && (A == C))
    # C can be always glued with the empty genus to obtain C
    return true
  elseif (A == zg) || (B == zg)
    # If A or B is empty but the other is not C, then there is no gluing
    return false
  end

  @req is_divisible_by(rank(B), p-1) "p-1 must divide the rank of B"

  lA = ngens(discriminant_group(A))
  lB = ngens(discriminant_group(B))
  if all(g -> abs(det(AperpB)) != p^(2*g)*abs(det(C)), 0:min(lA, lB, divexact(rank(B), p-1)))
    return false
  end

  # A+B and C must agree locally at every primes except p
  for q in filter(!=(p), union!(ZZRingElem[2], primes(AperpB), primes(C)))
    if local_symbol(AperpB, q) != local_symbol(C, q)
      return false
    end
  end

  g = divexact(valuation(div(det(AperpB), det(C)), p), 2)
  # If the determinants agree, A+B = C and, equivalently, they agree locally at p too.
  # Otherwise, their localizations at p must be rationally equal.
  if g == 0
    return local_symbol(AperpB, p) == local_symbol(C, p)
  elseif excess(local_symbol(AperpB, p)) != excess(local_symbol(C, p))
    return false
  end

  if !is_divisible_by(scale(AperpB), scale(C)) || !is_divisible_by(p*level(C), level(AperpB))
    return false
  end

  # At this point, if C is unimodular at p, the gluing condition is equivalent to having
  # an anti-isometry between the p-part of the (quadratic) discriminant forms of A and B
  qA = discriminant_group(A)
  qB = discriminant_group(B)
  if !is_divisible_by(numerator(det(C)), p)
    return is_anti_isometric_with_anti_isometry(primary_part(qA, p)[1], primary_part(qB, p)[1])[1]
  end

  l = valuation(level(C), p)
  Ap = local_symbol(A, p)
  Bp = local_symbol(B, p)
  a_max = symbol(Ap, l+1)[2]
  b_max = symbol(Bp, l+1)[2]
  # For the gluing, rho_{l+1}(A_p) and rho_{l+1}(B_p) are anti-isometric, so they must have the
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
      Ar = ZZLocalGenus(p, symbol(Ap)[1:end-1])
    else
      Ar = genus(matrix(ZZ, 0, 0, []), p)
    end

    if length(symbol(Bp)) > 1
      Br = ZZLocalGenus(p, symbol(Bp)[1:end-1])
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
  _Cp = deepcopy(symbol(Cp))
  for s in _Cp
    s[1] += 2
  end
  Cp = ZZLocalGenus(p, _Cp)

  if !represents(local_symbol(AperpB, p), Cp)
    return false
  end
  if !represents(C, AperpB)
    return false
  end

  qC = discriminant_group(C)
  special = (p == 2) && (_is_free(qA, p, l+1)) && (_is_free(qB, p, l+1)) && (_is_even(qC, p, l))
  rA = _rho_functor(qA, p, l+1; quad=special)
  rB = _rho_functor(qB, p, l+1; quad=special)
  return is_anti_isometric_with_anti_isometry(rA, rB)[1]
end

function is_admissible_triple(A::T, B::T, C::T, p::IntegerUnion) where T <: Union{ZZLat, ZZLatWithIsom}
  return is_admissible_triple(genus(A), genus(B), genus(C), p)
end

@doc raw"""
    admissible_triples(C::ZZGenus, p::Integer; pA::Int=-1,
                                               nA::Int=-1,
                                               pB::Int=-1,
                                               nB::Int=-1,
                                               b::Int=0)
                                               -> Vector{Tuple{ZZGenus, ZZGenus}}

Given a $\mathbb Z$-genus $C$ and a prime number $p$, return all tuples of
$\mathbb Z$-genera $(A, B)$ such that $(A, B, C)$ is $p$-admissible and
$B$ is of rank divisible by $p-1$.

One can choose the positive signatures for the genera $A$ and $B$ in output
respectively by setting `pA` and `pB` to the desired values. Similarly with
the negative signatures `nA` and `nB`.

If `b` is set to `0`, we allow in output the trivial pair, i.e. when $B$ is
the genus of rank 0 lattices. Otherwise, if `b` is set to `1`, the trivial
pair is discarded.

See also [`is_admissible_triple`](@ref).

# Examples
```jldoctest
julia> L = root_lattice(:A, 5);

julia> g = genus(L)
Genus symbol for integer lattices
Signatures: (5, 0, 0)
Local symbols:
  Local genus symbol at 2: 1^-4 2^1_7
  Local genus symbol at 3: 1^-4 3^1

julia> admissible_triples(g, 5)
2-element Vector{Tuple{ZZGenus, ZZGenus}}:
 (Genus symbol: II_(5, 0) 2^-1_3 3^1, Genus symbol: II_(0, 0))
 (Genus symbol: II_(1, 0) 2^1_7 3^1 5^1, Genus symbol: II_(4, 0) 5^1)

julia> admissible_triples(g, 2)
8-element Vector{Tuple{ZZGenus, ZZGenus}}:
 (Genus symbol: II_(5, 0) 2^-1_3 3^1, Genus symbol: II_(0, 0))
 (Genus symbol: II_(4, 0) 2^2_6 3^1, Genus symbol: II_(1, 0) 2^1_1)
 (Genus symbol: II_(3, 0) 2^-3_1 3^1, Genus symbol: II_(2, 0) 2^2_2)
 (Genus symbol: II_(3, 0) 2^3_3, Genus symbol: II_(2, 0) 2^-2 3^1)
 (Genus symbol: II_(2, 0) 2^-2 3^1, Genus symbol: II_(3, 0) 2^3_3)
 (Genus symbol: II_(2, 0) 2^2_2, Genus symbol: II_(3, 0) 2^-3_1 3^1)
 (Genus symbol: II_(1, 0) 2^1_1, Genus symbol: II_(4, 0) 2^2_6 3^1)
 (Genus symbol: II_(0, 0), Genus symbol: II_(5, 0) 2^-1_3 3^1)
```
"""
function admissible_triples(G::ZZGenus, p::IntegerUnion; pA::Int=-1, nA::Int=-1, pB::Int=-1, nB::Int=-1, b::Int=0)
  @req is_prime(p) "p must be a prime number"
  @req is_integral(G) "G must be a genus of integral lattices"
  @req b >= 0 "b must be non-negative"
  L = Tuple{ZZGenus, ZZGenus}[]
  rG = rank(G)
  sG = numerator(scale(G))
  lG = numerator(level(G))
  pG, nG = signature_pair(G)
  if pA >= 0
    pA > pG && return L
    if pB >= 0
      pA + pB != pG && return L
    else
      pB = pG - pA
    end
  elseif pB >= 0
    pB > pG && return L
    pA = pG - pB
  end
  if nA >= 0
    nA > nG && return L
    if nB >= 0
      nA + nB != nG && return L
    else
      nB = nG - nA
    end
  elseif nB >= 0
    nB > nG && return L
    nA = nG - nB
  end
  dG = numerator(det(G))
  even = iseven(G)
  for ep in b:div(rG, p-1)
    rp = (p-1)*ep
    if pB >= 0
      rp >= pB || continue
      if nB >= 0
        pB + nB != rp && continue
      end
    elseif nB >= 0
      rp >= nB || continue
    end
    r1 = rG - rp
    if pA >= 0
      r1 >= pA || continue
      if nA >= 0
        pA + nA != r1 && continue
      end
    elseif nA >= 0
      r1 >= nA || continue
    end
    m = Int(min(ep, r1))
    D = _find_D(dG, m, p)
    while !is_empty(D)
      d1, dp = pop!(D)
      L1 = _find_L(pG, nG, r1, d1, sG, lG, p, even; pos=pA, neg=nA)
      Lp = _find_L(pG, nG, rp, dp, sG, lG, p, even; pos=pB, neg=nB)
      for A in L1, B in Lp
        is_admissible_triple(A, B, G, p) && push!(L, (A, B))
      end
    end
  end
  return L
end

admissible_triples(L::T, p::IntegerUnion; pA::Int=-1, nA::Int=-1, pB::Int=-1, nB::Int=-1, b::Int=0) where T <: Union{ZZLat, ZZLatWithIsom} = admissible_triples(genus(L), p; pA, nA, pB, nB, b)

###############################################################################
#
# Representatives of lattices with isometry
#
###############################################################################

# we compute ideals of E/K whose absolute norm is equal to d
# The first function covers the general case and dispatch to the second function
# for `d` integral.

function _ideals_of_norm(E::Field, d::QQFieldElem)
  OE = maximal_order(E)
  if denominator(d) == 1
    return _ideals_of_norm(E, numerator(d))
  elseif numerator(d) == 1
    return Hecke.fractional_ideal_type(OE)[inv(I) for I in _ideals_of_norm(E, denominator(d))]
  else
    return Hecke.fractional_ideal_type(OE)[I*inv(J) for I in _ideals_of_norm(E, numerator(d)), J in _ideals_of_norm(E, denominator(d))]
  end
end

function _ideals_of_norm(E::Field, d::ZZRingElem)
  OE = maximal_order(E)
  isone(d) && return Hecke.fractional_ideal_type(OE)[fractional_ideal(OE, one(E))]
  @hassert :ZZLatWithIsom 1 E isa Hecke.RelSimpleNumField
  K = base_field(E)
  OK = maximal_order(K)
  DE = different(OE)
  ids = Hecke.fractional_ideal_type(OE)[]
  primes = Vector{ideal_type(OE)}[]
  for p in prime_divisors(d)
    v = valuation(d, p)
    pd = ideal_type(OK)[P[1] for P in prime_decomposition(OK, p)]
    for i in 1:length(pd)
      if !is_coprime(DE, ideal(OE, pd[i]))
        P = prime_decomposition(OE, pd[i])[1][1]
      else
        P = ideal(OE, pd[i])
      end
      nv = valuation(norm(P), pd[i])
      push!(primes, ideal_type(OE)[P^e for e in 0:divrem(v, nv)[1]])
    end
  end
  for I in Hecke.cartesian_product_iterator(primes)
    if prod(absolute_norm.(I)) != d
      continue
    end
    I = prod(I)
    @hassert :ZZLatWithIsom 1 absolute_norm(I) == d
    push!(ids, fractional_ideal(OE, I))
  end
  return ids
end

# Given a degree 2 extension of number fields E/K, return all
# the possible signatures dictionaries of any hermitian lattice over
# E/K of rank rk, and whose trace lattice has negative signature s2.
function _possible_signatures(s2::IntegerUnion, E::Field, rk::IntegerUnion)
  lb = iseven(s2) ? 0 : 1
  K = base_field(E)
  inf = Hecke.place_type(K)[p for p in real_places(K) if length(extend(p, E)) == 1]
  r = length(real_places(K)) - length(inf)
  s = length(inf)
  signs = Dict{Hecke.place_type(K), Int}[]
  perm = AllPerms(s)
  for l in lb:2:min(s2, rk*r)
    parts = Vector{Int}[]
    l = divexact(s2-l, 2)
    for v in AllParts(l)
      if any(>(rk), v)
        continue
      end
      if length(v) > s
        continue
      end
      while length(v) != s
        push!(v, 0)
      end
      push!(parts, copy(v))
      for vv in perm
        v2 = v[vv.d]
        push!(parts, v2)
      end
    end
    unique!(parts)
    for v in parts
      push!(signs, Dict(a => b for (a, b) in zip(inf, v)))
    end
  end
  return signs
end

function _action_on_genus(G::HermGenus, s::NumFieldHom)
  E = base_field(G)
  t = inv(s)
  lgs = local_symbols(G)
  si = signatures(G)
  si_new = empty(si)
  for (p, v) in si
    si_new[Hecke.induce_image(t, p)] = v
  end

  lgs_new = empty(lgs)
  for sym in lgs
    _s = sym.data
    # sym.data is always a list of triples, but we also need the norm valuations
    # whenever sym is dyadic and ramified. For those cases, we add it to our local symbol.
    if is_dyadic(sym) && is_ramified(sym)
      nn = sym.norm_val
      @assert length(nn) == length(_s)
      s = Tuple{Int, Int, Int, Int}[(_s[i][1], _s[i][2], _s[i][3], nn[i]) for i in 1:length(nn)]
      push!(lgs_new, genus(HermLat, E, t(sym.p), s))
    else
      push!(lgs_new, genus(HermLat, E, t(sym.p), _s))
    end
  end

  @hassert :ZZLatWithIsom 1 Hecke._check_global_genus(lgs_new, si_new)
  G_new = HermGenus(E, rank(G), lgs_new, si_new)
  return G_new
end

@doc raw"""
    representatives_of_hermitian_type(Lf::ZZLatWithIsom, m::Int = 1)
                                                       -> Vector{ZZLatWithIsom}

Given a lattice with isometry $(L, f)$ of finite hermitian type (i.e. the
minimal polynomial of $f$ is irreducible cyclotomic) and a positive integer $m$,
return a set of representatives of isomorphism classes of lattices with isometry
$(M, g)$ of finite hermitian type such that the type of $(M, g^m)$ is equal to the
type of $(L, f)$.

Note that in this case, the isometries $g$'s are of order $nm$.

See also [`type(::ZZLatWithIsom)`](@ref).

# Examples
```jldoctest
julia> L = root_lattice(:A, 2);

julia> Lf = integer_lattice_with_isometry(L);

julia> reps = representatives_of_hermitian_type(Lf, 6)
1-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 6

julia> is_of_hermitian_type(reps[1])
true
```
"""
function representatives_of_hermitian_type(Lf::ZZLatWithIsom, m::Int = 1, fix_root::Bool = false)
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  @req m >= 1 "m must be a positive integer"
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian type"

  n = order_of_isometry(Lf)
  @req is_finite(n) "Isometry must be of finite order"
  reps = representatives_of_hermitian_type(genus(Lf), cyclotomic_polynomial(n*m), fix_root)
  filter!(M -> is_of_same_type(M^m, Lf), reps)
  return reps
end

@doc raw"""
    representatives_of_hermitian_type(G::ZZGenus, m::Int; first::Bool=false)
    representatives_of_hermitian_type(L::ZZLat, m::Int; first::Bool=false)
                                                       -> Vector{ZZLatWithIsom}

Given a nonempty genus of integer lattices $G$, return a list of
representatives of isomorphism classes of pairs $(M, g)$ consisting of a lattice
$M$ in $G$ and $g \in O(M)$ is an isometry of minimal polynomial $\Phi_m(X)$,
the $m$th cyclotomic polynomial.

If $m = 1,2$, this goes back to enumerating $G$ as a genus of integer lattices.

One can also provide a representative $L$ of $G$ instead.

If `first` is set to `true`, only return the first representative computed.
"""
representatives_of_hermitian_type(::Union{ZZGenus, ZZLat}, ::Int, ::Bool)

representatives_of_hermitian_type(G::ZZGenus, m::Int, fix_root::Bool = false; first::Bool=false) = representatives_of_hermitian_type(G, cyclotomic_polynomial(m), fix_root; first)

representatives_of_hermitian_type(L::ZZLat, m::Int, fix_root::Bool = false; first::Bool=false) = representatives_of_hermitian_type(genus(L), cyclotomic_polynomial(m), fix_root; first)

@doc raw"""
    representatives_of_hermitian_type(G::ZZGenus, chi::Union{ZZPolyRingElem, QQPolyRingElem}; first::Bool=false)
    representatives_of_hermitian_type(L::ZZLat, chi::Union{ZZPolyRingElem, QQPolyRingElem}; first::Bool=false)
                                                                   -> Vector{ZZLatWithIsom}

Given a nonempty genus of integer lattices $G$ and a polynomial $chi$ irreducible
over $\mathbb Q$, such that the equation order of the associated number field is
maximal, return a list of representatives of isomorphism classes of pairs $(M, g)$
consisting of a lattice $M$ in $G$ and $g \in O(M)$ is an isometry of minimal polynomial
$chi$.

One can also provide a representative $L$ of $G$ instead.

If `first` is set to `true`, only return the first representative computed.
"""
representatives_of_hermitian_type(::Union{ZZLat, ZZGenus}, ::Union{ZZPolyRingElem, QQPolyRingElem}, ::Bool)

function representatives_of_hermitian_type(G::ZZGenus, chi::Union{ZZPolyRingElem, QQPolyRingElem}, fix_root::Bool = false; first::Bool=false)
  @req is_irreducible(chi) "Polynomial must be irreducible"
  @req is_integral(G) "For now G must be a genus symbol for integral lattices"

  rk = rank(G)
  d = abs(det(G))
  s1, _, s2 = signature_tuple(G)

  reps = ZZLatWithIsom[]

  rank(G) == 0 && return ZZLatWithIsom[integer_lattice_with_isometry(integer_lattice(; gram=matrix(QQ, 0, 0, QQFieldElem[])))]

  if degree(chi) == 1
    @hassert :ZZLatWithIsom 1 iszero(chi(1)*chi(-1))
    @vprintln :ZZLatWithIsom 1 "Order smaller than 3"
    f = iszero(chi(1)) ? identity_matrix(QQ, rk) : -identity_matrix(QQ, rk)
    repre = representatives(G)
    @vprintln :ZZLatWithIsom 1 "$(length(repre)) representative(s)"
    while !is_empty(repre)
      LL = pop!(repre)
      push!(reps, integer_lattice_with_isometry(LL, f; ambient_representation=false, check=false))
    end
    return reps
  end

  !iseven(degree(chi)) && return reps

  @vprintln :ZZLatWithIsom 1 "Order bigger than 3"
  ok, rk = divides(rk, degree(chi))
  ok || return reps

  R = parent(chi)
  list_cyc = Int[k for k in euler_phi_inv(degree(chi))]
  j = findfirst(k -> chi == cyclotomic_polynomial(k, R), list_cyc)
  if !isnothing(j)
    E, b = cyclotomic_field_as_cm_extension(list_cyc[j])
  else
    Etemp, btemp = number_field(chi; cached=false)
    @req is_maximal(equation_order(Etemp)) "For infinite isometries, the equation order of the associated number field must be maximal"
    K, a = number_field(minpoly(btemp + inv(btemp)), "a"; cached=false)
    Kt, t = K[:t]
    E, b = number_field(t^2-a*t+1, "b"; cached=false)
  end

  gene = Hecke.genus_herm_type(E)[]
  DEK = different(maximal_order(E))
  DK = different(base_ring(maximal_order(E)))
  DE = DK*maximal_order(E)*DEK

  @vprintln :ZZLatWithIsom 1 "We have the differents"

  ndE = d*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)
  isempty(detE) && return reps
  @vprintln :ZZLatWithIsom 1 "All possible ideal dets: $(length(detE))"

  signatures = _possible_signatures(s2, E, rk)
  isempty(signatures) && return reps
  @vprintln :ZZLatWithIsom 1 "All possible signatures: $(length(signatures))"

  for dd in detE, sign in signatures
    append!(gene, hermitian_genera(E, rk, sign, dd; min_scale=inv(DE), max_scale=numerator(dd)*DE))
  end
  unique!(gene)
  isempty(gene) && return reps

  # In the cyclotomic case, the Galois group of the fixed field K acts on the
  # set of genera by "change of fixed primitive root of unity".
  # On the bilinear level, this action corresponds to taking certain powers
  # of the isometry arising from the trace construction (given as multiplication
  # by a primitive root of unity).
  # If `fix_root == true`, then we consider the previous set of genera up to
  # this Galois action.
  if fix_root
    @req !isnothing(j) "fix_root argument available only for cyclotomic polynomials"
    K = base_field(E) # Totally real and Galois over QQ
    GKQ, j = automorphism_group(K)
    phi = inv(isomorphism(PermGroup, GKQ)) # Want a GAPGroup to create a gset
    omega = gset(domain(phi), (G, g) -> _action_on_genus(G, j(phi(inv(g)))), gene) # We have a mathematical left action to put into a right action -> double inverse
    gene = representative.(orbits(omega))
  end

  @vprintln :ZZLatWithIsom 1 "All possible genera: $(length(gene))"
  for g in gene
    if is_integral(G) && !is_integral(DE*scale(g))
      continue
    end
    if is_even(G) && !is_integral(DK*norm(g))
      continue
    end
    @v_do :ZZLatWithIsom 3 Base.show(stdout, MIME"text/plain"(), g)
    @vprintln :ZZLatWithIsom 1 ""

    H = representative(g)
    M, fM = trace_lattice_with_isometry(H)
    genus(M) != G && continue

    MfM = integer_lattice_with_isometry(M, fM; check=false)
    @hassert :ZZLatWithIsom 1 is_of_hermitian_type(MfM)
    first && return ZZLatWithIsom[MfM]

    gr = genus_representatives(H)
    for HH in gr
      M, fM = trace_lattice_with_isometry(HH)
      push!(reps, integer_lattice_with_isometry(M, fM; check=false))
    end
  end
  return reps
end

representatives_of_hermitian_type(L::ZZLat, chi::Union{ZZPolyRingElem, QQPolyRingElem}, fix_root::Bool = false; first::Bool=false) = representatives_of_hermitian_type(genus(L), chi, fix_root; first)

@doc raw"""
    splitting_of_hermitian_prime_power(Lf::ZZLatWithIsom, p::Int) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ of hermitian type with $f$ of order $q^e$
for some prime number $q$, and given another prime number $p \neq q$, return a
set of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

Note that $e$ can be `0`.

# Examples
```jldoctest
julia> L = root_lattice(:A, 2);

julia> f = matrix(QQ, 2, 2, [0 1; -1 -1]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> is_of_hermitian_type(Lf)
true

julia> reps = splitting_of_hermitian_prime_power(Lf, 2)
2-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 6
 Integer lattice with isometry of finite order 3

julia> all(is_of_hermitian_type, reps)
true

julia> is_of_same_type(Lf, reps[1]^2)
true

julia> is_of_same_type(Lf, reps[2]^2)
true
```
"""
function splitting_of_hermitian_prime_power(Lf::ZZLatWithIsom, p::IntegerUnion;
                                                               b::Int=0,
                                                               pA::Int=-1,
                                                               nA::Int=-1,
                                                               pB::Int=-1,
                                                               nB::Int=-1,
                                                               fix_root::Bool=false)
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  @req iseven(Lf) "Lattice must be even"
  @req is_prime(p) "p must be a prime number"
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian type"

  n = order_of_isometry(Lf)
  @req is_finite(n) "Isometry must be of finite order"

  ok, e, q = is_prime_power_with_data(n)

  @req ok || e == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = ZZLatWithIsom[]
  @vprintln :ZZLatWithIsom 1 "Compute admissible triples"
  atp = admissible_triples(Lf, p; pA, nA, pB, nB, b)
  @vprintln :ZZLatWithIsom 1 "$(length(atp)) admissible triple(s)"
  while !is_empty(atp)
    A, B = pop!(atp)
    RB = representatives_of_hermitian_type(B, p*q^e, fix_root)
    is_empty(RB) && continue
    RA = representatives_of_hermitian_type(A, q^e)
    is_empty(RA) && continue
    for L1 in RA, L2 in RB
      E = admissible_equivariant_primitive_extensions(L1, L2, Lf, p; check=false)
      append!(reps, E)
    end
  end
  return reps
end

@doc raw"""
    splitting_of_prime_power(Lf::ZZLatWithIsom, p::Int, b::Int = 0)
                                                       -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ with $f$ of order $q^e$ for some
prime number $q$, a prime number $p \neq q$ and an integer $b = 0, 1$, return
a set of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

If `b == 1`, return only the lattices with isometry $(M, g)$ where $g$ has
order $pq^e$.

Note that $e$ can be `0`.

# Examples
```jldoctest
julia> L = root_lattice(:A, 2);

julia> Lf = integer_lattice_with_isometry(L);

julia> splitting_of_prime_power(Lf, 2)
4-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 2
 Integer lattice with isometry of finite order 2
 Integer lattice with isometry of finite order 2
 Integer lattice with isometry of finite order 1

julia> splitting_of_prime_power(Lf, 3, 1)
1-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 3
```
"""
function splitting_of_prime_power(Lf::ZZLatWithIsom, p::IntegerUnion, b::Int = 0;
                                                                      p_inv::Int=-1,
                                                                      n_inv::Int=-1,
                                                                      pB::Int=-1,
                                                                      nB::Int=-1,
                                                                      fix_root::Bool=false)
  if rank(Lf) == 0
    (b == 0) && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  @req iseven(Lf) "Lattice must be even"
  @req is_prime(p) "p must be a prime number"
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"

  ord = order_of_isometry(Lf)
  @req ord isa Int "Order of isometry must be finite"

  ok, e, q = is_prime_power_with_data(ord)

  @req ok || e == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = ZZLatWithIsom[]

  if e == 0
    reps = splitting_of_hermitian_prime_power(Lf, p; b, pA=p_inv, nA=n_inv, pB, nB, fix_root)
    (b == 1) && filter!(M -> order_of_isometry(M) == p, reps)
    return reps
  end

  x = gen(Hecke.Globals.Qx)
  A0 = kernel_lattice(Lf, q^e)
  B0 = kernel_lattice(Lf, x^(q^(e-1))-1)
  RA = splitting_of_hermitian_prime_power(A0, p; pB, nB, fix_root)
  is_empty(RA) && return reps
  RB = splitting_of_prime_power(B0, p; p_inv, n_inv)
  is_empty(RB) && return reps
  for L1 in RA, L2 in RB
    n1 = order_of_isometry(L1)::Int
    n2 = order_of_isometry(L2)::Int
    b == 1 && !is_divisible_by(n1, p) && !is_divisible_by(n2, p) && continue
    E = admissible_equivariant_primitive_extensions(L2, L1, Lf, q, p; check=false)
    @hassert :ZZLatWithIsom 1 b == 0 || all(LL -> order_of_isometry(LL) == p*q^e, E)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_pure_mixed_prime_power(Lf::ZZLatWithIsom, p::Int)
                                                 -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ and a prime number $p$, such that
$\prod_{i=0}^e\Phi_{p^dq^i}(f)$ is trivial for some $d > 0$ and $e \geq 0$,
return a set of representatives of the isomorphism classes of lattices with
isometry $(M, g)$ such that the type of $(M, g^p)$ is equal to the type
of $(L, f)$.

Note that $e$ can be `0`, while $d$ has to be positive.
"""
function splitting_of_pure_mixed_prime_power(Lf::ZZLatWithIsom, _p::IntegerUnion;
                                                                p_inv::Int=-1,
                                                                n_inv::Int=-1,
                                                                pB::Int=-1,
                                                                nB::Int=-1,
                                                                fix_root::Bool=false)
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  @req iseven(Lf) "Lattice must be even"

  n = order_of_isometry(Lf)

  @req is_finite(n) "Isometry must be of finite order"

  p = typeof(n)(_p)
  @req is_prime(p) "p must be a prime number"

  pd = prime_divisors(n)

  @req 1 <= length(pd) <= 2 && p in pd "Order must be divisible by p and have at most 2 prime divisors"

  if length(pd) == 2
    q = pd[1] == p ? pd[2] : pd[1]
    d = valuation(n, p)::Int  # Weird ? Valuation is type stable and returns Int but it bugs here
    e = valuation(n, q)::Int
  else
    q = one(n)
    d = valuation(n, p)::Int
    e = zero(n)
  end

  phi = minimal_polynomial(Lf)
  chi = prod(cyclotomic_polynomial(p^d*q^i, parent(phi)) for i in 0:e; init=one(phi))

  @req is_divisible_by(chi, phi) "Minimal polynomial is not of the correct form"

  reps = ZZLatWithIsom[]

  if e == 0
    return representatives_of_hermitian_type(Lf, p, fix_root)
  end

  A0 = kernel_lattice(Lf, p^d*q^e)
  if pB >= 0 && signature_tuple(A0)[1] != pB
    return reps
  elseif nB >= 0 && signature_tuple(A0)[3] != nB
    return reps
  end
  bool, r = divides(phi, cyclotomic_polynomial(p^d*q^e, parent(phi)))
  @hassert :ZZLatWithIsom 1 bool

  B0 = kernel_lattice(Lf, r)
  RA = representatives_of_hermitian_type(A0, p, fix_root)
  is_empty(RA) && return reps
  RB = splitting_of_pure_mixed_prime_power(B0, p; p_inv, n_inv)
  is_empty(RB) && return reps
  for L1 in RA, L2 in RB
    E = admissible_equivariant_primitive_extensions(L2, L1, Lf, q, p; check=false)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_mixed_prime_power(Lf::ZZLatWithIsom, p::Int, b::Int = 1)
                                          -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ and a prime number $p$ such that
$f$ has order $p^dq^e$ for some prime number $q \neq p$, return a set
of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

If `b == 1`, return only the lattices with isometry $(M, g)$ where $g$ has
order $p^{d+1}q^e$.

Note that $d$ and $e$ can be both `0`.

# Examples
```jldoctest
julia> L = root_lattice(:E, 7);

julia> f = matrix(QQ, 7, 7, [ 1  1  2  1  0  0  1;
                             -1 -2 -3 -2 -1 -1 -1;
                              0  1  2  1  1  1  1;
                              0  0 -1 -1 -1 -1 -1;
                              1  1  2  2  2  1  1;
                              0  0 -1 -1 -1  0  0;
                              0  0  0  1  0  0  0]);

julia> Lf = integer_lattice_with_isometry(L, f)
Integer lattice of rank 7 and degree 7
  with isometry of finite order 6
  given by
  [ 1    1    2    1    0    0    1]
  [-1   -2   -3   -2   -1   -1   -1]
  [ 0    1    2    1    1    1    1]
  [ 0    0   -1   -1   -1   -1   -1]
  [ 1    1    2    2    2    1    1]
  [ 0    0   -1   -1   -1    0    0]
  [ 0    0    0    1    0    0    0]

julia> reps = splitting_of_mixed_prime_power(Lf, 2)
2-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 12
 Integer lattice with isometry of finite order 12

julia> all(LL -> is_of_same_type(Lf, LL^2), reps)
true
```
"""
function splitting_of_mixed_prime_power(Lf::ZZLatWithIsom, p::IntegerUnion, b::Int = 1;
                                                                            p_inv::Int=-1,
                                                                            n_inv::Int=-1,
                                                                            pB::Int=-1,
                                                                            nB::Int=-1,
                                                                            fix_root::Bool=false)
  if rank(Lf) == 0
    b == 0 && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  @req iseven(Lf) "Lattice must be even"
  n = order_of_isometry(Lf)

  @req is_finite(n) "Isometry must be of finite order"

  pd = prime_divisors(n)

  @req length(pd) <= 2 "Order must have at most 2 prime divisors"

  if !(p in pd)
    return splitting_of_prime_power(Lf, p, b; p_inv, n_inv, pB, nB, fix_root)
  end

  d = valuation(n, p)
  _, e, q = is_prime_power_with_data(divexact(n, p^d))

  reps = ZZLatWithIsom[]

  x = gen(parent(minimal_polynomial(Lf)))
  B0 = kernel_lattice(Lf, x^(divexact(n, p)) - 1)
  A0 = kernel_lattice(Lf, prod(cyclotomic_polynomial(p^d*q^i) for i in 0:e))
  RA = splitting_of_pure_mixed_prime_power(A0, p; pB, nB, fix_root)
  isempty(RA) && return reps
  RB = splitting_of_mixed_prime_power(B0, p, 0; p_inv, n_inv)
  is_empty(RB) && return reps
  for L1 in RA, L2 in RB
    E = admissible_equivariant_primitive_extensions(L2, L1, Lf, p; check=false)
    b == 1 && filter!(LL -> order_of_isometry(LL) == p^(d+1)*q^e, E)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    enumerate_classes_of_lattices_with_isometry(L::ZZLat, order::IntegerUnion)
                                                            -> Vector{ZZLatWithIsom}
    enumerate_classes_of_lattices_with_isometry(G::ZZGenus, order::IntegerUnion)
                                                            -> Vector{ZZLocalGenus}

Given an even integer lattice $L$, return representatives of isomorphism classes
of lattices with isometry $(M ,g)$ where $M$ is in the genus of $L$, and $g$ has order
`order`. Alternatively, one can input a given genus symbol $G$ for even integer
lattices as an input - the function first computes a representative of $G$.

Note that currently we support only orders which admit at most 2 prime divisors.
"""
function enumerate_classes_of_lattices_with_isometry(L::ZZLat, order::IntegerUnion)
  @req iseven(L) "Lattice must be even"
  @req is_finite(order) && order >= 1 "order must be positive and finite"
  if order == 1
    reps = representatives_of_hermitian_type(L, 1)
    return reps
  end
  pd = prime_divisors(order)
  @req length(pd) in [1, 2] "order = $(order) must have at most two prime divisors"
  if length(pd) == 1
    v = valuation(order, pd[1])
    reps = enumerate_prime_power_isometries(L, pd[1], v)
    return reps
  end
  p, q = sort!(pd)
  vp = valuation(order, p)
  vq = valuation(order, q)
  Lq = enumerate_prime_power_isometries(L, q, vq)
  reps = ZZLatWithIsom[]
  while !is_empty(Lq)
    N = pop!(Lq)
    append!(reps, _split_prime_power(N, p, vp))
  end
  @hassert :ZZLatWithIsom 6 all(N -> order_of_isometry(N) == order, reps)
  return reps
end

enumerate_classes_of_lattices_with_isometry(G::ZZGenus, order::IntegerUnion) =
                enumerate_classes_of_lattices_with_isometry(representative(G), order)

@doc raw"""
    enumerate_prime_power_isometries(L::ZZLat, q::IntegerUnion, vq::IntegerUnion) -> Vector{ZZLatWithIsom}
    enumerate_prime_power_isometries(G::ZZGenus, q::IntegerUnion, vq::IntegerUnion) -> Vector{ZZLatWithIsom}

Given a genus ``G`` of even integer lattices, or a representative ``L`` of ``G``,
return representatives of isomorphism classes of lattices with isometry ``(M, g)``
where ``M`` is a lattice in ``G`` and ``g`` has order $q^vq$.

``q`` must be a prime number and ``vq`` must be a positive integer.
"""
function enumerate_prime_power_isometries(L::ZZLat, q::IntegerUnion, vq::IntegerUnion)
  @hassert :ZZLatWithIsom 1 is_prime(q)
  @hassert :ZZLatWithIsom 1 vq >= 1
  Lq = splitting_of_prime_power(integer_lattice_with_isometry(L), q, 1)
  vq == 1 && return Lq
  rtypes = Dict[]
  reps = ZZLatWithIsom[]
  while !is_empty(Lq)
    N = popfirst!(Lq)
    any(t -> is_of_type(N, t), rtypes) && continue
    push!(rtypes, type(N))
    v = valuation(order_of_isometry(N), q)
    @hassert :ZZLatWithIsom 1 (1 <= v < vq)
    Nq = splitting_of_mixed_prime_power(N, q)
    @hassert :ZZLatWithIsom 1 all(NN -> valuation(order_of_isometry(NN), q) == v+1, Nq)
    if v == vq -1
      append!(reps, Nq)
    else
      append!(Lq, Nq)
    end
  end
  return reps
end

enumerate_prime_power_isometries(G::ZZGenus, q::IntegerUnion, vq::IntegerUnion) = enumerate_prime_power_isometries(representative(G), q, vq)

# `N` is lattice with isometry of order q^vq for some prime number q different
# from p. Computes representatives of isomorphism classes of lattice with
# isometry in the genus of `N`, of order p^vq*q^vq and whose q-type is the same
# as `N`.
function _split_prime_power(N::ZZLatWithIsom, p::IntegerUnion, vp::IntegerUnion)
  pd = prime_divisors(order_of_isometry(N))
  @hassert :ZZLatWithIsom 1 (length(pd) == 1 && !(p in pd))
  Np = splitting_of_prime_power(N, p, 1)
  vp == 1 && return Np
  rtypes = Dict[]
  reps = ZZLatWithIsom[]
  while !is_empty(Np)
    M = popfirst!(Np)
    any(t -> is_of_type(M, t), rtypes) && continue
    push!(rtypes, type(M))
    v = valuation(order_of_isometry(M), p)
    @hassert :ZZLatWithIsom 1 (1 <= v < vp)
    Mp = splitting_of_mixed_prime_power(M, p)
    @hassert :ZZLatWithIsom 1 all(MM -> valuation(order_of_isometry(MM), p) == v+1, Mp)
    if v == vp-1
      append!(reps, Mp)
    else
      append!(Np, Mp)
    end
  end
  return reps
end

###############################################################################
#
#  Testing functions
#
###############################################################################

function _get_isometry_prime_power!(D::Dict, L::ZZLat, p::IntegerUnion, j::IntegerUnion)
  if !haskey(D, p)
    Dp = splitting_of_prime_power(integer_lattice_with_isometry(L), p, 1)
    D[p] = Dp
  end
  for i in 2:j
    if !haskey(D, p^i)
      Dp = D[p^(i-1)]
      rtypes = Dict[]
      Dpi = ZZLatWithIsom[]
      for N in Dp
        any(t -> is_of_type(N, t), rtypes) && continue
        push!(rtypes, type(N))
        Np = splitting_of_mixed_prime_power(N, p)
        filter!(NN -> valuation(order_of_isometry(NN), p) == i, Np)
        append!(Dpi, Np)
      end
      D[p^i] = Dpi
    end
  end
  return nothing
end

function _get_isometry_composite!(D::Dict, n::IntegerUnion)
  p, q = sort(prime_divisors(n))
  i, j = valuation(n, p), valuation(n, q)
  for k in 1:i
    Dq = D[p^(k-1)*q^j]
    rtypes = Dict[]
    Dn = ZZLatWithIsom[]
    for N in Dq
      any(t -> is_of_type(N, t), rtypes) && continue
      push!(rtypes, type(N))
      Np = splitting_of_mixed_prime_power(N, p)
      filter!(NN -> order_of_isometry(NN) == p^k*q^j, Np)
      append!(Dn, Np)
    end
    D[p^k*q^j] = Dn
  end
  return nothing
end

function _test_isometry_enumeration(L::ZZLat, k::Int = 2*rank(L)^2)
  n = rank(L)
  ord = filter(m -> euler_phi(m) <= n && length(prime_divisors(m)) <= 2, 2:k)
  pds = unique!(reduce(vcat, prime_divisors.(ord)))
  vals = Int[maximum(Int[valuation(x, p) for x in ord]) for p in pds]
  D = Dict{Int, Vector{ZZLatWithIsom}}()
  D[1] = ZZLatWithIsom[integer_lattice_with_isometry(L)]
  for i in 1:length(vals)
    p = pds[i]
    j = vals[i]
    _get_isometry_prime_power!(D, L, p, j)
  end
  for n in ord
    is_prime_power_with_data(n)[1] && continue
    _get_isometry_composite!(D, n)
  end
  return D
end
