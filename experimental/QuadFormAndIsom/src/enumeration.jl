###############################################################################
#
# This is an import to Oscar of the methods written following the paper [BH23]
# on "Finite subgroups of automorphisms of K3 surfaces".
#
###############################################################################

###############################################################################
#
# Admissible triples
#
###############################################################################

# We collect the pairs `(d', d/d')` for all divisors `d'` of `d`.
function _tuples_divisors(d::T) where T <: IntegerUnion
  return Tuple{T, T}[(dd, abs(divexact(d, dd))) for dd in divisors(d)]
end

# This is line 8 of Algorithm 1, they correspond to the possible
# discriminant for the genera A and B to glue to fit in C. d is
# the determinant of C, m the maximal p-valuation of the gcd of
# d1 and dp.
function _find_D(d::T, m::Int, p::Int) where T <: IntegerUnion
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
function _find_L(
    pG::Int,
    nG::Int,
    r::Int,
    d::RationalUnion,
    s::ZZRingElem,
    l::ZZRingElem,
    p::Int,
    even::Bool = true;
    pos::AbstractVector{Int}=0:pG,
    neg::AbstractVector{Int}=0:nG
  )
  def = ZZGenus[genus(integer_lattice(; gram=matrix(QQ, 0, 0, [])))]
  if r == 0 && d == 1
    return def
  end
  gen = ZZGenus[]
  Ip = intersect(max(0, r-nG):min(pG, r), pos)
  In = intersect(max(0, r-pG):min(nG, r), neg)
  for s1 in Ip
    s2 = r-s1
    !(s2 in In) && continue
    L = integer_genera((s1, s2), d; even)
    filter!(G -> is_divisible_by(numerator(scale(G)), s), L)
    filter!(G -> is_divisible_by(p*l, numerator(level(G))), L)
    append!(gen, L)
  end
  return gen
end

@doc raw"""
    is_admissible_triple(
      A::ZZGenus,
      B::ZZGenus,
      C::ZZGenus,
      p::Int
    ) -> Bool

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
function is_admissible_triple(
    A::ZZGenus,
    B::ZZGenus,
    C::ZZGenus,
    p::Int
  )
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

function is_admissible_triple(
    A::T,
    B::T,
    C::T,
    p::IntegerUnion
  ) where T <: Union{ZZLat, ZZLatWithIsom}
  return is_admissible_triple(genus(A), genus(B), genus(C), p)
end

@doc raw"""
    admissible_triples(
      G::Union{ZZGenus, ZZLat, ZZLatWithIsom},
      p::Int;
      IrA::AbstractVector{Int}=0:rank(G),
      IpA::AbstractVector{Int}=0:rank(G),
      InA::AbstractVector{Int}=0:rank(G),
      IrB::AbstractVector{Int}=0:rank(G),
      IpB::AbstractVector{Int}=0:rank(G),
      InB::AbstractVector{Int}=0:rank(G),
      b::Int=0
    ) -> Vector{Tuple{ZZGenus, ZZGenus}}

Given a $\mathbb Z$-genus $G$ and a prime number $p$, return all tuples of
$\mathbb Z$-genera $(A, B)$ such that $(A, B, G)$ is $p$-admissible and
$B$ is of rank divisible by $p-1$.

One can choose the ranges `IrA`, `IpA` and `InA` in which to take the rank,
the positive signature and negative signature of `A` respectively. Similarly
for `B` by choosing `IrB`, `IpB` and `InB`.

If `b` is set to `0`, we allow in output the trivial pair, i.e. when $B$ is
the genus of rank 0 lattices. Otherwise, if `b` is set to `1`, the trivial
pair is discarded.

Alternatively, one can choose as input a lattice `L` or a lattice with
isometry `Lf`.

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
function admissible_triples(
    G::ZZGenus,
    p::Int;
    IrA::AbstractVector{Int}=0:rank(G),
    IpA::AbstractVector{Int}=0:rank(G),
    InA::AbstractVector{Int}=0:rank(G),
    IrB::AbstractVector{Int}=0:rank(G),
    IpB::AbstractVector{Int}=0:rank(G),
    InB::AbstractVector{Int}=0:rank(G),
    b::Int=0
  )
  @req is_prime(p) "p must be a prime number"
  @req is_integral(G) "G must be a genus of integral lattices"
  @req b >= 0 "b must be non-negative"
  atp = Tuple{ZZGenus, ZZGenus}[]
  rG = rank(G)
  sG = numerator(scale(G))
  lG = numerator(level(G))
  pG, nG = signature_pair(G)
  dG = numerator(det(G))
  even = iseven(G)
  intrB = intersect((p-1)*b:(p-1):rG, IrB)::AbstractVector{Int}
  for rB in intrB
    rA = rG - rB
    !(rA in IrA) && continue
    m = Int(min(div(rB, p-1), rA))
    D = _find_D(dG, m, p)
    while !is_empty(D)
      dA, dB = pop!(D)
      L1 = _find_L(pG, nG, rA, dA, sG, lG, p, even; pos=IpA, neg=InA)
      Lp = _find_L(pG, nG, rB, dB, sG, lG, p, even; pos=IpB, neg=InB)
      for A in L1, B in Lp
        is_admissible_triple(A, B, G, p) && push!(atp, (A, B))
      end
    end
  end
  return atp
end

function admissible_triples(
    L::T,
    p::IntegerUnion;
    IrA::AbstractVector{Int}=0:rank(G),
    IpA::AbstractVector{Int}=0:rank(G),
    InA::AbstractVector{Int}=0:rank(G),
    IrB::AbstractVector{Int}=0:rank(G),
    IpB::AbstractVector{Int}=0:rank(G),
    InB::AbstractVector{Int}=0:rank(G),
    b::Int=0
  ) where T <: Union{ZZLat, ZZLatWithIsom}
  return admissible_triples(genus(L), p; IrA, IpA, InA, IrB, IpB, InB, b)
end

###############################################################################
#
# Representatives of hermitian type
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
    representatives_of_hermitian_type(
      Lf::ZZLatWithIsom,
      m::Int = 1,
      fix_root::Int = -1
    )  -> Vector{ZZLatWithIsom}

Given a lattice with isometry $(L, f)$ of finite hermitian type, i.e. $f$ is of
finite order $n$ and its minimal polynomial is irreducible cyclotomic, and a
positive integer $m$, return a complete set of representatives for the
isomorphism classes of lattices with isometry $(M, g)$ of finite hermitian type
such that the type of $(M, g^m)$ is equal to the type of $(L, f)$.

Note that in this case, the isometries $g$'s are of order $nm$.

If the value of `fix_root` is exactly $nm$, then the function only returns
one generator for every conjugacy classes of finite cyclic groups generated
by an isometry $g$ as before.

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
function representatives_of_hermitian_type(
    Lf::ZZLatWithIsom,
    m::Int = 1,
    fix_root::Int = -1;
    cond::Vector{Int}=Int[-1, -1, -1],
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false
  )
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  @req m >= 1 "m must be a positive integer"
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian type"

  n = order_of_isometry(Lf)
  @req is_finite(n) "Isometry must be of finite order"
  reps = representatives_of_hermitian_type(genus(Lf), cyclotomic_polynomial(n*m), fix_root; cond, genusDB, root_test)
  filter!(M -> is_of_same_type(M^m, Lf), reps)
  return reps
end

@doc raw"""
    representatives_of_hermitian_type(
      G::Union{ZZGenus, ZZLat},
      m::Int,
      fix_root::Int = -1;
      first::Bool=false
    )  -> Vector{ZZLatWithIsom}

Given a nonempty genus of integer lattices $G$, return a list of
representatives of isomorphism classes of pairs $(M, g)$ consisting of a lattice
$M$ in $G$ and $g \in O(M)$ is an isometry of minimal polynomial $\Phi_m(X)$,
the $m$th cyclotomic polynomial.

If $m = 1,2$, this goes back to enumerating $G$ as a genus of integer lattices.

If the value of `fix_root` is exactly $m$, then the function only returns
one generator for every conjugacy classes of finite cyclic groups generated
by an isometry of minimal polynomial $\Phi_m(X)$.

One can also provide a representative $L$ of $G$ instead.

If `first` is set to `true`, only return the first representative computed.
"""
representatives_of_hermitian_type(::Union{ZZGenus, ZZLat}, ::Int, ::Int)

representatives_of_hermitian_type(
  G::ZZGenus,
  m::Int,
  fix_root::Int = -1;
  cond::Vector{Int}=Int[-1, -1, -1],
  first::Bool=false,
  genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
  root_test::Bool=false
 ) = representatives_of_hermitian_type(G, cyclotomic_polynomial(m), fix_root; cond, first, genusDB, root_test)

representatives_of_hermitian_type(
  L::ZZLat,
  m::Int,
  fix_root::Int = -1;
  cond::Vector{Int}=Int[-1, -1, -1],
  first::Bool=false,
  genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
  root_test::Bool=false
 ) = representatives_of_hermitian_type(genus(L), cyclotomic_polynomial(m), fix_root; cond, first, genusDB, root_test)

@doc raw"""
    representatives_of_hermitian_type(
      G::Union{ZZGenus, ZZLat},
      chi::Union{ZZPolyRingElem, QQPolyRingElem},
      fix_root::Int = -1;
      first::Bool=false
    ) -> Vector{ZZLatWithIsom}

Given a nonempty genus of integer lattices $G$ and a polynomial $chi$
irreducible over $\mathbb Q$, such that the equation order of the associated
number field is maximal, return a complete list of representatives for the
isomorphism classes of pairs $(M, g)$ consisting of a lattice $M$ in $G$ and
$g \in O(M)$ is an isometry of minimal polynomial $chi$.

One can also provide a representative $L$ of $G$ instead.

The value $n$ of `fix_root` matters only when $chi$ is cyclotomic.
In the case where $chi$ is the $n$th cyclotomic polynomial, the function
only returns only one generator for every conjugacy classes of
finite cyclic groups generated by an isometry of minimal polynomial $chi$.

If `first` is set to `true`, only return the first representative computed.
"""
representatives_of_hermitian_type(::Union{ZZLat, ZZGenus}, ::Union{ZZPolyRingElem, QQPolyRingElem}, ::Int)

function representatives_of_hermitian_type(
    G::ZZGenus,
    chi::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    cond::Vector{Int}=Int[-1, -1, -1],
    first::Bool=false,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false
  )
  @req is_irreducible(chi) "Polynomial must be irreducible"
  @req is_integral(G) "For now G must be a genus symbol for integral lattices"

  reps = ZZLatWithIsom[]
  rG, pG, nG = cond
  if rG >= 0 && rank(G) != rG
    return reps
  elseif pG >= 0 && signature_pair(G)[1] != pG
    return reps
  elseif nG >= 0 && signature_pair(G)[2] != nG
    return reps
  end

  rG = rank(G)
  dG = abs(det(G))
  pG, nG = signature_pair(G)

  if is_zero(rG)
    return ZZLatWithIsom[integer_lattice_with_isometry(integer_lattice(; gram=matrix(QQ, 0, 0, QQFieldElem[])))]
  end

  if degree(chi) == 1
    @hassert :ZZLatWithIsom 1 is_zero(chi(1)*chi(-1))
    @vprintln :ZZLatWithIsom 1 "Order smaller than 3"
    f = is_zero(chi(1)) ? identity_matrix(QQ, rG) : -identity_matrix(QQ, rG)
    repre = oscar_genus_representatives(G; genusDB, root_test)
    @vprintln :ZZLatWithIsom 1 "$(length(repre)) representative(s)"
    while !is_empty(repre)
      LL = pop!(repre)
      push!(reps, integer_lattice_with_isometry(LL, f; ambient_representation=false, check=false))
    end
    return reps
  end

  !iseven(degree(chi)) && return reps

  @vprintln :ZZLatWithIsom 1 "Order bigger than 3"
  ok, rk = divides(rG, degree(chi))
  ok || return reps

  R = parent(chi)
  is_cyclo, n = is_cyclotomic_polynomial_with_data(chi)
  if is_cyclo
    E, b = cyclotomic_field_as_cm_extension(n)
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

  ndE = dG*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)
  isempty(detE) && return reps
  @vprintln :ZZLatWithIsom 1 "All possible ideal dets: $(length(detE))"

  signatures = _possible_signatures(nG, E, rk)
  isempty(signatures) && return reps
  @vprintln :ZZLatWithIsom 1 "All possible signatures: $(length(signatures))"

  for dd in detE, sign in signatures
    append!(gene, hermitian_genera(E, rk, sign, dd; min_scale=inv(DE), max_scale=numerator(dd)*DE))
  end
  isempty(gene) && return reps
  unique!(gene)

  # In the cyclotomic case, the Galois group of the fixed field K acts on the
  # set of genera by "change of fixed primitive root of unity".
  # On the bilinear level, this action corresponds to taking certain powers
  # of the isometry arising from the trace construction (given as multiplication
  # by a primitive root of unity).
  # If `fix_root == n`, then we consider the previous set of genera up to
  # this Galois action.
  if is_cyclo && fix_root == n
    K = base_field(E) # Totally real and Galois over QQ
    GKQ, j = automorphism_group(K)
    phi = inv(isomorphism(PermGroup, GKQ)) # Want a GAPGroup to create a gset
    # We have a mathematical left action to put into a right action -> double inverse
    omega = gset(domain(phi), (G, g) -> _action_on_genus(G, j(phi(inv(g)))), gene)
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

function representatives_of_hermitian_type(
    L::ZZLat,
    chi::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    cond::Vector{Int}=Int[-1, -1, -1],
    first::Bool=false,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false
  )
  return representatives_of_hermitian_type(genus(L), chi, fix_root; eiglat_cond, first, genusDB, root_test)
end

###############################################################################
#
#  Splitting subprocedures
#
###############################################################################

@doc raw"""
    splitting_of_hermitian_type(
      Lf::ZZLatWithIsom,
      p::IntegerUnion,
      b::Int = 0;
      eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ of hermitian type with $f$ of
finite order $n$, and given a prime number $p$, return a complete set of
representatives for the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

The integer `b` can be set to `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*n$.

For every $(M, g)$ in output, one may decide on the rank, positive signature
and negative signature of the eigenlattices of $(M, g)$ using the keyword
argument `eiglat_cond`. It should consist of a dictionary where each key
is a divisor of the order of $g$, and the corresponding value is a tuple
`(r, p, n)` of integers.

# Examples
```jldoctest
julia> L = root_lattice(:A, 2);

julia> f = matrix(QQ, 2, 2, [0 1; -1 -1]);

julia> Lf = integer_lattice_with_isometry(L, f);

julia> is_of_hermitian_type(Lf)
true

julia> reps = splitting_of_hermitian_type(Lf, 2)
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
function splitting_of_hermitian_type(
    Lf::ZZLatWithIsom,
    p::IntegerUnion,
    b::Int = 0;
    eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=true
  )
  if rank(Lf) == 0
    (b == 0) && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @check begin
    @req iseven(Lf) "Lattice must be even"
    @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"
    @req is_finite(n) "Isometry must be of finite order"
    @req is_prime(p) "p must be a prime number"
  end

  k = p*n
  if iszero(mod(n, p))
    return representatives_of_hermitian_type(Lf, p, fix_root; cond=get(eiglat_cond, k, Int[-1, -1, -1]), genusDB, root_test)
  end

  reps = ZZLatWithIsom[]
  phi_n = euler_phi(n)
  phi_k = euler_phi(k)
  rL = rank(Lf)
  rA, pA, nA = get(eiglat_cond, n, Int[-1, -1, -1])
  intrA::AbstractVector{Int} = rA >= 0 ? Int[rA] : 0:phi_n:rL
  rB, pB, nB = get(eiglat_cond, k, Int[-1, -1, -1])
  intrB::AbstractVector{Int} = rB >= 0 ? Int[rB] : b*phi_k:phi_k:rL
  intrA = intersect(intrA, (rL).-(intrB))
  pL, _, nL = signature_tuple(Lf)
  if pA >= 0
    IpA = Int[pA]
    IpB = Int[pL - pA]
  elseif pB >= 0
    IpA = Int[pL - pB]
    IpB = Int[pB]
  end
  if nA >= 0
    InA = Int[nA]
    InB = Int[nL - nA]
  elseif nB >= 0
    InA = Int[nL - nB]
    InB = Int[nB]
  end
  for _rA in intrA
    _rB = rL - _rA
    if pA < 0 && pB < 0
      IpA = 0:max(_rA, pL)
      IpB = (pL).-(intrA)
    end
    if nA < 0 && nB < 0
      InA = 0:max(_rA, nL)
      InB = (nL).-(intnA)
    end
    atp = admissible_triples(Lf, p; IrA=[_rA], IpA, InA, IrB=[_rB], IpB, InB, b)
    for (A, B) in atp
      As = representatives_of_hermitian_type(A, n, fix_root; genusDB, root_test)
      isempty(As) && continue
      if root_test && is_definite(As[1])
        filter!(LA -> rank(LA) == 0 || minimum(LA) != 2, As)
      end
      Bs = representatives_of_hermitian_type(B, k, fix_root; genusDB, root_test)
      isempty(Bs) && continue
      if root_test && is_definite(Bs[1])
        filter!(LB -> rank(LB) == 0 || minimum(LB) != 2, Bs)
      end
      for LA in As, LB in Bs
        Es = admissible_equivariant_primitive_extensions(LA, LB, Lf, p)
        append!(reps, Es)
      end
    end
  end
  return reps
end

@doc raw"""
    splitting_of_prime_power(
      Lf::ZZLatWithIsom,
      p::Int,
      b::Int = 0;
      eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ with $f$ of order $n = q^e$ for
some prime number $q$ and a prime number $p \neq q$, return a complete set of
representatives for the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

The integer `b` can be set to `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*n$.

For every $(M, g)$ in output, one may decide on the rank, positive signature
and negative signature of the eigenlattices of $(M, g)$ using the keyword
argument `eiglat_cond`. It should consist of a dictionary where each key
is a divisor of the order of $g$, and the corresponding value is a tuple
`(r, p, n)` of integers.

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
function splitting_of_prime_power(
    Lf::ZZLatWithIsom,
    p::Int,
    b::Int = 0;
    eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=true
  )
  if rank(Lf) == 0
    (b == 0) && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @check begin
    @req iseven(Lf) "Lattice must be even"
    @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"
    @req is_finite(n) "Isometry must be of finite order"
    @req is_prime(p) "p must be a prime number"
  end

  if isone(n)
    return splitting_of_hermitian_type(Lf, p, b; eiglat_cond, fix_root, genusDB, root_test, check=false)
  end

  ok, e, q = is_prime_power_with_data(n)

  @req ok || e == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = ZZLatWithIsom[]

  x = gen(Hecke.Globals.Qx)
  A0 = kernel_lattice(Lf, x^(q^(e-1))-1)
  B0 = kernel_lattice(Lf, q^e)
  RB = splitting_of_hermitian_type(B0, p; eiglat_cond, fix_root, genusDB, root_test, check=false)
  is_empty(RB) && return reps
  RA = splitting_of_prime_power(A0, p; eiglat_cond, genusDB, root_test, check=false)
  is_empty(RA) && return reps
  for L1 in RA, L2 in RB
    n1 = order_of_isometry(L1)::Int
    n2 = order_of_isometry(L2)::Int
    b == 1 && !is_divisible_by(n1, p) && !is_divisible_by(n2, p) && continue
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, q, p; check=false)
    @hassert :ZZLatWithIsom 1 b == 0 || all(LL -> order_of_isometry(LL) == p*q^e, E)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_pure_mixed_prime_power(
      Lf::ZZLatWithIsom,
      p::Int;
      eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ with $f$ of order $n = p^d*q^e$
where $d > 0$ is positive, $e \geq 0$ is nonnegative and $q \neq p$ is a prime
number distinct from $p$, and such that $\prod_{i=0}^e\Phi_{p^dq^i}(f)$ is
trivial, return a complete set of representatives for the isomorphism classes
of lattices with isometry $(M, g)$ such that the type of $(M, g^p)$ is equal
to the type of $(L, f)$.

For every $(M, g)$ in output, one may decide on the rank, positive signature
and negative signature of the eigenlattices of $(M, g)$ using the keyword
argument `eiglat_cond`. It should consist of a dictionary where each key
is a divisor of the order of $g$, and the corresponding value is a tuple
`(r, p, n)` of integers.

Note that $e$ can be `0`, while $d$ has to be positive.
"""
function splitting_of_pure_mixed_prime_power(
    Lf::ZZLatWithIsom,
    p::Int;
    eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=true
  )
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  n = order_of_isometry(Lf)
  @check begin
    @req iseven(Lf) "Lattice must be even"
    @req is_finite(n) "Isometry must be of finite order"
    @req is_prime(_p) "p must be a prime number"
  end

  pd = prime_divisors(n)

  @req length(pd) <= 2 && p in pd "Order must be divisible by p and have at most 2 prime divisors"

  if length(pd) == 1
    return representatives_of_hermitian_type(Lf, p, fix_root; cond=get(eiglat_cond, p*n, Int[-1, -1, -1]), genusDB, root_test)
  end

  q = pd[1] == p ? pd[2] : pd[1]
  d = valuation(n, p)::Int  # Weird ? Valuation is type stable and returns Int but it bugs here
  e = valuation(n, q)::Int

  phi = minimal_polynomial(Lf)
  chi = prod(cyclotomic_polynomial(p^d*q^i, parent(phi)) for i in 0:e; init=one(phi))

  @req is_divisible_by(chi, phi) "Minimal polynomial is not of the correct form"

  reps = ZZLatWithIsom[]

  k = p^d*q^e
  bool, r = divides(phi, cyclotomic_polynomial(k, parent(phi)))
  @hassert :ZZLatWithIsom 1 bool

  A0 = kernel_lattice(Lf, r)
  B0 = kernel_lattice(Lf, k)
  RB = representatives_of_hermitian_type(B0, p, fix_root; cond=get(eiglat_cond, p*k, Int[-1, -1, -1]), genusDB, root_test)
  is_empty(RB) && return reps
  RA = splitting_of_pure_mixed_prime_power(A0, p; eiglat_cond, genusDB, root_test, check=false)
  is_empty(RA) && return reps
  for L1 in RA, L2 in RB
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, q, p; check=false)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_mixed_prime_power(
      Lf::ZZLatWithIsom,
      p::IntegerUnion,
      b::Int = 1;
      eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ and a prime number $p$ such that
$f$ has order $n=p^dq^e$ for some prime number $q \neq p$, return a set
of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

The integer `b` can be set to `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*n$.

For every $(M, g)$ in output, one may decide on the rank, positive signature
and negative signature of the eigenlattices of $(M, g)$ using the keyword
argument `eiglat_cond`. It should consist of a dictionary where each key
is a divisor of the order of $g$, and the corresponding value is a tuple
`(r, p, n)` of integers.

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
function splitting_of_mixed_prime_power(
    Lf::ZZLatWithIsom,
    p::Int,
    b::Int = 1;
    eiglat_cond::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}(),
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=true
  )
  if rank(Lf) == 0
    b == 0 && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @check begin
    @req iseven(Lf) "Lattice must be even"
    @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"
    @req is_finite(n) "Isometry must be of finite order"
    @req is_prime(p) "p must be a prime number"
  end

  pd = prime_divisors(n)

  @req length(pd) <= 2 "Order must have at most 2 prime divisors"

  if !(p in pd)
    return splitting_of_prime_power(Lf, p, b; eiglat_cond, fix_root, genusDB, root_test, check=false)
  end

  d = valuation(n, p)
  _, e, q = is_prime_power_with_data(divexact(n, p^d))

  reps = ZZLatWithIsom[]

  x = gen(parent(minimal_polynomial(Lf)))
  A0 = kernel_lattice(Lf, x^(divexact(n, p)) - 1)
  B0 = kernel_lattice(Lf, prod(cyclotomic_polynomial(p^d*q^i) for i in 0:e))
  RB = splitting_of_pure_mixed_prime_power(B0, p; eiglat_cond, fix_root, genusDB, root_test, check=false)
  isempty(RB) && return reps
  RA = splitting_of_mixed_prime_power(A0, p, 0; eiglat_cond, fix_root, genusDB, root_test, check=false)
  is_empty(RA) && return reps
  for L1 in RA, L2 in RB
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, p; check=false)
    b == 1 && filter!(LL -> order_of_isometry(LL) == p*n, E)
    append!(reps, E)
  end
  return reps
end

###############################################################################
#
#  Generic functions
#
###############################################################################

@doc raw"""
    splitting(
      Lf::ZZLatWithIsom,
      p::Int,
      b::Int = 0;
      char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
      min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
      rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[]
  ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ where $f$ is of finite order $n$,
and given a prime number $p$ return a complete set of representatives for the
isomorphism classes of lattices with isometry $(M, g)$ such that the type of
$(M, g^p)$ is equal to the type of $(L, f)$.

The integer `b` can be set to `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*n$.

For every $(M, g)$ in output, one may decide on the minimal polynomial (resp.
the characteristic polynomial) of $g$ by setting the keyword argument `min_poly`
(resp. `char_poly`) to the desired value. Moreover, one can also decide on the
rank, positive signature and negative signature of the eigenlattices of $(M, g)$
using the keyword arguments `rks`, `pos_sigs` and `neg_sigs` respectively. Each
list should consist of tuples $(m, i)$ where $m$ is a divisor of $p*n$ and $i$
is the value to be assigned for the given property (rank, positive/negative
signature) to the corresponding $\Phi_m$-eigenlattice.
"""
function splitting(
    Lf::ZZLatWithIsom,
    p::Int,
    b::Int = 0;
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    eiglat_cond::Dict{Int64, Vector{Int64}}=Dict{Int64, Vector{Int64}}(),
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=false
  )
  if rank(Lf) == 0
    b == 0 && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @check begin
    @req iseven(Lf) "Lattice must be even"
    @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"
    @req is_finite(n) "Isometry must be of finite order"
    @req is_prime(p) "p must be a prime number"
  end

  if isempty(eiglat_cond)
    eiglat_cond = _conditions_from_input(p*n, char_poly, min_poly, rks, pos_sigs, neg_sigs)
  end
  pds = prime_divisors(n)
  if (length(pds) <= 1) || (length(pds) == 2 && p in pds)
    return splitting_of_mixed_prime_power(Lf, p, b; eiglat_cond, fix_root, genusDB, root_test, check=false)
  end

  reps = ZZLatWithIsom[]
  ds = sort!(divisors(n))
  popfirst!(ds)

  I = invariant_lattice(Lf)
  Ns = splitting_of_hermitian_type(I, p; eiglat_cond, fix_root, genusDB, root_test, check=false)
  is_empty(Ns) && return reps

  x = gen(Hecke.Globals.Zx)
  chi = x - 1
  # TODO: implement a smart gluing procedure where we try to find the best
  # order for the gluings, instead of doing as currently
  for k in ds
    chi *= cyclotomic(k, x)
    M = kernel_lattice(Lf, k)
    Ms = splitting_of_hermitian_type(M, p; eiglat_cond, fix_root, genusDB, root_test, check=false)
    is_empty(Ms) && return reps
    Lq = kernel_lattice(Lf, chi)
    l = length(Ns)
    for i in 1:l
      N = popfirst!(Ns)
      for M in Ms
        ok, _Es = equivariant_primitive_extensions(N, M; q=first(discriminant_group(Lq)))
        !ok && continue
        Es = first.(_Es)
        filter!(T -> is_of_type(T^p, type(Lq)), Es)
        append!(Ns, Es)
      end
    end
    isempty(Ns) && return reps
  end
  return reps
end

@doc raw"""
    enumerate_classes_of_lattices_with_isometry(
      G::Union{ZZGenus, ZZLat},
      n::Int;
      char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
      min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
      rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[]
    ) -> Vector{ZZLatWithIsom}

Given an even integer lattice $L$, return a complete set of representatives
for the isomorphism classes of lattices with isometry $(M ,g)$ where $M$ is in
the genus of $L$, and $g$ is an isometry of $M$ of finite order $n$.
Alternatively, one can input a given genus symbol $G$ for even integer lattices
as an input --- the function first computes a representative of $G$.

For every $(M, g)$ in output, one may decide on the minimal polynomial (resp.
the characteristic polynomial) of $g$ by setting the keyword argument `min_poly`
(resp. `char_poly`) to the desired value. Moreover, one can also decide on the
rank, positive signature and negative signature of the eigenlattices of
$(M, g)$ using the keyword arguments `rks`, `pos_sigs` and `neg_sigs`
respectively. Each list should consist of tuples $(m, i)$ where $m$ is a
divisor of $p*n$ and $i$ is the value to be assigned for the given property
(rank, positive/negative signature) to the corresponding $\Phi_m$-eigenlattice.
"""
enumerate_classes_of_lattices_with_isometry(::Union{ZZGenus, ZZLat}, ::Int)

function enumerate_classes_of_lattices_with_isometry(
    L::ZZLat,
    n::Int;
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    fix_root::Int=-1
  )
  @req iseven(L) "Lattice must be even"
  @req is_finite(n) && n >= 1 "order must be positive and finite"

  eiglat_cond = _conditions_from_input(n, char_poly, min_poly, rks, pos_sigs, neg_sigs)
  @vprintln :ZZLatWithIsom 1 "Conditions computed"
  if n == 1
    reps = representatives_of_hermitian_type(L, 1; cond=get(eiglat_cond, 1, Int[-1, -1, -1]))
    return reps
  end

  Lq = ZZLatWithIsom[integer_lattice_with_isometry(L)]
  o = Int(1)
  pds = reverse!(sort!(prime_divisors(n)))
  for p in pds
    v = valuation(n, p)
    o *= p^v
    Lq = _split_prime_power!(Lq, p, v; eiglat_cond=_conditions_after_power(eiglat_cond, div(n, o)), fix_root=gcd(o, fix_root))
  end
  @hassert :ZZLatWithIsom 6 all(N -> order_of_isometry(N) == n, Lq)
  return Lq
end

function enumerate_classes_of_lattices_with_isometry(
    G::ZZGenus,
    n::Int;
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    fix_root::Int=-1,
  )
  return enumerate_classes_of_lattices_with_isometry(representative(G), n; char_poly, min_poly, rks, pos_sigs, neg_sigs, fix_root)
end

@doc raw"""
    _split_prime_power!(
      Np::Vector{ZZLatWithIsom},
      p::Int,
      v::Int;
      eiglat_cond::Dict{Int, Vector{Int}}=Dict{Int, Vector{Int}}()
    ) -> Vector{ZZLatWithIsom}

Given a list `Np` of even lattices equipped with an isometry of finite order,
return a complete set of representatives for the isomorphism classes of
lattices with isometry $(M, g)$ such that the type of $(M, g^{p^v})$ is the
same as the type of some element $(N, h)$ of `Np` and the order of $g$ is
exactly $p^v*n$ where $n$ is the order of $h$.

Not that the function empties the list `Np` in entry.

For every $(M, g)$ in output, one may decide on the rank, positive signature
and negative signature of the eigenlattices of $(M, g)$ using the keyword
argument `eiglat_cond`. It should consist of a dictionary where each key
is a divisor of the order of $g$, and the corresponding value is a tuple
`(r, p, n)` of integers.
"""
function _split_prime_power!(
    Np::Vector{ZZLatWithIsom},
    p::Int,
    v::Int;
    eiglat_cond::Dict{Int, Vector{Int}}=Dict{Int, Vector{Int}}(),
    fix_root::Int=-1
  )
  @hassert :ZZLatWithIsom 1 is_prime(p)
  @hassert :ZZLatWithIsom 1 v >= 1
  rtypes = Dict{Int, Vector{Dict}}()
  reps = ZZLatWithIsom[]
  while !is_empty(Np)
    M = popfirst!(Np)
    k = order_of_isometry(M)
    is_finite(k)
    if haskey(rtypes, k)
      any(t -> is_of_type(M, t), rtypes[k]) && continue
      push!(rtypes[k], type(M))
    else
      rtypes[k] = Dict[type(M)]
    end
    vp = valuation(k, p)
    @hassert :ZZLatWithIsom 1 (0 <= vp < v)
    q = p^(v-vp-1)
    Mp = splitting(M, p, 1; eiglat_cond=_conditions_after_power(eiglat_cond, q), check=false, fix_root=divexact(fix_root, gcd(fix_root, q)))
    @hassert :ZZLatWithIsom 1 all(MM -> valuation(order_of_isometry(MM), p) == vp+1, Mp)
    if vp == v-1
      append!(reps, Mp)
    else
      append!(Np, Mp)
    end
  end
  return reps
end

###############################################################################
#
#  Helper for managing conditions on isometries to enumerate
#
###############################################################################

function _from_cyclotomic_polynomial_to_dict(
    chi::Union{ZZPolyRingElem, QQPolyRingElem}
  )
  fac = factor(chi)
  V = Dict{Int, Int}()
  for (f, e) in fac
    ok, n = is_cyclotomic_polynomial_with_data(f)
    @req ok "Not a product of cyclotomic polynomials"
    V[n] = e
  end
  return V
end

function _conditions_from_input(
    m::Int,
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing},
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing},
    rks::Vector{Tuple{Int, Int}},
    pos_sigs::Vector{Tuple{Int, Int}},
    neg_sigs::Vector{Tuple{Int, Int}}
  )
  eiglat_cond = Dict{Int, Vector{Int}}()
  divs = divisors(m)
  if !isnothing(char_poly)
    V = _from_cyclotomic_polynomial_to_dict(char_poly)
    for (n, e) in V
      rn = euler_phi(n)*e
      jp = findfirst(a -> a[1] == n, pos_sigs)
      pn = isnothing(jp) ? -1 : pos_sigs[jp][2]
      jn = findfirst(a -> a[1] == n, neg_sigs)
      nn = isnothing(jp) ? -1 : neg_sigs[jp][2]
      eiglat_cond[n] = Int[rn, pn, nn]
    end
    for k in setdiff(divs, keys(eiglat_cond))
      eiglat_cond[k] = Int[0, 0, 0]
    end
  elseif !isnothing(min_poly)
    V = _from_cyclotomic_polynomial_to_dict(min_poly)
    for (n, e) in V
      jr = findfirst(a -> a[1] == n, rks)
      jp = findfirst(a -> a[1] == n, pos_sigs)
      pn = isnothing(jp) ? -1 : pos_sigs[jp][2]
      jn = findfirst(a -> a[1] == n, neg_sigs)
      nn = isnothing(jn) ? -1 : neg_sigs[jn][2]
      if !isnothing(jp) && !isnothing(jn)
        rn = pn + nn
      elseif !isnothing(jr)
        rn = rks[jr][2]
      else
        rn = -1
      end
      eiglat_cond[n] = Int[rn, pn, nn]
    end
    for k in setdiff(divs, keys(eiglat_cond))
      eiglat_cond[k] = Int[0, 0, 0]
    end
  else
    for k in divs
      jr = findfirst(a -> a[1] == n, rks)
      rn = isnothing(jr) ? -1 : rks[jr][2]
      jp = findfirst(a -> a[1] == n, pos_sigs)
      pn = isnothing(jp) ? -1 : pos_sigs[jp][2]
      jn = findfirst(a -> a[1] == n, neg_sigs)
      nn = isnothing(jn) ? -1 : neg_sigs[jn][2]
      eiglat_cond[k] = Int[rn, pn, nn]
    end
  end
  return eiglat_cond
end

function _conditions_after_power(
    eiglat_cond::Dict{Int, Vector{Int}},
    p::Int
  )
  isempty(eiglat_cond) && return eiglat_cond
  V = empty(eiglat_cond)
  for (n, t) in eiglat_cond
    m = div(n, gcd(p, n))
    if haskey(V, m)
      for i in 1:3
        if V[m][i] > -1 && t[i] > -1
          V[m][i] += t[i]
        else
          V[m][i] = -1
        end
      end
    else
      V[m] = max.([-1, -1, -1], t)
    end
  end
  return V
end

###############################################################################
#
#  Enhanced genus enumeration
#
###############################################################################

function oscar_genus_representatives(
  G::ZZGenus;
  genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
  root_test::Bool=false
)
  if !is_definite(G) || rank(G) <= 2
    return Hecke.representatives(G)
  end
  if !isnothing(genusDB)
    haskey(genusDB, G) && return genusDB[G]
  end
  r = rank(G)
  if root_test
    bn = Float64[0.5, 0.28868, 0.1847, 0.13127, 0.09987, 0.08112, 0.06981, 0.06326,
	       0.06007, 0.05953, 0.06136, 0.06559, 0.07253, 0.08278, 0.09735, 0.11774,
	       0.14624, 0.18629, 0.24308, 0.32454, 0.44289, 0.61722, 0.87767, 1.27241]
    if r <= 24 && abs(det(G)) < inv(bn[Int(r)])^2
      return ZZLat[]
    end
  end
  mm, l = enumerate_definite_genus(G; stop_after=1000)
  if !iszero(mm)
    inv_lat = Hecke.default_invariant_function(l[1])
    inv_dict = Dict{typeof(inv_lat), Vector{ZZLat}}(inv_lat => ZZLat[l[1]])
    for N in edg[2:end]
      inv_lat = Hecke.default_invariant_function(N)
      if haskey(inv_dict, inv_lat)
        push!(inv_dict[inv_lat], N)
      else
        inv_dict[inv_lat] = ZZLat[N]
      end
    end
    q = next_prime(last(Hecke.primes_up_to(r+1)))
    Lf = integer_lattice_with_isometry(l[1])
    pos = is_positive_definite(Lf)
    # Looking for certain lattices with isometry
    while !iszero(mm)
      d = denominator(mm)
      if isone(d)
        p = last(Hecke.primes_up_to(q-1))
      else
        p = maximum(prime_divisors(d))
      end
      q = p
      if p == 2
        interv = div(r, 2):-1:1
      else
        interv = reverse(p-1:p-1:r)
      end
      for k in interv
        if pos
          Ns = splitting_of_prime_power(Lf, Int(p), 1; eiglat_cond=Dict(1=>[r-k, r-k, 0], p=>[k, k, 0]), genusDB, root_test=false, check=false)
        else
          Ns = splitting_of_prime_power(Lf, Int(p), 1; eiglat_cond=Dict(1=>[r-k, 0, r-k], p=>[k, 0, k]), genusDB, root_test=false, check=false)
        end
        for Nf in Ns
          N = lll(lattice(Nf))
          invN = Hecke.default_invariant_function(N)
          if !haskey(inv_dict, invN)
            inv_dict[invN] = ZZLat[N]
	          push!(edg, N)
	          s = isometry_group_order(N)
	          sub!(mm, mm, 1//s)
          elseif all(M -> !is_isometric(N, M), inv_dict[invN])
            push!(inv_dict[invN], N)
            push!(edg, N)
	          s = isometry_group_order(N)
	          sub!(mm, mm, 1//s)
          end
          is_zero(mm) && break
        end
        is_zero(mm) && break
      end
    end
  end
  if !isnothing(genusDB)
    gesnuDB[G] = l
  end
  return l
end

###############################################################################
#
#  Functions for test coverage
#
###############################################################################

function _get_isometry!(
    D::Dict,
    n::Int
  )
  p = last(sort!(prime_divisors(n)))
  m = divexact(n, p)
  Dm = D[m]
  rtypes = Dict{Int, Vector{Dict}}()
  Dn = ZZLatWithIsom[]
  for N in Dm
    if haskey(rtypes, m)
      any(t -> is_of_type(N, t), rtypes[m]) && continue
      push!(rtypes[m], type(N))
    else
      rtypes[m] = Dict[type(N)]
    end
    Np = splitting(N, p, 1; check=false)
    append!(Dn, Np)
  end
  D[n] = Dn
  return nothing
end

function _test_isometry_enumeration(
    L::ZZLat,
    k::Int = 2*rank(L)^2
  )
  r = rank(L)
  ord = filter(m -> euler_phi(m) <= r, 2:k)
  pds = sort!(ord)
  popfirst!(ord)
  D = Dict{Int, Vector{ZZLatWithIsom}}()
  D[1] = ZZLatWithIsom[integer_lattice_with_isometry(L)]
  for n in ord
    _get_isometry!(D, n)
  end
  return D
end
