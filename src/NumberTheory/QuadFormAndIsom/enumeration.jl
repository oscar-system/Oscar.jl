##############################################################################
#
# This is an import to Oscar of the methods written following the paper [BH23]
# on "Finite subgroups of automorphisms of K3 surfaces".
#
# In this file, Algorithms 1--7 of [BH23] are implemented, together with extra
# methods generalising the previously mentioned algorithms. The code currently
# covers the enumeration of conjugacy classes of finite order isometries for
# even lattices. The most generic functions which are part of the UI are
# `splitting` and `enumerate_classes_of_lattices_with_isometry`.
#
###############################################################################

###############################################################################
#
# Admissible triples
#
###############################################################################

# Return all the pairs `(d', d/d')` where `d'` runs over the divisors of `d`.
function _tuples_divisors(d::T) where T <: IntegerUnion
  return Tuple{T, T}[(dd, abs(divexact(d, dd))) for dd in divisors(d)]
end

# This is line 8 of Algorithm 1 of [BH23]. Return all the pairs `(d1, dp)` of
# potential discriminants for the genera `A` and `B` in the admissible triples.
# Here `d` is the determinant of the third genus `C`, `p` is a prime number and
# `m` the maximal `p`-valuation of the gcd of `d1` and `dp`.
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

# This is line 10 of Algorithm 1 of [BH23]. We need the condition on the parity
# of `C` since subgenera of an even genus are even too.
# `r` is the rank of the subgenus, `d` its determinant, and `s` and `l` are
# respectively the scale and the level of `C`.
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
    neg::AbstractVector{Int}=0:nG,
  )
  if r == 0 && d == 1 # In this case, we return the trivial genus as default
    (!(0 in pos) || !(0 in neg)) && return ZZGenus[]
    return ZZGenus[genus(integer_lattice(; gram=matrix(QQ, 0, 0, [])))]
  end
  gen = ZZGenus[]
  Ip = intersect(max(0, r-nG):min(pG, r), pos) # Possible positive signature
  In = intersect(max(0, r-pG):min(nG, r), neg) # Possible negative signature
  for s1 in Ip
    s2 = r-s1
    !(s2 in In) && continue
    # We enumerate all potential genera with the given invariants
    append!(gen, integer_genera((s1, s2), d; min_scale=s, max_scale=p*l, even))
  end
  return gen
end

function _length_discriminant_group(g::ZZLocalGenus)
  return sum(i[2] for i in symbol(g) if i[1]>0;init = 0)
end 

function _length_discriminant_group(g::ZZGenus)
  return maximum(_length_discriminant_group(i) for i in local_symbols(g); init=0)
end 

function _is_anti_isometric_odd(a::Vector{Int}, b::Vector{Int}, is_minus_one_a_square::Bool)
  a[1]==b[1] || return false
  a[2]==b[2] || return false
  if iszero(mod(a[2], 2)) || is_minus_one_a_square
    return a[3]==b[3]
  end
  return a[3]!=b[3]  
end

function _is_anti_isometric_bilinear_2(a::Vector{Int}, b::Vector{Int})
  a[1]==b[1] || return false
  a[2]==b[2] || return false
  a[4]==b[4] || return false
end

function _is_anti_isometric_quadratic_2(a::Vector{Int}, b::Vector{Int})
  a[1]==b[1] || return false
  a[2]==b[2] || return false
  a[4] == b[4] || return false
  dets = kronecker_symbol(a[3],2)==kronecker_symbol(b[3],2)
  if a[4] == 0
    return dets
  end
  exc_a = a[5]+2*(1-kronecker_symbol(a[3],2))
  exc_b = -b[5]+2*(1-kronecker_symbol(b[3],2))
  return iszero(mod(exc_a-exc_b, 8) )
end

function _is_anti_isometric_bilinear(a::ZZLocalGenus, b::ZZLocalGenus, l::Int)
  al = symbol(a, l)
  bl = symbol(b, l)
  p = prime(a)
  if p == 2
    return _is_anti_isometric_bilinear_2(al, bl)
  end
  is_minus_one_a_square = isone(kronecker_symbol(ZZ(-1), p))
  return _is_anti_isometric_odd(al, bl, is_minus_one_a_square)
end

function _is_anti_isometric_quadratic(a::ZZLocalGenus, b::ZZLocalGenus, l::Int)
  al = symbol(a, l)
  bl = symbol(b, l)
  p = prime(a)
  if p == 2
    return _is_anti_isometric_quadratic_2(al, bl)
  end
  is_minus_one_a_square = isone(kronecker_symbol(ZZ(-1), p))
  return _is_anti_isometric_odd(al, bl, is_minus_one_a_square)
end 

@doc raw"""
    is_admissible_triple(
      A::ZZGenus,
      B::ZZGenus,
      C::ZZGenus,
      p::Int,
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

julia> multiplicative_order(f)
5

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
    p::Int,
  )
  AperpB = direct_sum(A, B)

  # Signatures conditions
  # A+B and C must have the same signatures
  (signature_tuple(AperpB) == signature_tuple(C)) || (return false)
  if ((rank(A) == 0) && (B == C)) || ((rank(B) == 0) && (A == C))
    # C can be always glued with the empty genus to obtain C
    return true
  elseif (rank(A) == 0) || (rank(B) == 0)
    # If A or B is empty but the other is not C, then there is no gluing
    return false
  end

  # Rank condition
  # The rank of B should be divisible by p-1
  !is_divisible_by(rank(B), p-1) && return false

  # Gluing condition
  # The determinant of C should divide the one of A+B, and the quotient should
  # be an even power of the prime p
  _q, _r = divrem(numerator(det(AperpB)), numerator(det(C)))
  !iszero(_r) && return false
  if isone(_q)
    g = 0
  else
    ok, _g, _p = is_prime_power_with_data(_q)
    !ok && return false
    _p != p && return false
    g, r = divrem(_g, 2)
    !iszero(r) && return false
  end

  # Condition (2) of Definition 4.13
  # Half the p-valuation of the quotient of the determinants should be
  # smaller than lA, lB and rank(B)/(p-1)
  lA = _length_discriminant_group(A)
  lB = _length_discriminant_group(B)
  if g > min(lA, lB, divexact(rank(B), p-1))
    return false
  end

  # Condition (1) of Definition 4.13
  # A+B and C must agree locally at all primes except p
  for q in filter(!=(p), union!(ZZRingElem[2], primes(AperpB), primes(C)))
    if local_symbol(AperpB, q) != local_symbol(C, q)
      return false
    end
  end


  # Gluing condition at p
  # If the determinants agree, A+B = C and, equivalently, they agree locally
  # at p.
  # Otherwise, they must be rationally equivalent locally at p.
  # Because their determinants differ by a square, Theorem 3 in Chapter 5,
  # section 5.1 (page 372) of Conway--Sloane tells that their localizations at
  # p are rationally equal if and only if the p-excess agree
  if g == 0
    return local_symbol(AperpB, p) == local_symbol(C, p)
  elseif excess(local_symbol(AperpB, p)) != excess(local_symbol(C, p))
    return false
  end

  # Condition (3) of Definition 4.13
  # The scale ideal of A+B should be contained in the scale ideal of C and
  # p times the level of C should be contained in the level of A+B
  # (here the level is the multiplicative inverse, in QQ, of the scale of the
  # dual)
  if !is_divisible_by(numerator(scale(AperpB)), numerator(scale(C)))
    return false
  elseif !is_divisible_by(p*numerator(level(C)), numerator(level(AperpB)))
    return false
  end


  # At this point, if C is unimodular at p, the gluing condition is equivalent to having
  # an anti-isometry between the p-part of the (quadratic) discriminant forms of A and B
  l = valuation(level(C), p)
  Ap = local_symbol(A, p)
  Bp = local_symbol(B, p)

  if iszero(valuation(det(C), p))
    # The lattice C_p is unimodular, so the level of A_p and B_p is at most p
    fl1 = _is_anti_isometric_quadratic(Ap, Bp, 1)
    return return fl1 
  end

  a_max = symbol(Ap, l+1)[2]
  b_max = symbol(Bp, l+1)[2]

  # Condition (4) of Definition 4.13
  # For the gluing, the rho functors rho_{l+1}(A_p) and rho_{l+1}(B_p) are
  # anti-isometric, so they must have the same order
  if a_max != b_max
    return false
  end

  # Condition (2) of Proposition 4.10
  # Since p^l*A^\vee/A is in the glue, its order is less than the order of the glue
  if g < a_max
    return false
  end

  # Ranks of the p-modular Jordan constituents
  a1 = symbol(Ap, 1)[2]
  b1 = symbol(Bp, 1)[2]
  # Sum of the ranks of the Jordan constituents of scale >= p^2
  a2 = rank(Ap) - a1 - symbol(Ap, 0)[2]
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
  Cp1 = ZZLocalGenus(p, _Cp)

  # Condition (5) of Definition 4.13
  if !represents(local_symbol(AperpB, p), Cp1)
    return false
  elseif !represents(C, AperpB)
    return false
  end

  # Condition (4) of Definition 4.13 
  _is_anti_isometric_bilinear(Ap, Bp, l+1) || return false 
    p != 2 && return true
  is_freeA = 0 == symbol(Ap, l)[4]
  is_freeB = 0 == symbol(Bp, l)[4]
  if is_freeA && is_freeB
    # condition 4'
    Cp_l_is_even = symbol(Cp, l)[4] == 0
    glues_to_even = _is_anti_isometric_quadratic(Ap, Bp, l+1)
    return Cp_l_is_even == glues_to_even
  end
  return true
end

function is_admissible_triple(
    A::T,
    B::T,
    C::T,
    p::Int,
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
      b::Int=0,
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

julia> length(admissible_triples(g, 5))
2

julia> length(admissible_triples(g, 2))
8
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
    b::Int=0,
  )
  @req is_prime(p) "p must be a prime number"
  @req is_integral(G) "G must be a genus of integral lattices"
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"

  atp = Tuple{ZZGenus, ZZGenus}[]
  # Collect some invariants of G
  rG = rank(G)
  sG = numerator(scale(G))
  lG = numerator(level(G)) # 1//lG is the scale of the dual genus
  pG, nG = signature_pair(G)
  dG = numerator(det(G))
  even = iseven(G) # If G is even, subgenera must also be even
  # Here we do a bit differently as in Algorithm 1 of [BH23]:
  # we directly iterate over the potential rank of B, which should be
  # a multiple of p-1
  intrB = intersect((p-1)*b:(p-1):rG, IrB)::AbstractVector{Int}
  for rB in intrB
    rA = rG - rB
    !(rA in IrA) && continue
    m = Int(min(div(rB, p-1), rA))
    D = _find_D(dG, m, p) # Set of tuples of potential determinants
                          # (here sign does not matter)
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
    p::Int;
    kwargs...,
  ) where T <: Union{ZZLat, ZZLatWithIsom}
  return admissible_triples(genus(L), p; kwargs...)
end

###############################################################################
#
#  Representatives of hermitian type
#
###############################################################################

@doc raw"""
    _ideals_of_norm(
      E::Field,
      d::Union{ZZRingElem, QQFieldElem},
    ) -> Vector{fractional_ideal_type(E)}

Given a degree 2 extension $E/K$ of number fields, with maximal order $O_E$,
return the (fractional) $O_E$-ideals of abolute norm $d$.

# Examples
```jldoctest
julia> E, b = cyclotomic_field_as_cm_extension(8; cached=false);

julia> absolute_norm.(Oscar._ideals_of_norm(E, QQ(81//2)))
1-element Vector{QQFieldElem}:
 81//2
```
"""
function _ideals_of_norm(E::Field, d::QQFieldElem)
  OE = maximal_order(E)
  if denominator(d) == 1 # Integral ideals
    return _ideals_of_norm(E, numerator(d))
  elseif numerator(d) == 1 # Inverse of integral ideals
    return Hecke.fractional_ideal_type(OE)[inv(I) for I in _ideals_of_norm(E, denominator(d))]
  else # Honest fractional ideals
    return Hecke.fractional_ideal_type(OE)[I*inv(J) for I in _ideals_of_norm(E, numerator(d)) for J in _ideals_of_norm(E, denominator(d))]
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
  # We look for ideals by working p-adically at each prime dividing
  # d, and we then multiply everything back together
  for p in prime_divisors(d)
    v = valuation(d, p)
    pd = ideal_type(OK)[P[1] for P in prime_decomposition(OK, p)]
    # We look for the largest O_E-ideal over q which is invariant
    # under the involution of E/K. If q is inert or split, then the
    # corresponding ideal is q*O_E. Otherwise, we q*O_E = P^2 for some
    # P prime, and P is actually that ideal.
    for q in pd
      if !is_coprime(DE, ideal(OE, q)) # That q ramifies in O_E
        P = prime_decomposition(OE, q)[1][1]
      else
        P = ideal(OE, q)
      end
      nv = valuation(norm(P), q)
      # Here we overcount, but it should not be a serious bottleneck.
      # A better way would be to only keep products of ideals whose absolute
      # norm is exactly p^v
      push!(primes, ideal_type(OE)[P^e for e in 0:divrem(v, nv)[1]])
    end
  end
  for _I in Hecke.cartesian_product_iterator(primes)
    if prod(absolute_norm.(_I)) != d
      continue
    end
    I = prod(_I)
    @hassert :ZZLatWithIsom 1 absolute_norm(I) == d
    push!(ids, fractional_ideal(OE, I))
  end
  return ids
end

# Given a degree 2 extension of number fields E/K, return all
# the possible signatures dictionaries of any hermitian lattice over
# `E/K` of rank `rk`, and whose trace lattice has negative signature `s2`.
# Note that `E/K` need not be CM.
function _possible_signatures(s2::Int, E::Field, rk::Int)
  # The negative signature of the trace lattice, which is s2, is the sum
  # of twice the sum of signatures at the real places of K which extend to
  # complex conjugate places in E, plus a remaining part at the other real
  # places of K.
  lb = iseven(s2) ? 0 : 1
  K = base_field(E)
  inf = Hecke.place_type(K)[p for p in real_places(K) if length(extend(p, E)) == 1]
  r = length(real_places(K)) - length(inf)
  s = length(inf)
  signs = Dict{Hecke.place_type(K), Int}[]
  perm = AllPerms(s)
  for _l in lb:2:min(s2, rk*r)
    parts = Vector{Int}[]
    l = divexact(s2-_l, 2) # s2 = 2*l + _l
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

# Given a genus `G` of hermitian lattices over a degree 2 extension of number
# field `E/K` and given a QQ-automorphism `s` of `K`, return the genus `G*s`
# obtained after letting `s` acts on the primes of `G` and on the real
# embeddings of `K`.
function _action_on_genus(G::HermGenus, s::NumFieldHom)
  E = base_field(G)
  t = inv(s)
  lgs = local_symbols(G)
  si = signatures(G)
  si_new = empty(si)
  for (p, v) in si # Action on signatures
    si_new[Hecke.induce_image(t, p)] = v
  end

  lgs_new = empty(lgs)
  for sym in lgs # Action on prime
    _s = sym.data
    # sym.data is always a list of triples, but we also need the norm valuations
    # whenever sym is dyadic and ramified. For those cases, we add it to our local symbol.
    if is_dyadic(sym) && is_ramified(sym)
      nn = sym.norm_val
      @hassert :ZZLatWithIsom 1 length(nn) == length(_s)
      s = Tuple{Int, Int, Int, Int}[(_s[i][1], _s[i][2], _s[i][3], nn[i]) for i in 1:length(nn)]
      push!(lgs_new, genus(HermLat, E, t(sym.p), s))
    else
      push!(lgs_new, genus(HermLat, E, t(sym.p), _s))
    end
  end

  @hassert :ZZLatWithIsom 1 Hecke._check_global_genus(lgs_new, si_new)
  G_new = HermGenus(E, rank(G), lgs_new, si_new) # = G*s
  return G_new
end

@doc raw"""
    representatives_of_hermitian_type(
      Lf::ZZLatWithIsom,
      m::Int = 1,
      fix_root::Int = -1;
      cond::Vector{Int}=Int[-1, -1, -1],
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    )  -> Vector{ZZLatWithIsom}

Given a lattice with isometry $(L, f)$ of finite hermitian type, i.e. $f$ is of
finite order $n$ and its minimal polynomial is irreducible cyclotomic, and a
positive integer $m$, return a complete set of representatives for the
isomorphism classes of lattices with isometry $(M, g)$ of finite hermitian type
such that the type of $(M, g^m)$ is equal to the type of $(L, f)$.

Note that in this case, the isometries $g$'s are of order $nm$.

If the value of `fix_root` is exactly $nm$, then the function returns only
one generator for every conjugacy classes of finite cyclic groups generated
by an isometry $g$ as before.

Note that if `Lf` is trivial, the algorithm returns `Lf` by default.

See also [`type(::ZZLatWithIsom)`](@ref).

!!! note "For the advanced users"
    When using this function in a larger algorithm, one can use some keyword
    arguments which can be carried along the computations.
    - by setting the values in `cond` to the desired values (in order: rank,
      positive signature, negative signature), one can control whether the
      lattice $L$ in input has the good rank or signatures;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists of lattices of maximum $-2$. In such a case, the enumeration is
      skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

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
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  @req m >= 1 "m must be a positive integer"
  @req is_of_hermitian_type(Lf) "Lf must be of hermitian type"

  n = order_of_isometry(Lf)
  @req is_finite(n) "Isometry must be of finite order"
  k = n*m
  _reps = representatives_of_hermitian_type(genus(Lf), cyclotomic_polynomial(k), fix_root; cond, genusDB, root_test, info_depth, discriminant_annihilator, _local)
  # We test the type condition
  if fix_root == k
    # In this case, we have fixed a generator for a cyclic group: we need
    # to find a generator which satisfies the type condition. If there is
    # one, we keep it, otherwise we discard the lattice with isometry
    Sk = Int[i for i in 1:k if isone(gcd(i, k))]
    reps = empty(_reps)
    while !isempty(_reps)
      M = pop!(_reps)
      j = findfirst(l -> is_of_same_type(M^(l*m), Lf), Sk)
      if isnothing(j)
        continue
      end
      l = Sk[j]
      push!(reps, M^l)
    end
  else
    reps = filter!(M -> is_of_same_type(M^m, Lf), _reps)
  end
  return reps
end

@doc raw"""
    representatives_of_hermitian_type(
      G::Union{ZZGenus, ZZLat},
      m::Int,
      fix_root::Int = -1;
      first::Bool=false,
      cond::Vector{Int}=Int[-1, -1, -1],
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_discriminant_annihilator(G),
    )  -> Vector{ZZLatWithIsom}

Given a nonempty genus of integer lattices $G$, return a list of
representatives of isomorphism classes of pairs $(M, g)$ consisting of a
lattice $M$ in $G$ and $g \in O(M)$ is an isometry of minimal polynomial
$\Phi_m(X)$, the $m$th cyclotomic polynomial.

If $m = 1, 2$, this goes back to enumerating $G$ as a genus of integer
lattices.

If the value of `fix_root` is exactly $m$, then the function returns only one
generator for every conjugacy classes of finite cyclic groups generated
by an isometry of minimal polynomial $\Phi_m(X)$.

One can also provide a representative $L$ of $G$ instead.

If `first` is set to `true`, only return the first representative computed.

Note that if `G` is trivial, the algorithm returns the trivial lattice with
isometry by default.

!!! note "For the advanced users"
    When using this function in a larger algorithm, one can use some
    keyword arguments which can be carried along the computations.
    - by setting the values in `cond` to the desired values (in order: rank,
      positive signature, negative signature), one can control whether the
      genus $G$ in input has the good rank or signatures;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists of lattices of maximum $-2$. In such a case, the enumeration is
      skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

# Examples
```jldoctest
julia> L = root_lattice(:A, 4)
Integer lattice of rank 4 and degree 4
with gram matrix
[ 2   -1    0    0]
[-1    2   -1    0]
[ 0   -1    2   -1]
[ 0    0   -1    2]

julia> reps = representatives_of_hermitian_type(L, 5)
1-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 5

julia> is_of_hermitian_type(reps[1])
true
```
"""
representatives_of_hermitian_type(::Union{ZZGenus, ZZLat}, ::Int, ::Int)

function representatives_of_hermitian_type(
    G::ZZGenus,
    m::Int,
    fix_root::Int = -1;
    kwargs...,
  )
  return representatives_of_hermitian_type(G, cyclotomic_polynomial(m), fix_root; kwargs...)
end

function representatives_of_hermitian_type(
    L::ZZLat,
    m::Int,
    fix_root::Int = -1;
    kwargs...,
  )
  return representatives_of_hermitian_type(genus(L), cyclotomic_polynomial(m), fix_root; kwargs...)
end

@doc raw"""
    representatives_of_hermitian_type(
      G::Union{ZZGenus, ZZLat},
      min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
      fix_root::Int = -1;
      first::Bool=false,
      cond::Vector{Int}=Int[-1, -1, -1],
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_discriminant_annihilator(G),
    ) -> Vector{ZZLatWithIsom}

Given a nonempty genus of integer lattices ``G`` and a polynomial `min_poly`
which is irreducible over $\mathbb Q$ and so that the equation order of the
associated number field is maximal, return a complete list of representatives
for the isomorphism classes of pairs $(M, g)$ consisting of a lattice $M$ in
$G$ and $g \in O(M)$ is an isometry of minimal polynomial `min_poly`

One can also provide a representative $L$ of $G$ instead.

The value $n$ of `fix_root` matters only when `min_poly` is cyclotomic.
In the case where `min_poly` is the $n$th cyclotomic polynomial, the function
returns only one generator for every conjugacy classes of finite cyclic
groups generated by an isometry of minimal polynomial `min_poly`

If `first` is set to `true`, only return the first representative computed.

Note that if `G` is trivial, the algorithm returns the trivial lattice with
isometry by default.

!!! note "For the advanced users"
    When using this function in a larger algorithm, one can use some
    keyword arguments which can be carried along the computations.
    - by setting the values in `cond` to the desired values (in order: rank,
      positive signature, negative signature), one can control whether the
      genus $G$ in input has the good rank or signatures;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists of lattices of maximum $-2$. In such a case, the enumeration is
      skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

# Examples
```jldoctest
julia> L = root_lattice(:A, 6)
Integer lattice of rank 6 and degree 6
with gram matrix
[ 2   -1    0    0    0    0]
[-1    2   -1    0    0    0]
[ 0   -1    2   -1    0    0]
[ 0    0   -1    2   -1    0]
[ 0    0    0   -1    2   -1]
[ 0    0    0    0   -1    2]

julia> chi = cyclotomic_polynomial(7)
x^6 + x^5 + x^4 + x^3 + x^2 + x + 1

julia> reps = representatives_of_hermitian_type(L, chi)
1-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 7

julia> is_of_hermitian_type(reps[1])
true
```
"""
representatives_of_hermitian_type(::Union{ZZLat, ZZGenus}, ::Union{ZZPolyRingElem, QQPolyRingElem}, ::Int)

function representatives_of_hermitian_type(
    G::ZZGenus,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    first::Bool=false,
    cond::Vector{Int}=Int[-1, -1, -1],
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_discriminant_annihilator(G),
    _local::Bool=false,
  )
  chi = min_poly
  @req is_irreducible(chi) "Polynomial must be irreducible"
  @req is_integral(G) "For now G must be a genus symbol for integral lattices"
  allow_info = get_verbosity_level(:ZZLatWithIsom) >= info_depth
  reps = ZZLatWithIsom[]
  # Only relevant in bigger algorithm where we need to control invariants
  # of eigenlattices
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

  # Default output
  if is_zero(rG)
    return ZZLatWithIsom[integer_lattice_with_isometry(integer_lattice(; gram=matrix(QQ, 0, 0, QQFieldElem[])))]
  end

  d_chi = degree(chi)
  # In that case the isometry is +- id, so we do not need hermitian genera.
  if isone(d_chi)
    !is_zero(chi(1)*chi(-1)) && return reps
    f = is_zero(chi(1)) ? identity_matrix(QQ, rG) : -identity_matrix(QQ, rG)
    # the discriminant annihilator is of local nature, we can test before the genus 
    T = representative(G)
    Tf = integer_lattice_with_isometry(T, f)
    _, f_on_disc = discriminant_group(Tf)
    if !is_annihilated(f_on_disc, discriminant_annihilator)
      return reps
    end
    # local conditions are okay, enumerate the genus
    allow_info && println("Enumerate Z-genus $G")
    repre = oscar_genus_representatives(G; genusDB, root_test, info_depth, max_lat=first ? 1 : inf, _local)
    allow_info && println("$(length(repre)) representative(s)")
    while !is_empty(repre)
      LL = pop!(repre)
      push!(reps, integer_lattice_with_isometry(LL, f; check=false, ambient_representation=false))
      first && return reps
    end
    return reps
  end

  # Polynomial must be symmetric
  if !iseven(d_chi)
    return reps
  elseif any(i -> coeff(chi, i) != coeff(chi, d_chi-i), 0:div(d_chi, 2))
    return reps
  end

  ok, rk = divides(rG, d_chi)
  ok || return reps

  # Detect if we have a finite order isometry
  R = parent(chi)
  is_cyclo, n = is_cyclotomic_polynomial_with_data(chi)
  if is_cyclo
    E, b = cyclotomic_field_as_cm_extension(n)
  else
    Etemp, btemp = number_field(chi; cached=false)
    @req is_maximal(equation_order(Etemp)) "For isometries of infinite order, the equation order of the associated number field must be maximal (for now)"
    K, a = number_field(minpoly(btemp + inv(btemp)), "a"; cached=false)
    Kt, t = K[:t]
    E, b = number_field(t^2-a*t+1, "b"; cached=false)
  end

  gene = Hecke.genus_herm_type(E)[]
  DEK = different(maximal_order(E))
  DK = different(base_ring(maximal_order(E)))
  DE = DK*maximal_order(E)*DEK

  ndE = dG*inv(QQ(absolute_norm(DE)))^rk
  detE = _ideals_of_norm(E, ndE)
  isempty(detE) && return reps

  signatures = _possible_signatures(nG, E, rk)
  isempty(signatures) && return reps

  discriminant_annihilatorE = begin
    if discriminant_annihilator isa ZZPolyRingElem
      [evaluate(discriminant_annihilator, b)]
    elseif discriminant_annihilator isa ZZMPolyRingElem
      [evaluate(discriminant_annihilator, [b])]
    else
     [evaluate(p, [b]) for p in gens(discriminant_annihilator)]
    end
  end

  discriminant_annihilatorOE = fractional_ideal(maximal_order(E), discriminant_annihilatorE)
  for dd in detE, sign in signatures
    append!(gene, hermitian_genera(E, rk, sign, dd; min_scale=inv(DE), max_scale=inv(DE)*discriminant_annihilatorOE + numerator(dd)*DE))
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

  allow_info &&  println("All possible hermitian genera $G  $min_poly: $(length(gene))")
  for g in gene
    if is_integral(G) && !is_integral(DE*scale(g))
      continue
    end
    if is_even(G) && !is_integral(DK*norm(g))
      continue
    end

    H = representative(g)
    M, fM = trace_lattice_with_isometry(H)
    genus(M) != G && continue

    MfM = integer_lattice_with_isometry(M, fM; check=false)
    @hassert :ZZLatWithIsom 1 is_of_hermitian_type(MfM)
    first && return ZZLatWithIsom[MfM]

    allow_info && println("Enumerate hermitian genus of rank $(rank(H))")
    if !_local
      gr = genus_representatives(H)
    else 
      gr = [H]
    end
    for HH in gr
      M, fM = trace_lattice_with_isometry(HH)
      push!(reps, integer_lattice_with_isometry(M, fM; check=false))
    end
  end
  @hassert :ZZLatWithIsom 1 all(Base.Fix2(is_annihilated_on_discriminant, discriminant_annihilator), reps)
  return reps
end

function representatives_of_hermitian_type(
    L::ZZLat,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    kwargs...,
  )
  return representatives_of_hermitian_type(genus(L), min_poly, fix_root; kwargs...)
end

@doc raw"""
    representatives_of_hermitian_type(
      signatures::Vector{Tuple{Int, Int}},
      determinants::AbstractVector{T},
      min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
      fix_root::Int = -1;
      even::Bool=true,
      first::Bool=false,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
    ) where T <: IntegerUnion -> Vector{ZZLatWithIsom}

Given a finite list of pairs of nonnegative integers `signatures`, a finite
list of integers `determinants` and a polynomial `min_poly` which is
irreducible over $\mathbb Q$ and such that the equation order of the associated
number field is maximal, return a complete list of representatives for the
isomorphism classes of pairs $(M, g)$ consisting of an integral
$\mathbb{Z}$-lattice $M$ of signatures in `signatures` and determinant in
`determinants` (up to sign) and $g \in O(M)$ is an isometry with minimal
polynomial `min_poly`.

If `even` is set to false, then the lattices in output are allowed to be odd.

The value $n$ of `fix_root` matters only when `min_poly` is cyclotomic.
In the case where `min_poly` is the $n$th cyclotomic polynomial, the function
returns only one generator for every conjugacy classes of finite cyclic
groups generated by an isometry of minimal polynomial `min_poly`.

If `first` is set to `true`, only return the first representative computed.

One can also decide to only input a single pair of signatures `signatures` or a
single value for `determinants` without having to create lists consisting of one
element.

!!! note "For the advanced users"
    One can use some keyword arguments for further features:
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists of lattices of maximum $-2$. In such a case, the enumeration is
      skipped.

# Examples
```jldoctest
julia> r = representatives_of_hermitian_type(Tuple{Int, Int}[(2,4), (2, 10)], Int[1, 3, 9, 27], cyclotomic_polynomial(9), 9)
4-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 9
 Integer lattice with isometry of finite order 9
 Integer lattice with isometry of finite order 9
 Integer lattice with isometry of finite order 9

julia> genus.(r)
4-element Vector{ZZGenus}:
 Genus symbol: II_(2, 10)
 Genus symbol: II_(2, 10) 3^-2
 Genus symbol: II_(2, 4) 3^1
 Genus symbol: II_(2, 4) 3^-3
```
"""
function representatives_of_hermitian_type(
    signatures::Vector{Tuple{Int, Int}},
    determinants::AbstractVector{T},
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    even::Bool=true,
    first::Bool=false,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    _local::Bool=false,
  ) where T <: Hecke.IntegerUnion
  chi = min_poly
  @req is_irreducible(chi) "Polynomial must be irreducible"
  @req all(>((0, 0)), signatures) "Signature tuples must consist of nonnegative integers with nonzero sum"
  allow_info = get_verbosity_level(:ZZLatWithIsom) >= info_depth
  reps = ZZLatWithIsom[]

  signs = unique(signatures)
  dets = unique(determinants)
  rGs = Dict{Int, Vector{Tuple{Int, Int}}}()
  for s in signs
    r = sum(s)
    if haskey(rGs, r)
      push!(rGs[r], s)
    else
      rGs[r] = typeof(s)[s]
    end
  end

  d_chi = degree(chi)
  if isone(d_chi)
    !is_zero(chi(1)*chi(-1)) && return reps
    for rG in keys(rGs)
      f = is_zero(chi(1)) ? identity_matrix(QQ, rG) : -identity_matrix(QQ, rG)
      int_gene = ZZGenus[]
      for sG in rGs[rG], dG in dets
        allow_info && println("Enumerate Z-genera $(sG)")
        append!(int_gene, integer_genera(sG, abs(dG); even))
        first && !isempty(int_gene) && break
      end
      repre = ZZLat[]
      for G in int_gene
        append!(repre, oscar_genus_representatives(G; genusDB, root_test, info_depth, max_lat=first ? 1 : inf, _local))
        first && !isempty(repre) && break
      end
      allow_info && println("$(length(repre)) representative(s)")
      while !is_empty(repre)
        LL = pop!(repre)
        push!(reps, integer_lattice_with_isometry(LL, f; check=false))
        first && return reps
      end
    end
    return reps
  end

  if !iseven(d_chi)
    return reps
  elseif any(i -> coeff(chi, i) != coeff(chi, d_chi-i), 0:div(d_chi, 2))
    return reps
  end

  rks = Dict{Int, Vector{Int}}()
  new_signs = empty(signs)
  for rG in keys(rGs)
    ok, rk = divides(rG, d_chi)
    if ok
      rks[rk] = Int[s[2] for s in rGs[rG]]
    end
  end

  R = parent(chi)
  is_cyclo, n = is_cyclotomic_polynomial_with_data(chi)
  if is_cyclo
    E, b = cyclotomic_field_as_cm_extension(n)
  else
    Etemp, btemp = number_field(chi; cached=false)
    @req is_maximal(equation_order(Etemp)) "For isometries of infinite order, the equation order of the associated number field must be maximal (for now)"
    K, a = number_field(minpoly(btemp + inv(btemp)), "a"; cached=false)
    Kt, t = K[:t]
    E, b = number_field(t^2-a*t+1, "b"; cached=false)
  end

  DEK = different(maximal_order(E))
  DK = different(base_ring(maximal_order(E)))
  DE = DK*maximal_order(E)*DEK

  for rk in keys(rks)
    gene = Hecke.genus_herm_type(E)[]
    ndEs = QQFieldElem[abs(dG)*inv(QQ(absolute_norm(DE)))^rk for dG in dets]
    unique!(ndEs)

    detE = Hecke.fractional_ideal_type(maximal_order(E))[]
    for ndE in ndEs
      append!(detE, _ideals_of_norm(E, ndE))
    end
    isempty(detE) && return reps

    hermitian_signatures = Dict{Hecke.place_type(base_field(E)), Int}[]
    for nG in rks[rk]
      append!(hermitian_signatures, _possible_signatures(nG, E, rk))
    end
    isempty(signatures) && return reps

    for dd in detE, sign in hermitian_signatures
      append!(gene, hermitian_genera(E, rk, sign, dd; min_scale=inv(DE), max_scale=numerator(dd)*DE))
    end
    isempty(gene) && return reps
    unique!(gene)

    if is_cyclo && fix_root == n
      K = base_field(E)
      GKQ, j = automorphism_group(K)
      phi = inv(isomorphism(PermGroup, GKQ))
      omega = gset(domain(phi), (G, g) -> _action_on_genus(G, j(phi(inv(g)))), gene)
      gene = representative.(orbits(omega))
    end

    allow_info && println("All possible hermitian genera: $(length(gene))")
    for g in gene
      if !is_integral(DE*scale(g))
        continue
      end
      if even && !is_integral(DK*norm(g))
        continue
      end

      H = representative(g)
      if !(volume(H) in detE)
        continue
      end
      M, fM = trace_lattice_with_isometry(H)

      MfM = integer_lattice_with_isometry(M, fM; check=false)
      @hassert :ZZLatWithIsom 1 is_of_hermitian_type(MfM)
      first && return ZZLatWithIsom[MfM]

      allow_info && println("Enumerate hermitian genus of rank $(rank(H))")
      if !_local
        gr = genus_representatives(H)
      else 
        gr = [H]
      end
      for HH in gr
        M, fM = trace_lattice_with_isometry(HH)
        push!(reps, integer_lattice_with_isometry(M, fM; check=false))
      end
    end
  end
  @hassert :ZZLatWithIsom 2 all(M -> signature_tuple(M)[[1,3]] in signatures, reps)
  @hassert :ZZLatWithIsom 2 all(M -> abs(det(M)) in abs.(dets), reps)
  return reps
end

function representatives_of_hermitian_type(
    signature::Tuple{Int, Int},
    determinant::Hecke.IntegerUnion,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    kwargs...,
  )
  return representatives_of_hermitian_type(
                                           Tuple{Int, Int}[signature],
                                           typeof(determinant)[determinant],
                                           min_poly,
                                           fix_root;
                                           kwargs...,
                                          )
end

function representatives_of_hermitian_type(
    signature::Tuple{Int, Int},
    determinants::AbstractVector{T},
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    kwargs...,
  ) where T <: Hecke.IntegerUnion
  return representatives_of_hermitian_type(
                                           Tuple{Int, Int}[signature],
                                           determinants,
                                           min_poly,
                                           fix_root;
                                           kwargs...,
                                          )
end

function representatives_of_hermitian_type(
    signatures::Vector{Tuple{Int, Int}},
    determinant::Hecke.IntegerUnion,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem},
    fix_root::Int = -1;
    kwargs...,
  )
  return representatives_of_hermitian_type(
                                           signatures,
                                           typeof(determinant)[determinant],
                                           min_poly,
                                           fix_root;
                                           kwargs...,
                                          )
end

###############################################################################
#
#  Splitting procedures from [BH23]
#
###############################################################################

@doc raw"""
    splitting_of_hermitian_type(
      Lf::ZZLatWithIsom,
      p::IntegerUnion,
      b::Int = 0;
      eiglat_cond::Vector{Dict{Int,Vector{Int}}}=[Dict{Int, Vector{Int}}()],
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      check::Bool=true,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ of hermitian type with $f$ of
finite order $m$, and given a prime number $p$, return a complete set of
representatives for the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.

The integer `b` can be set to be `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*m$.

For every $(M, g)$ in output, one may decide on the rank `r`, the positive
signature `p` and the negative signature `n` of the eigenlattices of $(M, g)$
using the keyword argument `eiglat_cond`. It should consist of a dictionary
where each key is a divisor of $p*m$, and the corresponding value is a tuple
`(r, p, n)` of integers. Any undetermined value can be set to a negative
number; for instance $(-1, 2, -4)$ means that the associated eigenlattice
must have positive signature 2, without restriction on its rank and
negative signature.

If the keyword argument `check` is set to `true`, the function tests whether
$(L, f)$ is of hermitian type.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

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
    eiglat_cond::Vector{Dict{Int, Vector{Int}}}=[Dict{Int, Vector{Int}}()],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=true,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"

  if rank(Lf) == 0
    (b == 0) && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @req iseven(Lf) "Lattice must be even"
  @req is_finite(n) "Isometry must be of finite order"
  @req is_prime(p) "p must be a prime number"
  @check begin
    @req is_of_hermitian_type(Lf) "Lattice with isometry must be of hermitian type"
  end

  reps = ZZLatWithIsom[]
  k = p*n
  # If p divides n, then any output is still of hermitian type since taking pth
  # power decreases the order of the isometry
  if iszero(mod(n, p))
    for cond in unique!([get(i, k, Int[-1, -1, -1]) for i in eiglat_cond])
      append!(reps,representatives_of_hermitian_type(Lf, p, fix_root; cond, genusDB, root_test, info_depth, discriminant_annihilator))
    end
    return reps
  end

  V = _from_cyclotomic_polynomial_to_dict(characteristic_polynomial(Lf))
  for i in keys(V)
    V[i] = euler_phi(i)*V[i]
  end
  
  # avoid overcounting
  eiglat_cond_trimmed = Vector{Dict{Int, Vector{Int}}}()
  for cond in eiglat_cond
    T = Dict{Int64,Vector{Int64}}()
    T[n] = get(cond, n, Int[-1, -1, -1])
    T[k] = get(cond, k, Int[-1, -1, -1])
    push!(eiglat_cond_trimmed, T)
  end 
  unique!(eiglat_cond_trimmed)


  for cond in eiglat_cond_trimmed
    condp = _conditions_after_power(cond, p)
    if !all(get(condp, i, [-1])[1] == V[i] || get(condp, i, [-1])[1] == -1 for i in keys(V))
      continue
    end
    tmp = __splitting_of_hermitian_type(Lf, p, b; eiglat_cond=cond, fix_root, genusDB, root_test, check, info_depth, discriminant_annihilator, _local)
    append!(reps, tmp)
  end
  return reps
end
  

function __splitting_of_hermitian_type(
    Lf::ZZLatWithIsom,
    p::IntegerUnion,
    b::Int = 0;
    eiglat_cond::Dict{Int, Vector{Int}}=Dict{Int, Vector{Int}}(),
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    check::Bool=true,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  n = order_of_isometry(Lf)
  k = p*n
  
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"
  if fix_root == k
    Sk = Int[i for i in 1:k if isone(gcd(i, k))]
  end
  reps = ZZLatWithIsom[]
  # If p does not divide n, then the characteristic polynomials of the
  # isometries in output are of the form \Phi_n^a*\Phi_k^b where a,b are
  # nonnegative. We can therefore use admissible triples to find the
  # potential genera of the corresponding eigenlattices (Lemma 4.15 [BH23])
  # First we do a bit a work on the potential conditions imposed to the
  # invariants of the corresponding eigenlattices

  # Ranks conditions
  phi_n = euler_phi(n)
  phi_k = euler_phi(k)
  rL = rank(Lf)

  rA, pA, nA = get(eiglat_cond, n, Int[-1, -1, -1])
  intrA::AbstractVector{Int} = rA >= 0 ? Int[rA] : 0:phi_n:rL
  rB, pB, nB = get(eiglat_cond, k, Int[-1, -1, -1])
  intrB::AbstractVector{Int} = rB >= 0 ? Int[rB] : b*phi_k:phi_k:rL
  # We put the rank conditions all in one range of values for the rank of
  # the \Phi_n-eigenlattice
  intrA = intersect(intrA, (rL).-(intrB))
  isempty(intrA) && return reps

  # Signatures conditions
  pL, _, nL = signature_tuple(Lf)
  if pA >= 0
    IpA = Int[pA]
    IpB = Int[pL - pA]
    if pB >= 0
      (pB == IpB[1]) || return reps
    end
  elseif pB >= 0
    IpA = Int[pL - pB]
    IpB = Int[pB]
  end
  if nA >= 0
    InA = Int[nA]
    InB = Int[nL - nA]
    if nB >= 0
      (nB == InB[1]) || return reps
    end
  elseif nB >= 0
    InA = Int[nL - nB]
    InB = Int[nB]
  end

  for _rA in intrA
    _rB = rL - _rA
    if pA < 0 && pB < 0
      IpA = 0:max(_rA, pL)
      IpB = (pL).-(IpA)
    end
    if nA < 0 && nB < 0
      InA = 0:max(_rA, nL)
      InB = (nL).-(InA)
    end
    # We follow part the same ideas as in Algorithm 4 of [BH23]
    atp = admissible_triples(Lf, p; IrA=[_rA], IpA, InA, IrB=[_rB], IpB, InB, b)
    for (A, B) in atp
      if root_test
        if rank(A) > 0 && _packing_density_test(A)
          continue
        elseif rank(B) > 0 && _packing_density_test(B)
          continue
        end
      end
      As = representatives_of_hermitian_type(A, n, fix_root; genusDB, info_depth, discriminant_annihilator=p*discriminant_annihilator, _local)
      if root_test && iszero(signature_tuple(A)[1]) # Remove lattices with (-2)-vectors
        filter!(LA -> rank(LA) == 0 || minimum(LA) != 2, As)
      end
      isempty(As) && continue
      Bs = representatives_of_hermitian_type(B, k, fix_root; genusDB, info_depth, discriminant_annihilator=p*discriminant_annihilator, _local)
      if root_test && iszero(signature_tuple(B)[1]) # Remove lattices with (-2)-vectors
        filter!(LB -> rank(LB) == 0 || minimum(LB) != 2, Bs)
      end
      isempty(Bs) && continue
      for LA in As, LB in Bs
        Es = admissible_equivariant_primitive_extensions(LA, LB, Lf, p; check=false, test_type=false, _local)
        if fix_root == k
          while !isempty(Es)
            M = pop!(Es)
            j = findfirst(l -> is_of_same_type(M^(l*p), Lf), Sk)
            if isnothing(j)
              continue
            end
            l = Sk[j]
            push!(reps, M^l)
          end
        else
          filter!(M -> is_of_same_type(M^p, Lf), Es)
          append!(reps, Es)
        end
        filter!(Base.Fix2(is_annihilated_on_discriminant, discriminant_annihilator), Es)
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
      eiglat_cond::Vector{Dict{Int,Vector{Int}}}=[Dict{Int, Vector{Int}}()],
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ with $f$ of order $m = q^e$ for
some prime number $q$ and a prime number $p \neq q$, return a complete set of
representatives for the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.
Note that $e$ can be `0`, in which case $f = id$ is trivial.

The integer `b` can be set to be `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*m$.

For every $(M, g)$ in output, one may decide on the rank `rM`, the positive
signature `pM` and the negative signature `nM` of the eigenlattices of $(M, g)$
using the keyword argument `eiglat_cond`. It should consist of a dictionary
where each key is a divisor of $p*m$, and the corresponding value is a tuple
`(rM, pM, nM)` of integers. Any undetermined value can be set to a negative
number; for instance $(-1, 2, -4)$ means that the associated eigenlattice
must have positive signature 2, without restriction on its rank and
negative signature.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

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
    eiglat_cond::Vector{Dict{Int, Vector{Int}}}=[Dict{Int, Vector{Int}}()],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"
  # Default output
  if rank(Lf) == 0
    (b == 0) && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end
  n = order_of_isometry(Lf)
  @req iseven(Lf) "Lattice must be even"
  @req is_finite(n) "Isometry must be of finite order"
  @req is_prime(p) "p must be a prime number"

  # In this case the pair (L, f) is of hermitian type so we can fallback to the
  # previous function.
  if isone(n)
    return splitting_of_hermitian_type(Lf, p, b; eiglat_cond, fix_root, genusDB, root_test, info_depth, check=false, discriminant_annihilator, _local)
  end

  ok, e, q = is_prime_power_with_data(n)

  @req ok || e == 0 "Order of isometry must be a prime power"
  @req p != q "Prime numbers must be distinct"

  reps = ZZLatWithIsom[]

  # We follow Algorithm 5 of [BH23]: there is a slight mistake in the pseudocode
  # though, `A_0` and `B_0` should be switched for calling the algorithm
  # `PrimitiveExtensions`.
  x = gen(Hecke.Globals.Qx)
  A0 = kernel_lattice(Lf, x^(q^(e-1))-1)
  B0 = kernel_lattice(Lf, q^e)
  # Compute this one first because it is faster to decide whether it is empty
  RB = splitting_of_hermitian_type(B0, p; eiglat_cond, fix_root, genusDB, root_test, info_depth, check=false, discriminant_annihilator=q*discriminant_annihilator, _local)
  is_empty(RB) && return reps
  # Recursive part of the function, with termination when the isometry of A0
  # is trivial
  RA = splitting_of_prime_power(A0, p; eiglat_cond, genusDB, root_test, info_depth=info_depth+1, discriminant_annihilator=q*discriminant_annihilator, _local)
  is_empty(RA) && return reps
  for L1 in RA, L2 in RB
    n1 = order_of_isometry(L1)::Int
    n2 = order_of_isometry(L2)::Int
    # Condition Line 11 of Algorithm 5 of [BH23]
    # We can decide before gluings whether the isometries in output will have
    # the correct order
    b == 1 && !is_divisible_by(n1, p) && !is_divisible_by(n2, p) && continue
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, q, p; check=false, _local)
    @hassert :ZZLatWithIsom 1 b == 0 || all(LL -> order_of_isometry(LL) == p*q^e, E)
    filter!(Base.Fix2(is_annihilated_on_discriminant, discriminant_annihilator), E)
    append!(reps, E)
  end
  return reps
end

@doc raw"""
    splitting_of_pure_mixed_prime_power(
      Lf::ZZLatWithIsom,
      p::Int;
      eiglat_cond::Vector{Dict{Int, Vector{Int}}}=[Dict{Int, Vector{Int}}()],
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ with $f$ of order $m = p^d*q^e$
where $d > 0$ is positive, $e \geq 0$ is nonnegative and $q \neq p$ is a prime
number distinct from $p$, and such that $\prod_{i=0}^e\Phi_{p^dq^i}(f)$ is
trivial, return a complete set of representatives for the isomorphism classes
of lattices with isometry $(M, g)$ such that the type of $(M, g^p)$ is equal
to the type of $(L, f)$.
Note that $e$ can be `0`, while $d$ has to be positive.

For every $(M, g)$ in output, one may decide on the rank `rM`, the positive
signature `pM` and the negative signature `nM` of the eigenlattices of $(M, g)$
using the keyword argument `eiglat_cond`. It should consist of a dictionary
where each key is a divisor of $p*m$, and the corresponding value is a tuple
`(rM, pM, nM)` of integers. Any undetermined value can be set to a negative
number; for instance $(-1, 2, -4)$ means that the associated eigenlattice
must have positive signature 2, without restriction on its rank and
negative signature.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

"""
function splitting_of_pure_mixed_prime_power(
    Lf::ZZLatWithIsom,
    p::Int;
    eiglat_cond::Vector{Dict{Int, Vector{Int}}}=[Dict{Int, Vector{Int}}()],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  rank(Lf) == 0 && return ZZLatWithIsom[Lf]

  n = order_of_isometry(Lf)

  @req iseven(Lf) "Lattice must be even"
  @req is_finite(n) "Isometry must be of finite order"
  @req is_prime(p) "p must be a prime number"

  pd = prime_divisors(n)

  reps = ZZLatWithIsom[]
  @req length(pd) <= 2 && p in pd "Order must be divisible by p and have at most 2 prime divisors"
  # In that case (L, f) is of hermitian type, so we can call the appropriate
  # function
  if length(pd) == 1
    for cond in unique!([get(i, p*n, Int[-1, -1, -1]) for i in eiglat_cond])
      append!(reps, representatives_of_hermitian_type(Lf, p, fix_root; cond, genusDB, root_test, info_depth, discriminant_annihilator, _local))
    end
    return reps
  end

  q = pd[1] == p ? pd[2] : pd[1]
  d = valuation(n, p)::Int  # Weird ? Valuation is type stable and returns Int but it bugs here
  e = valuation(n, q)::Int

  phi = minimal_polynomial(Lf)
  chi = prod(cyclotomic_polynomial(p^d*q^i, parent(phi)) for i in 0:e; init=one(phi))

  @req is_divisible_by(chi, phi) "Minimal polynomial is not of the correct form"

  bool, r = divides(phi, cyclotomic_polynomial(n, parent(phi)))
  @hassert :ZZLatWithIsom 1 bool

  # We follow Algorithm 6 of [BH23]: there is a slight mistake in the pseudocode
  # though, `A_0` and `B_0` should be switched for calling the algorithm
  # `PrimitiveExtensions`.
  A0 = kernel_lattice(Lf, r)
  B0 = kernel_lattice(Lf, n)
  # Compute this one first because it is faster to decide whether it is empty
  condB = unique!([get(i, p*n, Int[-1, -1, -1]) for i in eiglat_cond])
  RB = ZZLatWithIsom[]
  for cond in condB
    append!(RB, representatives_of_hermitian_type(B0, p, fix_root; cond=cond, genusDB, root_test, info_depth, discriminant_annihilator=q*discriminant_annihilator, _local))
  end
  is_empty(RB) && return reps
  RA = splitting_of_pure_mixed_prime_power(A0, p; eiglat_cond, genusDB, root_test, info_depth=info_depth+1, discriminant_annihilator=q*discriminant_annihilator, _local)
  is_empty(RA) && return reps
  for L1 in RA, L2 in RB
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, q, p; check=false, _local)
    filter!(Base.Fix2(is_annihilated_on_discriminant, discriminant_annihilator), E)
    append!(reps, E)
  end
  return reps
end
  

@doc raw"""
    splitting_of_mixed_prime_power(
      Lf::ZZLatWithIsom,
      p::IntegerUnion,
      b::Int = 1;
      eiglat_cond::Vector{Dict{Int,Vector{Int}}}=[Dict{Int, Vector{Int}}()],
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ and a prime number $p$ such that
$f$ has order $m=p^dq^e$ for some prime number $q \neq p$, return a set
of representatives of the isomorphism classes of lattices with isometry
$(M, g)$ such that the type of $(M, g^p)$ is equal to the type of $(L, f)$.
Note that $d$ and $e$ can be both `0`.

The integer `b` can be set to be `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*m$.

For every $(M, g)$ in output, one may decide on the rank `rM`, the positive
signature `pM` and the negative signature `nM` of the eigenlattices of $(M, g)$
using the keyword argument `eiglat_cond`. It should consist of a dictionary
where each key is a divisor of $p*m$, and the corresponding value is a tuple
`(rM, pM, nM)` of integers. Any undetermined value can be set to a negative
number; for instance $(-1, 2, -4)$ means that the associated eigenlattice
must have positive signature 2, without restriction on its rank and
negative signature.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

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
    eiglat_cond::Vector{Dict{Int, Vector{Int}}} = [Dict{Int, Vector{Int}}()],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"

  # Default output
  if rank(Lf) == 0
    b == 0 && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @req iseven(Lf) "Lattice must be even"
  @req is_finite(n) "Isometry must be of finite order"
  @req is_prime(p) "p must be a prime number"

  pd = prime_divisors(n)

  @req length(pd) <= 2 "Order must have at most 2 prime divisors"

  # In this case, the isometry f is of prime power order, so we can call
  # the appropriate function
  if !(p in pd)
    return splitting_of_prime_power(Lf, p, b; eiglat_cond, fix_root, genusDB, root_test, info_depth, discriminant_annihilator, _local)
  end

  d = valuation(n, p)
  _, e, q = is_prime_power_with_data(divexact(n, p^d))

  reps = ZZLatWithIsom[]

  x = gen(parent(minimal_polynomial(Lf)))
  # We follow Algorithm 7 of [BH23]: there is a slight mistake in the pseudocode
  # though, `A_0` and `B_0` should be switched for calling the algorithm
  # `PrimitiveExtensions`.
  A0 = kernel_lattice(Lf, x^(divexact(n, p)) - 1)
  B0 = kernel_lattice(Lf, prod(cyclotomic_polynomial(p^d*q^i) for i in 0:e))
  # Compute this one first because it is faster to decide whether it is empty
  RB = splitting_of_pure_mixed_prime_power(B0, p; eiglat_cond, fix_root, genusDB, root_test, info_depth, discriminant_annihilator=p*discriminant_annihilator, _local)
  isempty(RB) && return reps
  RA = splitting_of_mixed_prime_power(A0, p, 0; eiglat_cond, fix_root, genusDB, root_test, info_depth=info_depth+1, discriminant_annihilator=p*discriminant_annihilator, _local)
  is_empty(RA) && return reps
  for L1 in RA, L2 in RB
    E = admissible_equivariant_primitive_extensions(L1, L2, Lf, p; check=false, _local)
    b == 1 && filter!(LL -> order_of_isometry(LL) == p*n, E)
    filter!(Base.Fix2(is_annihilated_on_discriminant, discriminant_annihilator), E)
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
      char_polys::Vector=[nothing],
      min_polys::Vector=[nothing],
      rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      eiglat_cond::Vector{Dict{Int64, Vector{Int64}}}=Dict{Int64, Vector{Int64}}[],
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    ) -> Vector{ZZLatWithIsom}

Given an even lattice with isometry $(L, f)$ where $f$ is of finite order $m$,
and given a prime number $p$ return a complete set of representatives for the
isomorphism classes of lattices with isometry $(M, g)$ such that the type of
$(M, g^p)$ is equal to the type of $(L, f)$.

The integer `b` can be set to be `0` or `1`, depending on whether one allows
`Lf` to be among the outputs or not. For instance, setting `b = 1`
would enforce every $(M, g)$ in output to satisfy that $g$ has order $p*m$.

For every $(M, g)$ in output, one may decide on the minimal polynomial (resp.
the characteristic polynomial) of $g$ by setting the keyword argument `min_poly`
(resp. `char_poly`) to the desired value (or `char_polys` resp. `min_poly` 
to a list of possibilities). Moreover, one can also decide on the rank, 
positive signature and negative signature of the eigenlattices of $(M, g)$
using the keyword arguments `rks`, `pos_sigs` and `neg_sigs` respectively. Each
list should consist of tuples $(k, i)$ where $k$ is a divisor of $p*m$ and $i$
is the value to be assigned for the given property (rank, positive/negative
signature) to the corresponding $\Phi_k$-eigenlattice. All these conditions
will be first condensed in a dictionary keeping track, for each divisor $k$
of $p*m$ of the potential rank `rk`, positive signature `pk` and negative
signature `nk` of the corresponding $\Phi_k$-eigenlattice. The keys of such
dictionary are the divisors $k$, and the corresponding value is the vector
`[rk, pk, nk]`. Any undetermined value will be set automatically to $-1$ by
default.
If one already knows such a dictionary, one can choose it as input under the
keyword argument `eiglat_cond`.

!!! warning
    In the case where the order of the isometries in output has at most
    2 prime divisors, we rely on the machinery of [BH23](@cite). Otherwise,
    we use a naive approach which consists of decomposing $(L, f)$ into its
    irreducible eigenlattices, splitting each of them by `p` (using the
    function [`splitting_of_hermitian_type`](@ref)), and gluing back all
    blocks together.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

"""
function splitting(
    Lf::ZZLatWithIsom,
    p::Int,
    b::Int = 0;
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    char_polys::Vector=[nothing],
    min_polys::Vector=[nothing],
    rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    eiglat_cond::Vector{Dict{Int64, Vector{Int64}}}= Dict{Int64, Vector{Int64}}[],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_default_discriminant_annihilator(Lf),
    _local::Bool=false,
  )
  @req b == 0 || b == 1 "b must be an integer equal to 0 or 1"

  # Default output
  if rank(Lf) == 0
    (b == 0) && return ZZLatWithIsom[Lf]
    return ZZLatWithIsom[]
  end

  n = order_of_isometry(Lf)
  @req iseven(Lf) "Lattice must be even"
  @req is_finite(n) "Isometry must be of finite order"
  @req is_prime(p) "p must be a prime number"

  if min_poly === nothing
    if char_poly !== nothing
      min_poly = _radical(char_poly)
    end
  end
  
  # If the user does not already input a dictionary of conditions on the
  # eigenlattices, we create one based on the conditions imposed on the
  # characteristic/minimal polynomials, on the ranks and on the signatures
  if char_poly !== nothing
    char_polys = [char_poly]
  end 
  if min_poly !==nothing
    min_polys=[min_poly]
  end
  if isempty(eiglat_cond)
    eiglat_cond = [_conditions_from_input(p*n, char_p, min_p, rks, pos_sigs, neg_sigs) for char_p in char_polys for min_p in min_polys]
  end
  
  pds = prime_divisors(n)
  # If the order of the isometry f is a prime power, or a power of p times
  # another prime power, then we can call the machinery from [BH23].
  if (length(pds) <= 1) || (length(pds) == 2 && p in pds)
    return splitting_of_mixed_prime_power(Lf, p, b; eiglat_cond, fix_root, genusDB, root_test, info_depth, discriminant_annihilator, _local)
  end

  x = gen(Hecke.Globals.Zx)
  cyclo_in_minpoly = _from_cyclotomic_polynomial_to_dict(minpoly(Lf))
  if min_poly === nothing
    min_polyf = minimal_polynomial(Lf)
    annihi_poly = one(parent(x))
    for k in keys(cyclo_in_minpoly)
      if !is_divisible_by(k, p)
        annihi_poly *= cyclotomic(k, x)*cyclotomic(p*k, x)
      else
        annihi_poly *= cyclotomic(p*k, x)
      end
    end
  else
    annihi_poly = change_base_ring(ZZ, min_poly; parent=parent(x))
  end

  # The isometries in output will have at least 3 prime divisors, so we need
  # to change the approach. Our current approach is quite naive:
  # we split the initial lattice with isometry (L, f) into its irreducible
  # eigenlattices, we split each of them by p, and we glue back everything
  # together.
  ds = sort!(collect(keys(cyclo_in_minpoly)))
  k = popfirst!(ds)
  N = kernel_lattice(Lf, k)

  chi = cyclotomic(k, x)
  t = reduced_resultant(chi, remove(annihi_poly, chi)[2])
  Ns = splitting_of_hermitian_type(N, p; eiglat_cond, fix_root, genusDB, root_test, check=false, info_depth, discriminant_annihilator=t*discriminant_annihilator, _local)
  isempty(Ns) && return Ns

  for k in ds
    ck = cyclotomic(k, x)
    chi *= ck
    M = kernel_lattice(Lf, k)
    t = reduced_resultant(ck, remove(annihi_poly, ck)[2])
    Ms = splitting_of_hermitian_type(M, p; eiglat_cond, fix_root, genusDB, root_test, check=false, info_depth, discriminant_annihilator=t*discriminant_annihilator, _local)
    is_empty(Ms) && return Ms

    tt = annihi_poly
    for (c, ) in factor(chi)
      tt = remove(tt, c)[2]
    end

    tc = reduced_resultant(chi, tt)
    Lq = kernel_lattice(Lf, chi)
    l = length(Ns)
    for i in 1:l
      N = popfirst!(Ns)
      for M in Ms
        ok, _Es = equivariant_primitive_extensions(N, M; form_over=TorQuadModule[first(discriminant_group(Lq))], _local)
        !ok && continue
        Es = first.(_Es)
        filter!(T -> is_of_type(T^p, type(Lq)), Es)
        filter!(Base.Fix2(is_annihilated_on_discriminant, tc*discriminant_annihilator), Es)
        append!(Ns, Es)
      end
    end
    isempty(Ns) && return Ns
  end
  (b == 1) && filter!(N -> order_of_isometry(N) == p*n, Ns)
  return Ns
end

@doc raw"""
    enumerate_classes_of_lattices_with_isometry(
      G::Union{ZZGenus, ZZLat},
      m::Int;
      char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
      min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
      rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
      eiglat_cond::Vector{Dict{Int64, Vector{Int64}}}=Dict{Int64, Vector{Int64}}[],
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      keep_partial_result::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_discriminant_annihilator(L),
      _local::Bool=false
    ) -> Vector{ZZLatWithIsom}

Given an even integer lattice $L$, return a complete set of representatives
for the isomorphism classes of lattices with isometry $(M, g)$ where $M$ is in
the genus of $L$, and $g$ is an isometry of $M$ of finite order $m$.
Alternatively, one can input a given genus symbol $G$ for even integer lattices
as an input --- the function first computes a representative of $G$.

For every $(M, g)$ in output, one may decide on the minimal polynomial (resp.
the characteristic polynomial) of $g$ by setting the keyword argument `min_poly`
(resp. `char_poly`) to the desired value. Moreover, one can also decide on the
rank, positive signature and negative signature of the eigenlattices of $(M, g)$
using the keyword arguments `rks`, `pos_sigs` and `neg_sigs` respectively. Each
list should consist of tuples $(k, i)$ where $k$ is a divisor of $p*m$ and $i$
is the value to be assigned for the given property (rank, positive/negative
signature) to the corresponding $\Phi_k$-eigenlattice. All these conditions
will be first condensed in a dictionary keeping track, for each divisor $k$
of $p*m$ of the potential rank `rk`, positive signature `pk` and negative
signature `nk` of the corresponding $\Phi_k$-eigenlattice. The keys of such
dictionary are the divisors $k$, and the corresponding value is the vector
`[rk, pk, nk]`. Any undetermined value will be set automatically to $-1$ by
default.
If one already knows such a dictionary, one can choose it as input under the
keyword argument `eiglat_cond`.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - by setting `keep_partial_result` to `true`, the function returns all
      pairs of lattices of isometries which have been computed;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.
    - `_local` -- if `true` return only a single genus representative 
      of each genus of lattices with isometry matching the given conditions.
      This keyword argument is not considered as part of the interface and 
      may change in future versions.

# Examples
```jldoctest
julia> r = enumerate_classes_of_lattices_with_isometry(root_lattice(:A, 3), 4; rks=[(1, 0)])
1-element Vector{ZZLatWithIsom}:
 Integer lattice with isometry of finite order 4

julia> rank(invariant_lattice(r[1]))
0
```
"""
enumerate_classes_of_lattices_with_isometry(::Union{ZZGenus, ZZLat}, ::Int)

function enumerate_classes_of_lattices_with_isometry(
    L::ZZLat,
    m::Int;
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing}=nothing,
    char_polys::Vector=[nothing],
    min_polys::Vector=[nothing],
    rks::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    pos_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    neg_sigs::Vector{Tuple{Int, Int}}=Tuple{Int, Int}[],
    eiglat_cond::Vector{Dict{Int64, Vector{Int64}}} = Dict{Int64, Vector{Int64}}[],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    keep_partial_result::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}}=_discriminant_annihilator(L),
    _local::Bool=false,
  )
  @req iseven(L) "Lattice must be even"
  @req is_finite(m) && m >= 1 "Order must be positive and finite"
  allow_info = get_verbosity_level(:ZZLatWithIsom) >= info_depth
  # If the user does not already input a dictionary of conditions on the
  # eigenlattices, we create one based on the conditions imposed on the
  # characteristic/minimal polynomials, on the ranks and on the signatures
  if char_poly !== nothing
    char_polys = [char_poly]
  end 
  if min_poly !==nothing
    min_polys=[min_poly]
  end
  if length(eiglat_cond) == 0
    eiglat_cond = [_conditions_from_input(m, char_p, min_p, rks, pos_sigs, neg_sigs) for char_p in char_polys for min_p in min_polys]
  end

  allow_info && println("Conditions computed") 

  if m == 1
    reps = ZZLatWithIsom[]
    for cond in unique!([get(i, 1, Int[-1, -1, -1]) for i in eiglat_cond])
      append!(reps, representatives_of_hermitian_type(L, 1, fix_root; cond, genusDB, root_test, info_depth, discriminant_annihilator, _local))
    end
    return reps
  end

  Lq = ZZLatWithIsom[integer_lattice_with_isometry(L)]
  if keep_partial_result
    out = deepcopy(Lq)
  end
  o = Int(1)
  pds = reverse!(sort!(prime_divisors(m)))
  for p in pds
    v = valuation(m, p)
    o *= p^v
    eco = _conditions_after_power(eiglat_cond, div(m, o))
    disc_ann = _discriminant_annihilator_after_power(discriminant_annihilator, div(m, o))
    Lq = splitting_by_prime_power!(Lq, p, v; eiglat_cond=eco, fix_root=gcd(o, fix_root), genusDB, root_test, info_depth, discriminant_annihilator=disc_ann, _local)
    if keep_partial_result
      append!(out, Lq)
    end
  end
  @hassert :ZZLatWithIsom 6 all(N -> order_of_isometry(N) == m, Lq)
  if keep_partial_result
    return out
  else
    return Lq
  end
end

function enumerate_classes_of_lattices_with_isometry(
    G::ZZGenus,
    n::Int;
    kwargs...,
  )
  return enumerate_classes_of_lattices_with_isometry(representative(G), n; kwargs...)
end

@doc raw"""
    splitting_by_prime_power!(
      Np::Vector{ZZLatWithIsom},
      p::Int,
      v::Int;
      eiglat_cond::Vector{Dict{Int, Vector{Int}}}=Dict{Int, Vector{Int}}(),
      fix_root::Int=-1,
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
      discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}, Nothing}=nothing
    ) -> Vector{ZZLatWithIsom}

Given a list `Np` of even lattices equipped with an isometry of finite order,
return a complete set of representatives for the isomorphism classes of
lattices with isometry $(M, g)$ such that the type of $(M, g^{p^v})$ is the
same as the type of some element $(N, h)$ of `Np` and the order of $g$ is
exactly $p^v*m$ where $m$ is the order of $h$.

!!! warning
    The function empties the list `Np` in input.

For every $(M, g)$ in output, one may decide on the rank `rM`, the positive
signature `pM` and the negative signature `nM` of the eigenlattices of $(M, g)$
using the keyword argument `eiglat_cond`. It should consist of a dictionary
where each key is a divisor of $p*m$, and the corresponding value is a tuple
`(rM, pM, nM)` of integers. Any undetermined value can be set to a negative
number; for instance $(-1, 2, -4)$ means that the associated eigenlattice
must have positive signature 2, without restriction on its rank and
negative signature.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    - by setting the value of `fix_root` to a certain integer $k$, the function
      only computes one generator for every conjugacy class of finite cyclic
      groups when computing the corresponding $\Phi_k$-kernel sublattices;
    - if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    - if `root_test` is set to true, the genus enumeration algorithm determines
      whether any new genus of negative definite lattices to be enumerated
      consists only of lattices of maximum $-2$. In such a case, the
      enumeration is skipped;
    - `discriminant_annihilator` -- an ideal ``I  \leq \mathbb{Z}[x]`` such 
      that the annihilator of $D_g$, for every lattices with isometry `(M,g)`
      in output, contains ``I``. If ``I`` is prinicipal, one can also give a
      generator as input.

"""
function splitting_by_prime_power!(
    Np::Vector{ZZLatWithIsom},
    p::Int,
    v::Int;
    eiglat_cond::Vector{Dict{Int,Vector{Int}}}=Dict{Int, Vector{Int}}[],
    fix_root::Int=-1,
    genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
    root_test::Bool=false,
    info_depth::Int=1,
    discriminant_annihilator::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}, Nothing}=nothing,
    _local::Bool=false,
  )
  @req is_prime(p) "p must be a prime number"
  @req all(N -> is_finite(order_of_isometry(N)), Np) "Isometries must be of finite order"
  @req all(N -> is_even(N), Np) "Lattices must be even"
  @req v >= 1 "v must be a positive integer"
  if discriminant_annihilator === nothing
    P, x = polynomial_ring(ZZ,[:x]; cached=false)
    discriminant_annihilator = ideal(P, [x[1]^(p^v)-1])
  end

  # Keep track of the types of the isometries we have already done (since two
  # isometries with the same type will give rise to similar outputs)
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
    Mp = splitting(M, p, 1; eiglat_cond=_conditions_after_power(eiglat_cond, q), fix_root=divexact(fix_root, gcd(fix_root, q)), genusDB, root_test, info_depth, discriminant_annihilator=_discriminant_annihilator_after_power(discriminant_annihilator, q), _local)
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

# squarefree part of the polynomial
function _radical(p::Union{ZZPolyRingElem, QQPolyRingElem})
  return prod(a[1] for a in factor(p))
end

# smallest ideal whose pth power is contained in I
function _discriminant_annihilator_after_power(
  I::MPolyIdeal{ZZMPolyRingElem},
  p::Hecke.IntegerUnion,
)
  P = base_ring(I)
  x = gen(P, 1)
  phi = hom(P, P, [x^p])
  return preimage(phi, I)
end

function _discriminant_annihilator_after_power(
  f::Union{ZZPolyRingElem, ZZMPolyRingElem},
  p::Hecke.IntegerUnion,
)
  P, (x, ) = polynomial_ring(ZZ, [:x]; cached=false)
  I = ideal(f(x))
  phi = hom(P, P, [x^p])
  g = first(gens(preimage(phi, I)))
  return first(gens(parent(f)))(g)
end


# If chi = prod_{n}\Phi_n^{e_n} where n's and e_n's are positive integers,
# return the dictionary that associated e_n to n.
function _from_cyclotomic_polynomial_to_dict(
    chi::Union{ZZPolyRingElem, QQPolyRingElem},
  )
  fac = factor(chi)
  V = Dict{Int, Int}()
  for (f, en) in fac
    ok, n = is_cyclotomic_polynomial_with_data(f)
    @req ok "Not a product of cyclotomic polynomials"
    V[n] = en
  end
  return V
end

# Create a dictionary of restrictions on rank and signatures
# for the eigenlattices of an isometry of order m with possible
# characteristic polynomial `char_poly` or minimal polynomial
# `min_poly`, and with respect to the conditions imposed by the
# lists `ranks`, `pos_sigs` and `neg_sigs`.
function _conditions_from_input(
    m::Int,
    char_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing},
    min_poly::Union{ZZPolyRingElem, QQPolyRingElem, Nothing},
    rks::Vector{Tuple{Int, Int}},
    pos_sigs::Vector{Tuple{Int, Int}},
    neg_sigs::Vector{Tuple{Int, Int}},
  )
  eiglat_cond = Dict{Int, Vector{Int}}()
  divs = divisors(m)
  # We have fixed a characteristic polynomial, so we know the
  # nontrivial irreducible eigenlattices and their rank
  if !isnothing(char_poly)
    V = _from_cyclotomic_polynomial_to_dict(char_poly)
    for (n, e) in V
      rn = euler_phi(n)*e
      jp = findfirst(a -> a[1] == n, pos_sigs)
      pn = isnothing(jp) ? -1 : pos_sigs[jp][2]
      jn = findfirst(a -> a[1] == n, neg_sigs)
      nn = isnothing(jn) ? -1 : neg_sigs[jn][2]
      eiglat_cond[n] = Int[rn, pn, nn]
    end
    # For the other eigenlattices, we know that the rank is 0 so we
    # can enforce it
    for n in setdiff(divs, keys(eiglat_cond))
      eiglat_cond[n] = Int[0, 0, 0]
    end
  # We do not know the characteristic polynomial but the minimal one,
  # so we can still decide which irreducible eigenlattices are nontrivial,
  # but we do not know, a priori, their rank
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
    for n in setdiff(divs, keys(eiglat_cond))
      eiglat_cond[n] = Int[0, 0, 0]
    end
  else
    for n in divs
      jr = findfirst(a -> a[1] == n, rks)
      rn = isnothing(jr) ? -1 : rks[jr][2]
      jp = findfirst(a -> a[1] == n, pos_sigs)
      pn = isnothing(jp) ? -1 : pos_sigs[jp][2]
      jn = findfirst(a -> a[1] == n, neg_sigs)
      nn = isnothing(jn) ? -1 : neg_sigs[jn][2]
      eiglat_cond[n] = Int[rn, pn, nn]
    end
  end
  return eiglat_cond
end

# Given a dictionary of restrictions on rank and signatures
# for the irreducible eigenlattices of an isometry f,
# return the corresponding restrictions on the irreducible
# eigenlattices of f^p.
function _conditions_after_power(
    eiglat_cond::Dict{Int, Vector{Int}},
    p::Int,
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

function _conditions_after_power(
    eiglat_cond::Vector{Dict{Int, Vector{Int}}},
    p::Int,
  )
  return unique!([_conditions_after_power(i, p) for i in eiglat_cond])
end 

###############################################################################
#
#  Enhanced genus enumeration
#
###############################################################################
function _packing_density_test(G::ZZGenus)
  # The center of density of a definite lattice L is defined as
  # delta(L) = rho^n/sqrt(|det(L)|) where rho(L) = 1//2*sqrt(|min(L)|).
  # If there exists a lattice L in a genus G with absolute minimum at least 4
  # then rho(L) >= 1, and therefore
  #
  #              center_density_ub >= delta(L) >= 1/sqrt(|det(G))
  # Hence, if |det(G)| < 1//(bn[rank(G)])^2, all the lattices in G have
  # absolute minimum equal to 2.
  !iszero(signature_tuple(G)[1]) && return false
  
  RR = real_field()
  # Spherepacking bounds taken from Henry cohn
  # https://dspace.mit.edu/bitstream/handle/1721.1/153311/table.pdf?sequence=8
  _packing_density_upper_bound = ["1", "0.9068996821171089", "0.7404804896930610", "0.6361073321551329", "0.5126451306253027", "0.4103032818801865", "0.3211471056675559", "0.2536695079010480", "0.1911204152968963", "0.1434100871082547", "0.1067252934567631", "0.0797117710668987", "0.0601644380983860", "0.0450612211935181", "0.0337564432797899", "0.0249944093845237", "0.0184640903350649", "0.0134853404450862", "0.0098179551395438", "0.0071270536033763", "0.0051596603948176", "0.0037259419689206", "0.0026842798864291", "0.0019295743094039", "0.0013841907222857", "0.0009910238892216", "0.0007082297958617", "0.0005052542161057", "0.0003598581852089", "0.0002559028743732", "0.0001817083813917", "0.0001288432887595", "0.0000912356039023", "0.0000645221967438", "0.0000455743843107", "0.0000321530553313", "0.0000226586900106", "0.0000159506499105", "0.0000112168687009", "0.0000078801051697", "0.0000055306464395", "0.0000038780970907", "0.0000027169074727", "0.0000019017703144", "0.0000013300905665", "0.0000009295151556", "0.0000006490757338", "0.0000004529067791"]
  packing_density_upper_bound = [RR(i) for i in _packing_density_upper_bound]
  center_density_upper_bounds = [RR(i)*RR(pi)^(-n//2)*(gamma(QQ(n)//2+1,RR))  for (n,i) in enumerate(packing_density_upper_bound)]
  r = Int(rank(G))
  center_density_ub = center_density_upper_bounds[r]
  if r <= 48 && abs(det(G)) < inv(center_density_ub)^2
    # strict inequality saves us from rounding errors
    # in the computation of gamma
    return true
  end
  return false
end 

# Legacy
_roger_upper_bound_test(G::ZZGenus) =_packing_density_test(G)

@doc raw"""
    oscar_genus_representatives(
      G::ZZGenus,
      algorithm::Symbol = :default;
      rand_neigh::Int = 10,
      invariant_function::Function=Hecke.default_invariant_function,
      save_partial::Bool=false,
      save_path::Union{IO, String, Nothing}=nothing,
      stop_after::IntExt=1000,
      max_lat::IntExt=inf
      genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
      root_test::Bool=false,
    ) -> Vector{ZZLat}

Return a complete list of representatives for the isometry classes in the genus
`G`.

!!! note
    This is a complement of the Hecke implementation for large genera of
    definite lattices. If `G` cannot be enumerated using the usual
    `line_orbits` algorithm, we enumerate it using random search in the
    corresponding neighbour graph. If after a certain number of vain
    iterations (1000 by default) the genus `G` is still not enumerated,
    we complete the enumeration by looking for lattices of isometries `(L, f)`
    so that `L` lies in `G` and `f` is of prime order. This is done using the
    algorithms of [BH23](@cite). The algorithm is recursive as it calls itself
    to enumerate smaller genera of definite lattices during the procedure.

The second input `algorithm` gives the choice to which algorithm to use for the
initial enumeration using neighbours. We currently support two algorithms:
  * `:random` which finds new isometry classes by constructing neighbours from
    random isotropic lines;
  * `:orbit` which computes orbits of isotropic lines before constructing
    neighbours.
If `algorithm = :default`, the function chooses the most appropriate algorithm
depending on the rank and determinant of the genus to be enumerated.

There are possible extra optional arguments:
  * `rand_neigh::Int` (default = `10`) -> for random enumeration, how many
    random neighbours are computed at each iteration;
  * `invariant_function::Function` (default = `default_invariant_function`) ->
    a function to compute isometry invariants in order to avoid unnecessary
    isometry tests;
  * `save_partial::Bool` (default = `false`) -> whether one wants to save
    iteratively new isometry classes;
  * `save_path::String` (default = `nothing`) -> a path to a folder where
    to save new lattices in the case where `save_partial` is true;
  * `stop_after::IntExt` (default = `1000`) -> the inital enumeration algorithm
    stops after the specified amount of vain iterations without finding a new
    isometry class is reached;
  * `max_lat::IntExt` (default = `inf`) -> the algorithm stops after finding
    `max` isometry classes.

!!! warning
    The algorithm uses the mass by default, in order to use the codes of
    [BH23](@cite). To enumerate `G` without the mass formula, please use
    the Hecke function `enumerate_definite_genus`.

If `save_partial = true`, the lattices are stored in a compact way in a `.txt`
file. The storing only remembers the rank of a lattice, half of its Gram matrix
(which is enough to reconstruct the lattice as a standalone object) and the
order of the isometry group of the lattice if it has been computed.

The `default_invariant_function` currently computes:
  * the absolute length of a shortest vector in the given lattice
    (also known as [`minimum`](@ref));
  * an ordered list of tuples consisting of the decomposition of the root
    sublattice of the given lattice (see [`root_lattice_recognition`](@ref));
  * the kissing number of the given lattice, which is proportional to the
    number of vectors of shortest length;
  * the order of the isometry group of the given lattice.

!!! note "For the advanced users"
    When using this function, one can use some extra keyword arguments
    which are carried along the computation:
    * if available, one can use any database of genera of definite lattices
      using the keyword argument `genusDB` (which should be a dictionary whose
      keys are genus symbols, and the corresponding value is a list of lattices
      of this genus);
    * if `root_test` is set to true, the algorithm determines whether the genus
      `G` consists only of negative definite lattices of maximum $-2$ (which
      can sometimes be predicted using sphere packing conditions). In such a
      case, the enumeration is skipped.
"""
function oscar_genus_representatives(
  G::ZZGenus,
  algorithm::Symbol = :default;
  rand_neigh::Int=10,
  invariant_function::Function=Hecke.default_invariant_function,
  save_partial::Bool=false,
  save_path::Union{IO, String, Nothing}=nothing,
  stop_after::IntExt=1000,
  max_lat::IntExt=inf,
  genusDB::Union{Nothing, Dict{ZZGenus, Vector{ZZLat}}}=nothing,
  root_test::Bool=false,
  info_depth::Int=1,
  _local::Bool=false
)
  if _local
    return [representative(G)]
  end
  allow_info = get_verbosity_level(:ZZLatWithIsom) >= info_depth
  # We do not need anything new, Hecke can handle this perfectly
  if !is_definite(G) || rank(G) <= 2
    allow_info && println("Indefinite genus or of small rank")
    return Hecke.representatives(G)
  end

  # Maybe the genus `G` is already known in the datatabse genusDB
  if !isnothing(genusDB)
    if haskey(genusDB, G)
      return deepcopy(genusDB[G])
    end
    G2 = rescale(G, -1)
    if haskey(genusDB, G2)
      return ZZLat[rescale(LL, -1) for LL in genusDB[G2]]
    end
  end
  r = rank(G)

  # Here we use a sphere packing condition as used in Section 2.4 of
  # "Symplectic rigidity of O'Grady's tenfolds" by L. Giovenzana, Grossi,
  # Onorati and Veniani.
  if root_test && _packing_density_test(G)
    return ZZLat[]
  end
  # Enumerate G using Hecke. If the rank and deteterminant of G are reasonable,
  # it will call `line_orbits` computations and the list l will be complete.
  # Otherwise, we proceed by random search, and as soon as we reach a point
  # where after `stop_after` vain iterations we do not find any new isometry
  # class, we stop Kneser's algorithm and we start isometry enumeration instead.
  allow_info && println("Definite genus of rank bigger than 2")
  l = enumerate_definite_genus(G, algorithm; rand_neigh, invariant_function, save_partial, save_path, stop_after, max=max_lat)
  length(l) == max_lat && return l

  # Part of the mass of G which is missing
  mm = mass(G) - sum(1//automorphism_group_order(LL) for LL in l; init=QQ(0))

  # If `mm` is nonzero, we are missing some isometry classes
  if !iszero(mm)
    allow_info && println("Need to enumerate isometries")
    # Recollect a dictionary of invariants, which should be fast to compute
    # about the lattices already known (to ease comparison of lattices)
    inv_lat = invariant_function(l[1])
    inv_dict = Dict{typeof(inv_lat), Vector{ZZLat}}(inv_lat => ZZLat[l[1]])
    for N in l[2:end]
      inv_lat = invariant_function(N)
      if haskey(inv_dict, inv_lat)
        push!(inv_dict[inv_lat], N)
      else
        inv_dict[inv_lat] = ZZLat[N]
      end
    end
    # Setup a default prime number for looking for certain isometries of
    # lattices in G not already computed
    Lf = integer_lattice_with_isometry(l[1])
    pos = is_positive_definite(Lf)
    Ps = reverse!(Int.(Hecke.primes_up_to(r+1)))
    D = Dict{Int, AbstractVector{Int}}(p => p == 2 ? collect(div(r, 2, RoundDown):-1:1) : collect(r.-reverse(p-1:p-1:r)) for p in Ps)
    # Looking for certain lattices with isometry
    while !iszero(mm)
      d = denominator(mm)
      if isone(d) # Very unlikely, but still
        i = 1
      else
        # Wants to minimize the rank of genera to enumerate
        # So we look, among the primes dividing d, for which
        # one we haven't yet computed isometries with very
        # small rank for the invariant part. If several primes
        # have the same of smallest value, we keep the largest
        # of those primes to minimize the rank of the hermitian
        # genus to enumerate on the other side. For now this
        # seems to be a good optimization of this part of the
        # function
        Pd = filter(i -> iszero(mod(d, Ps[i])), 1:length(Ps))
        @hassert :ZZLatWithIsom 3 !isempty(Pd)
        i = first(Pd)
        for j in Pd[2:end]
          if first(D[Ps[j]]) < first(D[Ps[i]])
            i = j
          end
        end
      end
      p = Ps[i]
      k = popfirst!(D[p])
      allow_info && println("(k, p) = $((k, p))")
      if isempty(D[p])
        popat!(Ps, i)
      end
      # Take care of how to distribute the signatures between invariant
      # and coinvariant sublattices
      if pos
        atp = admissible_triples(Lf, p; IrA=Int[k], IpA=Int[k], InA=Int[0], IrB=Int[r-k], IpB=Int[r-k], InB=Int[0], b=1)
      else
        atp = admissible_triples(Lf, p; IrA=Int[k], IpA=Int[0], InA=Int[k], IrB=Int[r-k], IpB=Int[0], InB=Int[r-k], b=1)
      end
      allow_info && println("$(length(atp)) admissible triples")
      for (A, B) in atp
        Bs = representatives_of_hermitian_type(B, p; genusDB, info_depth=info_depth+1)
        isempty(Bs) && continue
        As = representatives_of_hermitian_type(A, 1; genusDB, info_depth=info_depth+1)
        isempty(As) && continue
        for LA in As, LB in Bs
          Ns = admissible_equivariant_primitive_extensions(LA, LB, Lf, p; check=false, _local)
          allow_info &&  println("$(length(Ns)) lattices to try")
          for Nf in Ns
            flag = false
            N = lll(lattice(Nf))
            invN = invariant_function(N)
            # If no other known lattices have the same invariants as N
            # then N is not isometric to any of them and we have found
            # a new isometry class
            if !haskey(inv_dict, invN)
              flag = true
              inv_dict[invN] = ZZLat[N]
              push!(l, N)
              s = automorphism_group_order(N)
              if save_partial
                Hecke.save_lattice(N, save_path)
              end
              sub!(mm, mm, 1//s)
            # Otherwise we compare N with every other lattices with the same
            # invariant as N
            elseif all(M -> !is_isometric(N, M), inv_dict[invN])
              flag = true
              push!(inv_dict[invN], N)
              push!(l, N)
              s = automorphism_group_order(N)
              if save_partial
                Hecke.save_lattice(N, save_path)
              end
              sub!(mm, mm, 1//s)
            end
            length(l) == max_lat && return l
            is_zero(mm) && break
            if flag && allow_info
              perc = Float64(mm//mass(G)) * 100
              println("Lattices: $(length(l)), Target mass: $(mass(G)). missing: $(mm) ($(perc)%)")
            end
          end
        end
      end
    end
  end
  # If we have a lattice database, then we add the new genus there to be
  # used later... or to update the global database on
  # https://github.com/StevellM/DefLatDB
  if !isnothing(genusDB)
    genusDB[G] = deepcopy(l)
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
    n::Int,
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
    Np = splitting(N, p, 1)
    append!(Dn, Np)
  end
  D[n] = Dn
  return nothing
end

function _test_isometry_enumeration(
    L::ZZLat,
    k::Int = 2*rank(L)^2,
  )
  r = rank(L)
  ord = filter(m -> euler_phi(m) <= r, 2:k)
  sort!(ord)
  D = Dict{Int, Vector{ZZLatWithIsom}}()
  D[1] = ZZLatWithIsom[integer_lattice_with_isometry(L)]
  for n in ord
    _get_isometry!(D, n)
  end
  return D
end

###############################################################################
#
#  Discriminant annihilators
#
###############################################################################

@doc raw"""
    is_annihilated(
      f::AutomorphismGroupElem{TorQuadModule},
      p::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}},
    ) -> Bool

Return whether the isometry `f` is annihilated by `p`, where `p` is
either a polynomial in $\mathbb{Z}[x]$ or an ideal of such.
"""
is_annhilated(::AutomorphismGroupElem{TorQuadModule}, ::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}})

function is_annihilated(
  f::AutomorphismGroupElem{TorQuadModule},
  p::ZZMPolyRingElem,
)
  P, y = polynomial_ring(ZZ, :x; cached=false)
  phi = hom(parent(p), P,  [y])
  pZ = phi(p)
  return iszero(pZ(abelian_group_homomorphism(hom(f))))
end

function is_annihilated(
  f::AutomorphismGroupElem{TorQuadModule},
  p::ZZPolyRingElem,
)
  return iszero(p(abelian_group_homomorphism(hom(f))))
end
function is_annihilated(
  f::AutomorphismGroupElem{TorQuadModule},
  I::MPolyIdeal{ZZMPolyRingElem}
)
  return all(is_annihilated(f, p) for p in gens(I))
end

@doc raw"""
    is_annihilated_on_discriminant(
      Lf::ZZLatWithIsom,
      p::Union{ZZPolyRingElem, ZZMPolyRingElem, MPolyIdeal{ZZMPolyRingElem}},
    ) -> Bool

Given a lattice with isometry `(L, f)`,  whether the isometry $D_f$ induced
by `f` on the discriminant group of `L`  is annihilated by `p`, where
`p` is either a polynomial in $\mathbb{Z}[x]$ or an ideal of such.
"""
is_annhilated_on_discriminant

function is_annihilated_on_discriminant(
  Lf::ZZLatWithIsom,
  I::MPolyIdeal,
)
  _, Df= discriminant_group(Lf)
  return is_annihilated(Df, I)
end

function is_annihilated_on_discriminant(
  Lf::ZZLatWithIsom,
  p::Union{ZZPolyRingElem, ZZMPolyRingElem},
)
  _, Df= discriminant_group(Lf)
  return is_annihilated(Df, p)
end

_default_discriminant_annihilator(Lf::ZZLatWithIsom) = _annihilator(discriminant_group(Lf)[1])

_discriminant_annihilator(L::Union{ZZLat, ZZGenus}; parent=polynomial_ring(ZZ,[:x]; cached=false)[1]) = _annihilator(discriminant_group(L); parent)

_discriminant_annihilator(Lf::ZZLatWithIsom; parent=polynomial_ring(ZZ,[:x]; cached=false)[1]) = _annihilator(discriminant_group(Lf)[2]; parent)

function _annihilator(
  T::TorQuadModule;
  parent=polynomial_ring(ZZ, [:x]; cached=false)[1],
)
  P = parent
  d = exponent(T)
  return ideal(P, [P(d)])
end

raw"""
  _annihilator(
  f::AutomorphismGroupElem{TorQuadModule}; parent=polynomial_ring(ZZ, [:x]; cached=false)[1])

Return some ideal ``I \subseteq \ZZ[x]`` such that 
``g(f)=0`` for all ``g \in I``. 
  
!warning
  ``I`` is merely contained in the annihilator but may not be equal to it.
"""
function _annihilator(
  f::AutomorphismGroupElem{TorQuadModule};
  parent=polynomial_ring(ZZ, [:x]; cached=false)[1],
)
  D = domain(f)
  J = _annihilator(D; parent)
  P = base_ring(J)
  m = minpoly(Hecke.Globals.Zx, matrix(f))
  x = gen(parent,1)
  m = m(x)
  J = J + ideal(parent, m) + ideal(parent, x^order(f) - 1)
  # a bound for the discriminant annihilator .... but I don't think that it is everything.
  return J
end

function _divisors(d::T) where T <: RingElem
  F = factor(d)
  B = collect(F)
  D = [0:j for (_, j) in B]

  res = T[]
  for b in Iterators.product(D...)
    J = prod(B[i][1]^b[i] for i in 1:length(b);init=one(d))
    push!(res, J)
  end
  return res
end

function *(x::Hecke.IntegerUnion, I::MPolyIdeal{ZZMPolyRingElem})
  return ideal(base_ring(I), [x*i for i in gens(I)])
end

######################################################################################
#
# Legacy interface ... to be deprecated at some point
#
######################################################################################
splitting_by_prime_power!(args...;eiglat_cond::Dict{Int,Vector{Int}}, kwargs...) = splitting_by_prime_power!(args...;eiglat_cond=[eiglat_cond], kwargs...)

splitting_of_hermitian_type(args...; eiglat_cond::Dict{Int, Vector{Int}}, kwargs...) = splitting_of_hermitian_type(args...;eiglat_cond=[eiglat_cond], kwargs...)
  
splitting_of_prime_power(args...;eiglat_cond::Dict{Int,Vector{Int}}, kwargs...) = splitting_of_prime_power(args...;eiglat_cond=[eiglat_cond], kwargs...)

splitting_of_pure_mixed_prime_power(args...;eiglat_cond::Dict{Int,Vector{Int}}, kwargs...) = splitting_of_pure_mixed_prime_power(args...;eiglat_cond=[eiglat_cond], kwargs...)

splitting_of_mixed_prime_power(args...;eiglat_cond::Dict{Int,Vector{Int}}, kwargs...) = splitting_of_mixed_prime_power(args...;eiglat_cond=[eiglat_cond], kwargs...)

splitting(args...;eiglat_cond::Dict{Int,Vector{Int}}, kwargs...) = splitting(args...;eiglat_cond=[eiglat_cond], kwargs...)

enumerate_classes_of_lattices_with_isometry(args...;eiglat_cond::Dict{Int,Vector{Int}}, kwargs...) = enumerate_classes_of_lattices_with_isometry(args...;eiglat_cond=[eiglat_cond], kwargs...)
