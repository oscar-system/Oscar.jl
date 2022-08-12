export admissible_triples, is_admissble_triple

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
function _find_D(d::fmpz, m::Integer, p::Integer)
  @req is_prime(p) "p must be a prime number"
  @req d != 0 "The discriminant of a non-degenerate form is non-zero"

  # If m == 0, there are no conditions on the gcd of d1 and dp
  if m == 0
    return _tuples_divisors(d)
  end
  
  D = Tuple{Int64,Int64}[]
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
function _find_L(r::Integer, d, s::fmpz, l::fmpz, p::Integer, even = true)
  L = ZGenus[]
  for (s1,s2) in [(s,t) for s=0:r for t=0:r if s+t==r]
    gen = genera((s1,s2), d, even=even)
    filter!(G -> divides(scale(G), s)[1], gen)
    filter!(G -> divides(p*l, level(G))[1], gen)
    append!(L, gen)
  end
  return L
end

primes(G::ZGenus) = prime.(local_symbols(G))

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
  if !divides(det(C), p)[1]
    qA = gram_matrix_quadratic(normal_form(primary_part(discriminant_group(A), p)[1])[1])
    qB = gram_matrix_quadratic(normal_form(rescale(primary_part(discriminant_group(B), p)[1], -1))[1])
    return qA == qB
  end

  l = valuation(level(C), p)
  Ap = symbol(local_symbol(A, p))
  Bp = symbol(local_symbol(B,p))
  a_max = sum(Int64[s[2] for s in Ap if s[1] == l+1])
  b_max = sum(Int64[s[2] for s in Bp if s[1] == l+1])
  # For the glueing, rho_{l+1}(A) and rho_{l+1}(B) are anti-isometric, so they must have the
  # same order
  if a_max != b_max
    return false
  end
  
  # Since p^l*A^\vee/A is in the glue, its order is less than the order of the glue
  if g < a_max
    return false
  end

  a1 = sum(Int64[s[2] for s in Ap if s[1] == 1])
  a2 = sum(Int64[s[2] for s in Ap if s[1] >= 2])
  b1 = sum(Int64[s[2] for s in Bp if s[1] == 1])
  b2 = sum(Int64[s[2] for s in Bp if s[1] >= 2])

  if l != 0
    ker_min = max(g-a1, g-b1, a_max)
  else
    ker_min = max(g-a1, g-b1)
  end
  ker_max = min(a2+div(a1,2), b2+div(b1,2), g)

  if ker_max < ker_min
    return false
  end

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
    s1 = symbol(Ap)[s3[1]+1]
    s2 = symbol(Bp)[s3[1]+1]
    if s1[3] != s2[3]
      return false
    end
  end

  for s in Cp
    s[1] += 2
  end
  Cp = ZpGenus(p, Cp)
  
  if !represented(local_symbol(AperpB, p), Cp)
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

admissible_triples(L::ZLat, p::Integer) = admissible_triple(genus(L), p)


function primitive_extensions(A::ZLat, B::ZLat, C::ZLat, p::Integer)
  @req is_admissible_triple(A, B, C, p) "(A, B, C) must be p-admissble"
  L = []
  g = valuation(divexact(det(A)*det(B), det(C)), p)
  

end
