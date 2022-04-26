########################################################################
#
# Localizations of finitely generated modules over multivariate 
# polynomial rings.
#
# This file implements the required functionality for the 
# Posur-interface; see [1].
#
# [1] Posur: Linear systems over localizations of rings, arXiv:1709.08180v2


# compute the syzygies of a matrix
function syz(A::MatrixType) where {T<:MPolyElem, MatrixType<:MatrixElem{T}}
  R = base_ring(A)
  m = nrows(A)
  n = ncols(A)
  F = FreeMod(R, m)
  G = FreeMod(R, n)
  f = hom(F, G, A)
  K, inc = kernel(f)
  return matrix(inc)
end

function has_nonepmty_intersection(U::MPolyPowersOfElement, I::MPolyIdeal)
  R = ambient_ring(U)
  R == base_ring(I) || error("the multiplicative set and the ideal must be defined over the same ring")

  d = prod(denominators(U))
  inradical(d, I) || return false, zero(R), zero(MatrixSpace(R, 1, ngens(I)))
  (k, f) = Oscar._minimal_power_such_that(d, (x->x in I))
  return true, f, coordinates(f, I)
end

function has_nonepmty_intersection(U::MPolyComplementOfPrimeIdeal, I::MPolyIdeal)
  R = ambient_ring(U)
  R == base_ring(I) || error("the multiplicative set and the ideal must be defined over the same ring")
  P = prime_ideal(U)
  candidates = [(f, i) for (f, i) in zip(gens(I), 1:ngens(I)) if !(f in P)]
  length(candidates) == 0 && return false, zero(R), zero(MatrixSpace(R, 1, ngens(I)))
  d = maximum([total_degree(f) for (f, i) in candidates])
  (g, j) = candidates[1]
  for (h, k) in candidates
    if total_degree(h) < total_degree(g) 
      (g, j) = (h, k)
    end
  end
  A = zero(MatrixSpace(R, 1, ngens(I)))
  A[1, j] = 1
  return true, g, A
end

function _minimal_power_such_that(f::RingElemType, P::PropertyType) where {RingElemType<:RingElem, PropertyType}
  P(one(parent(f))) && return (0, one(f))
  P(f) && return (1, f)
  f_powers = [(1,f)]

  while !P(last(f_powers)[2])
    push!(f_powers, (last(f_powers)[1]*2, last(f_powers)[2]^2))
  end
  upper = pop!(f_powers)
  lower = pop!(f_powers)
  while upper[1]!=lower[1]+1
    middle = pop!(f_powers)
    middle = (lower[1]+middle[1], lower[2]*middle[2])
    if P(middle[2])
      upper = middle
    else
      lower = middle
    end
  end
  return upper
end

