################################################################################
# q-analogs
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Mar 2023: Imported into OSCAR by UT
# Apr 2025: Cleaned up and finalized for OSCAR by MH
################################################################################

export q_integer, q_factorial, q_binomial

################################################################################
# Quantum integers
################################################################################
@doc raw"""
    q_integer(n::IntegerUnion, q::RingElement)
    q_integer(n::IntegerUnion)

Return the ``q``-integer $[n]_q$ which is defined as $\frac{q^n-1}{q-1}$
when ``q-1`` is invertible.

For general ring elements `q`, we use the following identities to compute
$[n]_q$: if `n` is non-negative, then $[n]_q = \sum_{i=0}^{n-1} q^i$. To
handle negative values `n` we use the identity $[n]_q = -q^{n} [-n]_q$. Thus for
negative `n` we require `q` to be invertible.

Note that for ``q=1`` we obtain $[n]_1 = n$ hence the ``q``-integers are
"deformations" of the usual integers. For details about these objects see
[Con00](@cite) or [KC02](@cite).

If `q` is omitted then it defaults to the generator of a Laurent polynomial
ring over the integers.

# Examples
```jldoctest
julia> q_integer(3)
q^2 + q + 1

julia> q_integer(-3)
-q^-1 - q^-2 - q^-3

julia> q_integer(3,2)
7

julia> q_integer(-3,2)
ERROR: DomainError with -3:
Cannot raise an integer x to a negative power -3.
[...]

julia> q_integer(-3,2//1)
-7//8

julia> K,i = cyclotomic_field(4, "i");

julia> q_integer(3, i)
i
```
"""
function q_integer(n::IntegerUnion, q::RingElement)
  R = parent(q)
  isone(q) && return R(n)
  n < 0 && return -q^n * q_integer(-n, q)

  z = zero(R)
  qi = one(R)
  for i = 0:n-1
    # at this point qi = q^i
    z = add!(z, qi)
    qi = mul!(qi, q)
  end
  return z
end

function q_integer(n::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  n >= 0 && return _laurent_poly(R, [one(ZZ) for _ in 1:n], 0)
  return _laurent_poly(R, [-one(ZZ) for _ in 1:-n], Int(n))
end

# Return the sum of `c[j] * q^(e + j - 1)`, with linear cost in `length(c)`;
# a sum of monomials is quadratic because `+` reallocates the dense storage.
function _laurent_poly(R::LaurentPolyRing, c::Vector{ZZRingElem}, e::Int)
  z = zero(R)

  # ascending degrees: the lowest one fixes the offset, the others only extend
  for (j, a) in enumerate(c)
    iszero(a) && continue
    z = set_coefficient!(z, e + j - 1, a)
  end
  return z
end


################################################################################
# Quantum factorials
################################################################################
@doc raw"""
    q_factorial(n::IntegerUnion, q::RingElement)
    q_factorial(n::IntegerUnion)

Return the ``q``-factorial $[n]_q!$ for a non-negative integer `n` and an
element ``q`` of a ring ``R`` which is defined as $[1]_q \cdots [n]_q$.

Note that for ``q=1`` we obtain $[n]_1! = n!$ hence the ``q``-factorial is a
"deformation" of the usual factorial. For details about these objects see
[Con00](@cite) or [KC02](@cite).

# Examples
```jldoctest
julia> q_factorial(3)
q^3 + 2*q^2 + 2*q + 1

julia> q_factorial(3,2)
21

julia> K,i = cyclotomic_field(4, "i");

julia> q_factorial(3, i)
i - 1
```
"""
function q_factorial(n::IntegerUnion, q::RingElement)
  @req n >= 0 "n >= 0 required"

  R = parent(q)
  isone(q) && return R(factorial(n))

  z = one(R)
  qi = one(R)
  qint = zero(R)
  for _ in 1:n
    # [i]_q = [i-1]_q + q^(i-1)
    qint = add!(qint, qi)
    z = mul!(z, z, qint)
    qi = mul!(qi, qi, q)
  end
  return z
end

function q_factorial(n::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return q_factorial(n,q)
end


################################################################################
# Quantum binomials
################################################################################

@doc raw"""
    q_binomial(n::IntegerUnion, k::IntegerUnion, q::RingElement)
    q_binomial(n::IntegerUnion, k::IntegerUnion)

Return the ``q``-binomial $\binom{n}{k}_q$ for an integer `n` and a
non-negative integer `k` which is defined as
```math
\binom{n}{k}_q ≔ \frac{[n]_q!}{[k]_q! [n-k]_q!} = \frac{[n]_q [n-1]_q \cdots [n-k+1]_q}{[k]_q!}
```

Note that the first expression is only defined for ``n`` greater or equal to
``k`` since the ``q``-factorials are only defined for non-negative integers,
but the second expression is well-defined for all integers `n` and is used for
the implementation.

Note that for ``q=1`` we obtain $\binom{n}{k}_1 = \binom{n}{k}$ hence the
``q``-binomial coefficient is a "deformation" of the usual binomial
coefficient. For details about these objects see [Con00](@cite) or
[KC02](@cite).

# Examples
```jldoctest
julia> q_binomial(4,2)
q^4 + q^3 + 2*q^2 + q + 1

julia> q_binomial(19,5,-1)
36

julia> K,i = cyclotomic_field(4);

julia> q_binomial(17,10,i)
0
```

# Extended help

In [Con00](@cite) it is shown that
```math
\binom{n}{k}_q = \sum_{i=0}^{n-k} q^i \binom{i+k-1}{k-1}_q
```
if $n ≥ k > 0$. Since $\binom{n}{0}_q = 1$ for all `n` and $\binom{n}{k}_q = 0$
if $0 ≤ n < k$ it follows inductively that
$\binom{n}{k}_q ∈ ℤ[q]$ if $n ≥ 0$.

For all ``n ∈ ℤ`` we have the relation
```math
\binom{n}{k}_q = (-1)^k q^{-k(k-1)/2+kn} \binom{k-n-1}{k}_q
```
which shows that $\binom{n}{k}_q ∈ ℤ[q^{-1}]$ if $n < 0$. In particular,
$\binom{n}{k}_q ∈ ℤ[q,q^{-1}]$ for all `n`.
Now, for an element ``q`` of a ring ``R`` we define ``\binom{n}{k}_q`` as the
specialization of ``\binom{n}{k}_q`` in ``q``, where ``q`` is
assumed to be invertible in ``R`` if ``n < 0``.
"""
function q_binomial(n::IntegerUnion, k::IntegerUnion, q::RingElement)
  @req k >= 0 "k >= 0 required"

  R = parent(q)
  isone(q) && return R(binomial(n,k))
  k == 0 && return one(R)
  0 <= n < k && return zero(R)

  n, k = Int(n), Int(k)
  n >= 0 && return _q_binomial(n, k, q)
  z = q^_q_binomial_negative_shift(n, k) * _q_binomial(k-n-1, k, q)
  return isodd(k) ? neg!(z) : z
end

function q_binomial(n::IntegerUnion, k::IntegerUnion)
  @req k >= 0 "k >= 0 required"

  R,q = laurent_polynomial_ring(ZZ, "q")
  k == 0 && return one(R)
  0 <= n < k && return zero(R)

  n, k = Int(n), Int(k)
  n >= 0 && return _laurent_poly(R, _q_binomial_coefficients(n, k), 0)
  c = _q_binomial_coefficients(k-n-1, k)
  if isodd(k)
    foreach(neg!, c)
  end
  return _laurent_poly(R, c, _q_binomial_negative_shift(n, k))
end

# exponent `e` in the relation `[n, k]_q = (-1)^k q^e [k-n-1, k]_q`
_q_binomial_negative_shift(n::Int, k::Int) = k*n - div(k*(k-1), 2)

# `[n, k]_q` for `0 <= k <= n` via the q-Pascal rule on `G(a, b) = [a+b, a]_q`:
#
#   G(a, b) = G(a-1, b) + q^a G(a, b-1),   G(a, 0) = G(0, b) = 1
#
# This costs `k(n-k)` multiplications and additions, without divisions, so it
# works over every commutative ring.
function _q_binomial(n::Int, k::Int, q::RingElement)
  R = parent(q)
  k = min(k, n - k)
  qpow = [q^a for a in 1:k]

  # c[a+1] holds G(a, b) for the current b
  c = [one(R) for _ in 0:k]
  for _ in 1:n-k
    for a in 1:k
      c[a+1] = mul!(c[a+1], c[a+1], qpow[a])
      c[a+1] = add!(c[a+1], c[a])
    end
  end
  return c[k+1]
end

# Coefficients of `[n, k]_q` in ascending degree for `0 <= k <= n`, from
#
#   [m, i]_q = [m-1, i-1]_q (1 - q^m) / (1 - q^i),   m = n-k+i
#
# in `O(k^2 (n-k))` integer additions, a factor `n-k` less than `_q_binomial`
# with polynomial arithmetic.
function _q_binomial_coefficients(n::Int, k::Int)
  k = min(k, n - k)
  deg = k*(n-k)
  c = [zero(ZZ) for _ in 0:deg+n]
  c[1] = one(ZZ)

  # c[j+1] is the coefficient of q^j; d the current degree
  d = 0
  for i in 1:k
    m = n - k + i
    for j in d+m:-1:m
      c[j+1] = sub!(c[j+1], c[j-m+1])
    end
    d += m - i

    # exact division: P = (1 - q^i) Q gives Q_j = P_j + Q_(j-i); what remains
    # above degree d is the negated top of Q
    for j in i:d
      c[j+1] = add!(c[j+1], c[j-i+1])
    end
    for j in d+1:d+i
      c[j+1] = zero!(c[j+1])
    end
  end
  return resize!(c, deg + 1)
end
