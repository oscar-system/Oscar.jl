################################################################################
# Quantum analogs
#
# Copyright (C) 2020 Ulrich Thiel, ulthiel.com/math
#
# Mar 2023: Imported into OSCAR by UT
#
# Future features could include:
# - Symmetric quantum numbers
# - Two-colored quantum numbers
################################################################################

export quantum_integer, quantum_factorial, quantum_binomial

################################################################################
# Quantum integers
################################################################################
@doc raw"""
    quantum_integer(n::IntegerUnion, q::RingElement)
    quantum_integer(n::IntegerUnion)

Return the quantum integer $[n]_q$ which is defined as $\frac{q^n-1}{q-1}$
when ``q-1`` is invertible.

For general ring elements `q`, we use the following identities to compute
$[n]_q$: if `n` is non-negative, then $[n]_q = \sum_{i=0}^{n-1} q^i$. To
handle negative values `n` we use the identity $[n]_q = -q^{n} [-n]_q$. Thus for
negative `n` we require `q` to be invertible.

Note that for ``q=1`` we obtain $[n]_1 = n$ hence the quantum integers are
"deformations" of the usual integers. For details about these objects see
[Con00](@cite) or [KC02](@cite).

If `q` is omitted then it defaults to the generator of a Laurent polynomial
ring over the integers.

# Examples
```jldoctest
julia> quantum_integer(3)
q^2 + q + 1

julia> quantum_integer(-3)
-q^-1 - q^-2 - q^-3

julia> quantum_integer(3,2)
7

julia> quantum_integer(-3,2)
ERROR: DomainError with -3:
Cannot raise an integer x to a negative power -3.

julia> quantum_integer(-3,2//1)
-7//8

julia> K,i = cyclotomic_field(4, "i");

julia> quantum_integer(3, i)
i
```
"""
function quantum_integer(n::IntegerUnion, q::RingElement)
  R = parent(q)
  isone(q) && return R(n)
  n < 0 && return -q^n * quantum_integer(-n, q)

  z = zero(R)
  qi = one(R)
  for i = 0:n-1
    # at this point qi = q^i
    z = add!(z, qi)
    qi = mul!(qi, q)
  end
  return z
end

function quantum_integer(n::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return quantum_integer(n,q)
end


################################################################################
# Quantum factorials
################################################################################
@doc raw"""
    quantum_factorial(n::IntegerUnion, q::RingElement)
    quantum_factorial(n::IntegerUnion)

Return the quantum factorial $[n]_q!$ for a non-negative integer `n` and an
element ``q`` of a ring ``R`` which is defined as $[1]_q \cdots [n]_q$.

Note that for ``q=1`` we obtain $[n]_1! = n!$ hence the quantum factorial is a
"deformation" of the usual factorial. For details about these objects see
[Con00](@cite) or [KC02](@cite).

# Examples
```jldoctest
julia> quantum_factorial(3)
q^3 + 2*q^2 + 2*q + 1

julia> quantum_factorial(3,2)
21

julia> K,i = cyclotomic_field(4, "i");

julia> quantum_factorial(3, i)
i - 1
```
"""
function quantum_factorial(n::IntegerUnion, q::RingElement)
  @req n >= 0 "n >= 0 required"

  R = parent(q)
  isone(q) && return R(factorial(n))

  z = one(R)
  for i = 1:n
    z *= quantum_integer(i,q)
  end
  return z
end

function quantum_factorial(n::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return quantum_factorial(n,q)
end


################################################################################
# Quantum binomials
################################################################################

@doc raw"""
    quantum_binomial(n::IntegerUnion, k::IntegerUnion, q::RingElement)
    quantum_binomial(n::IntegerUnion, k::IntegerUnion)

Return the quantum binomial $\binom{n}{k}_q$ for an integer `n` and a
non-negative integer `k` which is defined as
```math
\binom{n}{k}_q ≔ \frac{[n]_q!}{[k]_q! [n-k]_q!} = \frac{[n]_q [n-1]_q \cdots [n-k+1]_q}{[k]_q!}
```

Note that the first expression is only defined for ``n ≥ k`` since the quantum factorials
are only defined for non-negative integers, but the second expression is well-defined for
all integers `n` and is used for the implementation.

Note that for ``q=1`` we obtain $\binom{n}{k}_1 = \binom{n}{k}$ hence the
quantum binomial coefficient is a "deformation" of the usual binomial
coefficient. For details about these objects see [Con00](@cite) or
[KC02](@cite).

# Examples
```jldoctest
julia> quantum_binomial(4,2)
q^4 + q^3 + 2*q^2 + q + 1

julia> quantum_binomial(19,5,-1)
36

julia> K,i = cyclotomic_field(4);

julia> quantum_binomial(17,10,i)
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
function quantum_binomial(n::IntegerUnion, k::IntegerUnion, q::RingElement)
  @req k >= 0 "k >= 0 required"

  R = parent(q)
  isone(q) && return R(binomial(n,k))
  k == 0 && return one(R)
  k == 1 && return quantum_integer(n,q)
  n < 0 && return (-1)^k * q^(div(-k*(k-1),2) + k*n) * quantum_binomial(k-n-1,k,q)
  n < k && return zero(R)

  z = zero(R)
  for i = 0:n-k
    z += q^i * quantum_binomial(i+k-1,k-1,q)
  end
  return z
end

function quantum_binomial(n::IntegerUnion, k::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return quantum_binomial(n,k,q)
end
