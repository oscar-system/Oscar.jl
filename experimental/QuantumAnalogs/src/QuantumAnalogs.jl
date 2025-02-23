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
    quantum_integer(n::IntegerUnion, q::RingElem)
    quantum_integer(n::IntegerUnion, q::Integer)
    quantum_integer(n::IntegerUnion)

Let ``n ∈ ℤ`` and let ``ℚ(q)`` be the fraction field of the polynomial ring ``ℤ[q]`` in
one variable ``q``. The **quantum integer** ``[n]_q ∈ ℚ(q)`` of ``n`` is defined as
```math
[n]_q ≔ \frac{q^n-1}{q-1} \;.
```
We have
```math
[n]_q = \sum_{i=0}^{n-1} q^i ∈ ℤ[q] \quad \text{if } n ≥ 0
```
and
```math
[n]_q = -q^{n} [-n]_q \quad \text{for any } n ∈ ℤ \;,
```
hence
```math
[n]_q = - \sum_{i=0}^{-n-1} q^{n+i} ∈ ℤ[q^{-1}] \quad \text{ if } n < 0 \;.
```
This shows in particular that actually
```math
[n]_q ∈ ℤ[q,q^{-1}] ⊂ ℚ(q) \quad \text{ for any } n ∈ ℤ \;.
```
Now, for an element ``q`` of a ring ``R`` we define ``[n]_q ∈ R`` as the specialization of
``[n]_q`` in ``q`` using the two equations above—assuming that ``q`` is invertible in ``R``
if ``n<0``. Note that for ``q=1`` we obtain
```math
[n]_1 = n \quad \text{for any } n ∈ ℤ \;,
```
so the quantum integers are "deformations" of the usual integers.

# Functions
* `quantum_integer(n::IntegerUnion,q::RingElem)` returns ``[n]_q`` as an element of ``R``,
  where ``R`` is the parent ring of ``q``.
* `quantum_integer(n::IntegerUnion,q::Integer)` returns ``[n]_q``. Here, if ``n >= 0`` or
  ``q = ± 1``, then ``q`` is considered as an element of ``ℤ``, otherwise it is taken as an
  element of ``ℚ``.
* `quantum_integer(n::IntegerUnion)` returns ``[n]_q`` as an element of ``ℤ[q^{-1}]``.

# Examples
```jldoctest
julia> quantum_integer(3)
q^2 + q + 1

julia> quantum_integer(-3)
-q^-1 - q^-2 - q^-3

julia> quantum_integer(3,2)
7

julia> quantum_integer(-3,2)
-7//8

julia> K,i = cyclotomic_field(4, "i");

julia> quantum_integer(3, i)
i
```

# References
1. [Con00](@cite)
2. [KC02](@cite)
"""
function quantum_integer(n::IntegerUnion, q::RingElem)
  R = parent(q)
  isone(q) && return R(n)

  z = zero(R)
  if n >= 0
    for i = 0:n-1
      z += q^i
    end
  else
    for i = 0:-n-1
      z -= q^(n+i)
    end
  end
  return z
end

function quantum_integer(n::IntegerUnion, q::Integer)
  # this method is not type stable, but that is acceptable as it is only
  # intended for convenience during interactive use
  if n >= 0 || q == 1 || q == -1
    return quantum_integer(n,ZZ(q))
  else
    return quantum_integer(n,QQ(q))
  end
end

function quantum_integer(n::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return quantum_integer(n,q)
end


################################################################################
# Quantum factorials
################################################################################
@doc raw"""
    quantum_factorial(n::IntegerUnion, q::RingElem)
    quantum_factorial(n::IntegerUnion, q::Integer)
    quantum_factorial(n::IntegerUnion)

For a non-negative integer ``n`` and an element ``q`` of a ring ``R`` the **quantum
factorial** ``[n]_q! ∈ R`` is defined as
```math
[n]_q! ≔ [1]_q ⋅ … ⋅ [n]_q ∈ R \;.
```
Note that for ``q=1`` we obtain
```math
[n]_1! = n! \quad \text{ for all } n ∈ ℤ \;,
```
hence the quantum factorial is a "deformation" of the usual factorial.

# Functions

The functions `quantum_factorial` work analogously to [`quantum_integer`](@ref).

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
function quantum_factorial(n::IntegerUnion, q::RingElem)
  @req n >= 0 "n >= 0 required"

  R = parent(q)
  isone(q) && return R(factorial(n))

  z = one(R)
  for i = 1:n
    z *= quantum_integer(i,q)
  end
  return z
end

function quantum_factorial(n::IntegerUnion, q::Integer)
  return quantum_factorial(n,ZZ(q))
end

function quantum_factorial(n::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return quantum_factorial(n,q)
end


################################################################################
# Quantum binomials
################################################################################

@doc raw"""
    quantum_binomial(n::IntegerUnion, q::RingElem)
    quantum_binomial(n::IntegerUnion, k::IntegerUnion, q::Integer)
    quantum_binomial(n::IntegerUnion, k::IntegerUnion)

Let ``k`` be a non-negative integer and let ``n ∈ ℤ``. The **quantum binomial**
``\binom{n}{k}_q \in ℚ(q)`` is defined as
```math
\binom{n}{k}_q ≔ \frac{[n]_q!}{[k]_q! [n-k]_q!} = \frac{[n]_q [n-1]_q⋅ … ⋅ [n-k+1]_q}{[k]_q!}
```
Note that the first expression is only defined for ``n ≥ k`` since the quantum factorials
are only defined for non-negative integers, but the second expression is well-defined for
all ``n ∈ ℤ`` and is used for the definition. In [Con00](@cite) it is shown that
```math
\binom{n}{k}_q = \sum_{i=0}^{n-k} q^i \binom{i+k-1}{k-1}_q \quad \text{if } n ≥ k > 0 \;.
```
Since
```math
\binom{n}{0}_q = 1 \quad \text{for all } n ∈ ℤ
```
and
```math
\binom{n}{k}_q = 0 \quad \text{if } 0 ≤ n < k \;,
```
it follows inductively that
```math
\binom{n}{k}_q ∈ ℤ[q] \quad \text{if } n ≥ 0 \;.
```
For all ``n ∈ ℤ`` we have the relation
```math
\binom{n}{k}_q = (-1)^k q^{-k(k-1)/2+kn} \binom{k-n-1}{k}_q \;,
```
which shows that
```math
\binom{n}{k}_q ∈ ℤ[q^{-1}] \quad \text{if } n < 0 \;.
```
In particular,
```math
\binom{n}{k}_q ∈ ℤ[q,q^{-1}] \quad \text{for all } n ∈ ℤ \;.
```
Now, for an element ``q`` of a ring ``R`` we define ``\binom{n}{k}_q`` as the
specialization of ``\binom{n}{k}_{\mathbf{q}}`` in ``q``, where ``q`` is
assumed to be invertible in ``R`` if
``n < 0``.

Note that for ``q=1`` we obtain
```math
\binom{n}{k}_1 = \binom{n}{k} \;,
```
hence the quantum binomial coefficient is a "deformation" of the usual binomial coefficient.

# Functions

The functions `quantum_binomial` work analogously to [`quantum_integer`](@ref).

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

# References
1. [Con00](@cite)
"""
function quantum_binomial(n::IntegerUnion, k::IntegerUnion, q::RingElem)
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

function quantum_binomial(n::IntegerUnion, k::IntegerUnion, q::Integer)
  # this method is not type stable, but that is acceptable as it is only
  # intended for convenience during interactive use
  if n > 0
    return quantum_binomial(n,k,ZZ(q))
  else
    return quantum_binomial(n,k,QQ(q))
  end
end

function quantum_binomial(n::IntegerUnion, k::IntegerUnion)
  R,q = laurent_polynomial_ring(ZZ, "q")
  return quantum_binomial(n,k,q)
end
