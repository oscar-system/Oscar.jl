```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["integer.md"]
```

# Integers

An important design decision in Oscar.jl is to use Julia as the user language
by default. This means that integers typed at the
[REPL](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop)
are [Julia integers](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/). However, for performance reasons, Oscar has its own integer format.

Julia has a number of different integer types, but the two that are most
relevant here are `Base.Int` and `Base.BigInt`. All the Julia integer types
belong to `Base.Integer`.

The `Base.Int` type is for machine integers which are highly efficient, but
can only represent integers up to a certain hardware defined size before
wrapping around.

The `Base.BigInt` type is backed by GMP multiprecision integers and can
represent integers whose size is usually only limited by available memory.

Oscar currently only has one integer type, `Oscar.fmpz` which for performance
reasons scales internally from machine integers to GMP multiprecision integers.
The Oscar integer type belongs to `Oscar.Integer`.

This situation is illustrated in the following diagram.

![alt text](img/IntegerTypes.svg)

In the documentation below, we always use `Base.Integer` for a Julia integer
and `Oscar.Integer` for an Oscar integer. Some functions accept only machine
integers for certain arguments; in such cases, we refer to `Base.Int`.

## The ring of integers

Every object in Oscar representing a mathematical element has a parent. This is
an object encoding information about where that element belongs.

The parent of an Oscar integer is the ring of integers `ZZ`.

```@repl oscar
ZZ
```

### Integer constructors

Oscar integers are created using `ZZ`:

```@repl oscar
ZZ(2)^100
ZZ(618970019642690137449562111)
```
One can also construct the integer ``0`` with the empty constructor:

```@repl oscar
ZZ()
```

The following special constructors are also provided:

* `zero(ZZ)` : the integer 0
* `one(ZZ)` : the integer 1

Note that `ZZ` is not a Julia type, but the above methods of constructing
Oscar integers are similar to the way that Julia integer types can be used to
construct Julia integers.

```@repl oscar
Int(123)
BigInt(123456343567843598776327698374259876295438725)
zero(BigInt)
one(Int)
```

### Limitations

Oscar integers have the same limitations as [GMP](https://gmplib.org/)
multiprecision integers, namely that they are limited by the available memory
on the machine and in any case to signed integers whose absolute value does not
exceed ``2^{37}`` bits.

!!! note
    The Julia `Int` type is either a 32 or 64 bit integer, depending on the
    machine architecture (usually 64 bits on most modern machines). The range of
    values is machine dependent, but can be found by typing `typemin(Int)` and
    `typemax(Int)` in Julia.

## Julia integers in Oscar functions

For convenience, all basic arithmetic and exact division functions in Oscar
also accept Julia integers. For example:

```@repl oscar
divexact(ZZ(234), 2)
```

In this example, `2` is a Julia integer but is still valid in the
call to the Oscar function `divexact`.

## Predicates

* `iszero(n::Oscar.Integer) -> Bool`
* `isone(n::Oscar.Integer) -> Bool`
* `isunit(n::Oscar.Integer) -> Bool`
* `issquare(n::Oscar.Integer) -> Bool`
* `isprime(n::Oscar.Integer) -> Bool`
* `isprobable_prime(n::Oscar.Integer) -> Bool`

The `isprime` predicate will prove primality, whereas `isprobable_prime` may
declare a composite number to be prime with very low probability.

Negative numbers, ``0`` and ``1`` are not considered prime by `isprime` and
`isprobable_prime`.

```@repl oscar
isone(ZZ(1))
isunit(ZZ(-1))
issquare(ZZ(16))
isprobable_prime(ZZ(23))
```

## Properties

* `sign(n::Oscar.Integer) -> Base.Int`

Returns the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```@repl oscar
sign(ZZ(23))
sign(ZZ(0))
sign(ZZ(-1))
```

## Basic arithmetic

Oscar provides the basic arithmetic operations `+`, `-` and `*` and comparison
operators `==`, `!=`, `<`, `<=`, `>`, `>=`, including mixed operations between
Julia and Oscar integers. It also provides division and powering as described
below.

### Division in Oscar

Oscar distinguishes a number of different kinds of division:

* [Exact division](@ref integer_exact_division) (`divexact`)
* [Euclidean division](@ref integer_euclidean_division) (`div`, `rem`, `divrem` and `mod`)
* Construction of fractions (`a//b`)
* Floating point division (`a/b`)
* [Divisibility testing](@ref integer_divisibility_testing) (`divides`)

These choices have been made for maximum parsimony with the Julia language.

!!! note
    It is a common error to enter `1/2` for the fraction 'one half' in Julia.
    This expression is reserved for floating point division. Instead, the
    double slash operator `//` should be used for fractions.

### [Exact Division](@id integer_exact_division)

* `divexact(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer`

One or both arguments may be Julia integers, however if they both are, the
result will be a Julia integer not an Oscar integer.

The result of the exact division of two integers will always be another
integer. Exact division raises an exception if the division is not exact, or if
division by zero is attempted.

```@repl oscar
divexact(ZZ(6), ZZ(3))
divexact(ZZ(6), ZZ(0))
divexact(ZZ(6), ZZ(5))
```

### Powering

* `^(a::Oscar.Integer, b::Base.Int) -> Oscar.Integer`

Powering can be accomplished naturally using the special caret infix
operator:

```@repl oscar
ZZ(37)^37
ZZ(1)^(-2)
```

!!! note
    An exception will be raised if an integer other than ``-1`` or ``1`` is
    raised to a negative exponent.

The following is allowed for convenience.

```@repl oscar
ZZ(0)^0
```

!!! note
    In Julia, `2^64` will return zero, as the Julia integer ``2`` is a machine
    integer. In Oscar, the expression `ZZ(2)^64` will return the expected
    result, just as the Julia equivalent `BigInt(2)^64` does.


## [Euclidean division](@id integer_euclidean_division)

The ring of integers is a Euclidean domain and Oscar provides Euclidean
division through the functions `divrem`, `div` and `rem`.

* `divrem(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer` : division with remainder
* `div(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer` : quotient only
* `rem(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer` : remainder only

One or both arguments may be Julia integers, however if they both are, the
result will be a Julia integer not an Oscar integer.

Integer Euclidean division of ``a`` by ``b`` computes a quotient and
remainder such that
```@math
a = qb + r
```
where $|r| < |b|$. 

When ``r \neq 0`` the sign of ``r`` will be the same as the sign of ``a``.

* `mod(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer` : remainder only

One or both arguments may be Julia integers, however if they both are, the
result will be a Julia integer not an Oscar integer.

Returns a remainder which is either ``0`` or the same sign as ``b``. Thus, if
``b > 0`` then `mod(a, b)` will be in the range ``[0, b)``. There is no
function implemented to compute the corresponding quotient.

remainder | division   | sign             | rounding
----------|------------|------------------|---------------------
rem       | div/divrem | same as dividend | towards zero
mod       |            | same as divisor  | towards ``-\infty``

```@repl oscar
q, r = divrem(ZZ(5), ZZ(3))
q = div(ZZ(7), ZZ(2))
r = mod(ZZ(4), ZZ(3))
```

All three functions raise an exception if the modulus ``b`` is zero.

!!! note
    The rem function does not provide a minimal set of representatives, e.g.
    `rem(-2, 3) = -2` but `rem(1, 3) = 1`.

## [Divisibility testing](@id integer_divisibility_testing)

* `divides(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer`

In Oscar, we say that ``b`` divides ``a`` if there exists ``c`` in the same
ring such that ``a = bc``.

The call `divides(a, b)` returns a tuple `(flag, q)` where `flag` is either
`true` if `b` divides `a` in which case `q` will be a quotient, or `flag` is
`false` if `b` does not divide `a` in which case `q` will be an integer whose
value is not defined.
 
```@repl oscar
divides(ZZ(6), ZZ(3))
divides(ZZ(5), ZZ(2))
```

Note that for convenience we define:

```@repl oscar
divides(ZZ(0), ZZ(0))
```

## Greatest common divisor

* `gcd(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer`

One or both arguments may be Julia integers, however if they both are, the
result will be a Julia integer not an Oscar integer.

Returns the greatest common divisor of its inputs, which is by definition the
largest integer dividing the two inputs, unless both inputs are zero in which
case it returns zero. The result will always be non-negative and will only be
zero if both inputs are zero.

```@repl oscar
gcd(ZZ(34), ZZ(17))
gcd(ZZ(3), ZZ(0))
```

* `lcm(a::Oscar.Integer, b::Oscar.Integer) -> Oscar.Integer`

One or both arguments may be Julia integers, however if they both are, the
result will be a Julia integer not an Oscar integer.

Returns the least common multiple of ``a`` and ``b``. This is the least
positive multiple of ``a`` and ``b``, unless ``a = 0`` or ``b = 0``
which case we define the least common multiple to be zero.

```@repl oscar
lcm(ZZ(6), ZZ(21))
lcm(ZZ(0), ZZ(0))
```

!!! note
    The identity ``\gcd(m, n)\mathrm{lcm}(m, n) = mn`` does not hold for the
    definition that Oscar uses, unless both ``m`` and ``n`` are the same sign
    or one of them is zero.

## Roots

Julia and Oscar distinguish two kinds of square root:

* Integer square root (`isqrt`)
* Floating point square root (`sqrt`)

We describe only the first of these here.

* `isqrt(n::Oscar.Integer) -> Oscar.Integer`

Returns the floor of the square root of its argument, i.e. the largest integer
whose square does not exceed its input. An exception is raised if a negative
input is passed.

```@repl oscar
isqrt(ZZ(16))
isqrt(ZZ(0))
isqrt(ZZ(5))
isqrt(ZZ(-3))
```

* `isqrtrem(n::Oscar.Integer) -> (Oscar.Integer, Oscar.Integer)`

Returns the tuple `(s, r)` such that ``s`` is equal to `isqrt(n)` and
``n = s^2 + r``.

```@repl oscar
isqrtrem(ZZ(16))
isqrtrem(ZZ(5))
```

* `root(a::Oscar.Integer, n::Base.Int) -> Oscar.Integer`

Returns the value ``r`` of largest absolute value such that ``r^n \leq a``.
When ``a`` is a perfect ``n``-th power, the return value will be an ``n``-th
root of ``a``.

When ``n`` is even, the non-negative root is always returned. An exception is
raised if ``n \leq 0`` or if ``n`` is even and ``a < 0``.

```@repl oscar
root(ZZ(16), 4)
root(ZZ(5), 2)
root(ZZ(-5), 3)
root(ZZ(0), 4)
root(ZZ(-5), 2)
root(ZZ(12), -2)
```

## Conversions

* `Int(n::Oscar.Integer) -> Base.Int`
* `BigInt(n::Oscar.Integer) -> Base.BigInt`

Convert the Oscar integer to the respective Julia integer.

```@repl oscar
n = ZZ(123)
Int(n)
BigInt(n)
```

In the case of `Int`, if the Oscar integer is too large to fit, an exception
is raised.

```@repl oscar
Int(ZZ(12348732648732648763274868732687324))
```

* `fits(::Type{Int}, n::Oscar.Integer) -> Bool`

Returns `true` if the Oscar integer will fit in a `Base.Int`.

```@repl oscar
fits(Int, ZZ(123))
fits(Int, ZZ(12348732648732648763274868732687324))
```

## Factorisation

* `factor(n::Oscar.Integer) -> Fac`

Returns a factorisation of the given integer. The return value is a special
factorisation struct which can be manipulated using the functions below.

```@repl oscar
factor(ZZ(-6000361807272228723606))
factor(ZZ(0))
```

* `unit(F::Fac) -> Oscar.Integer`

```@repl oscar
F = factor(ZZ(-12))
unit(F)
```

### Factorisation are iterable

Once created, a factorisation is iterable:

```@repl oscar
F = factor(ZZ(-60))
for (p, e) in F; println("$p^$e"); end
```

The pairs `(p, e)` in a factorisation represent the prime power factors
``p^e`` of the non-unit part of the factorisation. They can be placed in an
array using `collect`:

```@repl oscar
F = factor(ZZ(-60))
collect(F)
```

### Accessing exponents in a factorisation

One can also determine whether a given prime is in the non-unit part of a
factorisation and if so return its exponent. If the exponent of a prime that
is not in a factorisation is requested, an exception is raised.

For convenience, a `Base.Int` can be used instead of an Oscar integer for this
functionality.

```@repl oscar
F = factor(ZZ(-60))
5 in F
ZZ(3) in F
7 in F
F[3]
F[ZZ(7)]
```

## Combinatorial functions

* `Oscar.factorial(n::Base.Int) -> Oscar.Integer`

Returns the factorial of ``n``, i.e. ``n!``. An exception is raised if
``n < 0``. We define ``0! = 1``.

!!! note
    The function `factorial` is already defined in Julia, but returns a Julia
    `Base.Int`, which overflows when the result is too large. To disambiguate
    the Oscar version of the function it is accessed via `Oscar.factorial`.

* `rising_factorial(x::Oscar.Integer, n::Base.Int) -> Oscar.Integer`

Returns ``x(x + 1)(x + 2)\ldots(x + n - 1)``. An exception is raised if
``n < 0``. We define `rising_factorial(x, 0)` to be ``1``.

* `primorial(n::Base.Int) -> Oscar.Integer`

Returns the ``n``-th primorial number ``P(n)``, i.e. the product of all primes
less than or equal to ``n``. An exception is raised if ``n < 0``. We define
``P(0) = P(1) = 1``.

* `bell(n::Base.Int) -> Oscar.Integer`

Returns the ``n``-th Bell number ``B(n)``, i.e. the number of ways of
partitioning a set of ``n`` elements. An exception is raised if ``n < 0``.

* `Oscar.binomial(n::Base.Int, k::Base.Int) -> Oscar.Integer`

Returns the binomial coefficient ``\frac{n!}{k!(n - k)!}``. If ``n, k < 0`` or
``k > n`` we return zero.

!!! note
    Julia already defines the `binomial` function,  which returns a `Base.Int`
    that may overflow when the result is too large. To disambiguate the Oscar
    version of the function it is accessed via `Oscar.binomial`.

* `number_of_partitions(n::Oscar.Integer) -> Oscar.Integer`
* `number_of_partitions(n::Base.Integer) -> Oscar.Integer`

Returns the number of integer partitions ``p(n)`` of ``n``, i.e. the number
of distinct ways to write ``n`` as a sum of positive integers. Note that
``p(0) = 1``, as the empty sum is counted. For ``n < 0`` we return zero.

## Number theoretic functionality

* `fibonacci(n::Base.Int) -> Oscar.Integer`

Returns the ``n``-th Fibonacci number ``F(n)``, defined by the recurrence
relation ``F(1) = 1``, ``F(2) = 1`` and ``F(n) = F(n - 1) + F(n - 2)`` for
``n \geq 3``. For convenience we define ``F(0) = 0``. An exception is raised
if ``n < 0``.

* `moebius_mu(n::Oscar.Integer) -> Base.Int`
* `moebius_mu(n::Base.Integer) -> Base.Int` 

Return the Moebius function ``\mu(n)``, which is defined to be ``0`` if
``n`` is not squarefree and otherwise is defined to be ``+1`` or ``-1`` if
``n`` has an even or odd number of prime factors, respectively. Alternatively,
``\mu(n)`` can be defined to be the sum of the primitive ``n``-th roots of
unity. An exception is raised if ``n < 0``.

* `jacobi_symbol(m::Oscar.Integer, n::Oscar.Integer) -> Base.Int`
* `jacobi_symbol(m::Base.Integer, n::Base.Integer) -> Base.Int`

Return the Jacobi symbol ``\left(\frac{m}{n}\right)``, which is defined for
integers ``m`` and odd positive integers ``n``. If the factorisation of ``n``
is ``n = p_1^{i_1}p_2^{i_2}\ldots p_r^{i_r}`` then we define
```math
\left(\frac{m}{n}\right) = \left(\frac{m}{p_1}\right)^{i_1}\left(\frac{m}{p_2}\right)^{i_2}\ldots \left(\frac{m}{p_r}\right)^{i_r}
```
where ``\left(\frac{m}{p}\right)`` on the right hand side is the Legendre
symbol, which is defined for an odd prime number ``p`` to be ``0`` if ``p``
divides ``m`` and otherwise ``+1`` or ``-1`` depending on whether ``m`` is
a square modulo ``p`` or not. An exception is raised if ``n`` is even or
not positive.

* `divisor_sigma(m::Oscar.Integer, n::Oscar.Integer) -> Oscar.Integer`
* `divisor_sigma(m::Base.Integer, n::Base.Integer) -> Oscar.Integer`

Return the sum of the ``n``-th powers of the divisors of ``m``
```math
\sigma(m, n) = \sum_{d\;|\;m} d^n.
```
We define ``\sigma(0, n) = 0`` for all ``n``. If ``n < 0`` we raise an
exception.

* `euler_phi(n::Oscar.Integer) -> Oscar.Integer`
* `euler_phi(n::Base.Integer) -> Oscar.Integer`

Return the Euler totient function ``\varphi(n)``, i.e. the number of positive
integers ``1 \leq x \leq n`` which are coprime to ``n``. Note that
``\varphi(1) = 1`` and ``\varphi(0) = 0``. We raise an exception if ``n < 0``.




