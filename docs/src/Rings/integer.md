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
relevant here are `Int` and `BigInt`. All the Julia integer types
belong to `Integer`.

The `Int` type is for machine integers which are highly efficient, but
can only represent integers up to a certain size, and most basic arithmetic
operations are performed unchecked, that is, they can [silently overflow](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Overflow-behavior). The
`Int` type is the type of literal input such as `12`, and should be used for
loop control flow, array indices, and other situations where the overflow can
be provably avoided.

The `BigInt` type is backed by GMP multiprecision integers and can represent
integers whose size is usually only limited by available memory. While the
`BigInt` type avoids overflow problems, it can be relatively slow in the
`Int` range.

Oscar currently has the integer type `fmpz`, which for performance
reasons scales internally from machine integers to GMP multiprecision integers.

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

* `zero(ZZ)`
* `one(ZZ)`

```@repl oscar
zero(ZZ)
one(ZZ)
```

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
also accept Julia integers. If all of the arguments to an Oscar function are
julia integers, the resulting integers should be julia integers. However, once
at least one of the arguments is an `fmpz`, the function will generally behave
as if all integer arguments were promoted to the type `fmpz`, and the integers
in the return generally should also be of type `fmpz`. For example:

```@repl oscar
divexact(ZZ(234), 2)
typeof(gcd(4, 6))
typeof(gcdx(4, 6))
typeof(gcd(4, ZZ(6)))
typeof(gcdx(4, ZZ(6)))
typeof(jacobi_symbol(ZZ(2), ZZ(3)))
```

In the first example, `2` is a Julia integer but is still valid in the
call to the Oscar function `divexact`. In the last example, the exceptional
function `jacobi_symbol` returns an `Int` as this will always be able to hold
the three possible return values of `-1`, `0`, or `1`.

## Predicates

* `iszero(n::fmpz) -> Bool`
* `isone(n::fmpz) -> Bool`
* `is_unit(n::fmpz) -> Bool`
* `isodd(n::fmpz) -> Bool`
* `iseven(n::fmpz) -> Bool`
* `is_square(n::fmpz) -> Bool`
* `is_prime(n::fmpz) -> Bool`
* `is_probable_prime(n::fmpz) -> Bool`

The `is_prime` predicate will prove primality, whereas `is_probable_prime` may
declare a composite number to be prime with very low probability.

Negative numbers, ``0`` and ``1`` are not considered prime by `is_prime` and
`is_probable_prime`.

```@repl oscar
isone(ZZ(1))
is_unit(ZZ(-1))
is_square(ZZ(16))
is_probable_prime(ZZ(23))
```

## Properties

* `sign(n::fmpz) -> fmpz`

Return the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```@repl oscar
sign(ZZ(23))
sign(ZZ(0))
sign(ZZ(-1))
```

* `abs(n::fmpz) -> fmpz`

Return the absolute value of ``n``, i.e. ``n`` if ``n \geq 0`` and ``-n``
otherwise


```@repl oscar
abs(ZZ(-3))
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

* `divexact(a::fmpz, b::fmpz) -> fmpz`

Return the quotient of ``a`` by ``b``. The result of the exact division of two
integers will always be another integer. Exact division raises an exception if
the division is not exact, or if division by zero is attempted.

```@repl oscar
divexact(ZZ(6), ZZ(3))
divexact(ZZ(6), ZZ(0))
divexact(ZZ(6), ZZ(5))
divexact(ZZ(6), 2)
```

### Powering

* `^(a::fmpz, b::Int) -> fmpz`

Return the result of powering ``a`` by ``b``.

```@repl oscar
ZZ(37)^37
ZZ(1)^(-2)
```

!!! note
    An exception will be raised if an integer other than ``-1`` or ``1`` is
    raised to a negative exponent.

!!! note
    In Julia `2^-2` is called a literal power. The value returned is a
    floating point value. To get behaviour that agrees with Oscar, one can
    write `2^Int(-2)`.

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

Integer Euclidean division of ``a`` by ``b`` computes a quotient and
remainder such that
```@math
a = qb + r
```
with ``|r| < |b|``.

### Division with remainder

* `divrem(a::fmpz, b::fmpz) -> (fmpz, fmpz)` : division with remainder
* `div(a::fmpz, b::fmpz) -> fmpz` : quotient only
* `rem(a::fmpz, b::fmpz) -> fmpz` : remainder only

Both `rem` and `divrem` compute the remainder ``r`` such that when ``r \neq 0``
the sign of ``r`` is the same as the sign of ``a``.

All three functions raise an exception if the modulus ``b`` is zero.

```@repl oscar
divrem(ZZ(5), ZZ(3))
div(ZZ(7), ZZ(2))
rem(ZZ(4), ZZ(3))
div(ZZ(2), ZZ(0))
```

!!! note
    The rem function does not provide a minimal set of representatives, e.g.
    `rem(-2, 3) = -2` but `rem(1, 3) = 1`.

## Modular arithmetic

### Modular reduction

* `mod(a::fmpz, b::fmpz) -> fmpz` : remainder only

The `mod` function computes a remainder ``r`` such that when ``r \neq 0`` the
sign of ``r`` is the same as the sign of ``b``. Thus, if ``b > 0`` then
`mod(a, b)` will be in the range ``[0, b)``. An exception is raised if the
modulus ``b`` is zero. This is summarised in the following table.

remainder | division   | sign             | rounding
----------|------------|------------------|---------------------
rem       | div/divrem | same as dividend | towards zero
mod       |            | same as divisor  | towards ``-\infty``

There is no function implemented to compute the quotient corresponding to
the remainder given by `mod`.

```@repl oscar
mod(ZZ(4), ZZ(3))
mod(ZZ(2), ZZ(0)) 
```

## [Divisibility testing](@id integer_divisibility_testing)

* `divides(a::fmpz, b::fmpz) -> (Bool, fmpz)`

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

### Greatest common divisor

* `gcd(a::fmpz, b::fmpz) -> fmpz`

Return the greatest common divisor of its inputs, which is by definition the
largest integer dividing the two inputs, unless both inputs are zero in which
case it returns zero. The result will always be non-negative and will only be
zero if both inputs are zero.

```@repl oscar
gcd(ZZ(34), ZZ(17))
gcd(ZZ(3), ZZ(0))
```

### Extended GCD

* `gcdx(a::fmpz, b::fmpz) -> (fmpz, fmpz, fmpz)`

Return a tuple ``(g, s, t)`` such that ``g`` is the greatest common divisor of
``a`` and ``b`` and ``g = as + bt``. Normally ``s`` and ``t`` are chosen so
that ``|s| < |b|/(2g)`` and ``|t| < |a|/(2g)``, where this uniquely defines
``s`` and ``t``. The following cases are handled specially:
1. if ``|a| = |b|`` then ``t = b/|b|``
2. if ``b = 0`` or ``|b| = 2g`` then ``s = a/|a|``
3. if ``a = 0`` or ``|a| = 2g`` then ``t = b/|b|``

### Least common multiple

* `lcm(a::fmpz, b::fmpz) -> fmpz`

Return the least common multiple of ``a`` and ``b``. This is the least
positive multiple of ``a`` and ``b``, unless ``a = 0`` or ``b = 0``
which case we define the least common multiple to be zero.

```@repl oscar
lcm(ZZ(6), ZZ(21))
lcm(ZZ(0), ZZ(0))
```

## Roots

### Square roots

Julia and Oscar distinguish two kinds of square root:

* Integer square root (`isqrt`)
* Floating point square root (`sqrt`)

We describe only the first of these here.

* `isqrt(n::fmpz) -> fmpz`

Return the floor of the square root of its argument, i.e. the largest integer
whose square does not exceed its input. An exception is raised if a negative
input is passed.

```@repl oscar
isqrt(ZZ(16))
isqrt(ZZ(0))
isqrt(ZZ(5))
isqrt(ZZ(-3))
```

* `isqrtrem(n::fmpz) -> (fmpz, fmpz)`

Return the tuple `(s, r)` such that ``s`` is equal to `isqrt(n)` and
``n = s^2 + r``.

```@repl oscar
isqrtrem(ZZ(16))
isqrtrem(ZZ(5))
```

### General roots

* `root(a::fmpz, n::Int) -> fmpz`

Return the value ``r`` of largest absolute value such that ``r^n \leq a``.
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

* `Int(n::fmpz) -> Int`
* `BigInt(n::fmpz) -> BigInt`

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

* `fits(::Type{Int}, n::fmpz) -> Bool`

Return `true` if the Oscar integer will fit in an `Int`.

```@repl oscar
fits(Int, ZZ(123))
fits(Int, ZZ(12348732648732648763274868732687324))
```

## Factorisation

* `factor(n::fmpz) -> Fac{fmpz}`

Return a factorisation of the given integer. The return value is a special
factorisation struct which can be manipulated using the functions below.

```@repl oscar
factor(ZZ(-6000361807272228723606))
factor(ZZ(0))
```

* `unit(F::Fac) -> fmpz`

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

For convenience, a `Int` can be used instead of an Oscar integer for this
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

!!! note
    The functions in this section that take `Int` arguments will return an
    `Int`, which may overflow or throw an error. Use the `fmpz` versions if
    this is not the desired behaviour.

### Factorial

* `factorial(n::fmpz) -> fmpz`

Return the factorial of ``n``, i.e. ``n!``. An exception is raised if
``n < 0``. We define ``0! = 1``.

* `rising_factorial(x::Int, n::Int) -> Int`
* `rising_factorial(x::fmpz, n::Int) -> fmpz`
* `rising_factorial(x::fmpz, n::fmpz) -> fmpz`

Return ``x(x + 1)(x + 2)\ldots(x + n - 1)``. An exception is raised if
``n < 0``. We define `rising_factorial(x, 0)` to be ``1``.

```@repl oscar
factorial(ZZ(30))
rising_factorial(ZZ(-30), 3)
```

### Primorial

* `primorial(n::Int) -> Int`
* `primorial(n::fmpz) -> fmpz`

Return the primorial ``P(n)``, i.e. the product of all primes less than or
equal to ``n``. An exception is raised if ``n < 0``. We define
``P(0) = P(1) = 1``.

```@repl oscar
primorial(ZZ(100))
```

### Bell numbers

* `bell(n::Int) -> Int`
* `bell(n::fmpz) -> fmpz`

Return the ``n``-th Bell number ``B(n)``, i.e. the number of ways of
partitioning a set of ``n`` elements. An exception is raised if ``n < 0``.

```@repl oscar
bell(ZZ(20))
```

### Binomial coefficients

* `binomial(n::fmpz, k::fmpz) -> fmpz`

Return the binomial coefficient ``\frac{n (n-1) \cdots (n-k+1)}{k!}`` for
``k \ge 0`` and returns `0` for `k < 0`.

!!! note
    Julia already defines the `binomial` function for `Int`, which throws an
    error on overflow.

```@repl oscar
binomial(ZZ(72), ZZ(15))
```

### Integer partitions

* `number_of_partitions(n::Int) -> Int`
* `number_of_partitions(n::fmpz) -> fmpz`

Return the number of integer partitions ``p(n)`` of ``n``, i.e. the number
of distinct ways to write ``n`` as a sum of positive integers. Note that
``p(0) = 1``, as the empty sum is counted. For ``n < 0`` we return zero.

```@repl oscar
number_of_partitions(ZZ(10^6))
```

### Fibonacci sequence

* `fibonacci(n::Int) -> Int`
* `fibonacci(n::fmpz) -> fmpz`

Return the ``n``-th Fibonacci number ``F(n)``, defined by the recurrence
relation ``F(1) = 1``, ``F(2) = 1`` and ``F(n) = F(n - 1) + F(n - 2)`` for
``n \geq 3``. We define ``F(0) = 0`` and for ``n > 0`` we have
``F(-n) = (-1)^{n+1}F(n)``.

```@repl oscar
fibonacci(ZZ(100))
fibonacci(-2)
```

## Number theoretic functionality

!!! note
    The functions in this section that take `Int` arguments will return a
    `Int`, which may overflow or throw an error. Use the `fmpz` versions if
    this is not the desired behaviour.

### Moebius mu function

* `moebius_mu(n::Int) -> Int`
* `moebius_mu(n::fmpz) -> Int` 

Return the Moebius function ``\mu(n)``, which is defined to be ``0`` if
``n`` is not squarefree and otherwise is defined to be ``+1`` or ``-1`` if
``n`` has an even or odd number of prime factors, respectively. Alternatively,
``\mu(n)`` can be defined to be the sum of the primitive ``n``-th roots of
unity. An exception is raised if ``n \leq 0``.

```@repl oscar
moebius_mu(30)
```

### Jacobi symbols

* `jacobi_symbol(m::Int, n::Int) -> Int`
* `jacobi_symbol(m::fmpz, n::fmpz) -> Int`

Return the Jacobi symbol ``\left(\frac{m}{n}\right)``, which is defined for
integers ``m`` and odd, positive integers ``n``. If the factorisation of ``n``
is ``n = p_1^{i_1}p_2^{i_2}\ldots p_r^{i_r}`` then we define
```math
\left(\frac{m}{n}\right) = \left(\frac{m}{p_1}\right)^{i_1}\left(\frac{m}{p_2}\right)^{i_2}\ldots \left(\frac{m}{p_r}\right)^{i_r}
```
where ``\left(\frac{m}{p}\right)`` on the right hand side is the Legendre
symbol, which is defined for an odd prime number ``p`` to be ``0`` if ``p``
divides ``m`` and otherwise ``+1`` or ``-1`` depending on whether ``m`` is
a square modulo ``p`` or not. An exception is raised if ``n`` is even or if
``n \leq 0``.

```@repl oscar
jacobi_symbol(3, 37)
```

### Sigma function

* `divisor_sigma(m::Int, n::Int) -> Int`
* `divisor_sigma(m::fmpz, n::Int) -> fmpz`
* `divisor_sigma(m::fmpz, n::fmpz) -> fmpz`

Return the sum of the ``n``-th powers of the divisors of ``m``
```math
\sigma(m, n) = \sum_{d\;|\;m} d^n.
```
If ``m \leq 0`` or ``n < 0`` we raise an exception.

```@repl oscar
divisor_sigma(60, 5)
```

### Euler totient function

* `euler_phi(n::Int) -> Int`
* `euler_phi(n::fmpz) -> fmpz`

Return the Euler totient function ``\varphi(n)``, i.e. the number of positive
integers ``1 \leq x \leq n`` which are coprime to ``n``. Note that
``\varphi(1) = 1``. We raise an exception if ``n \leq 0``.

```@repl oscar
euler_phi(200)
```

