```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Integers

An important design decision in Oscar.jl is to use Julia as the user language
by default. This means that integers typed at the
[REPL](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop)
are [Julia integers](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/). However, for performance reasons, OSCAR has its own integer format.

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

OSCAR currently has the integer type `ZZRingElem`, which for performance
reasons scales internally from machine integers to GMP multiprecision integers.

## The ring of integers

Every object in OSCAR representing a mathematical element has a parent. This is
an object encoding information about where that element belongs.

The parent of an OSCAR integer is the ring of integers `ZZ`.

```jldoctest
julia> ZZ
Integer Ring

```

### Integer constructors

OSCAR integers are created using `ZZ`:

```jldoctest
julia> ZZ(2)^100
1267650600228229401496703205376

julia> ZZ(618970019642690137449562111)
618970019642690137449562111

```
One can also construct the integer ``0`` with the empty constructor:

```jldoctest
julia> ZZ()
0

```

The following special constructors are also provided:

* `zero(ZZ)`
* `one(ZZ)`

```jldoctest
julia> zero(ZZ)
0

julia> one(ZZ)
1

```

Note that `ZZ` is not a Julia type, but the above methods of constructing
OSCAR integers are similar to the way that Julia integer types can be used to
construct Julia integers.

```jldoctest
julia> Int(123)
123

julia> BigInt(123456343567843598776327698374259876295438725)
123456343567843598776327698374259876295438725

julia> zero(BigInt)
0

julia> one(Int)
1

```

### Limitations

OSCAR integers have the same limitations as [GMP](https://gmplib.org/)
multiprecision integers, namely that they are limited by the available memory
on the machine and in any case to signed integers whose absolute value does not
exceed ``2^{37}`` bits.

!!! note
    The Julia `Int` type is either a 32 or 64 bit integer, depending on the
    machine architecture (usually 64 bits on most modern machines). The range of
    values is machine dependent, but can be found by typing `typemin(Int)` and
    `typemax(Int)` in Julia.

## Julia integers in OSCAR functions

For convenience, all basic arithmetic and exact division functions in OSCAR
also accept Julia integers. If all of the arguments to an OSCAR function are
julia integers, the resulting integers should be julia integers. However, once
at least one of the arguments is an `ZZRingElem`, the function will generally behave
as if all integer arguments were promoted to the type `ZZRingElem`, and the integers
in the return generally should also be of type `ZZRingElem`. For example:

```jldoctest
julia> divexact(ZZ(234), 2)
117

julia> typeof(gcd(4, 6))
Int64

julia> typeof(gcdx(4, 6))
Tuple{Int64, Int64, Int64}

julia> typeof(gcd(4, ZZ(6)))
ZZRingElem

julia> typeof(gcdx(4, ZZ(6)))
Tuple{ZZRingElem, ZZRingElem, ZZRingElem}

julia> typeof(jacobi_symbol(ZZ(2), ZZ(3)))
Int64

```

In the first example, `2` is a Julia integer but is still valid in the
call to the OSCAR function `divexact`. In the last example, the exceptional
function `jacobi_symbol` returns an `Int` as this will always be able to hold
the three possible return values of `-1`, `0`, or `1`.

## Predicates

* `iszero(n::ZZRingElem) -> Bool`
* `isone(n::ZZRingElem) -> Bool`
* `is_unit(n::ZZRingElem) -> Bool`
* `isodd(n::ZZRingElem) -> Bool`
* `iseven(n::ZZRingElem) -> Bool`
* `is_square(n::ZZRingElem) -> Bool`
* `is_prime(n::ZZRingElem) -> Bool`
* `is_probable_prime(n::ZZRingElem) -> Bool`

The `is_prime` predicate will prove primality, whereas `is_probable_prime` may
declare a composite number to be prime with very low probability.

Negative numbers, ``0`` and ``1`` are not considered prime by `is_prime` and
`is_probable_prime`.

```jldoctest
julia> isone(ZZ(1))
true

julia> is_unit(ZZ(-1))
true

julia> is_square(ZZ(16))
true

julia> is_probable_prime(ZZ(23))
true

```

## Properties

* `sign(n::ZZRingElem) -> ZZRingElem`

Return the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```jldoctest
julia> sign(ZZ(23))
1

julia> sign(ZZ(0))
0

julia> sign(ZZ(-1))
-1

```

* `abs(n::ZZRingElem) -> ZZRingElem`

Return the absolute value of ``n``, i.e. ``n`` if ``n \geq 0`` and ``-n``
otherwise


```jldoctest
julia> abs(ZZ(-3))
3

```

## Basic arithmetic

OSCAR provides the basic arithmetic operations `+`, `-` and `*` and comparison
operators `==`, `!=`, `<`, `<=`, `>`, `>=`, including mixed operations between
Julia and OSCAR integers. It also provides division and powering as described
below.

### Division in OSCAR

OSCAR distinguishes a number of different kinds of division:

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

* `divexact(a::ZZRingElem, b::ZZRingElem) -> ZZRingElem`

Return the quotient of ``a`` by ``b``. The result of the exact division of two
integers will always be another integer. Exact division raises an exception if
the division is not exact, or if division by zero is attempted.

```jldoctest
julia> divexact(ZZ(6), ZZ(3))
2

julia> divexact(ZZ(6), ZZ(0))
ERROR: DivideError: integer division error

julia> divexact(ZZ(6), ZZ(5))
ERROR: ArgumentError: Not an exact division

julia> divexact(ZZ(6), 2)
3

```

### Powering

* `^(a::ZZRingElem, b::Int) -> ZZRingElem`

Return the result of powering ``a`` by ``b``.

```jldoctest
julia> ZZ(37)^37
10555134955777783414078330085995832946127396083370199442517

julia> ZZ(1)^(-2)
1

```

!!! note
    An exception will be raised if an integer other than ``-1`` or ``1`` is
    raised to a negative exponent.

!!! note
    In Julia `2^-2` is called a literal power. The value returned is a
    floating point value. To get behaviour that agrees with OSCAR, one can
    write `2^Int(-2)`.

The following is allowed for convenience.

```jldoctest
julia> ZZ(0)^0
1

```

!!! note
    In Julia, `2^64` will return zero, as the Julia integer ``2`` is a machine
    integer. In OSCAR, the expression `ZZ(2)^64` will return the expected
    result, just as the Julia equivalent `BigInt(2)^64` does.


## [Euclidean division](@id integer_euclidean_division)

The ring of integers is a Euclidean domain and OSCAR provides Euclidean
division through the functions `divrem`, `div` and `rem`.

Integer Euclidean division of ``a`` by ``b`` computes a quotient and
remainder such that
```@math
a = qb + r
```
with ``|r| < |b|``.

### Division with remainder

* `divrem(a::ZZRingElem, b::ZZRingElem) -> (ZZRingElem, ZZRingElem)` : division with remainder
* `div(a::ZZRingElem, b::ZZRingElem) -> ZZRingElem` : quotient only
* `rem(a::ZZRingElem, b::ZZRingElem) -> ZZRingElem` : remainder only

Both `rem` and `divrem` compute the remainder ``r`` such that when ``r \neq 0``
the sign of ``r`` is the same as the sign of ``a``.

All three functions raise an exception if the modulus ``b`` is zero.

```jldoctest
julia> divrem(ZZ(5), ZZ(3))
(1, 2)

julia> div(ZZ(7), ZZ(2))
3

julia> rem(ZZ(4), ZZ(3))
1

julia> div(ZZ(2), ZZ(0))
ERROR: DivideError: integer division error

```

!!! note
    The rem function does not provide a minimal set of representatives, e.g.
    `rem(-2, 3) = -2` but `rem(1, 3) = 1`.

## Modular arithmetic

### Modular reduction

* `mod(a::ZZRingElem, b::ZZRingElem) -> ZZRingElem` : remainder only

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

```jldoctest
julia> mod(ZZ(4), ZZ(3))
1

julia> mod(ZZ(2), ZZ(0))
ERROR: DivideError: integer division error

```

## [Divisibility testing](@id integer_divisibility_testing)

* `divides(a::ZZRingElem, b::ZZRingElem) -> (Bool, ZZRingElem)`

In OSCAR, we say that ``b`` divides ``a`` if there exists ``c`` in the same
ring such that ``a = bc``.

The call `divides(a, b)` returns a tuple `(flag, q)` where `flag` is either
`true` if `b` divides `a` in which case `q` will be a quotient, or `flag` is
`false` if `b` does not divide `a` in which case `q` will be an integer whose
value is not defined.
 
```jldoctest
julia> divides(ZZ(6), ZZ(3))
(true, 2)

julia> divides(ZZ(5), ZZ(2))
(false, 0)

```

Note that for convenience we define:

```jldoctest
julia> divides(ZZ(0), ZZ(0))
(true, 0)

```

## Greatest common divisor

### Greatest common divisor

* `gcd(a::ZZRingElem, b::ZZRingElem) -> ZZRingElem`

Return the greatest common divisor of its inputs, which is by definition the
largest integer dividing the two inputs, unless both inputs are zero in which
case it returns zero. The result will always be non-negative and will only be
zero if both inputs are zero.

```jldoctest
julia> gcd(ZZ(34), ZZ(17))
17

julia> gcd(ZZ(3), ZZ(0))
3

```

### Extended GCD

* `gcdx(a::ZZRingElem, b::ZZRingElem) -> (ZZRingElem, ZZRingElem, ZZRingElem)`

Return a tuple ``(g, s, t)`` such that ``g`` is the greatest common divisor of
``a`` and ``b`` and ``g = as + bt``. Normally ``s`` and ``t`` are chosen so
that ``|s| < |b|/(2g)`` and ``|t| < |a|/(2g)``, where this uniquely defines
``s`` and ``t``. The following cases are handled specially:
1. if ``|a| = |b|`` then ``t = b/|b|``
2. if ``b = 0`` or ``|b| = 2g`` then ``s = a/|a|``
3. if ``a = 0`` or ``|a| = 2g`` then ``t = b/|b|``

### Least common multiple

* `lcm(a::ZZRingElem, b::ZZRingElem) -> ZZRingElem`

Return the least common multiple of ``a`` and ``b``. This is the least
positive multiple of ``a`` and ``b``, unless ``a = 0`` or ``b = 0``
which case we define the least common multiple to be zero.

```jldoctest
julia> lcm(ZZ(6), ZZ(21))
42

julia> lcm(ZZ(0), ZZ(0))
0

```

## Roots

### Square roots

Julia and OSCAR distinguish two kinds of square root:

* Integer square root (`isqrt`)
* Floating point square root (`sqrt`)

We describe only the first of these here.

* `isqrt(n::ZZRingElem) -> ZZRingElem`

Return the floor of the square root of its argument, i.e. the largest integer
whose square does not exceed its input. An exception is raised if a negative
input is passed.

```jldoctest
julia> isqrt(ZZ(16))
4

julia> isqrt(ZZ(0))
0

julia> isqrt(ZZ(5))
2

julia> isqrt(ZZ(-3))
ERROR: DomainError with -3:
Argument must be non-negative

```

* `isqrtrem(n::ZZRingElem) -> (ZZRingElem, ZZRingElem)`

Return the tuple `(s, r)` such that ``s`` is equal to `isqrt(n)` and
``n = s^2 + r``.

```jldoctest
julia> isqrtrem(ZZ(16))
(4, 0)

julia> isqrtrem(ZZ(5))
(2, 1)

```

### General roots

* `root(a::ZZRingElem, n::Int) -> ZZRingElem`

Return an ``n``-th root of ``a`` or throw an error if it does not exist.

When ``n`` is even, the non-negative root is always returned. An exception is
raised if ``n \leq 0`` or if ``n`` is even and ``a < 0``.

```jldoctest
julia> root(ZZ(16), 4)
2

julia> root(ZZ(-5), 2)
ERROR: DomainError with (-5, 2):
Argument `x` must be positive if exponent `n` is even

julia> root(ZZ(12), -2)
ERROR: DomainError with -2:
Exponent must be positive
```

## Conversions

* `Int(n::ZZRingElem) -> Int`
* `BigInt(n::ZZRingElem) -> BigInt`

Convert the OSCAR integer to the respective Julia integer.

```jldoctest
julia> n = ZZ(123)
123

julia> Int(n)
123

julia> BigInt(n)
123

```

In the case of `Int`, if the OSCAR integer is too large to fit, an exception
is raised.

```jldoctest
julia> Int(ZZ(12348732648732648763274868732687324))
ERROR: InexactError: convert(Int64, 12348732648732648763274868732687324)

```

* `fits(::Type{Int}, n::ZZRingElem) -> Bool`

Return `true` if the OSCAR integer will fit in an `Int`.

```jldoctest
julia> fits(Int, ZZ(123))
true

julia> fits(Int, ZZ(12348732648732648763274868732687324))
false

```

## Factorisation

* `factor(n::ZZRingElem) -> Fac{ZZRingElem}`

Return a factorisation of the given integer. The return value is a special
factorisation struct which can be manipulated using the functions below.

```jldoctest
julia> factor(ZZ(-6000361807272228723606))
-1 * 2 * 229^3 * 43669^3 * 3

julia> factor(ZZ(0))
ERROR: ArgumentError: Argument is not non-zero

```

* `unit(F::Fac) -> ZZRingElem`

```jldoctest
julia> F = factor(ZZ(-12))
-1 * 2^2 * 3

julia> unit(F)
-1

```

### Factorisation are iterable

Once created, a factorisation is iterable:

```jldoctest
julia> F = factor(ZZ(-60))
-1 * 5 * 2^2 * 3

julia> for (p, e) in F; println("$p^$e"); end
5^1
2^2
3^1

```

The pairs `(p, e)` in a factorisation represent the prime power factors
``p^e`` of the non-unit part of the factorisation. They can be placed in an
array using `collect`:

```jldoctest
julia> F = factor(ZZ(-60))
-1 * 5 * 2^2 * 3

julia> collect(F)
3-element Vector{Pair{ZZRingElem, Int64}}:
 5 => 1
 2 => 2
 3 => 1

```

### Accessing exponents in a factorisation

One can also determine whether a given prime is in the non-unit part of a
factorisation and if so return its exponent. If the exponent of a prime that
is not in a factorisation is requested, an exception is raised.

For convenience, a `Int` can be used instead of an OSCAR integer for this
functionality.

```jldoctest
julia> F = factor(ZZ(-60))
-1 * 5 * 2^2 * 3

julia> 5 in F
true

julia> ZZ(3) in F
true

julia> 7 in F
false

julia> F[3]
1

julia> F[ZZ(7)]
ERROR: 7 is not a factor of -1 * 5 * 2^2 * 3

```

## Combinatorial functions

!!! note
    The functions in this section that take `Int` arguments will return an
    `Int`, which may overflow or throw an error. Use the `ZZRingElem` versions if
    this is not the desired behaviour.

### Factorial

* `factorial(n::ZZRingElem) -> ZZRingElem`

Return the factorial of ``n``, i.e. ``n!``. An exception is raised if
``n < 0``. We define ``0! = 1``.

* `rising_factorial(x::Int, n::Int) -> Int`
* `rising_factorial(x::ZZRingElem, n::Int) -> ZZRingElem`
* `rising_factorial(x::ZZRingElem, n::ZZRingElem) -> ZZRingElem`

Return ``x(x + 1)(x + 2)\ldots(x + n - 1)``. An exception is raised if
``n < 0``. We define `rising_factorial(x, 0)` to be ``1``.

```jldoctest
julia> factorial(ZZ(30))
265252859812191058636308480000000

julia> rising_factorial(ZZ(-30), 3)
-24360

```

### Primorial

* `primorial(n::Int) -> Int`
* `primorial(n::ZZRingElem) -> ZZRingElem`

Return the primorial ``P(n)``, i.e. the product of all primes less than or
equal to ``n``. An exception is raised if ``n < 0``. We define
``P(0) = P(1) = 1``.

```jldoctest
julia> primorial(ZZ(100))
2305567963945518424753102147331756070

```

### Bell numbers

* `bell(n::Int) -> Int`
* `bell(n::ZZRingElem) -> ZZRingElem`

Return the ``n``-th Bell number ``B(n)``, i.e. the number of ways of
partitioning a set of ``n`` elements. An exception is raised if ``n < 0``.

```jldoctest
julia> bell(ZZ(20))
51724158235372

```

### Binomial coefficients

* `binomial(n::ZZRingElem, k::ZZRingElem) -> ZZRingElem`

Return the binomial coefficient ``\frac{n (n-1) \cdots (n-k+1)}{k!}`` for
``k \ge 0`` and returns `0` for `k < 0`.

!!! note
    Julia already defines the `binomial` function for `Int`, which throws an
    error on overflow.

```jldoctest
julia> binomial(ZZ(72), ZZ(15))
1155454041309504

```

### Integer partitions

* `number_of_partitions(n::Int) -> Int`
* `number_of_partitions(n::ZZRingElem) -> ZZRingElem`

Return the number of integer partitions ``p(n)`` of ``n``, i.e. the number
of distinct ways to write ``n`` as a sum of positive integers. Note that
``p(0) = 1``, as the empty sum is counted. For ``n < 0`` we return zero.

```jldoctest
julia> number_of_partitions(ZZ(10^6))
1471684986358223398631004760609895943484030484439142125334612747351666117418918618276330148873983597555842015374130600288095929387347128232270327849578001932784396072064228659048713020170971840761025676479860846908142829356706929785991290519899445490672219997823452874982974022288229850136767566294781887494687879003824699988197729200632068668735996662273816798266213482417208446631027428001918132198177180646511234542595026728424452592296781193448139994664730105742564359154794989181485285351370551399476719981691459022015599101959601417474075715430750022184895815209339012481734469448319323280150665384042994054179587751761294916248142479998802936507195257074485047571662771763903391442495113823298195263008336489826045837712202455304996382144601028531832004519046591968302787537418118486000612016852593542741980215046267245473237321845833427512524227465399130174076941280847400831542217999286071108336303316298289102444649696805395416791875480010852636774022023128467646919775022348562520747741843343657801534130704761975530375169707999287040285677841619347472368171772154046664303121315630003467104673818

```

### Fibonacci sequence

* `fibonacci(n::Int) -> Int`
* `fibonacci(n::ZZRingElem) -> ZZRingElem`

Return the ``n``-th Fibonacci number ``F(n)``, defined by the recurrence
relation ``F(1) = 1``, ``F(2) = 1`` and ``F(n) = F(n - 1) + F(n - 2)`` for
``n \geq 3``. We define ``F(0) = 0`` and for ``n > 0`` we have
``F(-n) = (-1)^{n+1}F(n)``.

```jldoctest
julia> fibonacci(ZZ(100))
354224848179261915075

julia> fibonacci(-2)
-1

```

## Number theoretic functionality

!!! note
    The functions in this section that take `Int` arguments will return a
    `Int`, which may overflow or throw an error. Use the `ZZRingElem` versions if
    this is not the desired behaviour.

### Moebius mu function

* `moebius_mu(n::Int) -> Int`
* `moebius_mu(n::ZZRingElem) -> Int` 

Return the Moebius function ``\mu(n)``, which is defined to be ``0`` if
``n`` is not squarefree and otherwise is defined to be ``+1`` or ``-1`` if
``n`` has an even or odd number of prime factors, respectively. Alternatively,
``\mu(n)`` can be defined to be the sum of the primitive ``n``-th roots of
unity. An exception is raised if ``n \leq 0``.

```jldoctest
julia> moebius_mu(30)
-1

```

### Jacobi symbols

* `jacobi_symbol(m::Int, n::Int) -> Int`
* `jacobi_symbol(m::ZZRingElem, n::ZZRingElem) -> Int`

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

```jldoctest
julia> jacobi_symbol(3, 37)
1

```

### Sigma function

* `divisor_sigma(m::Int, n::Int) -> Int`
* `divisor_sigma(m::ZZRingElem, n::Int) -> ZZRingElem`
* `divisor_sigma(m::ZZRingElem, n::ZZRingElem) -> ZZRingElem`

Return the sum of the ``n``-th powers of the divisors of ``m``
```math
\sigma(m, n) = \sum_{d\;|\;m} d^n.
```
If ``m \leq 0`` or ``n < 0`` we raise an exception.

```jldoctest
julia> divisor_sigma(60, 5)
806220408

```

### Euler totient function

* `euler_phi(n::Int) -> Int`
* `euler_phi(n::ZZRingElem) -> ZZRingElem`

Return the Euler totient function ``\varphi(n)``, i.e. the number of positive
integers ``1 \leq x \leq n`` which are coprime to ``n``. Note that
``\varphi(1) = 1``. We raise an exception if ``n \leq 0``.

```jldoctest
julia> euler_phi(200)
80

```

