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

In the following, unless stated otherwise, when we refer to integers we mean
Oscar integers; when we refer to an `Int` we mean the Julia `Int`.

## Constructors

Oscar integers are created using the `ZZ` constructor.

```@repl oscar
ZZ(2)^100
ZZ(618970019642690137449562111)
```
One can also construct the zero integer with the empty constructor:

```@repl oscar
ZZ()
```

The following special constructors are also provided:

* `zero(ZZ)` : the integer 0
* `one(ZZ)` : the integer 1

### Limitations

Oscar integers have the same limitations as [GMP](https://gmplib.org/)
multiprecision integers, namely that they are limited by the available memory
on the machine and in any case to signed integers whose absolute value does not
exceed ``2^{37}`` bits.

!!! note
    The Julia 'Int' type is either a 32 or 64 bit integer, depending on the
    machine architecture (usually 64 bits on most modern machines). The range of
    values is machine dependent, but can be found by typing 'typemin(Int)' and
    'typemax(Int)' in Julia.

## Julia integers in Oscar functions

For convenience, basic arithmetic and exact division functions in Oscar also
accept Julia integers. For example:

```@repl oscar
divexact(ZZ(234), 2)
```

In this example, `2` is a Julia integer but is still valid in the
call to the Oscar function `divexact`.

## Predicates and properties

The following predicates are provided, which return `true` or `false`:

* `iszero(n)` : ``n = 0``
* `isone(n)` : ``n = 1``
* `isunit(n)` : ``n = \pm 1``

The following properties can also be computed:

* `sign(n)` returns the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0`` or ``0``
  otherwise. The return value is a Julia `Int`.

```@repl oscar
sign(ZZ(23))
sign(ZZ(0))
sign(ZZ(-1))
```

Every object in Oscar representing a mathematical element has a parent. This is
an object encoding information about where that element belongs.

The parent of an Oscar integer is the ring of integers `ZZ`.

```@repl oscar
R = parent(ZZ(2))
R == ZZ
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
    It is a common error to enter '1/2' for the fraction 'one half' in Julia.
    This expression is reserved for floating point division. Instead, the
    double slash operator '//' should be used for fractions.

### [Exact Division](@id integer_exact_division)

Exact division is carried out using the `divexact` function.

The result of the exact division of two integers will always be another
integer. Exact division raises an exception if the division is not exact, or if
division by zero is attempted.

```@repl oscar
divexact(ZZ(6), ZZ(3))
divexact(ZZ(6), ZZ(0))
divexact(ZZ(6), ZZ(5))
```

### Powering

Powering of integers is performed using the caret operator `^`. The exponent
can be any Julia `Int`.

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
    In Julia, '2^64' will return 0, as the Julia integer 2 is a machine word.
    In Oscar, the expression 'ZZ(2)^64' will return the expected result.


## [Euclidean division](@id integer_euclidean_division)

The ring of integers is a Euclidean domain and Oscar provides Euclidean
division.

In a Euclidean domain in Oscar the `divrem` function returns both quotient
and remainder, `div` returns just the quotient and `rem` returns just the
remainder.

For integers, Euclidean division of ``a`` by ``n`` computes a quotient and
remainder such that
```@math
a = qn + r
```
where $|r| < |n|$. For conformity with Julia, when ``r \neq 0`` the sign of
``r`` will be the same as the sign of ``a``.

If one instead wants Euclidean remainder with ``r`` and ``n`` having the same
sign, one can use `mod`. Then if ``n > 0`` the remainder will be non-negative.

remainder | division   | sign             | rounding
----------|------------|------------------|---------------------
rem       | div/divrem | same as dividend | towards zero
mod       |            | same as divisor  | towards ``-\infty``

```@repl oscar
q, r = divrem(ZZ(5), ZZ(3))
q = div(ZZ(7), ZZ(2))
r = mod(ZZ(4), ZZ(3))
```

All three functions raise an exception if the modulus ``m`` is zero.

!!! note
    The rem function does not provide a minimal set of representatives, e.g.
    rem(-2, 3) = -2 but rem(1, 3) = 1.

## [Divisibility testing](@id integer_divisibility_testing)

Divisibility testing is performed using the `divides` function.

In Oscar, we say that ``b`` divides ``a`` if there exists ``c`` in the same
ring such that ``a = bc``.

The call `divides(a, b)` returns a tuple `(flag, q)` where `flag` is either
`true` if `b` divides `a` in which case `q` will be a quotient, or `flag` is
`false` if `b` does not divide `a` and `q` will be an integer whose value is
not defined.
 
```@repl oscar
divides(ZZ(6), ZZ(3))
divides(ZZ(5), ZZ(2))
```

Note that for convenience we define:

```@repl oscar
divides(ZZ(0), ZZ(0))
```

## GCD and LCM

The `gcd` function returns the greatest common divisor of its inputs, which is
by definition the largest integer dividing the two inputs. The result will
always be non-negative and will only be zero if both inputs are zero.

```@repl oscar
gcd(ZZ(34), ZZ(17))
gcd(ZZ(3), ZZ(0))
```

The `lcm` function returns the least positive multiple of its inputs, unless
one or more of its inputs is zero, in which case it returns zero.

```@repl oscar
lcm(ZZ(6), ZZ(21))
lcm(ZZ(0), ZZ(0))
```

!!! note
    The identity ``\gcd(m, n)\lcm(m, n) = mn`` does not hold for the definition
    that Oscar uses, unless both ``m`` and ``n`` are the same sign or one of
    them is zero.

## Conversions

Oscar integers can be converted to Julia `Int`'s and `BigInt`'s in the usual
Julia way:

```@repl oscar
n = ZZ(123)
Int(n)
BigInt(n)
```

If the Oscar integer is too large to fit in an `Int`, an exception is raised.

```@repl oscar
Int(ZZ(12348732648732648763274868732687324))
```

The `fits` function can be used to determine whether an Oscar integer will fit
in a Julia `Int`:

```@repl oscar
fits(Int, ZZ(12348732648732648763274868732687324))
```

