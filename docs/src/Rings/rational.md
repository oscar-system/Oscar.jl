```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["rational.md"]
```

# Rationals

Fractions are created in Julia with the double slash operator `//`. If a
fraction is created from Julia integers, a Julia fraction results, and if either
the numerator or denominator is an Oscar integer of type `fmpz`, an Oscar
fraction of type `fmpq` results.

Julia has its own parameterised type `Rational{T}` for its own fractions, where
`T` is the integer type of the numerator and denominator, e.g.
`Rational{Int}` and `Rational{BigInt}`. Unlike with `Int`, all of the basic
arithmetic operations on Julia's `Rational{Int}` are checked for overflow in
the numerator and denominator.

## The field of rationals

The parent of an Oscar rational number is the field of rationals. It can be
constructed from the ring of integers `ZZ` using the `FractionField`
constructor.

For convenience, `QQ` is already defined to be the field of rational numbers.

```@repl oscar
S = FractionField(ZZ)
QQ
```

### Integer constructors

Oscar rationals can be created using `QQ`. Two arguments can be passed to
specify numerator and denominator. If a single argument is passed, the
denominator is set to `1`.

For convenience, `QQ` also accepts Julia integers and rationals, but will
always construct an Oscar rational.

Naturally, Julia's double slash operator can also be used to construct
fractions. However, unlike `QQ`, the double slash operator only constructs an
Oscar rational if either the numerator or denominator is an Oscar integer.

An exception is raised if a fraction is constructed with denominator zero.

```@repl oscar
QQ(1, 2)
QQ(5)
ZZ(3)//5
1//ZZ(7)
QQ(2//3)
ZZ(3)//0
```
One can also construct the rational number ``0`` with the empty constructor:

```@repl oscar
QQ()
```

The following special constructors are also provided:

* `zero(QQ)`
* `one(QQ)`

```@repl oscar
zero(QQ)
one(QQ)
```

## Predicates

* `iszero(n::fmpq) -> Bool`
* `isone(n::fmpq) -> Bool`
* `is_unit(n::fmpq) -> Bool`

The `is_unit` function will return `true` iff ``n \neq 0``.

```@repl oscar
iszero(QQ())
isone(one(QQ))
is_unit(QQ(-2, 3))
```

## Properties

* `numerator(n::fmpq) -> fmpz`
* `denominator(n::fmpq) -> fmpz`

Return the numerator and denominator respectively, of $n$.

* `sign(n::fmpq) -> fmpq`

Return the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```@repl oscar
sign(QQ(2, 3))
sign(QQ())
sign(QQ(-1))
```

* `abs(n::fmpq) -> fmpq`

Return the absolute value of ``n``, i.e. ``n`` if ``n \geq 0`` and ``-n``
otherwise.


```@repl oscar
abs(QQ(-3, 2))
```

* `height(n::fmpq) -> fmpz`

Return the maximum of the absolute values of the numerator and denominator of
$n$.

```@repl oscar
height(QQ(324987329, -8372492324))
```

* `floor(n::fmpq) -> fmpq`

Return the greatest integer $m$ (as a rational number) such that $m \leq n$.

* `ceil(n::fmpq) -> fmpq`

Return the least integer $m$ (as a rational number) such that $m \geq n$.

```@repl oscar
floor(QQ(-2, 3))
ceil(QQ(7, 2))
typeof(ans)
ceil(QQ(5))
```
* `floor(fmpz, n::fmpq) -> fmpz`

Return the greatest integer $m$ such that $m \leq n$.

* `ceil(fmpz, n::fmpq) -> fmpz`

Return the least integer $m$ such that $m \geq n$.

```@repl oscar
floor(fmpz, QQ(-2, 3))
ceil(fmpz, QQ(7, 2))
typeof(ans)
ceil(fmpz, QQ(5))
```

## Basic arithmetic

Oscar provides the basic arithmetic operations `+`, `-` and `*` and comparison
operators `==`, `!=`, `<`, `<=`, `>`, `>=`, including mixed operations between
Julia and Oscar rationals and integers.

### [Exact Division]

* `divexact(a::fmpq, b::fmpq) -> fmpq`
* `divexact(a::fmpq, b::Union{fmpz,Base.Integer,Base.Rational}) -> fmpq`
* `divexact(a::Union{fmpz,Base.Integer,Base.Rational}, b::fmpq) -> fmpq`

Return the quotient of ``a`` by ``b``. Exact division raises an exception if
division by zero is attempted.

```@repl oscar
divexact(QQ(2, 3), QQ(3, 5))
divexact(QQ(1, 3), ZZ(0))
divexact(QQ(3, 4), ZZ(5))
divexact(ZZ(6), QQ(2, 3))
divexact(QQ(1, 3), 5)
```

### Powering

* `^(a::fmpq, b::Int) -> fmpq`

Return the result of powering ``a`` by ``b``.

```@repl oscar
QQ(5, 7)^32
QQ(1, 2)^(-2)
```

The following is allowed for convenience.

```@repl oscar
QQ(0)^0
```

!!! note
    In Julia, the rational number ``0//1`` when raised to a negative power
    returns ``1//0`` to indicate that the value is undefined. Oscar raises
    an exception.

```@repl oscar
QQ(0)^-2
```

* `is_power(a::fmpq, b::Int) -> Bool, fmpq`

Test if ``a`` is an ``n``-th power. If so, return ```true``` and the root,
```false``` and any rational otherwise.

* `is_power(a::fmpq) -> Int, fmpq`

Find the largest ``n`` such that ``a`` is an ``n``-th power. Return ``n`` and the root.

* `root(a::fmpq, b::Int) -> fmpq`

Compute an ``n``-th root of ``a``, raises an error if ``a`` is not an ``n``-th power.

```@repl oscar
is_power(QQ(8), 3)
is_power(QQ(8), 2)
is_power(QQ(9//16))
root(QQ(25//9), 2)
```

