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
fraction is created from Julia integers, a Julia fraction results and if either
the numerator or denominator is an Oscar integer, an Oscar fraction results.

Oscar has just a single rational fraction type `Oscar.fmpq`. The Oscar fraction
type belongs to `Oscar.Rational`.

Julia has a parameterised type `Base.Rational{T}` for its own fractions, where
`T` is the integer type the fractions are constructed from, e.g.
`Base.Rational{Int}` and `Base.Rational{BigInt}`.

Note that conveniently, both of these Julia rational types belong to
`Base.Rational`. Therefore in the description below we refer to
`Oscar.Rational` for Oscar rational numbers and `Base.Rational` for any of the
Julia rational number types.

The situation is described by the following diagram.

![alt text](../img/RationalTypes.svg)

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

* `iszero(n::Oscar.Rational) -> Bool`
* `isone(n::Oscar.Rational) -> Bool`
* `isunit(n::Oscar.Rational) -> Bool`

The `isunit` function will return `true` iff ``n \neq 0``.

```@repl oscar
iszero(QQ())
isone(one(QQ))
isunit(QQ(-2, 3))
```

## Properties

* `numerator(n::Oscar.Rational) -> Oscar.Integer`
* `denominator(n::Oscar.Rational) -> Oscar.Integer`

Return the numerator and denominator respectively, of $n$.

* `sign(n::Oscar.Rational) -> Oscar.Rational`

Returns the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```@repl oscar
sign(QQ(2, 3))
sign(QQ())
sign(QQ(-1))
```

* `abs(n::Oscar.Rational) -> Oscar.Rational`

Return the absolute value of ``n``, i.e. ``n`` if ``n \geq 0`` and ``-n``
otherwise.


```@repl oscar
abs(QQ(-3, 2))
```

* `height(n::Oscar.Rational) -> Oscar.Integer`

Return the maximum of the absolute values of the numerator and denominator of
$n$.

```@repl oscar
height(QQ(324987329, -8372492324))
```

* `floor(n::Oscar.Rational) -> Oscar.Rational`

Return the greatest integer $m$ such that $m \leq n$.

* `ceil(n::Oscar.Rational) -> Oscar.Rational`

Return the least integer $m$ such that $m \geq n$.

```@repl oscar
floor(QQ(-2, 3))
ceil(QQ(7, 2))
ceil(QQ(5))
```

## Basic arithmetic

Oscar provides the basic arithmetic operations `+`, `-` and `*` and comparison
operators `==`, `!=`, `<`, `<=`, `>`, `>=`, including mixed operations between
Julia and Oscar rationals and integers.

### [Exact Division]

* `divexact(a::Oscar.Rational, b::Oscar.Rational) -> Oscar.Rational`
* `divexact(a::Oscar.Rational, b::Oscar.Integer) -> Oscar.Rational`
* `divexact(a::Oscar.Integer, b::Oscar.Rational) -> Oscar.Rational`

In the first signature, one of the arguments may be a Julia rational and in the
other two signatures the integers may be Julia integers.

Returns the quotient of ``a`` by ``b``. Exact division raises an exception if
division by zero is attempted.

```@repl oscar
divexact(QQ(2, 3), QQ(3, 5))
divexact(QQ(1, 3), ZZ(0))
divexact(QQ(3, 4), ZZ(5))
divexact(ZZ(6), QQ(2, 3))
divexact(QQ(1, 3), 5)
```

### Powering

* `^(a::Oscar.Rational, b::Base.Int) -> Oscar.Rational`

Powering can be accomplished naturally using the special caret infix
operator:

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
    returns ``1/0`` to indicate that the value is undefined. Oscar follows
    this convention.

```@repl oscar
QQ(0)^-2
```

