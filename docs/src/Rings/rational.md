```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Rationals

Fractions are created in Julia with the double slash operator `//`. If a
fraction is created from Julia integers, a Julia fraction results, and if either
the numerator or denominator is an OSCAR integer of type `ZZRingElem`, an OSCAR
fraction of type `QQFieldElem` results.

Julia has its own parameterised type `Rational{T}` for its own fractions, where
`T` is the integer type of the numerator and denominator, e.g.
`Rational{Int}` and `Rational{BigInt}`. Unlike with `Int`, all of the basic
arithmetic operations on Julia's `Rational{Int}` are checked for overflow in
the numerator and denominator.

## The field of rationals

The parent of an OSCAR rational number is the field of rationals. It can be
constructed from the ring of integers `ZZ` using the `fraction_field`
constructor.

For convenience, `QQ` is already defined to be the field of rational numbers.

```jldoctest
julia> S = fraction_field(ZZ)
Rational field

julia> QQ
Rational field

```

### Integer constructors

OSCAR rationals can be created using `QQ`. Two arguments can be passed to
specify numerator and denominator. If a single argument is passed, the
denominator is set to `1`.

For convenience, `QQ` also accepts Julia integers and rationals, but will
always construct an OSCAR rational.

Naturally, Julia's double slash operator can also be used to construct
fractions. However, unlike `QQ`, the double slash operator only constructs an
OSCAR rational if either the numerator or denominator is an OSCAR integer.

An exception is raised if a fraction is constructed with denominator zero.

```jldoctest
julia> QQ(1, 2)
1//2

julia> QQ(5)
5

julia> ZZ(3)//5
3//5

julia> 1//ZZ(7)
1//7

julia> QQ(2//3)
2//3

julia> ZZ(3)//0
ERROR: DivideError: integer division error
[...]

```
One can also construct the rational number ``0`` with the empty constructor:

```jldoctest
julia> QQ()
0

```

The following special constructors are also provided:

* `zero(QQ)`
* `one(QQ)`

```jldoctest
julia> zero(QQ)
0

julia> one(QQ)
1

```

## Predicates

* `iszero(n::QQFieldElem) -> Bool`
* `isone(n::QQFieldElem) -> Bool`
* `is_unit(n::QQFieldElem) -> Bool`

The `is_unit` function will return `true` iff ``n \neq 0``.

```jldoctest
julia> iszero(QQ())
true

julia> isone(one(QQ))
true

julia> is_unit(QQ(-2, 3))
true

```

## Properties

* `numerator(n::QQFieldElem) -> ZZRingElem`
* `denominator(n::QQFieldElem) -> ZZRingElem`

Return the numerator and denominator respectively, of $n$.

* `sign(n::QQFieldElem) -> QQFieldElem`

Return the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```jldoctest
julia> sign(QQ(2, 3))
1

julia> sign(QQ())
0

julia> sign(QQ(-1))
-1

```

* `abs(n::QQFieldElem) -> QQFieldElem`

Return the absolute value of ``n``, i.e. ``n`` if ``n \geq 0`` and ``-n``
otherwise.


```jldoctest
julia> abs(QQ(-3, 2))
3//2

```

* `height(n::QQFieldElem) -> ZZRingElem`

Return the maximum of the absolute values of the numerator and denominator of
$n$.

```jldoctest
julia> height(QQ(324987329, -8372492324))
8372492324

```

* `floor(n::QQFieldElem) -> QQFieldElem`

Return the greatest integer $m$ (as a rational number) such that $m \leq n$.

* `ceil(n::QQFieldElem) -> QQFieldElem`

Return the least integer $m$ (as a rational number) such that $m \geq n$.

```jldoctest
julia> floor(QQ(-2, 3))
-1

julia> ceil(QQ(7, 2))
4

julia> typeof(ans)
QQFieldElem

julia> ceil(QQ(5))
5

```
* `floor(ZZRingElem, n::QQFieldElem) -> ZZRingElem`

Return the greatest integer $m$ such that $m \leq n$.

* `ceil(ZZRingElem, n::QQFieldElem) -> ZZRingElem`

Return the least integer $m$ such that $m \geq n$.

```jldoctest
julia> floor(ZZRingElem, QQ(-2, 3))
-1

julia> ceil(ZZRingElem, QQ(7, 2))
4

julia> typeof(ans)
ZZRingElem

julia> ceil(ZZRingElem, QQ(5))
5

```

## Basic arithmetic

OSCAR provides the basic arithmetic operations `+`, `-` and `*` and comparison
operators `==`, `!=`, `<`, `<=`, `>`, `>=`, including mixed operations between
Julia and OSCAR rationals and integers.

### [Exact Division]

* `divexact(a::QQFieldElem, b::QQFieldElem) -> QQFieldElem`
* `divexact(a::QQFieldElem, b::Union{ZZRingElem,Base.Integer,Base.Rational}) -> QQFieldElem`
* `divexact(a::Union{ZZRingElem,Base.Integer,Base.Rational}, b::QQFieldElem) -> QQFieldElem`

Return the quotient of ``a`` by ``b``. Exact division raises an exception if
division by zero is attempted.

```jldoctest
julia> divexact(QQ(2, 3), QQ(3, 5))
10//9

julia> divexact(QQ(1, 3), ZZ(0))
ERROR: DivideError: integer division error
[...]

julia> divexact(QQ(3, 4), ZZ(5))
3//20

julia> divexact(ZZ(6), QQ(2, 3))
9

julia> divexact(QQ(1, 3), 5)
1//15

```

### Powering

* `^(a::QQFieldElem, b::Int) -> QQFieldElem`

Return the result of powering ``a`` by ``b``.

```jldoctest
julia> QQ(5, 7)^32
23283064365386962890625//1104427674243920646305299201

julia> QQ(1, 2)^(-2)
4

```

The following is allowed for convenience.

```jldoctest
julia> QQ(0)^0
1

```

!!! note
    In Julia, the rational number ``0//1`` when raised to a negative power
    returns ``1//0`` to indicate that the value is undefined. OSCAR raises
    an exception.

```jldoctest
julia> QQ(0)^-2
ERROR: DivideError: integer division error
[...]

```

* `is_power(a::QQFieldElem, b::Int) -> Bool, QQFieldElem`

Test if ``a`` is an ``n``-th power. If so, return `true` and the root,
`false` and any rational otherwise.

* `is_perfect_power_with_data(a::QQFieldElem) -> Int, QQFieldElem`

Find the largest ``n`` such that ``a`` is an ``n``-th power. Return ``n`` and the root.

* `root(a::QQFieldElem, b::Int) -> QQFieldElem`

Compute an ``n``-th root of ``a``, raises an error if ``a`` is not an ``n``-th power.

```jldoctest
julia> is_power(QQ(8), 3)
(true, 2)

julia> is_power(QQ(8), 2)
(false, 8)

julia> is_perfect_power_with_data(QQ(9//16))
(2, 3//4)

julia> root(QQ(25//9), 2)
5//3

```

