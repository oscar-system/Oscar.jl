```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["affine.md"]
```

# Affine Rings

general bla explaining what this is...

## Quotient Rings

```@repl oscar
R, x = PolynomialRing(QQ, :x=>1:3)
Q, mQ = quo(R, ideal([x[1]^2-2, x[2]^2-3, x[3]^2-5]))
```

```@docs
quo(R::MPolyRing, I::MPolyIdeal)
```

### Elements

Elements are created by moving polynomials into the ring:

```@repl oscar
R, x = PolynomialRing(QQ, :x=>1:3)
Q, mQ = quo(R, ideal([x[1]^2-2, x[2]^2-3, x[3]^2-5]))
Q(x[2]+x[3])
Q(x[2])
ans^2
ans == 3
```

As shown: elements are not automatically simplified.

```@docs
simplify(::MPolyQuoElem)
Oscar.simplify!(::MPolyQuoElem)
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

## Predicates

* `iszero(n::Oscar.Integer) -> Bool`
* `isone(n::Oscar.Integer) -> Bool`
* `isprobable_prime(n::Oscar.Integer) -> Bool`


## Properties

* `sign(n::Oscar.Integer) -> Oscar.Integer`

Returns the sign of `n`, i.e. ``n/|n|`` if ``n \neq 0``, or ``0`` otherwise.

```@repl oscar
sign(ZZ(23))
sign(ZZ(0))
sign(ZZ(-1))
```

## Basic arithmetic

Oscar provides the basic arithmetic operations `+`, `-` and `*` and comparison
operators `==`, `!=`, including mixed operations between
Julia and Oscar integers. It also provides division and powering as described
below.
