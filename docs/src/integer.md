```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Integers

An important design decision in Oscar.jl is to use Julia as the user language
by default. This means that integers typed at the
[REPL](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop)
are [Julia integers](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/).

### Constructors

For performance reasons, Oscar has its own integer format. These are
entered using the `ZZ` constructor.

```@repl oscar
a = ZZ(2)^100
```

### Julia integers in Oscar functions

For convenience, many Oscar functions also accept Julia integers as inputs by
converting them to Oscar integers. For example:

```@repl oscar
divexact(ZZ(234), 2)
```

In this example, `2` is a Julia integer but is still valid in the
call to `divexact` because the first argument is an Oscar integer.

In general, Oscar can only automatically convert from Julia integers to Oscar
integers if they are combined with other Oscar objects or passed to Oscar
functions.

### Limitations

Oscar integers have the same limitations as [GMP](https://gmplib.org/)
multiprecision integers, namely that they are limited by the available memory
on the machine and in any case to signed integers whose absolute value does not
exceed ``2^{37}`` binary bits.

In the following, unless stated otherwise, when we refer to integers, we mean
Oscar integers. When we refer to an `Int` we mean the Julia `Int`.

!!! note
    The Julia 'Int' type is either a 32 or 64 bit integer, depending on the
    machine architecture (usually 64 bits on most modern machines). The range of
    values is machine dependent, but can be found by typing 'typemin(Int)' and
    'typemax(Int)' in Julia.

## Basic arithmetic

Oscar provides the basic arithmetic operations `+`, `-` and `*`, and division
and powering as described below.

### Division in Oscar

Oscar distinguishes a number of different kinds of division:

* Exact division (`divexact`)
* Euclidean division (`div`, `divrem`, `mod`, `rem`)
* Construction of fractions (`a//b`)
* Floating point division (`a/b`)
* Divisibility testing (`divides`)

These choices have been made for maximum parsimony with the Julia language.

!!! note
    It is a common error to enter '1/2' for the fraction 'one half' in Julia.
    In the Julia language, this expression is reserved for floating point
    division. Instead, the double slash operator is used for fractions.

### Exact Division

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
    An exception will be raised if an integer outside the range [-1, 1] is
    raised to a negative exponent.

The following is allowed for convenience.

```@repl oscar
ZZ(0)^0
```

!!! note
    In Julia, '2^64' will return 0, as the Julia integer 2 is a machine word.
    In Oscar, the expression 'ZZ(2)^64' will return the expected result.

