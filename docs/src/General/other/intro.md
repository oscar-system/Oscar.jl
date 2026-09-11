# Notes for users of other computer algebra systems

OSCAR is written in Julia, so an OSCAR session is just a Julia session.
Many of the differences described here are therefore differences in
behaviour and syntax between Julia and the CAS you are used to.
This page collects differences that affect users of all such systems.
Notes for users of specific systems follow on separate pages:

- [Notes for GAP users](@ref)
- [Notes for Magma users](@ref)
- [Notes for SageMath users](@ref)
- [Notes for Singular users](@ref)
- [Notes for polymake users](@ref)

!!! note "Help wanted"
    These pages are far from complete.
    If you are missing something here, or if a difference between OSCAR and
    the system you know took you a while to understand,
    please tell us about it, for example by opening an
    [issue on GitHub](https://github.com/oscar-system/Oscar.jl/issues)
    or on [Slack](https://oscar-system.org/slack).
    Contributions to these pages are very welcome.

## [Integers and rational numbers](@id other_integers)

An integer literal such as `2` is a 64 bit machine integer in Julia.
Arithmetic with machine integers silently wraps around on overflow,
and some functions throw an error instead.
```jldoctest
julia> 2^100
0

julia> factorial(21)
ERROR: OverflowError: 21 is too large to look up in the table; consider using `factorial(big(21))` instead
[...]
```
Use OSCAR's integers (elements of `ZZ`) whenever the numbers may get large.
```jldoctest
julia> ZZ(2)^100
1267650600228229401496703205376

julia> factorial(ZZ(21))
51090942171709440000
```
Integers returned by OSCAR functions are of this type already, only literals
that you type yourself are affected.

The quotient of two integers is a floating point number in Julia,
and `//` creates rational numbers.
OSCAR makes a subtle but important
[distinction between `/` and `//`](@ref subtle_distinction_for_rings).
```jldoctest
julia> 3/4
0.75

julia> QQ(3, 4)
3//4
```

## Every object has a parent

Each element of an algebraic structure knows its parent,
and `parent(x)` returns it.
Operations that involve several elements usually require that the
parents coincide.
Some conversions happen automatically, for example integers and
rational numbers can be combined with elements of most rings:
for any ring ``R`` there is exactly one ring homomorphism from the
integers into ``R``, and it extends uniquely to a partial map from the
rational numbers, defined wherever the denominator is invertible in
``R``. This is what `R(5)` and `R(7//2)` compute, and such an automatic
conversion is called a coercion.
In other cases, an element has to be moved into the required structure
explicitly, by calling the parent like a function.
```jldoctest
julia> ZZ(7) + QQ(3, 2)
17//2

julia> F = GF(7);

julia> x = F(3)
3

julia> parent(x)
Prime field of characteristic 7

julia> x + 5
1
```

## Many constructors return more than one object

Functions that construct a substructure or a quotient usually also
return the map that connects it with the original structure,
and functions that construct rings with generators also return the
generators.
In such cases, the result is a tuple, and the individual parts of it
can be assigned to variables in one go.
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> G = symmetric_group(4);

julia> H, emb = sub(G, [cperm(G, [1, 2, 3])]);

julia> H
Permutation group of degree 4

julia> emb
Group homomorphism
  from permutation group of degree 4
  to symmetric group of degree 4
```
Ignore parts that you do not need by assigning them to `_`,
or pick a part by indexing.
```jldoctest
julia> G = symmetric_group(4);

julia> Z, _ = center(G);

julia> A = derived_subgroup(G)[1]
Alternating group of degree 4
```
The same applies for example to `quo`, `kernel`, `image`, `number_field`,
and `residue_ring`.

## Semicolons and output

In an interactive session, Julia prints the value of each entered
expression, and a trailing semicolon suppresses this output.

Inside a Julia function, nothing is printed unless you call `println` or
use the `@show` macro.

## Indexing starts at 1

Like in GAP, Magma, and Singular, and unlike in Python (and hence in
SageMath) and polymake, the first entry of a list has index 1.
Ranges include both endpoints: `1:5` consists of the numbers 1 to 5.

## Names of functions

OSCAR functions have `snake_case` names, see the section on
[Naming conventions](@ref).
A predicate has a name that starts with `is_`, for example `is_prime`,
and the name of a function that modifies one of its arguments ends with
an exclamation mark, for example `push!`.

Functions are not attached to objects, that is, one writes `order(G)`
rather than `G.order()`.
Instead, Julia chooses the method to be called according to the types of
*all* arguments (this is called multiple dispatch),
and the same function name is used for analogous operations on different
kinds of objects, for example `order` for groups and group elements,
`degree` for polynomials, field extensions, and permutation groups,
and `length` for lists and other collections.

Optional arguments are usually given by name, after a semicolon,
for example `kernel(M; side = :right)`.

Some hints for finding the function you are looking for:

- Guessing the `snake_case` version of the name that you know from the
  system you are used to often works.
- Type the beginning of a name and hit the tab key to see all completions.
- `?name` shows the documentation of the function `name`,
  and `methods(name)` shows for which types of arguments it is defined.
- `methodswith(typeof(x); supertypes = true)` lists all functions that
  have a method for objects of the same type as `x`.
- `apropos("text")` lists all functions whose documentation mentions
  the given text.
- The search field at the top of this manual searches the whole manual.

## Interactive sessions

- When an error occurs, or when you hit ctrl-C, Julia returns to the prompt.
  There is no break loop, in contrast to GAP and Magma.
- Measure the runtime of a computation with the `@time` macro:
  `@time f(x)`.
- Read a file with Julia code with `include("file.jl")`.
- Save OSCAR objects to a file with `save("file.mrdi", x)` and read
  them back with `load("file.mrdi")`, see [Serialization](@ref).
