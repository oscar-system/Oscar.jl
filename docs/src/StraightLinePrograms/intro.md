```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```


# Introduction

This is the documentation of the straight-line programs (SLP) implementation by
[Rafael Fourquet](https://github.com/rfourquet). Originally this was supposed
to become a separate Julia module, however it has now been incorporated into
the OSCAR core.

The main SLP type is `SLProgram`, to which other types can "compile" (or
"transpile"). The easiest way to create an `SLProgram` is to combine
"generators":

```jldoctest
julia> using Oscar;

julia> using Oscar.StraightLinePrograms; const SLP = Oscar.StraightLinePrograms;

julia> x, y, z = SLP.gens(SLProgram, 3)
3-element Vector{SLProgram{Union{}}}:
 x
 y
 z

julia> p = (x*y^2 + 1.3*z)^-1
#1 = ^   y  2  ==>  y^2
#2 = *   x #1  ==>  (x*y^2)
#3 = * 1.3  z  ==>  (1.3*z)
#4 = +  #2 #3  ==>  ((x*y^2) + (1.3*z))
#5 = ^  #4 -1  ==>  ((x*y^2) + (1.3*z))^-1
return: #5
```

On the right side of the above output is the representation of the computation
so far. It's done via another SLP type (tentatively) called `Lazy` which
represent "formulas" as trees:

```jldoctest
julia> using Oscar;

julia> using Oscar.StraightLinePrograms; const SLP = Oscar.StraightLinePrograms;

julia> x, y, z = SLP.gens(SLProgram, 3);

julia> p = (x*y^2 + 1.3*z)^-1;

julia> X, Y, Z = SLP.gens(Lazy, 3)
3-element Vector{Lazy}:
 x
 y
 z

julia> q = (X*Y^2 + 1.3*Z)^-1
((x*y^2) + (1.3*z))^-1

julia> f = SLP.evaluate(p, [X, Y, Z])
((x*y^2) + (1.3*z))^-1

julia> SLP.evaluate(f, [X, Y, Z]) == f
true

julia> SLP.evaluate(p, Any[x, y, z]) == p
true

julia> dump(q) # q::Lazy is a tree
Oscar.StraightLinePrograms.Lazy
  x: Oscar.StraightLinePrograms.Exp
    p: Oscar.StraightLinePrograms.Plus
      xs: Array{Oscar.StraightLinePrograms.LazyRec}((2,))
        1: Oscar.StraightLinePrograms.Times
          xs: Array{Oscar.StraightLinePrograms.LazyRec}((2,))
            1: Oscar.StraightLinePrograms.Input
              n: Int64 1
            2: Oscar.StraightLinePrograms.Exp
              p: Oscar.StraightLinePrograms.Input
                n: Int64 2
              e: Int64 2
        2: Oscar.StraightLinePrograms.Times
          xs: Array{Oscar.StraightLinePrograms.LazyRec}((2,))
            1: Oscar.StraightLinePrograms.Const{Float64}
              c: Float64 1.3
            2: Oscar.StraightLinePrograms.Input
              n: Int64 3
    e: Int64 -1
  gens: Array{Symbol}((3,))
    1: Symbol x
    2: Symbol y
    3: Symbol z
```

Evaluation of SLPs is done via `evaluate`, which can take a vector of
anything which supports the operations used in the SLP (e.g. `*`, `+` and `^`
in this example; `-` is also supported but division not yet).
Note that currently, the `eltype` of the input vector for `SLProgram`
must be a supertype of any intermediate computation (so it's always safe to
pass a `Vector{Any}`).

```jldoctest
julia> using Oscar;

julia> using Oscar.StraightLinePrograms; const SLP = Oscar.StraightLinePrograms;

julia> x, y, z = SLP.gens(SLProgram, 3);

julia> p = (x*y^2 + 1.3*z)^-1;

julia> X, Y, Z = SLP.gens(Lazy, 3);


julia> SLP.evaluate(p, [2.0, 3.0, 5.0])
0.04081632653061224

julia> SLP.evaluate(X*Y^2, ['a', 'b'])
"abb"
```

# Returning multiple values

There is a low-level interface to return multiple values from an `SLProgram`;
for example, to return the second and last intermediate values from `p`
above, we would "assign" these values to positions `#1` and `#2`,
delete all other positions (via the "keep" operation), and return the
resulting array (the one used for intermediate computations):

```jldoctest
julia> using Oscar;

julia> using Oscar.StraightLinePrograms; const SLP = Oscar.StraightLinePrograms;

julia> x, y, z = SLP.gens(SLProgram, 3);

julia> p = (x*y^2 + 1.3*z)^-1;

julia> X, Y, Z = SLP.gens(Lazy, 3);


julia> SLP.pushop!(p, SLP.assign, SLP.Arg(2), SLP.Arg(1))
       SLP.pushop!(p, SLP.assign, SLP.Arg(5), SLP.Arg(2))
       SLP.pushop!(p, SLP.keep, SLP.Arg(2))
       SLP.setmultireturn!(p)
#1 = ^   y  2  ==>  y^2
#2 = *   x #1  ==>  (x*y^2)
#3 = * 1.3  z  ==>  (1.3*z)
#4 = +  #2 #3  ==>  ((x*y^2) + (1.3*z))
#5 = ^  #4 -1  ==>  ((x*y^2) + (1.3*z))^-1
#1 =    #2     ==>  (x*y^2)
#2 =    #5     ==>  ((x*y^2) + (1.3*z))^-1
keep: #1..#2
return: [#1, #2]

julia> SLP.evaluate(p, [X, Y, Z])
list([(x*y^2), ((x*y^2) + (1.3*z))^-1])
```

## Straight line decisions

A "decision" is a special operation which allows to stop prematurely the
execution of the program if a condition is `false`, and the program returns
`true` if no condition failed.
Currently, the interface is modeled after GAP's SLPs and defaults to testing
the `AbstractAlgebra.order` of an element. More specifically,
`test(prg, n::Integer)` tests whether the order of the result of `prg` is
equal to `n`, and `dec1 & dec2` chains two programs with a short-circuiting
behavior:

```julia
julia> p1 = SLP.test(x*y^2, 2)
#1 = ^ y  2  ==>  y^2
#2 = * x #1  ==>  (xy^2)
test: order(#2) == 2 || return false
return: true

julia> p2 = SLP.test(y, 4)
test: order(y) == 4 || return false
return: true

julia> p1 & p2
#1 = ^ y  2  ==>  y^2
#2 = * x #1  ==>  (xy^2)
test: order(#2) == 2 || return false
test: order(y) == 4 || return false
return: true

julia> SLP.evaluate(p1 & p2, [X, Y])
test((xy^2), 2) & test(y, 4)

julia> using AbstractAlgebra; perm1, perm2 = perm"(1, 4)", perm"(1, 3, 4, 2)";

julia> SLP.evaluate(p1 & p2, [perm1, perm2])
true

julia> SLP.evaluate(p1 & p2, [perm2, perm1])
false
```


