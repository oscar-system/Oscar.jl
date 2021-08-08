[![Build Status](https://travis-ci.org/rfourquet/StraightLinePrograms.jl.svg?branch=master)](https://travis-ci.org/rfourquet/StraightLinePrograms.jl)

# StraightLinePrograms

This is a WIP implementation of straight-line programs (SLP)
This is part of the [Oscar](https://oscar.computeralgebra.de/) project.


## Introduction

The main SLP type is `SLProgram`, to which other types can "compile" (or
"transpile"). The easiest way to create an `SLProgram` is to combine
"generators":

```julia
julia> using StraightLinePrograms; const SL = StraightLinePrograms

julia> x, y, z = gens(SLProgram, 3)
3-element Vector{SLProgram{Union{}}}:
 x
 y
 z

julia> p = (x*y^2 + 1.3*z)^-1
#1 = ^   y  2  ==>  y^2
#2 = *   x #1  ==>  (xy^2)
#3 = * 1.3  z  ==>  (1.3z)
#4 = +  #2 #3  ==>  ((xy^2) + (1.3z))
#5 = ^  #4 -1  ==>  ((xy^2) + (1.3z))^-1
return: #5
```

On the right side of the above output is the representation of the computation
so far. It's done via another SLP type (tentatively) called `Free` which
represent "formulas" as trees:

```julia
julia> X, Y, Z = gens(Free, 3)
3-element Vector{Free}:
 x
 y
 z

julia> q = (X*Y^2 + 1.3*Z)^-1
((xy^2) + (1.3z))^-1

julia> f = evaluate(p, [X, Y, Z])
((xy^2) + (1.3z))^-1

julia> evaluate(f, [X, Y, Z]) == f
true

julia> evaluate(p, Any[x, y, z]) == p
true

julia> dump(q) # q::Free is a tree
Free
  x: StraightLinePrograms.Exp
    p: StraightLinePrograms.Plus
      xs: Array{StraightLinePrograms.Lazy}((2,))
        1: StraightLinePrograms.Times
          xs: Array{StraightLinePrograms.Lazy}((2,))
            1: StraightLinePrograms.Input
              n: Int64 1
            2: StraightLinePrograms.Exp
              p: StraightLinePrograms.Input
                n: Int64 2
              e: Int64 2
        2: StraightLinePrograms.Times
          xs: Array{StraightLinePrograms.Lazy}((2,))
            1: StraightLinePrograms.Const{Float64}
              c: Float64 1.3
            2: StraightLinePrograms.Input
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

```julia
julia> evaluate(p, [2.0, 3.0, 5.0])
0.04081632653061224

julia> evaluate(X*Y^2, ['a', 'b'])
"abb"
```

## Returning multiple values

There is a low-level interface to return multiple values from an `SLProgram`;
for example, to return the second and last intermediate values from `p`
above, we would "assign" these values to positions `#1` and `#2`,
delete all other positions (via the "keep" operation), and return the
resulting array (the one used for intermediate computations):

```julia
julia> SL.pushop!(p, SL.assign, SL.Arg(2), SL.Arg(1))
       SL.pushop!(p, SL.assign, SL.Arg(5), SL.Arg(2))
       SL.pushop!(p, SL.keep, SL.Arg(2))
       SL.setmultireturn!(p)
#1 = ^   y  2  ==>  y^2
#2 = *   x #1  ==>  (xy^2)
#3 = * 1.3  z  ==>  (1.3z)
#4 = +  #2 #3  ==>  ((xy^2) + (1.3z))
#5 = ^  #4 -1  ==>  ((xy^2) + (1.3z))^-1
#1 =    #2     ==>  (xy^2)
#2 =    #5     ==>  ((xy^2) + (1.3z))^-1
keep: #1..#2
return: [#1, #2]

julia> evaluate(p, [X, Y, Z])
2-element Vector{Free}:
 (xy^2)
 ((xy^2) + (1.3z))^-1
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
julia> p1 = SL.test(x*y^2, 2)
#1 = ^ y  2  ==>  y^2
#2 = * x #1  ==>  (xy^2)
test: order(#2) == 2 || return false
return: true

julia> p2 = SL.test(y, 4)
test: order(y) == 4 || return false
return: true

julia> p1 & p2
#1 = ^ y  2  ==>  y^2
#2 = * x #1  ==>  (xy^2)
test: order(#2) == 2 || return false
test: order(y) == 4 || return false
return: true

julia> evaluate(p1 & p2, [X, Y])
test((xy^2), 2) & test(y, 4)

julia> using AbstractAlgebra; perm1, perm2 = perm"(1, 4)", perm"(1, 3, 4, 2)";

julia> evaluate(p1 & p2, [perm1, perm2])
true

julia> evaluate(p1 & p2, [perm2, perm1])
false
```

## GAP's SLPs

There are two other available SLP types: `GAPSLProgram` and `AtlasSLProgram`,
and related `GAPSLDecision` and `AtlasSLDecision`, which are constructed
similarly as in GAP:

```julia
julia> prg = GAPSLProgram( [ [1,2,2,3], [3,-1] ], 2 )
# input:
r = [ g1, g2 ]
# program:
r[3] = r[1]^2*r[2]^3
r[4] = r[3]^-1
# return value:
r[4]

julia> evaluate(prg, [perm1, perm2])
(1,3,4,2)

julia> evaluate(prg, [x, y])
#1 = ^  x  2  ==>  x^2
#2 = ^  y  3  ==>  y^3
#3 = * #1 #2  ==>  (x^2y^3)
#4 = ^ #3 -1  ==>  (x^2y^3)^-1
return: #4

julia> SLProgram(prg) # direct compilation (with room for optimizations obviously)
#1 =    x     ==>  x
#2 =    y     ==>  y
#3 = ^ #1  2  ==>  x^2
#4 = ^ #2  3  ==>  y^3
#5 = * #3 #4  ==>  (x^2y^3)
#3 =   #5     ==>  (x^2y^3)
keep: #1..#3
#4 = ^ #3 -1  ==>  (x^2y^3)^-1
keep: #1..#4
return: #4

julia> GAPSLProgram( [ [2,3], [ [3,1,1,4], [1,2,3,1] ] ], 2 )
# input:
r = [ g1, g2 ]
# program:
r[3] = r[2]^3
# return values:
[ r[3]*r[1]^4, r[1]^2*r[3] ]

julia> GAPSLDecision([ [ [ 1, 1, 2, 1 ], 3 ], [ "Order", 1, 2 ], [ "Order", 2, 3 ], [ "Order", 3, 5 ] ] )
# input:
r = [ g1, g2 ]
# program:
r[3] = r[1]*r[2]
order( r[1] ) == 2 || return false
order( r[2] ) == 3 || return false
order( r[3] ) == 5 || return false
# return value:
true

julia> SLProgram(ans)
#1 =    x     ==>  x
#2 =    y     ==>  y
#3 = * #1 #2  ==>  (xy)
keep: #1..#3
test: order(#1) == 2 || return false
test: order(#2) == 3 || return false
test: order(#3) == 5 || return false
return: true

julia> d = AtlasSLDecision("inp 2\nchor 1 2\nchor 2 3\nmu 1 2 3\nchor 3 5")
inp 2
chor 1 2
chor 2 3
mu 1 2 3
chor 3 5

376> evaluate(d, [perm1, perm2])
false

julia> GAPSLDecision(d)
# input:
r = [ g1, g2 ]
# program:
order( r[1] ) == 2 || return false
order( r[2] ) == 3 || return false
r[3] = r[1]*r[2]
order( r[3] ) == 5 || return false
# return value:
true

julia> SLProgram(d)
#1 =    x     ==>  x
#2 =    y     ==>  y
test: order(#1) == 2 || return false
test: order(#2) == 3 || return false
#3 = * #1 #2  ==>  (xy)
keep: #1..#3
test: order(#3) == 5 || return false
return: true
```

## AbstractAlgebra's polynomial interface

<details>
<summary>
This is the initial API of SLPs which hasn't been updated in a while
and might not work as-is with the current state of the package.
</summary>

Currently, SLPs have a polynomial interface (`SLPoly`).

## Examples

```julia
julia> using AbstractAlgebra, StraightLinePrograms, BenchmarkTools;

julia> S = SLPolyRing(zz, [:x, :y]); x, y = gens(S)
2-element Vector{SLPoly{Int64,SLPolyRing{Int64,AbstractAlgebra.Integers{Int64}}}}:
 x
 y

julia> p = 3 + 2x * y^2 # each line of the SLP is shown with current value
  #1 = * 2 x    ==>     (2x)
  #2 = ^ y 2    ==>     y^2
  #3 = * #1 #2  ==>     (2xy^2)
  #4 = + 3 #3   ==>     (3 + (2xy^2))

julia> p.cs # constants used in the program
2-element Vector{Int64}:
 3
 2

julia> p.lines # each line is a UInt64 encoding an opcode and 2 arguments
 Line(0x0500000028000001)
 Line(0x8780000020000002)
 Line(0x0500000030000004)
 Line(0x0300000010000005)

julia> evaluate(p, [2, 3])
39

julia> p2 = (p*(x*y))^6
#1 = *  2  x  ==>  (2x)
#2 = ^  y  2  ==>  y^2
#3 = * #1 #2  ==>  (2xy^2)
#4 = +  3 #3  ==>  (3 + (2xy^2))
#5 = *  x  y  ==>  (xy)
#6 = * #4 #5  ==>  ((3 + (2xy^2))xy)
#7 = ^ #6  6  ==>  ((3 + (2xy^2))xy)^6

julia> R, (x1, y1) = PolynomialRing(zz, ["x", "y"]); R
Multivariate Polynomial Ring in x, y over Integers

julia> q = convert(R, p2)
64*x^12*y^18+576*x^11*y^16+2160*x^10*y^14+4320*x^9*y^12+4860*x^8*y^10+2916*x^7*y^8+729*x^6*y^6

julia> v = [3, 5]; @btime evaluate($q, $v)
  32.101 μs (634 allocations: 45.45 KiB)
-1458502820125772303

julia> @btime evaluate($p2, $v)
  270.341 ns (4 allocations: 352 bytes)
-1458502820125772303

julia> res = Int[]; @btime StraightLinePrograms.evaluate!($res, $p2, $v)
  171.013 ns (0 allocations: 0 bytes)
-1458502820125772303

julia> res # intermediate computations (first 2 elements are constants)
9-element Vector{Int64}:
                    3
                    2
                    6
                   25
                  150
                  153
                   15
                 2295
 -1458502820125772303

julia> f2 = StraightLinePrograms.compile!(p2) # compile to machine code
#3 (generic function with 1 method)

julia> @btime evaluate($p2, $v)
  31.153 ns (1 allocation: 16 bytes)
-1458502820125772303

julia> @btime $f2($v) # use a function barrier for last bit of efficiency
 7.980 ns (0 allocations: 0 bytes)
-1458502820125772303

julia> q
64*x^12*y^18+576*x^11*y^16+2160*x^10*y^14+4320*x^9*y^12+4860*x^8*y^10+2916*x^7*y^8+729*x^6*y^6

julia> p3 = convert(S, q) # convert back q::Mpoly to an SLPoly
 #1 = ^    x   6  ==>  x^6
 #2 = ^    x   7  ==>  x^7
 #3 = ^    x   8  ==>  x^8
 #4 = ^    x   9  ==>  x^9
 #5 = ^    x  10  ==>  x^10
 #6 = ^    x  11  ==>  x^11
 #7 = ^    x  12  ==>  x^12
 #8 = ^    y   6  ==>  y^6
 #9 = ^    y   8  ==>  y^8
#10 = ^    y  10  ==>  y^10
#11 = ^    y  12  ==>  y^12
#12 = ^    y  14  ==>  y^14
#13 = ^    y  16  ==>  y^16
#14 = ^    y  18  ==>  y^18
#15 = *   64  #7  ==>  (64x^12)
#16 = *  #15 #14  ==>  (64x^12y^18)
#17 = *  576  #6  ==>  (576x^11)
#18 = *  #17 #13  ==>  (576x^11y^16)
#19 = +  #16 #18  ==>  ((64x^12y^18) + (576x^11y^16))
#20 = * 2160  #5  ==>  (2160x^10)
#21 = *  #20 #12  ==>  (2160x^10y^14)
#22 = +  #19 #21  ==>  ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14))
#23 = * 4320  #4  ==>  (4320x^9)
#24 = *  #23 #11  ==>  (4320x^9y^12)
#25 = +  #22 #24  ==>  ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12))
#26 = * 4860  #3  ==>  (4860x^8)
#27 = *  #26 #10  ==>  (4860x^8y^10)
#28 = +  #25 #27  ==>  ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12) + (4860x^8y^10))
#29 = * 2916  #2  ==>  (2916x^7)
#30 = *  #29  #9  ==>  (2916x^7y^8)
#31 = +  #28 #30  ==>  ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12) + (4860x^8y^10) + (2916x^7y^8))
#32 = *  729  #1  ==>  (729x^6)
#33 = *  #32  #8  ==>  (729x^6y^6)
#34 = +  #31 #33  ==>  ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12) + (4860x^8y^10) + (2916x^7y^8) + (729x^6y^6))

julia> @btime evaluate($p3, $v)
  699.465 ns (5 allocations: 1008 bytes)
-1458502820125772303

julia> @time f3 = StraightLinePrograms.compile!(p3);
  0.002830 seconds (1.40 k allocations: 90.930 KiB)

julia> @btime $f3($v)
 80.229 ns (0 allocations: 0 bytes)
-1458502820125772303

julia> p4 = convert(S, q, limit_exp=true); # use different encoding

julia> @btime evaluate($p4, $v)
  766.864 ns (5 allocations: 1008 bytes)
-1458502820125772303

julia> @time f4 = StraightLinePrograms.compile!(p4);
  0.002731 seconds (1.74 k allocations: 108.676 KiB)

julia> @btime $f4($v)
  11.781 ns (0 allocations: 0 bytes)
-1458502820125772303
```
</details>
