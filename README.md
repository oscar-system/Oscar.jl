# StraightLinePrograms

This is a buggy WIP implementation of straight-line programs (SLP).
This is part of the [Oscar](https://oscar.computeralgebra.de/) project.

Currently, SLPs have a polynomial interface (`SLPoly`), and the simplest way
to construct them is through `LazyPoly`s.

## Examples

```julia
julia> using AbstractAlgebra, StraightLinePrograms, BenchmarkTools;

julia> L = LazyPolyRing(zz); x, y = L(:x), L(:y)
(x, y)

julia> S = SLPolyRing(zz, [:x, :y]);

julia> p = S(3 + 2x * y^2) # each line of the SLP is shown with current value
  #1 = * 2 x    ==>     (2x)
  #2 = ^ y 2    ==>     y^2
  #3 = * #1 #2  ==>     (2xy^2)
  #4 = + 3 #3   ==>     (3 + (2xy^2))

julia> p.cs # constants used in the program
2-element Array{Int64,1}:
 3
 2

julia> p.lines # each line is a UInt64 encoding an opcode and 2 arguments
4-element Array{UInt64,1}:
 0x0500000028000001
 0x8780000020000002
 0x0500000030000004
 0x0300000010000005

julia> evaluate(p, [2, 3])
39

julia> p2 = (p*S(x*y))^6
  #1 = * 2 x    ==>     (2x)
  #2 = ^ y 2    ==>     y^2
  #3 = * #1 #2  ==>     (2xy^2)
  #4 = + 3 #3   ==>     (3 + (2xy^2))
  #5 = * x y    ==>     (xy)
  #6 = * #4 #5  ==>     ((3 + (2xy^2))xy)
  #7 = ^ #6 6   ==>     ((3 + (2xy^2))xy)^6

julia> R, (x1, y1) = PolynomialRing(zz, ["x", "y"]); R
Multivariate Polynomial Ring in x, y over Integers

julia> q = convert(R, p2)
64*x^12*y^18+576*x^11*y^16+2160*x^10*y^14+4320*x^9*y^12+4860*x^8*y^10+2916*x^7*y^8+729*x^6*y^6

julia> v = [3, 5]; @btime evaluate($q, $v)
  32.101 μs (634 allocations: 45.45 KiB)
-1458502820125772303

julia> @btime evaluate($p2, $v)
 415.709 ns (6 allocations: 384 bytes)
-1458502820125772303

julia> res = Int[]; @btime StraightLinePrograms.evaluate!($res, $p2, $v)
  307.940 ns (2 allocations: 32 bytes)
-1458502820125772303

julia> res # intermediate computations (first 2 elements are constants)
9-element Array{Int64,1}:
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
  #1 = ^ x 6    ==>     x^6
  #2 = ^ x 7    ==>     x^7
  #3 = ^ x 8    ==>     x^8
  #4 = ^ x 9    ==>     x^9
  #5 = ^ x 10   ==>     x^10
  #6 = ^ x 11   ==>     x^11
  #7 = ^ x 12   ==>     x^12
  #8 = ^ y 6    ==>     y^6
  #9 = ^ y 8    ==>     y^8
 #10 = ^ y 10   ==>     y^10
 #11 = ^ y 12   ==>     y^12
 #12 = ^ y 14   ==>     y^14
 #13 = ^ y 16   ==>     y^16
 #14 = ^ y 18   ==>     y^18
 #15 = * 64 #7  ==>     (64x^12)
 #16 = * #15 #14        ==>     (64x^12y^18)
 #17 = * 576 #6 ==>     (576x^11)
 #18 = * #17 #13        ==>     (576x^11y^16)
 #19 = + #16 #18        ==>     ((64x^12y^18) + (576x^11y^16))
 #20 = * 2160 #5        ==>     (2160x^10)
 #21 = * #20 #12        ==>     (2160x^10y^14)
 #22 = + #19 #21        ==>     ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14))
 #23 = * 4320 #4        ==>     (4320x^9)
 #24 = * #23 #11        ==>     (4320x^9y^12)
 #25 = + #22 #24        ==>     ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12))
 #26 = * 4860 #3        ==>     (4860x^8)
 #27 = * #26 #10        ==>     (4860x^8y^10)
 #28 = + #25 #27        ==>     ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12) + (4860x^8y^10))
 #29 = * 2916 #2        ==>     (2916x^7)
 #30 = * #29 #9 ==>     (2916x^7y^8)
 #31 = + #28 #30        ==>     ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12) + (4860x^8y^10) + (2916x^7y^8))
 #32 = * 729 #1 ==>     (729x^6)
 #33 = * #32 #8 ==>     (729x^6y^6)
 #34 = + #31 #33        ==>     ((64x^12y^18) + (576x^11y^16) + (2160x^10y^14) + (4320x^9y^12) + (4860x^8y^10) + (2916x^7y^8) + (729x^6y^6))

julia> @btime evaluate($p3, $v)
  1.452 μs (39 allocations: 1.52 KiB)
-1458502820125772303

julia> @time f3 = StraightLinePrograms.compile!(p3);
  0.003005 seconds (1.40 k allocations: 90.930 KiB)

julia> @btime $f3($v)
 80.570 ns (0 allocations: 0 bytes)
-1458502820125772303
```
