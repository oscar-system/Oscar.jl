# Oscar.jl

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ga-img]][ga-url] [![][codecov-img]][codecov-url] |


Welcome to the OSCAR project, a visionary new computer algebra system
which combines the capabilities of four cornerstone systems: GAP,
Polymake, Antic and Singular.

## Installation

OSCAR requires Julia 1.3 or newer. In principle it can be installed and used
like any other Julia package; doing so will take a couple of minutes:

```julia
julia> using Pkg
julia> Pkg.add("Oscar")
julia> using Oscar
```

However, some of Oscar's components have additional requirements.
For more detailed information, please consult the [installation
instructions](https://oscar.computeralgebra.de/install/) on our website.

## Examples of usage

```julia
julia> using Oscar
...
 -----    -----    -----      -      -----   
|     |  |     |  |     |    | |    |     |  
|     |  |        |         |   |   |     |  
|     |   -----   |        |     |  |-----   
|     |        |  |        |-----|  |   |    
|     |  |     |  |     |  |     |  |    |   
 -----    -----    -----   -     -  -     -  

...combining (and extending) GAP, Hecke, Nemo, Polymake and Singular
Version 0.5.1-DEV ... 
 ... which comes with absolutely no warranty whatsoever
Type: '?Oscar' for more information
(c) 2019-2021 by The Oscar Development Team


julia> k, a = quadratic_field(-5)
(Number field over Rational Field with defining polynomial x^2+5, sqrt(-5))

julia> zk = maximal_order(k)
Maximal order of Number field over Rational Field with defining polynomial x^2+5
with basis nf_elem[1, sqrt(-5)]

julia> factorisations(zk(6))
2-element Vector{Fac{NfAbsOrdElem{AnticNumberField,nf_elem}}}:
 -1 * (2) * (-3)
 -1 * (sqrt(-5)+1) * (sqrt(-5)-1)

julia> Qx, x = PolynomialRing(QQ, :x=>1:2)
(Multivariate Polynomial Ring in x1, x2 over Rational Field, fmpq_mpoly[x1, x2])

julia> R = grade(Qx, [1,2])
Multivariate Polynomial Ring in x1, x2 over Rational Field graded by 
        x1 -> [1]
        x2 -> [2]

julia> f = R(x[1]^2+x[2])
x1^2 + x2
julia> degree(f)
graded by [2]

julia> F = FreeModule(R, 1)
Free module of rank 1 over R, graded as R^1([0])

julia> s = sub(F, [f*F[1]])
Subquotient by Array of length 1
1 -> (x1^2 + x2)*e[1]

a> mH(H[1])
Map with following data
Domain:
=======
s
Codomain:
=========
Subquotient of Array of length 1
1 -> (1)*e[1]
 by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

julia> H, mH = hom(s, quo(F, s))
(hom of (s, Subquotient of Array of length 1
1 -> (1)*e[1]
 by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

), Map from
H to Set of all homomorphisms from Subquotient by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

 to Subquotient of Array of length 1
1 -> (1)*e[1]
 by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

 defined by a julia-function with inverse
)

julia> D = decoration(H)
GrpAb: Z

julia> homogeneous_component(H, D[0])
(H_[0] of dim 2, Map from
H_[0] of dim 2 to H defined by a julia-function with inverse
)
```

Of course, the cornerstones are also available directly:

```julia
julia> C = Polymake.polytope.cube(3);

julia> C.F_VECTOR
pm::Vector<pm::Integer>
8 12 6

julia> RP2 = Polymake.topaz.real_projective_plane();

julia> RP2.HOMOLOGY
PropertyValue wrapping pm::Array<polymake::topaz::HomologyGroup<pm::Integer>>
({} 0)
({(2 1)} 0)
({} 0)

```

## Funding

The development of this Julia package is supported by the Deutsche
Forschungsgemeinschaft DFG within the
[Collaborative Research Center TRR 195](https://www.computeralgebra.de/sfb/).

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://oscar-system.github.io/Oscar.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://oscar-system.github.io/Oscar.jl/stable/

[ga-img]: https://github.com/oscar-system/Oscar.jl/workflows/Run%20tests/badge.svg
[ga-url]: https://github.com/oscar-system/Oscar.jl/actions?query=workflow%3A%22Run+tests%22

[codecov-img]: https://codecov.io/gh/oscar-system/Oscar.jl/branch/master/graph/badge.svg?branch=master
[codecov-url]: https://codecov.io/gh/oscar-system/Oscar.jl

[coveralls-img]: https://coveralls.io/repos/github/oscar-system/Oscar.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/oscar-system/Oscar.jl?branch=master
