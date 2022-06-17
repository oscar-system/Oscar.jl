# Oscar.jl

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ga-img]][ga-url] [![][codecov-img]][codecov-url] |


Welcome to the OSCAR project, a visionary new computer algebra system
which combines the capabilities of four cornerstone systems: GAP,
Polymake, Antic and Singular.

## Installation

OSCAR requires Julia 1.6 or newer. In principle it can be installed and used
like any other Julia package; doing so will take a couple of minutes:

```julia
julia> using Pkg
julia> Pkg.add("Oscar")
julia> using Oscar
```

However, some of Oscar's components have additional requirements.
For more detailed information, please consult the [installation
instructions](https://oscar.computeralgebra.de/install/) on our website.

## Contributing to OSCAR

Please read the [introduction for new developers](https://oscar-system.github.io/Oscar.jl/dev/DeveloperDocumentation/new_developers/)
in the OSCAR manual to learn more on how to contribute to OSCAR.

## Examples of usage

```julia
julia> using Oscar
 -----    -----    -----      -      -----
|     |  |     |  |     |    | |    |     |
|     |  |        |         |   |   |     |
|     |   -----   |        |     |  |-----
|     |        |  |        |-----|  |   |
|     |  |     |  |     |  |     |  |    |
 -----    -----    -----   -     -  -     -

...combining (and extending) ANTIC, GAP, Polymake and Singular
Version 0.10.0 ...
... which comes with absolutely no warranty whatsoever
Type: '?Oscar' for more information
(c) 2019-2022 by The Oscar Development Team

julia> k, a = quadratic_field(-5)
(Imaginary quadratic field defined by x^2 + 5, sqrt(-5))

julia> zk = maximal_order(k)
Maximal order of Imaginary quadratic field defined by x^2 + 5
with basis nf_elem[1, sqrt(-5)]

julia> factorisations(zk(6))
2-element Vector{Fac{NfOrdElem}}:
 -1 * (sqrt(-5) + 1) * (sqrt(-5) - 1)
 -1 * 2 * -3

julia> Qx, x = PolynomialRing(QQ, [:x1,:x2])
(Multivariate Polynomial Ring in x1, x2 over Rational Field, fmpq_mpoly[x1, x2])

julia> R = grade(Qx, [1,2])[1]
Multivariate Polynomial Ring in x1, x2 over Rational Field graded by
  x1 -> [1]
  x2 -> [2]

julia> f = R(x[1]^2+x[2])
x1^2 + x2

julia> degree(f)
graded by [2]

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> s = sub(F, [f*F[1]])
Submodule with 1 generator
1 -> (x1^2 + x2)*e[1]
represented as subquotient with no relations.

julia> H, mH = hom(s, quo(F, s))
(hom of (s, Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> (x1^2 + x2)*e[1]), Map from
H to Set of all homomorphisms from Submodule with 1 generator
1 -> (x1^2 + x2)*e[1]
represented as subquotient with no relations. to Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> (x1^2 + x2)*e[1] defined by a julia-function with inverse)

julia> mH(H[1])
Map with following data
Domain:
=======
s
Codomain:
=========
Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> (x1^2 + x2)*e[1]
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

## Citing OSCAR

If you have used OSCAR in the preparation of a paper please cite it as described below:

    [OSCAR]
        OSCAR -- Open Source Computer Algebra Research system, Version 0.10.0 The OSCAR Team, 2022. (https://oscar.computeralgebra.de)
    [OSCAR-book]
        Christian Eder, Wolfram Decker, Claus Fieker, Max Horn, Michael Joswig, The OSCAR book, 2024.

If you are using BibTeX, you can use the following BibTeX entries:

    @misc{OSCAR,
      key          = {OSCAR},
      organization = {The OSCAR Team},
      title        = {OSCAR -- Open Source Computer Algebra Research system,
                      Version 0.10.0},
      year         = {2022},
      url          = {https://oscar.computeralgebra.de},
      }

    @Book{OSCAR-book,
      editor = {Eder, Christian and Decker, Wolfram and Fieker, Claus and Horn, Max and Joswig, Michael},
      title = {The OSCAR book},
      year = {2024},
    }

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
