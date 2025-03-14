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

```
julia> using Pkg
julia> Pkg.add("Oscar")
julia> using Oscar
```

However, some of OSCAR's components have additional requirements.
For more detailed information, please consult the [installation
instructions](https://www.oscar-system.org/install/) on our website.

## Contributing to OSCAR

Please read the [introduction for new developers](https://docs.oscar-system.org/dev/DeveloperDocumentation/new_developers/)
in the OSCAR manual to learn more on how to contribute to OSCAR.

## Examples of usage

```
julia> using Oscar
  ___   ____   ____    _    ____
 / _ \ / ___| / ___|  / \  |  _ \   |  Combining ANTIC, GAP, Polymake, Singular
| | | |\___ \| |     / _ \ | |_) |  |  Type "?Oscar" for more information
| |_| | ___) | |___ / ___ \|  _ <   |  Manual: https://docs.oscar-system.org
 \___/ |____/ \____/_/   \_\_| \_\  |  Version 1.3.1

julia> k, a = quadratic_field(-5)
(Imaginary quadratic field defined by x^2 + 5, sqrt(-5))

julia> zk = maximal_order(k)
Maximal order of Imaginary quadratic field defined by x^2 + 5
with basis AbsSimpleNumFieldElem[1, sqrt(-5)]

julia> factorizations(zk(6))
2-element Vector{Fac{AbsSimpleNumFieldOrderElem}}:
 -1 * -3 * 2
 -1 * (-sqrt(-5) - 1) * (-sqrt(-5) + 1)

julia> Qx, x = polynomial_ring(QQ, [:x1,:x2])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x1, x2])

julia> R = grade(Qx, [1,2])[1]
Multivariate polynomial ring in 2 variables over QQ graded by
  x1 -> [1]
  x2 -> [2]

julia> f = R(x[1]^2+x[2])
x1^2 + x2

julia> degree(f)
[2]

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> s = sub(F, [f*F[1]])[1]
Submodule with 1 generator
  1: (x1^2 + x2)*e[1]
represented as subquotient with no relations

julia> H, mH = hom(s, quo(F, s)[1])
(hom of (s, Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 1 generator
  1: (x1^2 + x2)*e[1]), Map: H -> set of all homomorphisms from s to subquotient of submodule with 1 generator
  1: e[1]
by submodule with 1 generator
  1: (x1^2 + x2)*e[1])

julia> mH(H[1])
Module homomorphism
  from s
  to subquotient of submodule with 1 generator
    1: e[1]
  by submodule with 1 generator
    1: (x1^2 + x2)*e[1]
```

Of course, the cornerstones are also available directly. For example:

```
julia> C = Polymake.polytope.cube(3);

julia> C.F_VECTOR
pm::Vector<pm::Integer>
8 12 6

julia> RP2 = Polymake.topaz.real_projective_plane();

julia> RP2.HOMOLOGY
pm::Array<topaz::HomologyGroup<pm::Integer> >
({} 0)
({(2 1)} 0)
({} 0)
```

## Citing OSCAR

If you have used OSCAR in the preparation of a paper please cite it as described below:

    [OSCAR]
        OSCAR -- Open Source Computer Algebra Research system, Version 1.3.1,
        The OSCAR Team, 2025. (https://www.oscar-system.org)
    [OSCAR-book]
        Wolfram Decker, Christian Eder, Claus Fieker, Max Horn, Michael Joswig, eds.
        The Computer Algebra System OSCAR: Algorithms and Examples,
        Algorithms and Computation in Mathematics, Springer, 2025.

If you are using BibTeX, you can use the following BibTeX entries:

    @misc{OSCAR,
      key          = {OSCAR},
      organization = {The OSCAR Team},
      title        = {OSCAR -- Open Source Computer Algebra Research system,
                      Version 1.3.1},
      year         = {2025},
      url          = {https://www.oscar-system.org},
      }

    @book{OSCAR-book,
      editor = {Decker, Wolfram and Eder, Christian and Fieker, Claus and Horn, Max and Joswig, Michael},
      title = {The {C}omputer {A}lgebra {S}ystem {OSCAR}: {A}lgorithms and {E}xamples},
      year = {2025},
      publisher = {Springer},
      series = {Algorithms and {C}omputation in {M}athematics},
      volume = {32},
      edition = {1},
      url = {https://link.springer.com/book/9783031621260},
      issn = {1431-1550},
      doi = {10.1007/978-3-031-62127-7},
    }

## Funding

The development of this Julia package is supported by the
German Research Foundation (DFG) within the
[Collaborative Research Center TRR 195](https://www.computeralgebra.de/sfb/).

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://docs.oscar-system.org/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://docs.oscar-system.org/stable/

[ga-img]: https://github.com/oscar-system/Oscar.jl/actions/workflows/CI.yml/badge.svg?branch=master&event=push
[ga-url]: https://github.com/oscar-system/Oscar.jl/actions?query=workflow%3A%22Run+tests%22

[codecov-img]: https://codecov.io/gh/oscar-system/Oscar.jl/branch/master/graph/badge.svg?branch=master
[codecov-url]: https://codecov.io/gh/oscar-system/Oscar.jl
