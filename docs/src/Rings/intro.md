```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Introduction](@id rings)

The rings part of OSCAR provides functionality for handling
various kinds of rings: 
- the ring of integers
- polynomial rings (univariate and multivariate, see [Generic univariate polynomial types](@ref) and [Generic sparse distributed multivariable polynomial types](@ref)),
- orders in number fields
- series rings


## [Subtle but important: `/` vs `//`](@id subtle_distinction_for_rings)

OSCAR distinguishes a number of different kinds of division. In particular, the operators `/` and `//` have distinct meanings.

- `x / y` performs exact division within the current algebraic structure, if defined (and provides a convenient shorthand for `divexact`).
- `x // y` constructs a formal quotient, placing the result in a fraction-field parent.

Example:

```julia
julia> R, x = polynomial_ring(QQ, :x)
(Univariate polynomial ring in x over QQ, x)

julia> f = 1+x
x + 1

julia> parent(f//2)
Fraction field
  of univariate polynomial ring in x over QQ

julia> parent(f/2)
Univariate polynomial ring in x over QQ
```

!!! note
    The above behavior applies to OSCAR types. For plain Julia numbers, `/`
    denotes floating-point division. Indeed, it is a common error to enter
    `1/2` for the fraction 'one half' in Julia. See [here](@ref division_of_integers_in_OSCAR)
    for further details on the division of integers in OSCAR.

If the parent is already a field (which means it is canonically isomorphic to its field of fractions), `//` coincides with exact division:

```julia
julia> F101 = GF(101)
Finite field of characteristic 101

julia> j = F101(9)
9

julia> parent(j//j)
Finite field of characteristic 101

julia> parent(j) == parent(j//j)
true
```

If the parent is not an integral domain, a field of fractions does not exist in a strict mathematical sense. Nevertheless, `//` may still construct formal fractions. However, computations with such objects may fail. Use with care!

The following example illustrates such failures for ``\mathbb{Z}/100 \mathbb{Z}`` (which is not an integral domain, since ``100 = 2^2 \times 5^2``, hence has zero divisors).

```julia
julia> RR,_ = residue_ring(ZZ, 100)
(Integers modulo 100, Map: ZZ -> ZZ/(100))

julia> a = RR(9)
9

julia> parent(a)
Integers modulo 100

julia> W = parent(a//a)
Fraction field
  of integers modulo 100

julia> new_j = 10*W(1)
10

julia> inv_new_j = 1/new_j
1//10

julia> inv_new_j^2
Error showing value of type AbstractAlgebra.Generic.FracFieldElem{zzModRingElem}:
ERROR: Impossible inverse in Integers modulo 100

julia> weird = (inv_new_j)^2;

julia> is_zero(1//weird)
true
```

In contrast, running the above code over ``\mathbb{Z}/101\mathbb{Z}`` completes successfully.



## Contact

Please direct questions about this part of OSCAR to the following people:
* [Claus Fieker](https://math.rptu.de/en/wgs/agag/people/head/fieker),
* [Tommy Hofmann](https://www.thofma.com/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
