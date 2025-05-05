```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# [Introduction](@id commutative_algebra)

The commutative algebra part of OSCAR provides functionality for dealing with

- multivariate polynomial rings and their ideals,
- quotients of multivariate polynomial rings modulo ideals and ideals of such quotients,
- localizations of the above rings and ideals of such localizations, and 
- modules over all rings above.

We use *affine algebra* as a synonym for *quotient of a multivariate polynomial ring modulo an ideal*.

Fundamental to computational commutative algebra is the concept of *standard bases*. Each such basis
is defined relative to a *monomial ordering*. If this ordering is a well-ordering, a standard basis is also called
a *Gröbner basis*. We refer to the corresponding [section](@ref gb_fields) in this chapter for details.

!!! note
    Each multivariate polynomial ring in OSCAR comes equipped with a monomial ordering according to which the
    polynomials are stored and displayed. Independently of this ordering, standard bases can be computed with respect
    to any monomial ordering: The `groebner_basis` and `standard_basis` functions provided by OSCAR allow one to
    specify the desired monomial `ordering` as a key word argument. Typically, however, the user does not have
    to worry about Gröbner (standard) bases: The functions discussed in this chapter compute such bases behind
    the scenes when needed. Once computed, each such basis is cached for later reuse.

!!! note
    Our main focus in this chapter is on Gröbner (standard) basis methods for computations in multivariate polynomial rings over fields (exact fields supported by OSCAR). Where not indicated otherwise, such methods apply to polynomial rings over $\mathbb Z$, too. Similarly for polynomial rings over rings of type $\mathbb Z/ m\mathbb Z$. Note, however, that computing radicals and primary decompositions requires in addition methods for sqarefree decomposition and polynomial factorization, respectively. In the case of coefficient rings (fields) for which such methods are not implemented, an error will be thrown.
	
!!! note
    In Oscar, it is possible to equip multivariate polynomial rings with gradings by finitely presented groups. 
    Most functions related to multivariate polynomial rings discussed in this chapter apply to both the ungraded and graded case.
	However, for simplicity of the presentation, in this documentation, the functions are often only illustrated by examples with
    focus on the former case, but work similarly for homogeneous ideals and graded modules in the latter case.

General textbooks offering details on theory and algorithms include: 
- [GP08](@cite)
- [DL06](@cite)
- [DP13](@cite)


## Tutorials

We encourage you to take a look at the tutorials on commutative algebra in
OSCAR, which can be found [here](https://www.oscar-system.org/tutorials/CommutativeAlgebra/).


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Wolfram Decker](https://math.rptu.de/en/wgs/agag/people/head/decker).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
