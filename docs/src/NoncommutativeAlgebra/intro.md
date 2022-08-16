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

The noncommutative algebra part of OSCAR provides functionality for dealing with finitely
presented associative algebras. Every such algebra is represented as the quotient of a free
associative algebra modulo a two-sided ideal of relations. Currently, only algebras over
fields are supported.


!!! note
    Most functions discussed in this chapter rely on noncommutative Gröbner bases.
    They are implemented for algebras over exact fields supported by OSCAR
    and execute corresponding Singular functionality.

In contrast to the case of multivariate polynomial rings considered in the commutative
algebra part, the algebras dicussed here need not be Noetherian. From a practical point
of view, Gröbner bases need not be finite, even if we start from a finitely generated ideal.
The situation is different for the large class of G-algebras whose computational study will
be a main topic in what follows. In fact, G-algebras enjoy structural properties reminiscent
of those of multivariate polynomial rings. In analogy to the commutative algebra part, where
we use Gröbner bases in multivariate polynomial rings to study more generally affine algebras,
that is, quotients of multivariate polynomial rings modulo ideals, we now also study GR-algebras,
that is, quotients of G-algebras modulo two-sided ideals. Note that the class of G-algebras includes
the Weyl algebras together with quite a variety of other important algebras. Examples of GR-algebras
are the Clifford algebras (in particular, the exterior algebras).

General textbooks offering details on the theory and algorithms presented in this chapter include:
- [GP08](@cite)
- [DL06](@cite)
- [BGV03](@cite)
