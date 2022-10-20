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

The commutative algebra part of OSCAR provides functionality for handling

- ideals of multivariate polynomial rings,
- quotients of multivariate polynomial rings modulo ideals, as well as ideals of such quotients, and 
- modules over the above rings.

In describing this functionality, we will refer to quotients of multivariate polynomial rings also as *affine algebras*.

!!! note
    Most functions discussed in this chapter rely on Gröbner basis techniques. They either execute corresponding Singular functionality, or are written as pieces of OSCAR code which rely on Singular for Gröbner basis computations.The functions are  implemented for multivariate polynomial rings over fields (exact fields supported by OSCAR) and, if not indicated otherwise, for multivariate polynomial rings over the integers.

!!! note
    The workhorse behind the concept of Gröbner bases is Buchberger's algorithm for computing Gröbner bases. For both the concept and the algorithm a convenient way of
    ordering the monomials appearing in multivariate polynomials and, thus, to distinguish leading terms of such polynomials is needed.  Note that the choice of such a monomial
	ordering has a crucial impact on the performance of the algorithm and the resulting Gröbner basis.

!!! note
    The lexicograpical monomial ordering `lex` specifies the default way of storing and displaying multivariate polynomials in OSCAR (terms are sorted in descending order).
    The other orderings which can be attached to a multivariate polynomial ring are the degree lexicographical ordering `deglex` and the degree reverse lexicographical
	ordering`degrevlex`. Independently of the attached orderings, Gröbner bases can be computed with respect to any monomial ordering. See the sections on
	monomial orderings and Gröbner bases.

!!! note
     OSCAR provides functionality for equipping multivariate polynomial rings with gradings. These gradings descend to quotients of multivariate polynomial rings modulo homogeneous ideals.
     A large majority of the functions discussed in what follows apply to both the ungraded and graded case. For simplicity of the presentation in this documentation, however, such functions
	 are often only illustrated by examples with focus on the former case, but work similarly for homogeneous ideals and graded modules in the latter case.

General textbooks offering details on theory and algorithms include: 
- [GP08](@cite)
- [DL06](@cite)
- [DP13](@cite)

