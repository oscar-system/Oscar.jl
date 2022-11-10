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

The commutative algebra part of OSCAR provides functionality for dealing with

- multivariate polynomial rings and their ideals,
- quotients of multivariate polynomial rings modulo ideals, as well as ideals of such quotients,
- localizations of the above rings and their ideals, and 
- modules over all rings above.

In describing this functionality, we refer to quotients of multivariate polynomial rings also as *affine algebras*.

!!! note
    Fundamental to the algorithmic treatment of multivariate polynomial rings and their quotients is the concept of *Gröbner bases*.
    Buchberger's algorithm for computing such bases makes use of reduction steps which rely on multivariate division with remainder.
	Here, the natural idea is to extend Euclidean division with remainder, allowing more than one divisor. For this (in fact already for the
	very definition of Gröbner bases), we need consistent ways of distinguishing leading terms (leading monomials) of multivariate
	polynomials.  That is, a consistent concept of *monomial orderings* is required. The choice of ordering has a crucial
	impact on the performance of Buchberger's algorithm and on the resulting Gröbner basis.
	
!!! note
    The monomial orderings used in the context of Gröbner bases are all well-orderings. This guarantees the termination of the
    extended Euclidean division algorithm. In order to deal with localizations of multivariate polynomial rings and their quotients, however, we have to
	drop the well-ordering assumption. In this more general context, we use the name *standard basis* instead of Gröbner basis. Computationally,
	the termination of Buchberger's algorithm is guaranteed by using a variant of the multivariate division algorithm due to Mora.
    
!!! note
    The lexicograpical monomial ordering specifies the default way of storing and displaying multivariate polynomials in OSCAR (terms are sorted in descending order).
    The other orderings which can be attached to a multivariate polynomial ring are the degree lexicographical ordering  and the degree reverse lexicographical ordering.
	Independently of the orderings attached to the rings, Gröbner and standard bases can be computed with respect to any monomial ordering. See the corresponding
	sections below for details.

!!! note
     OSCAR provides functionality for equipping multivariate polynomial rings with gradings by finitely presented abelian groups. These gradings descend to quotients of multivariate polynomial rings modulo homogeneous ideals.
     A large majority of the functions discussed in what follows apply to both the ungraded and graded case. For simplicity of the presentation in this documentation, however, they
	 are often only illustrated by examples with focus on the former case, but work similarly for homogeneous ideals and graded modules in the latter case.

!!! note
    The main focus of the OSCAR functions discussed in what follows is on polynomial rings over fields (exact fields supported by OSCAR). Where not indicated otherwise,
    the functions also apply to polynomial rings over $\mathbb Z$. 
    
General textbooks offering details on theory and algorithms include: 
- [GP08](@cite)
- [DL06](@cite)
- [DP13](@cite)

