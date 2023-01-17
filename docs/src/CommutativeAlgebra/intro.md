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
- quotients of multivariate polynomial rings modulo ideals and ideals of such quotients,
- localizations of the above rings and ideals of such localizations, and 
- modules over all rings above.

We use *affine algebra* as a synonym for *quotient ring of a multivariate polynomial ring modulo an ideal*.

!!! note
    Fundamental to most functions described below is the concept of *standard bases*. Each such basis
    is defined relative to a *monomial ordering*. If this ordering is a well-ordering, a standard basis is also called
    a *Gröbner basis*.

!!! note
    Every multivariate polynomial ring in OSCAR comes equipped with a monomial ordering according to which the
    polynomials are stored and displayed. Independently of this ordering, standard bases can be computed with respect
	to any monomial ordering: The `groebner_basis` and `standard_basis` functions provided by OSCAR allow us to
	specify the desired monomial `ordering` as a key word argument.
	
!!! note
	OSCAR provides functionality for equipping multivariate polynomial rings with gradings by finitely presented groups. 
    Most functions discussed below apply to both ungraded and graded polynomial rings.
	For simplicity of the presentation in this documentation, however, they are often only illustrated by examples with
	focus on the former case, but work similarly for homogeneous ideals and graded modules in the latter case.

!!! note
    Our main focus in this chapter is on multivariate polynomial rings over fields (exact fields supported by OSCAR). Where not indicated
    otherwise, the presented functions also apply to polynomial rings over $\mathbb Z$. 
    
General textbooks offering details on theory and algorithms include: 
- [GP08](@cite)
- [DL06](@cite)
- [DP13](@cite)

