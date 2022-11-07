```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["para_rational_curves.md"]
```

# Rational Parametrizations of Rational Plane Curves

!!! note
    In this section, $C$ will denote a complex projective plane curve, defined by an absolutely irreducible,
    homogeneous polynomial in three variables, with coefficients in $\mathbb Q$. Moreover, we will write $n = \deg C$.

Recall that the curve $C$ is *rational* if it is birationally equivalent to the projective line $\mathbb P^1(\mathbb C)$.
In other words, there exists a *rational parametrization* of $C$, that is, a birational map $\mathbb P^1(\mathbb C)\dashrightarrow C$.
Note that such a parametrization is given by three homogeneous polynomials of the same degree in the homogeneous coordinates on
$\mathbb P^1(\mathbb C)$.

!!! note
    The curve $C$ is rational iff its geometric genus is zero.

Based on work of Max Noether on adjoint curves, Hilbert und Hurwitz showed that if
$C$ is rational, then there is a birational map $C \dashrightarrow D$ defined over $\mathbb Q$ such
that $D = \mathbb P^1(\mathbb C)$ if $n$ is odd, and $D\subset\mathbb P^2(\mathbb C)$ is a conic if $n$ is even.

!!! note
    If a conic $D$ contains a rational point, then there exists a parametrization of $D$ defined over $\mathbb Q$;
    otherwise, there exists a parametrization of $D$ defined over a quadratic field extension of $\mathbb Q$.

The approach of Hilbert und Hurwitz is constructive and allows one, in principle, to find rational parametrizations.
The resulting algorithm is not very practical, however, as the approach asks to compute adjoint curves repeatedly,
at each of a number of reduction steps.

The algorithm implemented in OSCAR relies on reduction steps of a different type and requires the computation of adjoint
curves only once. Its individual steps are interesting in their own right:

 - Assure that the curve  $C$ is rational by checking that its geometric genus is zero;
 - compute a basis of the adjoint curves of $C$ of degree ${n-2}$; each such basis defines a birational map $C \dashrightarrow C_{n-2},$
    where $C_{n-2}$ is a rational normal curve in $\mathbb P^{n-2}(\mathbb C)$;
 - the anticanonical linear system on $C_{n-2}$ defines a birational map $C_{n-2}\dashrightarrow C_{n-4}$, where $C_{n-4}$ is a rational normal curve in in $\mathbb P^{n-4}(\mathbb C)$;
 - iterate the previous step to obtain a birational map  $C_{n-2} \dashrightarrow \dots \dashrightarrow D$,
    where $D = \mathbb P^1(\mathbb C)$ if $n$ is odd, and $D\subset\mathbb P^2(\mathbb C)$ is a conic if $n$ is even;
 - invert the birational map  $C \dashrightarrow C_{n-2} \dashrightarrow \dots \dashrightarrow D$; 
 - if $n$ is even, compute a parametrization of the conic $D$ and compose it with the inverted map above.

!!! note
    The defining property of an adjoint curve is that it passes with “sufficiently high” multiplicity through the singularities of $C$.
    There are several concepts of making this precise. For each such concept, there is a corresponding  *adjoint ideal* of $C$,
	namely the homogeneous ideal formed by the defining polynomials of the adjoint curves. In OSCAR, we follow
	the concept of Gorenstein which leads to the largest possible adjoint ideal.
	
See [Bhm99](@cite) and [BDLP17](@cite) for details and further references.
 
## Creating Projective Plane Curves

The data structures for algebraic curves in OSCAR are still under development
and subject to change. Here is the current constructor for projective plane curves:

```@docs
ProjPlaneCurve(f::MPolyElem{T}) where {T <: FieldElem}
```

## The Genus of a Plane Curve

```@docs
 geometric_genus(C::ProjectivePlaneCurve{T}) where T <: FieldElem
```

## Adjoint Ideals of Plane Curves

```@docs
adjoint_ideal(C::ProjPlaneCurve{fmpq})
```

## Rational Points on Conics

```@docs
rational_point_conic(D::ProjPlaneCurve{fmpq})
```
## Parametrizing Rational Plane Curves

```@docs
 parametrization_plane_curve(C::ProjPlaneCurve{fmpq})
```






