```@meta
CurrentModule = Oscar
```
# Projective Plane Curves
```@docs
ProjectivePlaneCurve
```

Projective plane curves are modeled in Oscar as projective
algebraic sets. See `AbsProjectiveAlgebraicSet`(@ref).
In addition to the methods for algebraic sets
the following methods special to plane curves are available.

```@docs
defining_equation(C::ProjectivePlaneCurve{S,MPolyQuoRing{T}}) where {S,T}
degree(C::ProjectivePlaneCurve)
common_components(C::S, D::S) where {S<:ProjectivePlaneCurve}
multiplicity(C::ProjectivePlaneCurve, P::AbsProjectiveRationalPoint)
tangent_lines(C::ProjectivePlaneCurve, P::AbsProjectiveRationalPoint)
intersection_multiplicity(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve
is_transverse_intersection(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve
arithmetic_genus(C::ProjectivePlaneCurve)
geometric_genus(C::AbsProjectiveCurve)
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




## Adjoint Ideals of Plane Curves

```@docs
adjoint_ideal(C::ProjectivePlaneCurve{QQField})
```

## Rational Points on Conics

```@docs
rational_point_conic(D::ProjectivePlaneCurve{QQField})
```
## Parametrizing Rational Plane Curves

```@docs
parametrization(C::ProjectivePlaneCurve{QQField})
```


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Janko Böhm](https://www.mathematik.uni-kl.de/~boehm/),
* [Wolfram Decker](https://www.mathematik.uni-kl.de/en/agag/people/head/prof-dr-wolfram-decker/seite).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
