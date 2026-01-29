```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Adjunction Process for Surfaces

A surface in this section is a smooth projective surface over $\mathbb C$.

Blowing up a surface in a point means to replace the point by an exceptional curve. Each such curve $E$ is a smooth, rational
curve with self-intersection number $E^{2}=-1$. We speak of  a *$(-1)$-curve*. A surface is *minimal* if it contains no
$(-1)$-curves. That is, the surface cannot be obtained by blowing up a point on another surface. A surface $X_{\text{min}}$
is called a *minimal model* of a surface $X$ if $X_{\text{min}}$ is minimal and $X$ can be obtained from $X_{\text{min}}$
by repeatedly blowing up a point. Each surface $X$ has a minimal model which is unique if $X$ has non-negative Kodaira dimension.
The Enriques-Kodaira classification classifies surfaces according to their minimal models. See [BHPV-D-V04](@cite) for more on this.

Given a surface, we may apply the *adjunction process* of Van de Ven and Sommese [SV-D-V87](@cite) to discover a minimal model.
To describe this process, consider a surface $X \subset \mathbb P^{n}$ of codimension $c$. Let $S$ and $S_{X}$
denote the homogeneous coordinate rings of $\mathbb P^{n}$ and $X$, respectively. Consider $\omega_{X}=\text{Ext}^{c}_{S}(S_{X},S(-n-1)),$
the graded *dualizing module* of $S_{X}$. A basis of the graded piece $(\omega_{X})_{{1}}$ corresponds to the linear system $|K_X +H|$, where $K_X$ is a canonical
divisor on $X$ and $H$ is the hyperplane class. Except for some exceptional cases, this linear system defines a birational morphism $\varphi_{|K_X+H|}\colon X \to X'$
onto another smooth projective surface $X'$ such that  $\varphi_{|K_X+H|}$ blows down precisely all $(-1)$-lines on $X$.
As shown by Van de Ven and Sommese, in the exceptional cases,

- ``X`` is a linearly or quadratically embedded $ \mathbb P^{2}$ or $X$ is ruled by lines, in which case $|K_X+H| = \emptyset$,
- ``X`` is an anti-canonically embedded del Pezzo surface, in which case $\varphi_{|K_X+H|}$ maps $X$ to a point,
- ``X`` is a conic bundle, in which case  $\varphi_{|K_{X}+H|}\colon X \to B$ maps $X$ to a curve $B$ such that the fibers of $\varphi_{|K_{X}+H|}$ are the conics, or
- ``X`` is a surface in one of four explicit families identified by Sommese and Van de Ven, and $\varphi_{|K_X+H|}\colon X \to X'$ is not birational, but finite to one.

If we are not in one of these cases, a $(-1)$-conic $C$ in $X$ is mapped to a $(-1)$-line in $X'$ since $(K_X+H)\;\!. \;\! C=-1+2=1$.
Thus, the *adjunction process*, which  consists of applying the *adjunction maps* $\varphi_{|K_X+H|}$, $\varphi_{|K_{X'}+H'|}$ and so on, yields finitely many surfaces
$X \rightarrow X^{\prime} \rightarrow X^{\prime\prime} \rightarrow \dots$ which are called the *adjoint surfaces* of $X$. The last adjoint surface is either minimal or belongs to one of the
exceptional cases. In particular, if  $X$ has non-negative Kodaira dimension, the adjunction process yields the uniquely determined minimal model of $X$.

!!! note
    If $X$ is rational, the last adjoint surface is either $\mathbb P^{2}$, the Veronese surface, a Hirzebruch surface, a Del Pezzo surface, a conic bundle, or one of the four explicit families identified by Sommese and Van de Ven.

!!! note
    In explicit computations, we consider surfaces which are defined by polynomial equations with coefficients in a subfield of $\mathbb C$ which can be handled by OSCAR.

!!! note
    The surfaces in the examples below are taken from the OSCAR data base of nongeneral type surfaces in $\mathbb P^4$. To ease subsequent computations, the surfaces in the data base where constructed over finite fields. Note, however, that the recipes used in the constructions also work in characteristic zero. So all computations can be confirmed in characteristic zero, although this may be time consuming.

## Adjunction Process

What we describe here goes back to joint work of Wolfram Decker and Frank-Olaf Schreyer. See [DES93](@cite),  [DS00](@cite).

```@docs
adjunction_process(X::AbsProjectiveVariety, steps::Int=0)
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Wolfram Decker](https://math.rptu.de/en/wgs/agag/people/head/decker).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

