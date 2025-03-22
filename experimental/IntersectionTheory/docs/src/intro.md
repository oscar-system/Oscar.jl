```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

In this chapter, we
- introduce OSCAR tools which support computations in intersection theory, and
- give examples which illustrate how intersection theory is used to solve problems from enumerative geometry.

!!! note
    Making use of OSCAR, a first version of what we present here was written by Jeiao Song as a julia package.
    This package was "heavily inspired by the Macaulay2 package Schubert2 and the Sage library Chow. Some
    functionalities from [the Sage library] Schubert3 are also implemented."

!!! note
    Schubert2 was written by Daniel R. Grayson, Michael E. Stillman, Stein A. Strømme, David Eisenbud, and Charley Crissman
    while Chow is due to Manfred Lehn and Christoph Sorger. Schubert3  as well as the Singular library schubert.lib were
    written by Dang Tuan Hiep. The basis for all this work, including ours, is the Maple package Schubert written
    by Sheldon Katz and Stein A. Strømme.

Throughout the chapter, the varieties we consider are smooth projective varieties over the complex numbers.

!!! note
    The Chow ring of a variety `X` is the group of cycles of `X` modulo an equivalence relation,
    together with the intersection pairing which defines the multiplication of the ring. Here,
    in contrast to most textbooks, we consider numerical equivalence classes of cycles rather than
    rational equivalence classes.

Our approach is abstract in the sense that we do not work with concrete varieties; that is,
our varieties are not given by equations. Instead, we represent a variety by specifying its
dimension together with its Chow ring and, possibly, further data. We refer to such such
a collection of data as an *abstract variety*, and to results obtained from manipulating
the data as results which apply to all (smooth projective complex) varieties sharing the data. 

Of particular interest is the tangent bundle of a variety (recall that the Todd class of the
tangent bundle enters the Hirzebruch-Riemann-Roch formula). As with any other vector bundle,
the tangent bundle is represented as a collection of data referred to as an *abstract vector bundle*.
The main data here is the Chern character polynomial of the vector bundle.

In the same spirit, we introduce  *abstract variety maps*. Their key ingredient is the pullback 
morphism between the Chow rings.

General textbooks offering details on theory and algorithms include: 
- [EH16](@cite)
- [Ful98](@cite)

For computations in the Chow rings of abstract flag bundles see 
- [GSS22](@cite).


## Tutorials

We encourage you to take a look at the tutorials on intersection theory in OSCAR,
which can be found [here](https://www.oscar-system.org/tutorials/IntersectionTheory/).


## Contact

Please direct questions about this part of OSCAR to the following
people:
* [Pieter Belmans](https://pbelmans.ncag.info/),
* [Wolfram Decker](https://math.rptu.de/en/wgs/agag/people/head/decker),
* [Tommy Hofmann](https://www.thofma.com/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
