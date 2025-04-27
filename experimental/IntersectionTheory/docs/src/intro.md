```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

In this chapter, we
- introduce OSCAR tools which support computations in intersection theory, and
- give examples which illustrate how intersection theory is used to solve problems from enumerative geometry.
The varieties we are interested in are smooth projective varieties over the complex numbers.

!!! note
    Making use of OSCAR, a first version of what we present here was written by Jeiao Song as a `julia` package.
    This package was "heavily inspired by the Macaulay2 package `Schubert2` and the `sage` library `Chow`.
	Some functionalities from [the `sage` library] `Schubert3` are also implemented."
    The authors of `Schubert2` are Daniel R. Grayson, Michael E. Stillman, Stein A. Strømme, David Eisenbud, and Charley Crissman
    while `Chow` is due to Manfred Lehn and Christoph Sorger. `Schubert3` as well as the `Singular` library `schubert.lib` is due
    to Dang Tuan Hiep. All this work, including ours, is inspired by the `Maple` package `Schubert` written
    by Sheldon Katz and Stein A. Strømme.

The starting point of the original `Schubert` package was the problem of enumerating twisted cubic curves
on a general quintic hypersuface in $\mathbb P^4$, see [ES02](@cite). We quote from that paper:

> One way to approach enumerative problems is to find a suitable complete parameter space for the objects that one wants to count, and express the locus of objects satisfying given conditions as a certain zero-cycle on the parameter space. For this method to yield an explicit numerical answer, one needs in particular to be able to evaluate the degree of a given zerodimensional cycle class [(integration)]. This is possible in principle whenever the numerical intersection ring (cycles modulo numerical equivalence) of the parameter space is known, say in terms of generators and relations.


!!! note
    Following the authors of `Schubert`, we work with cycles modulo numerical equivalence rather than
    rational equivalence. Nevertheless, abusing our notation, we refer to the resulting intersection rings as
	Chow rings. These rings are graded by the codimension of cycles.

As in `Schubert`, our approach is abstract in the sense that we do not work with explicit
varieties given by equations. Instead, we represent a variety by specifying its
dimension together with its Chow ring and, possibly, further data. We refer to such
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

For the Chow rings of abstract flag bundles see 
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
