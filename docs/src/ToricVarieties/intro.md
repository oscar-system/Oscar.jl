```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["intro.md"]
```


# Introduction

## Content

The toric geometry part of OSCAR comprises algorithms addressing normal toric varieties
and objects from commutative algebra and polyhedral geometry derived thereof. In particular,
we provide support for the following:
- torus-invariant divisor (classes),
- line bundles,
- line bundle cohomology via `cohomCalg` (cf. [cohomCalg:Implementation](@cite)),
- vanishing sets of line bundle cohomology (cf. `Appendix B` of [Bie18](@cite)),
- cohomology ring and cohomology classes,
- Chow ring, algebraic cycles and intersection theory,
- elementary support for closed subvarieties.


## Status

This project is work-in-progress.


## Tutorial

We provde a [tutorial for toric geometry in OSCAR ](https://nbviewer.org/github/oscar-system/oscar-website/blob/gh-pages/tutorials/ToricGeometryInOSCAR.ipynb).


## Long term goals

We follow [CLS11](@cite). Our long term goals include the following:
- Ensure that one can perform all computations of `Appendix B` in [CLS11](@cite).
- Provide support for coherent sheaves and their sheaf cohomologies. In particular, the existing algorithms in [ToricVarieties_project](https://github.com/homalg-project/ToricVarieties_project) (based on [Bie18](@cite)) should eventually be available in OSCAR.


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Martin Bies](https://martinbies.github.io/),
* [Lars Kastner](https://lkastner.github.io/).
