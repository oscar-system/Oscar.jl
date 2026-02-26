```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Blow-ups

The *blow-up* of a smooth variety $X$ along a smooth subvariety $Z$ replaces $Z$ by the
exceptional divisor $E = \mathbb P(N_{Z/X})$, the projectivization of the normal bundle of
$Z$ in $X$. The resulting variety $\widetilde{X}$ comes equipped with a map $\pi\colon \widetilde{X}\to X$
that is an isomorphism away from $E$.

In OSCAR, a blow-up is constructed from an inclusion map $i\colon Z \hookrightarrow X$ (an
`AbstractVarietyMap` with `inclusion = true`). The function `blowup` returns the map
$\pi\colon \widetilde{X}\to X$. The exceptional divisor $E$ and the
tautological line bundle $\mathcal{O}_E(-1)$ are accessible via the attributes of the
returned map.

```@docs
blowup(i::AbstractVarietyMap; symbol::String = "e")
```

```@docs
blowup_points(X::AbstractVariety, n::Int; symbol::String = "e")
```
