# JToric -- toric geometry in Julia

```@contents
```

```@meta
CurrentModule = JToric
```

## Goal

The goal of the package `JToric` is to make computations on toric geometry accessible in the [OSCAR computer algebra system](https://oscar.computeralgebra.de/). As such, we follow [the book](https://www.ams.org/publications/authors/books/postpub/gsm-124) `Toric Varieties` by `David A. Cox`, `John B. Little` and `Henry K. Schenck`]. Currently, this project is work-in-progress. Our goal is to ensure that one can perform all computations of `Appendix B`.


### Generalities

```@docs
version
```

## Toric Varieties

### Constructors

```@docs
NormalToricVariety
```

```@docs
projective_space
```

### Conversion among GAP and Polymake toric varieties

```@docs
ntv_gap2polymake
ntv_polymake2gap
```

### Properties of toric varieties

```@docs
is_normal_variety
is_affine
is_projective
is_smooth
is_complete
has_torusfactor
is_orbifold
is_simplicial
is_isomorphic_to_projective_space
is_direct_product_of_projective_spaces
```
