```@meta
CurrentModule = Oscar
```

```@contents
Pages = ["ToricLineBundles.md"]
```


# Toric Line Bundles


## Constructors

### Generic constructors

```@docs
ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{fmpz})
ToricLineBundle(v::AbstractNormalToricVariety, c::Vector{Int})
ToricLineBundle(v::AbstractNormalToricVariety, d::ToricDivisor)
```

### Tensor products

```@docs
Base.:*(l1::ToricLineBundle, l2::ToricLineBundle)
Base.:^(l::ToricLineBundle, p::fmpz)
Base.:inv(l::ToricLineBundle)
```

### Equality

```@docs
Base.:(==)(l1::ToricLineBundle, l2::ToricLineBundle)
```


## Properties

```@docs
istrivial(l::ToricLineBundle)
is_basepoint_free(l::ToricLineBundle)
isample(l::ToricLineBundle)
is_very_ample(l::ToricLineBundle)
```


## Attributes

```@docs
degree(l::ToricLineBundle)
divisor_class(l::ToricLineBundle)
toric_divisor(l::ToricLineBundle)
toric_variety(l::ToricLineBundle)
```

## Method

We use [cohomCalg](https://github.com/BenjaminJurke/cohomCalg)
to compute the dimension of line bundle cohomologies. This is achieved with the following methods.

```@docs
all_cohomologies(l::ToricLineBundle)
cohomology(l::ToricLineBundle, i::Int)
```

In addition, we can also compute a basis of the global sections
from characters.

```@docs
basis_of_global_sections(l::ToricLineBundle)
```
