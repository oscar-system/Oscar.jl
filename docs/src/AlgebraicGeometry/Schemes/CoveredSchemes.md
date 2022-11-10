```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["CoveredSchemes.md"]
```

# Covered schemes

Oscar supports modeling abstract schemes by means of a covering by affine charts. 

## Types
The abstract type for these is:
```@docs
    AbsCoveredScheme{BaseRingType}
```
The basic concrete instance of an `AbsCoveredScheme` is:
```@docs
    CoveredScheme{BaseRingType}
```

## Constructors
You can manually construct a `CoveredScheme` from a `Covering` using 
```@docs
    CoveredScheme(C::Covering)
```
In most cases, however, you may wish for the computer to provide you with a ready-made 
`Covering` and use a more high-level constructor, such as, for instance, 
```@docs
    covered_scheme(P::ProjectiveScheme)
```

## Attributes
To access the affine charts of a `CoveredScheme` $X$ use 
```@docs
    affine_charts(X::AbsCoveredScheme)
```
Other attributes are the `base_ring` over which the scheme is defined and 
```@docs
    default_covering(X::AbsCoveredScheme)
    coverings(X::AbsCoveredScheme)
```

## Properties
An `AbsCoveredScheme` may have different properties such as 
```
    is_empty(X::AbsCoveredScheme)
    is_smooth(X::AbsCoveredScheme)
```


    
