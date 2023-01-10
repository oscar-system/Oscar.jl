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
```

## Properties
An `AbsCoveredScheme` may have different properties such as 
```
    is_empty(X::AbsCoveredScheme)
    is_smooth(X::AbsCoveredScheme)
```

## The modeling of covered schemes and their expected behaviour 

Any `AbsCoveredScheme` may posess several `Covering`s. This is necessary for 
several reasons; for instance, a morphism $f : X \to Y$ between `AbsCoveredScheme`s 
will in general only be given on affine patches on a refinement of the `default_covering` of `X`.
The list of available `Covering`s can be obtained using 
```@docs
    coverings(X::AbsCoveredScheme)
```
Every `AbsCoveredScheme` $X$ has to be modeled using one original `default_covering` $C$, simply 
to gather the data necessary to fully describe $X$. The `affine_charts` of $X$ return the 
patches of this covering. For any refinement $D < C$, we require the following to hold: 
Every element $U$ of the `affine_charts` of $D$ is either 

  * directly an element of the `affine_charts` of $C$;
  * a `PrincipalOpenSubset` with some ancestor in the `affine_charts` of $C$; 
  * a `SimplifiedSpec` with some original in the `affine_charts` of $C$.

In all these cases, the affine subsets in the refinements form a tree and thus remember 
their origins and ambient spaces. In particular, affine patches and also their glueings can be recycled 
and reused in different coverings and the latter should be merely seen as lists pointing 
to the objects involved. 

    
