```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["CoveredSchemes.md"]
```

# Coverings

`Covering`s are the backbone data structure for `CoveredScheme`s in Oscar. 
```@docs
    Covering
```

## Constructors
```@docs
    Covering(patches::Vector{<:AbsSpec})
    disjoint_union(C1::Covering, C2::Covering)
```

## Attributes
```@docs
    affine_charts(C::Covering)
    glueings(C::Covering)
```

## Methods
```@docs
    add_glueing!(C::Covering, G::AbsGlueing)
```

# Glueings

Glueings are used to identify open subsets $U \subset X$ and $V \subset Y$ 
of affine schemes along an isomorphism $f \colon U \leftrightarrow V \colon g$. 

## Types 
The abstract type of any such glueing is 
```@docs
    AbsGlueing
```
The available concrete types are 
```@docs
    Glueing
    SimpleGlueing
```

## Constructors
```@docs
    Glueing(X::AbsSpec, Y::AbsSpec, f::SchemeMor, g::SchemeMor)
```

## Attributes
```@docs
    patches(G::AbsGlueing)
    glueing_domains(G::AbsGlueing)
    glueing_morphisms(G::AbsGlueing)
    inverse(G::AbsGlueing)
```

## Methods
```@docs
    compose(G::AbsGlueing, H::AbsGlueing)
    maximal_extension(G::Glueing)
    restrict(G::AbsGlueing, f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)
```



