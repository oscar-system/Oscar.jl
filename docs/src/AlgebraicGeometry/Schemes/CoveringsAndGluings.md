```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# Coverings

`Covering`s are the backbone data structure for `CoveredScheme`s in Oscar. 
```@docs
Covering
```

## Constructors
```@docs
Covering(patches::Vector{<:AbsAffineScheme})
disjoint_union(C1::Covering, C2::Covering)
```

## Attributes
```@docs
affine_charts(C::Covering)
gluings(C::Covering)
```

## Methods
```@docs
add_gluing!(C::Covering, G::AbsGluing)
```

# Gluings

Gluings are used to identify open subsets $U \subset X$ and $V \subset Y$ 
of affine schemes along an isomorphism $f \colon U \leftrightarrow V \colon g$. 

## Types 
The abstract type of any such gluing is 
```@docs
AbsGluing
```
The available concrete types are 
```@docs
Gluing
SimpleGluing
```

## Constructors
```@docs
Gluing(X::AbsAffineScheme, Y::AbsAffineScheme, f::SchemeMor, g::SchemeMor)
```

## Attributes
```@docs
patches(G::AbsGluing)
gluing_domains(G::AbsGluing)
gluing_morphisms(G::AbsGluing)
inverse(G::AbsGluing)
```

## Methods
```@docs
compose(G::AbsGluing, H::AbsGluing)
maximal_extension(G::Gluing)
restrict(G::AbsGluing, f::AbsAffineSchemeMor, g::AbsAffineSchemeMor; check::Bool=true)
```



