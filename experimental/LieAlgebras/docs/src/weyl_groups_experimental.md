```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Weyl groups (experimental features)

This page is an addition to the documentation of [Weyl groups](@ref) with the additional experimental features.


## Conversion to other group types

For many computations, it may be suitable to have a `WeylGroup` as a different kind of group object, to e.g. use functionality that is only available for that other type.

The conversion functions come in pairs: one only creates an isomorphic group object, the other also computes the isomorphism.

```@docs
fp_group(::WeylGroup)
isomorphism(::Type{FPGroup}, ::WeylGroup)
```
