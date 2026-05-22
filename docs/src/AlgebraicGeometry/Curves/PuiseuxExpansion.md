```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Puiseux expansion 
For a (possibly singular) plane curve through the origin over `QQ` we have
```@docs
puiseux_expansion(
    f::MPolyRingElem{T}, 
    max_ord::Int=10;
    precision::Int=max_ord
  ) where {T <: QQFieldElem}
```
