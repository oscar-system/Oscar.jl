```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Algebraic closure of finite prime fields

It is sometimes useful to consider various finite fields in a fixed
characteristic at the same time, together with natural embeddings
between these fields.
The fields returned by [`abelian_closure`](@ref) are intended for that
purpose.

```@docs
algebraic_closure(F::T) where T <: FinField
ext_of_degree
```
