```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["solving.md"]
```

# Multivariate Polynomial System Solvers

```@docs
real_solutions(I::MPolyIdeal, <keyword arguments>)
rational_solutions(I::Ideal{T} where T <: MPolyElem, <keyword arguments>)
```

