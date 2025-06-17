```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

This represents the field of rational functions over $\mathbb{Q}$. By default the indeterminate is labelled `q`.

```@docs
quantum_field()
quantum_field(::Symbol)
```

```@docs
q_integer(::Int, ::QuantumFieldElem)
q_factorial(::Int, ::QuantumFieldElem)
q_binomial(::Int, ::Int, ::QuantumFieldElem)
```
