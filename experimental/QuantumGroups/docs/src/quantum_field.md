```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

The field of rational functions over ``\mathbb{Q}`` and its elements are represented by the types `QuantumField` and `QuantumFieldElem`, respectively. By default the indeterminate is labelled `q`.

The q-analogs for `QuantumFieldElem`s are defined as follows:
```math
  \begin{aligned}
    [n] &= q^{n-1} + q^{n-3} + \ldots + q^{-n+1} \\
    [n]! &= [n] [n-1] \cdots [1] \\
    \begin{bmatrix} n \\ k \end{bmatrix} &= \frac{[n]!}{[k]! [n-k]!}
  \end{aligned}
```

# Constructors

```@docs
quantum_field()
quantum_field(::Symbol)
```

# Q-analogs

```@docs
q_integer(::Int, ::QuantumFieldElem)
q_factorial(::Int, ::QuantumFieldElem)
q_binomial(::Int, ::Int, ::QuantumFieldElem)
```
