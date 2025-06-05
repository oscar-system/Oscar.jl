```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Constructors
```@docs
quantum_group(::Symbol, ::Int)
```

# Generators
```@docs
gens(::QuantumGroup)
negative_chevalley_gens(::QuantumGroup)
```

# Canonical basis
```@docs
canonical_basis_elem(::QuantumGroup, ::Vector{Int})
canonical_basis_expansion(::QuantumGroupElem)
```
