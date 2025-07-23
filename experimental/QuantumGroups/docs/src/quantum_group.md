```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Constructors
```@docs
quantum_group(::Symbol, ::Int)
quantum_group(::QuantumField, ::Symbol, ::Int)
```

# Generators
```@docs
gens(::QuantumGroup)
negative_chevalley_gens(::QuantumGroup)
```

# Monomials
```@docs
monomial(::QuantumGroup, ::Vector{Int} ::Vector{Int})
```

# Canonical basis
```@docs
canonical_basis_elem(::QuantumGroup, ::Vector{Int})
canonical_basis_expansion(::QuantumGroupElem)
```
