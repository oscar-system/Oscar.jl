```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Compositions

A **weak composition** of a non-negative integer $n$ is a sequence $\lambda_1,\dots,\lambda_k$ of non-negative integers $\lambda_i$ such that $n = \lambda_1 + \dots + \lambda_k$.
A **composition** of $n$ is a weak composition consisting of *positive* integers.
The $\lambda_i$ are called the **parts** of the (weak) composition.

```@docs
weak_composition
composition
```

## Generating and counting

### Unrestricted compositions
```@docs
compositions(::Oscar.IntegerUnion)
number_of_compositions(::Oscar.IntegerUnion)
```
Note that an integer $n$ has infinitely many weak compositions as one may always append zeros to the end of a given weak composition.
Without restrictions on the number of parts, we can hence only generate compositions, but not weak compositions.

### Restricted compositions
```@docs
compositions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
number_of_compositions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
```
### Restricted weak compositions
```@docs
weak_compositions
number_of_weak_compositions(::Oscar.IntegerUnion, ::Oscar.IntegerUnion)
```

### Ascending compositions
```@docs
ascending_compositions
```
The number of ascending compositions of $n$ coincides with the number of partitions of $n$.
