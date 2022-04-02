```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

The commutative algebra part of OSCAR provides functionality for handling

- ideals of multivariate polynomial rings,
- quotients of multivariate polynomial rings modulo ideals, as well as ideals of such quotients, and 
- modules over the above rings.

In describing this functionality, we will refer to quotients of multivariate polynomial rings also as *affine algebras*.

!!! note
    Most functions discussed here rely on Gröbner basis techniques. They either execute corresponding Singular functionality, or are written as pieces of OSCAR code which rely on Singular for Gröbner basis computations.
    The functions are  implemented for multivariate polynomial rings over fields (exact fields supported by OSCAR) and, if not indicated otherwise, for multivariate polynomial rings over the integers.

!!! note
     Multivariate polynomial rings may or may not carry a distinguished (multi)grading. If not mentioned otherwise, the discussed functions apply to both the ungraded and graded case.
     For simplicity of the presentation, such functions are often only illustrated by examples with focus on the former case, but work similarly for homogeneous ideals and graded modules in the former case.

General textbooks offering details on theory and algorithms include: 
- [GP08](@cite)
- [DL06](@cite)
- [DP13](@cite)

