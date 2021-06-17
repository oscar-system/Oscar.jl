```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca.md"]
```

# Introduction

The commutative algebra part of OSCAR provides functionality for handling
- ideals and modules over multivariated polynomial rings,
- quotient rings of such rings, with particular emphasis on affine $K$-algebras.

It is under development with regard to providing both the functionality and the documentation. 

!!! note
    - The functions discussed here essentially rely on Gröbner basis techniques. If not indicated otherwise, they are  implemented for multivariate polynomial rings over fields and over the integers.
    - Multivariate polynomial rings may or may not carry a distinguished grading. For simplicity of the presentation, functions are often only illustrated by examples with focus on the latter case, but work similarly for homogeneous ideals and graded modules in the former case.


General textbooks offering details on theory and algorithms include: 
- [GP07](@cite)
- [DL06](@cite)
- [DP13](@cite)

