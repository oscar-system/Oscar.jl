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
- ideals and modules over multivariate polynomial rings,
- quotient rings of such rings, with particular emphasis on affine $K$-algebras.

It is under development with regard to providing both the functionality and the documentation. 

!!! note
    - Most of the functions discussed here rely on Gröbner basis techniques. They are  implemented for multivariate polynomial rings over fields (exact fields supported by OSCAR) and, if not indicated otherwise, for multivariate polynomial rings over the integers.
    - Multivariate polynomial rings may or may not carry a distinguished grading. For simplicity of the presentation, functions are often only illustrated by examples with focus on the latter case, but work similarly for homogeneous ideals and graded modules in the former case.


General textbooks offering details on theory and algorithms include: 
- [GP08](@cite)
- [DL06](@cite)
- [DP13](@cite)

