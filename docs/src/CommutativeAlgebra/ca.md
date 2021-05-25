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
- ideals and modules over polynomial rings,
- quotient rings of polynomial rings, with particular emphasis on affine $K$-algebras.

It is under development with regard to providing
both the functionality and the documentation. 

!!! note
    - The functions discussed here essentially rely on Gröbner basis techniques. If not indicated otherwise, they are  implemented for both polynomial rings over fields and polynomial rings over the integers.
    - We distinguish between multivariate polynomial rings and multivariate polynomial rings which carry a distinguished grading. For simplicity of the presentation, functions are often only described with focus on the former case, but work similarly for homogeneous ideals and graded modules in the latter case.


General textbooks offering details on theory and algorithms include: 
- [GP07](@cite)
- [DL06](@cite)
- [DP13](@cite)

