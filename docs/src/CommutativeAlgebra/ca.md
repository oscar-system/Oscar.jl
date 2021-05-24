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

The commutative algebra part of OSCAR provides functionality for
ideals and modules over multivariate polynomial rings or quotient
rings of these. This part is under development with regard to providing
both the functionality and the documentation. At current state, the
documentation serves in particular as a guide on what to do and
will be constantly updated.

General textbooks offering details on theory and algorithms include: 
- [GP07](@cite)
- [DL06](@cite)
- [DP13](@cite)

!!! note
    Most of the commutative algebra functionality provided by OSCAR
    relies on Gröbner basis methods. If not indicated otherwise, these
    methods are implemented for both polynomial rings over fields
	and  polynomial rings over the integers.
    ungraded or graded, GB, fields or ZZ
