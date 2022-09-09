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

Working over a field $K$, our focus in this chapter is on noncommutative Gröbner bases and their application
to the computational study of finitely presented associative $K$-algebras. At present state, OSCAR offers
- a comprehensive toolkit for dealing with PBW-algebras and their quotients modulo two-sided ideals,
- functionality for computing and applying (partial) two-sided Gröbner bases in free associative algebras on fiitely many letters.

!!! note
    In contrast to the general case of finitely presented associative algebras, (left, right, two-sided) ideals in PBW-algebras
    admit finite  (left, right, two-sided) Gröbner bases. In particular, PBW-algebras are Noetherian.

!!! note
    The class of PBW-algebras includes the Weyl algebras. Algebras which arise as quotients of PBW-algebras
    include the Clifford algebras (in particular, the exterior algebras).

The textbooks
- [GP08](@cite)
- [DL06](@cite)
- [BGV03](@cite)
and the thesis
- [Lev05](@cite)
offer details on theory and algorithms.
