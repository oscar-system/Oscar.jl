```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["quotients.md"]
```

# Quotients of PBW-Algebras

In analogy to the affine algebras section in the commutative algebra chapter, we describe OSCAR
functionality for dealing with quotients of PBW-algebras modulo two-sided ideals.

**Example.** ``\;`` The *$n$-th exterior algebra over $K$* is the quotient of the PBW-algebra

$A=K \langle e_1,\dots, e_n \mid e_ie_j = - e_je_i \ \text { for }\ i\neq j\rangle$

modulo the ideal

$\langle e_1^2,\dots, e_n^2\rangle.$

