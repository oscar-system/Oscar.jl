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

Building on the previous section,  our focus is now on GR-algebras. We recall their definition:

**Definition.** ``\;`` A *GR-algebra* is the quotient $A/J$ of a $G$-algebra $A$ modulo a two-sided ideal $J \subset A$.


**Example.** ``\;`` The *$n$-th exterior algebra over $K$* is the quotient of the $G$-algebra

$A=K \langle e_1,\dots, e_n \mid e_ie_j = - e_je_i \ \text { for }\ i\neq j\rangle$

modulo the ideal

$\langle e_1^2,\dots, e_n^2\rangle.$

In analogy to the affine algebras section in the commutative algebra chapter, we describe OSCAR
functionality for dealing with GR-algebras.
