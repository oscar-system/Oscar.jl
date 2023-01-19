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

Our focus in this section is on finitely presented modules over rings from the following list:
- multivariate polynomial rings (OSCAR type `MPolyRing`),
- quotients of multivariate polynomial rings  (OSCAR type `MPolyQuo`), and
- localizations of the above rings (OSCAR types `MPolyLocalizedRing`, `MPolyQuoLocalizedRing`).
Hence, the word *module* will always refer to a finitely presented $R$-module, where $R$ is
one of the above types of a ring. The most general way of implementing such a module in OSCAR is
that of a *subquotient module*, that is, as a submodule of a quotient of a free $R$-module (of
finite rank). Offering a sparse way of implementing free $R$-modules, the OSCAR type `FreeMod` 
provides the basis for implementing all modules discussed here.

In what follows, after introducing types and constructors for free modules and subquotient modules,
together with some related functionality, we discuss various operations on modules, with particular emphasis
on concepts from homological algebra.
