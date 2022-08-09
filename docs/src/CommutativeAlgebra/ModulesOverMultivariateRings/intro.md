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

In this section, the term module will refer to a finitely presented module over a multivariate polynomial ring.
The most general way of implementing such a module in OSCAR is that of a *subquotient module*, that is,
as a submodule of a quotient of a free module. After introducing OSCAR types for free modules and
subquotient modules, together with some related functionality, we will discuss various operations
on modules, with particular emphasis on concepts from homological algebra.
