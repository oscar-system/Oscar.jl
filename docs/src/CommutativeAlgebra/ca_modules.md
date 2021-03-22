```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca_modules.md"]
```

## Modules over Multivariate Polynomial Rings

Modules over a multivariate polynomial ring are implemented as `subquotients`.
That is, they are  submodules of quotients of free modules. Explicitly, the subquotient
$M$ with generators $A$ and relations $B$ over the ring $R$ is defined to be the module
$M = (\im A + \im B)/\im B$, where $A: R^m\rightarrow R^p$ and $B: R^m\rightarrow R^p$
are two matrices representing maps of free $R$-modules with the same codomain.

In all commands, be sure that ideals and rings are considered as modules. Are free modules modules?

Other rings? E.g. Euclidean rings
