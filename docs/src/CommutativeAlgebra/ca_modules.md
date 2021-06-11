```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["ca_modules.md"]
```

# Modules Over Polynomial Rings

Modules over a multivariate polynomial ring are implemented as `subquotients`.
That is, they are  submodules of quotients of free modules. Explicitly, the subquotient
$M$ with generators $A$ and relations $B$ over the ring $R$ is defined to be the module
$M = (\text{im } A + \text{im }  B)/\text{im }  B$, where $A: R^m\rightarrow R^p$ and $B: R^m\rightarrow R^p$
are two matrices representing maps of free $R$-modules with the same codomain.

!!! note
    Functionality and docu are under construction.


